from collections import Counter
import numpy as np
import pickle
import sympy
import utils
import polynomials

LOG = utils.LOG
POINCARE_POLY_VAR = sympy.Symbol('t')


class Differential(object):
    """Differential of a DGA element.
    """

    def __init__(self, expression, coeff_mod=2):
        self.expression = expression
        self.coeff_mod = coeff_mod

    def __repr__(self):
        return str(self.expression)

    def is_linear(self):
        return polynomials.is_linear(self.expression)

    def linearize(self, subs):
        """Return linearized differential"""
        return self.bilinearize(subs, subs)

    def bilinearize(self, subs_1, subs_2):
        """Return bilinearized differential using substitutions"""
        poly = sympy.sympify(self.expression)
        output = 0
        if not poly.is_number:
            # override use of defaultdict
            coeff_dict = dict(poly.as_coefficients_dict())
            output_summands = []
            for monom, coeff in coeff_dict.items():
                if monom.is_symbol:
                    output_summands.append(coeff * monom)
                elif not monom.is_number:
                    args = monom.args
                    for i in range(len(args)):
                        factor_1 = utils.prod([s.subs(subs_1) for s in args[:i]])
                        linear_term = args[i]
                        factor_2 = utils.prod([s.subs(subs_2) for s in args[i+1:]])
                        summand = coeff * factor_1 * linear_term * factor_2
                        output_summands.append(summand)
            output = sum(output_summands)
        return Differential(expression=output, coeff_mod=self.coeff_mod)


class Matrix(object):
    """Attributes:

    lazy_ref: If False, compute row echelon form upon initialization
    values: np.array
    n_rows:
    n_cols:
    coeff_mod: Consider as matrix with coeffs modulus this number
    ref: Row echelon form
    ref_q: Matrix such that ref_q * ref = self
    ref_q_inv: Inverse of ref_q
    rank_im: Dimension of the image
    rank_ker: Dimension of the kernel
    """

    def __init__(self, values, coeff_mod=0, lazy_ref=False):
        self.lazy_ref = lazy_ref
        self.values = np.array(values)

        wrong_shape = len(self.values.shape) != 2
        if wrong_shape:
            self._init_wrong_shape()
            return

        self.n_rows = self.values.shape[0]
        self.n_cols = self.values.shape[1]
        self.coeff_mod = coeff_mod
        if coeff_mod != 0:
            self.values %= self.coeff_mod
        self.ref = None
        self.ref_q = None
        self.ref_q_inv = None
        self.rank_im = None
        self.rank_ker = None
        if not lazy_ref:
            self.set_row_echelon()

    def __repr__(self):
        return str(self.values)

    def _init_wrong_shape(self):
        if self.values.shape[0] == 0:
            self.n_rows = 0
            self.n_cols = 0
            self.ref = self
            self.ref_q = self
            self.ref_q_inv = self
            self.rank_im = 0
            self.rank_ker = 0
            return
        else:
            raise ValueError(f'Matrix has unfriendly shape {self.values.shape}')

    def multiply(self, other):
        """Multiply self with another Matrix. Must have appropriate dimensions and have the same coeff modulus.

        :param other: Other instance of Matrix
        :return: Matrix
        """
        if not isinstance(other, Matrix):
            raise ValueError(f"Trying to multiply matrix with {other.__class__}")
        if not self.coeff_mod == other.coeff_mod:
            raise ValueError(f"Trying to multiply matrix having coeff_modulus={self.coeff_mod} "
                             f"with other having coeff_modulus={other.coeff_mod}")
        values = np.matmul(self.values, other.values)
        return Matrix(values=values, coeff_mod=self.coeff_mod, lazy_ref=self.lazy_ref)

    def set_row_echelon(self):
        """Set ref, ref_q, ref_q_inv, rank_im, rank_ker attributes for self.

        :return: None
        """
        ref = self.values.copy()
        q = np.identity(self.n_rows)
        q_inv = np.identity(self.n_rows)
        k = 0
        l = 0
        while k < self.n_rows:
            while l < self.n_cols and not (np.any(ref[k:, l])):
                l += 1
            if l == self.n_cols:
                break
            ref, q, q_inv = self._row_reduce(ref, q, q_inv, k, l)
            k += 1
        self.ref = Matrix(ref, coeff_mod=self.coeff_mod, lazy_ref=True)
        self.ref_q = Matrix(q, coeff_mod=self.coeff_mod, lazy_ref=True)
        self.ref_q_inv = Matrix(q_inv, coeff_mod=self.coeff_mod, lazy_ref=True)
        self.rank_im = k
        self.rank_ker = self.n_cols - k

    def is_square(self):
        return self.n_rows == self.n_cols

    def _row_reduce(self, ref, q, q_inv, k, l):
        while np.any(ref[k + 1:, l]):
            ref, q, q_inv = self._row_prepare(ref, q, q_inv, k, l)
            ref, q, q_inv = self._partial_row_reduce(ref, q, q_inv, k, l)
        return ref, q, q_inv

    def _row_prepare(self, ref, q, q_inv, k, l):
        (a, i) = self._smallest_nonzero_index(ref[:, l], k)
        ref[[i, k], :] = ref[[k, i], :]
        q_inv[[i, k], :] = q_inv[[k, i], :]  # row swap
        q[:, [i, k]] = q[:, [k, i]]  # column swap
        return ref, q, q_inv

    @staticmethod
    def _smallest_nonzero_index(v, k):
        # TODO: This pattern is bad
        try:
            alpha = min(abs(v[k:][np.nonzero(v[k:])]))
            i = min(i for i in range(k, len(v)) if abs(v[i]) == alpha)
            return alpha, i
        except:
            return 0, np.nan  # minNonZero causes this sometimes

    def _partial_row_reduce(self, ref, q, q_inv, k, l):
        for i in range(k + 1, q.shape[0]):
            c = (ref[i, l] // ref[k, l])
            if self.coeff_mod != 0:
                c %= self.coeff_mod
            # row add i,k,-q
            ref[i] += (-c * ref[k])
            q_inv[i] += (-c * q_inv[k])  # row add
            q[:, k] += (c * q[:, i])  # column add (note i,k are switched)
            if self.coeff_mod != 0:
                ref[i] %= self.coeff_mod
                q_inv[i] %= self.coeff_mod
                q[:, k] %= self.coeff_mod
        return ref, q, q_inv


class LinearMap(object):
    """Interface to matrices for linear polynomial expressions"""

    def __init__(self, coeff_dict, range_symbols, coeff_mod=2):
        """
        :param coeff_dict: Dict of the form {symbol: sympy expression}
        :param coeff_mod: Coefficient modulus to use
        """
        self.coeff_dict = coeff_dict
        self.range_symbols = range_symbols
        self.coeff_mod = coeff_mod
        self._set_dims()
        self._set_input_vars()
        self._set_matrix()

    def __repr__(self):
        return str(self.coeff_dict)

    def _set_dims(self):
        self.dim_dom = len(self.coeff_dict.keys())
        self.dim_range = len(self.range_symbols)

    def _set_input_vars(self):
        """Encode input and output symbols as `onehot` vectors so that they could be turned into a matrix"""
        self.input_vars = list(self.coeff_dict.keys())

    def _set_matrix(self):
        # The dtype should be essential to avoid float arithmetic errors
        values = np.zeros((len(self.range_symbols), len(self.input_vars)), dtype=np.integer)
        for i in range(len(self.input_vars)):
            for j in range(len(self.range_symbols)):
                temp_subs = {s: 0 if s != self.range_symbols[j] else 1 for s in self.range_symbols}
                value = sympy.sympify(self.coeff_dict[self.input_vars[i]]).subs(temp_subs)
                values[j][i] = value
        self.matrix = Matrix(values=values, coeff_mod=self.coeff_mod)


class DGBase(object):
    """Base class for differential graded objects"""

    def __init__(self, gradings, differentials, filtration_levels=None, coeff_mod=0, grading_mod=0):
        """
        :param symbols: Dictionary mapping generator -> grading.
        :param differentials: Dictionary mapping generator -> polynomial.
        :param coeff_mod: m when using Z/mZ coeffs. Must be zero or prime.
        :param grading_mod: m when using Z/mZ grading.
        """
        self.gradings = gradings
        self.symbols = self.gradings.keys()
        self.differentials = differentials
        self.filtration_levels = filtration_levels
        self.coeff_mod = coeff_mod
        self.grading_mod = grading_mod
        self._verify_init_args()
        self._correct_gradings()
        self._correct_filtration_levels()

    def reduced(self, ceoff_mod, grading_mod):
        if self.coeff_mod % ceoff_mod != 0:
            raise ValueError(f"New coefficient modulus {ceoff_mod} doesn't divide old modulus {self.coeff_mod}")
        if self.grading_mod % grading_mod != 0:
            raise ValueError(f"New grading modulus {grading_mod} doesn't divide old modulus {self.grading_mod}")
        return self.__init__(
            gradings=self.gradings,
            differentials=self.differentials,
            coeff_mod=ceoff_mod,
            grading_mod=grading_mod
        )

    def _correct_filtration_levels(self):
        if self.filtration_levels is None:
            self.filtration_levels = {k: 1 for k in self.gradings.keys()}
            self.max_filtration_level = 1
            return
        self.max_filtration_level = max(self.filtration_levels.values())

    def _verify_init_args(self):
        if self.coeff_mod != 0:
            if not sympy.isprime(self.coeff_mod):
                raise ValueError(f"Coefficient modulus {self.coeff_mod} is not 0 or prime")
        if self.symbols != self.differentials.keys():
            raise ValueError("generators don't match in symbols and keys")
        for g, d in self.differentials.items():
            if not isinstance(d, Differential):
                raise ValueError(f"Differential for {g} is not instance of class Differential.")

    def _correct_gradings(self):
        if self.grading_mod != 0:
            corrected_gradings = {g: self.gradings[g] % self.grading_mod for g in self.gradings.keys()}
            self.gradings = corrected_gradings


class ChainComplex(DGBase):

    def __init__(self, gradings, differentials, filtration_levels=None, coeff_mod=0, grading_mod=0):
        super(ChainComplex, self).__init__(
            gradings=gradings,
            differentials=differentials,
            filtration_levels=filtration_levels,
            coeff_mod=coeff_mod,
            grading_mod=grading_mod
        )
        self._verify_linear_differentials()
        self._set_linear_maps()
        self._set_poincare_poly()

    def _verify_linear_differentials(self):
        for k, v in self.differentials.items():
            if not v.is_linear():
                raise ValueError(f"Trying to instantiate ChainComplex with non-linear differential {k}: {v}")

    def _set_linear_maps(self):
        # We have to use a dictionary because the degrees may be non-positive
        self.linear_maps = dict()
        for i in self.gradings.values():
            i_min_1 = i - 1
            if self.grading_mod != 0:
                i_min_1 %= self.grading_mod
            range_symbols = [s for s in self.symbols if self.gradings[s] == i_min_1]
            domain_symbols = [s for s in self.symbols if self.gradings[s] == i]
            diffs = {s: self.differentials[s].expression for s in domain_symbols}
            self.linear_maps[i] = LinearMap(
                coeff_dict=diffs,
                range_symbols=range_symbols,
                coeff_mod=self.coeff_mod
            )

    def _set_poincare_poly(self):
        output_dict = dict()
        output_dual_dict = dict()
        for i in self.linear_maps.keys():
            i_min_1 = i - 1
            if self.grading_mod != 0:
                i_min_1 %= self.grading_mod
            linear_map = self.linear_maps[i]
            rank_im = linear_map.matrix.rank_im
            rank_im_dual = rank_im
            rank_ker = linear_map.matrix.rank_ker
            rank_ker_dual = linear_map.dim_range - rank_im_dual
            # Sum up bettis
            if i in output_dict.keys():
                output_dict[i] += rank_ker
            else:
                output_dict[i] = rank_ker
            if i_min_1 in output_dict.keys():
                output_dict[i_min_1] -= rank_im
            else:
                output_dict[i_min_1] = -rank_im
            # Sum up dual bettis
            if i_min_1 in output_dual_dict.keys():
                output_dual_dict[i_min_1] += rank_ker_dual
            else:
                output_dual_dict[i_min_1] = rank_ker_dual
            if i in output_dual_dict.keys():
                output_dual_dict[i] -= rank_im_dual
            else:
                output_dual_dict[i] = -rank_im_dual
        # poincare poly
        output = 0
        for i in output_dict.keys():
            output += output_dict[i] * (POINCARE_POLY_VAR ** i)
        # dual poly
        output_dual = 0
        for i in output_dual_dict.keys():
            output_dual += output_dual_dict[i] * (POINCARE_POLY_VAR ** i)
        self.poincare_poly = output
        self.poincare_poly_dual = output_dual


class DGA(DGBase):

    def __init__(self, gradings, differentials, filtration_levels=None,
                 coeff_mod=0, grading_mod=0, lazy_aug_data=False, lazy_augs=False, lazy_bilin=False):
        super(DGA, self).__init__(
            gradings=gradings,
            differentials=differentials,
            filtration_levels=filtration_levels,
            coeff_mod=coeff_mod,
            grading_mod=grading_mod
        )
        self.augmentations = None
        self.augmentations_compressed = None
        self.n_augs = None
        self.n_augs_compressed = None
        self.bilin_polys = dict()
        self.bilin_polys_dual = dict()
        self.lin_poly_list = None
        self.bilin_poly_list = None
        self.lazy_aug_data = lazy_aug_data
        self.lazy_augs = lazy_augs
        self.lazy_bilin = lazy_bilin
        if lazy_augs and not lazy_bilin:
            raise ValueError("DGA cannot autocompute bilinearized polys without autocomputing augmentations")

        # Variables needed to setup augs
        self._aug_vars_set = False
        self.aug_comm_symbols_to_symbols = None
        self.aug_symbols_to_comm_symbols = None
        self.aug_analysis_table = None
        self.aug_polys = None
        if not lazy_aug_data:
            self.set_aug_data()
        if not lazy_augs:
            self.set_augmentations()
            if not lazy_bilin:
                self.set_all_bilin()

    def pickle(self, file_name, compress=True):
        storage = dict()
        # Required for init
        storage['gradings'] = self.gradings
        storage['differentials'] = self.differentials
        storage['filtration_levels'] = self.filtration_levels
        storage['bilin_polys'] = self.bilin_polys
        storage['lin_poly_list'] = self.lin_poly_list
        storage['bilin_poly_list'] = self.bilin_poly_list
        storage['coeff_mod'] = self.coeff_mod
        storage['grading_mod'] = self.grading_mod
        # Data we want to set up
        storage['_aug_vars_set'] = self._aug_vars_set
        storage['aug_comm_symbols_to_symbols'] = self.aug_comm_symbols_to_symbols
        storage['aug_symbols_to_comm_symbols'] = self.aug_symbols_to_comm_symbols
        storage['aug_analysis_table'] = self.aug_analysis_table
        storage['aug_polys'] = self.aug_polys
        storage['augmentations_compressed'] = self.augmentations_compressed
        if compress:
            storage['augmentations'] = None
        else:
            storage['augmentations'] = self.augmentations
        with open(file_name, 'wb') as f:
            pickle.dump(storage, f)

    @classmethod
    def from_pickle(cls, file_name, decompress=False):
        with open(file_name, 'rb') as f:
            storage = pickle.load(f)
        dga = DGA(
            gradings=storage['gradings'],
            differentials=storage['differentials'],
            filtration_levels=storage['filtration_levels'],
            coeff_mod=storage['coeff_mod'],
            grading_mod=storage['grading_mod'],
            lazy_augs=True,
            lazy_bilin=True,
            lazy_aug_data=True
        )
        dga._aug_vars_set = storage.get('_aug_vars_set')
        dga.aug_comm_symbols_to_symbols = storage.get('aug_comm_symbols_to_symbols')
        dga.aug_symbols_to_comm_symbols = storage.get('aug_symbols_to_comm_symbols')
        dga.aug_analysis_table = storage.get('aug_analysis_table')
        dga.aug_polys = storage.get('aug_polys')
        dga.bilin_polys = storage.get('bilin_polys')
        dga.lin_poly_list = storage.get('lin_poly_list')
        dga.bilin_poly_list = storage.get('bilin_poly_list')
        dga.augmentations_compressed = storage.get('augmentations_compressed')
        dga.augmentations = storage.get('augmentations')
        if decompress and dga.augmentations is None:
            dga.decompress_augmentations()
        return dga

    def get_verbose_subs_from_aug(self, aug):
        """Extend aug to all generators. This means setting non-zero graded generators to zero.

        :param aug: dict
        :return: dict
        """
        for s in self.symbols:
            if s not in aug.keys():
                aug[s] = 0
        return aug

    def get_bilin_differential(self, aug_1, aug_2):
        subs_1 = self.get_verbose_subs_from_aug(aug_1)
        subs_2 = self.get_verbose_subs_from_aug(aug_2)
        return {
            g: v.bilinearize(subs_1, subs_2)
            for g, v in self.differentials.items()
        }

    @utils.log_start_stop
    def set_all_bilin(self):
        # How frequently to log?
        bilin_counter = 0
        for i in range(self.n_augs):
            for j in range(self.n_augs):
                aug_1 = self.augmentations[i]
                aug_2 = self.augmentations[j]
                bilin_diff = self.get_bilin_differential(aug_1, aug_2)
                cx = ChainComplex(
                    gradings=self.gradings,
                    differentials=bilin_diff,
                    grading_mod=self.grading_mod,
                    coeff_mod=self.coeff_mod
                )
                self.bilin_polys[(i,j)] = cx.poincare_poly
                self.bilin_polys_dual[(i,j)] = cx.poincare_poly_dual
                bilin_counter += 1
                LOG.info(f"Computed {bilin_counter} of {self.n_augs ** 2} bilinearized Poincare polys so far")
        self.lin_poly_list = sorted(
            utils.unique_elements([v for k, v in self.bilin_polys.items() if k[0] == k[1]]),
            key=str
        )
        self.bilin_poly_list = sorted(
            utils.unique_elements(self.bilin_polys.values()),
            key=str
        )

    def are_homotopy_equivalent(self, aug_1, aug_2):
        """The augmentations are homotopy equivalent iff the induced map on bilinearized homology to the base field
        given by aug = (aug_1 - aug_2) is zero. How do we figure this out using only generators of the underlying vector
        spaces of generators? Since aug(del(x)) always vanishes, it is already zero on Image(del). Therefore aug == 0
        on homology iff it vanishes on ker del. Another way to check is if the lin homologies agree and if the bilin
        homologies agree with those of the lin homologies.

        This is actually annoying to do in sympy and easy to do in sage :(
        However, we can do a trick... Represent ker as the intersection of the varieties determined by the row vectors
        of a marix representation of del. Say these linear equations are p_i. Since aug != 0 iff 1 in Im(aug), we see
        that aug restricted to ker is non-zero iff there is a solution to [p_1, ..., p_k, aug + 1]. This can be worked
        out with the equation solver using existence_only=True.

        It would actually save time also to use transitivity to avoid checking all of the homotopy equiv classes.

        :param aug_1:
        :param aug_2:
        :return:
        """
        raise NotImplementedError()

    @utils.log_start_stop
    def set_aug_data(self):
        """Set variables needed to compute augmentations.

        :return: None
        """
        zero_graded_symbols = [g for g in self.symbols if self.gradings[g] == 0]
        LOG.info(f"DGA has {len(zero_graded_symbols)} generators with 0 grading")
        self.aug_symbols_to_comm_symbols = {g: sympy.Symbol(str(g), commutative=True) for g in zero_graded_symbols}
        self.aug_comm_symbols_to_symbols = {v: k for k, v in self.aug_symbols_to_comm_symbols.items()}
        # Only need to mod out by differentials of degree 1 elements
        d_expressions = {k: v.expression for k, v in self.differentials.items() if self.gradings[k] == 1}
        LOG.info(f"Differentials required to compute augs: {len(d_expressions)}")
        # Set all non-zero graded elements to zero.
        zero_substitutions = {g: 0 for g in self.symbols if self.gradings[g] != 0}
        # Make polynomials in commutative variables of the deg=0 generators.
        # This is required to apply zero_set which only works with commuting variables!
        LOG.info(f"Eliminating degree-non-zero terms from polys")
        aug_polys = {
            k: sympy.sympify(v).subs(zero_substitutions)
            for k, v in d_expressions.items()
        }
        LOG.info(f"Replacing non-commutative symbols with commutative")
        aug_polys = {
            k: v.subs(self.aug_symbols_to_comm_symbols)
            for k, v in aug_polys.items()
        }
        aug_polys = utils.unique_elements(aug_polys.values())
        if 0 in aug_polys:
            aug_polys.remove(0)
        LOG.info(f"Polynomials required to compute augs after simplification: {len(aug_polys)}")
        self.aug_polys = aug_polys
        # Create table to help inform batch size
        # TODO: This is replicated in batch_zero_set. Should make into a function
        poly_to_symbols = {p: polynomials.poly_symbols(p) for p in aug_polys}
        symbol_freq = Counter()
        for symbol_set in poly_to_symbols.values():
            symbol_freq.update(list(symbol_set))
        symbol_freq = dict(symbol_freq)
        for sym in self.aug_comm_symbols_to_symbols.keys():
            if sym not in symbol_freq.keys():
                symbol_freq[sym] = 0
        symbol_table = [
            {
                "symbol": k,
                "freq": v
            }
            for k, v in symbol_freq.items()
        ]
        symbol_table.sort(key=lambda d: d["freq"], reverse=True)
        running_symbol_set = set()
        for i in range(len(symbol_table)):
            running_symbol_set.add(symbol_table[i]["symbol"])
            running_polys = [k for k, v in poly_to_symbols.items() if v.issubset(running_symbol_set)]
            symbol_table[i]["n_cum_polys"] = len(running_polys)
            symbol_table[i]["poly_to_sym_ratio"] = 1.0*len(running_polys) / (i+1)
        self.aug_analysis_table = symbol_table
        LOG.info(f"Aug analysis table: \n" + "\n".join([str(r) for r in symbol_table]))
        self._aug_vars_set = True

    @utils.log_start_stop
    def set_augmentations(self, batch_size=None, filtered=False, decompress=True):
        # this pattern is bad, because we allow 0 coeff_mod upon instance creation and then this method always runs
        if self.coeff_mod == 0:
            raise ValueError("We cannot search for augmentations over ZZ. It's too difficult :(")
        if filtered and batch_size is not None:
            raise ValueError("Cannot do batch and filtered simultaneously")
        if not self._aug_vars_set:
            self.set_aug_data()
        zero_graded_symbols = list(self.aug_symbols_to_comm_symbols.keys())
        comm_symbols = list(self.aug_symbols_to_comm_symbols.values())
        if len(zero_graded_symbols) == 0:
            comm_augmentations = []
        else:
            if batch_size is not None:
                comm_augmentations = polynomials.batch_zero_set(
                    polys=self.aug_polys, symbols=comm_symbols, modulus=self.coeff_mod, batch_size=batch_size)
            elif filtered:
                comm_filtration_levels = {
                    v: self.filtration_levels[k]
                    for k, v in self.aug_symbols_to_comm_symbols.items()
                }
                comm_augmentations = polynomials.filtered_zero_set(
                    polys=self.aug_polys, symbols=comm_symbols,
                    filtration_levels=comm_filtration_levels, modulus=self.coeff_mod)
            else:
                comm_augmentations = polynomials.zero_set(
                    polys=self.aug_polys, symbols=comm_symbols, modulus=self.coeff_mod)
        self.augmentations_compressed = comm_augmentations
        self.n_augs_compressed = len(self.augmentations_compressed)
        LOG.info(f"Found {self.n_augs_compressed} compressed augmentations of DGA")
        if decompress:
            self.decompress_augmentations()
        else:
            self._set_n_augs_from_compressed()
            LOG.info(f"Found {self.n_augs} (uncompressed) augmentations of DGA")

    def get_decompressed_augmentations(self, start_i=None, end_i=None, fill_na=None):
        if fill_na is not None:
            LOG.info(f"Working with a subset of augmentations using default value {fill_na}")
        if start_i is None:
            start_i = 0
        if end_i is None:
            end_i = len(self.augmentations_compressed)
        augs_compressed = self.augmentations_compressed[start_i:end_i]
        augmentations = []
        comm_augmentations = polynomials.expand_zero_set_from_unset(
            augs_compressed, modulus=self.coeff_mod, fill_na=fill_na
        )
        for aug in comm_augmentations:
            augmentations.append({self.aug_comm_symbols_to_symbols[k]: v for k, v in aug.items()})
        return augmentations

    @utils.log_start_stop
    def decompress_augmentations(self, fill_na=None):
        # comm_augs will be a list of dicts whose keys are the commutative symbols.
        # We need to switch them back!
        self.augmentations = self.get_decompressed_augmentations(fill_na=fill_na)
        self.n_augs = len(self.augmentations)
        LOG.info(f"Found {self.n_augs} augmentations of DGA")

    def _set_n_augs_from_compressed(self):
        count = 0
        for aug_c in self.augmentations_compressed:
            n_unset_vars = len([k for k, v in aug_c.items() if v is polynomials.UNSET_VAR])
            count += self.coeff_mod ** n_unset_vars
        self.n_augs = count
