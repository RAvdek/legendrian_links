from collections import Counter
from math import prod
import numpy as np
import sympy
import utils
import polynomials

LOG = utils.LOG


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
        """Return bilinearized differential should use expression.args to capture monomials and then .args again to get
        individual terms"""
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
                        factor_1 = prod([s.subs(subs_1) for s in args[:i]])
                        linear_term = args[i]
                        factor_2 = prod([s.subs(subs_2) for s in args[i+1:]])
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
    ref_q: Matrix q such that q*ref = self
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
        self._set_input_vars()
        self._set_matrix()

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

    def __init__(self, gradings, differentials, coeff_mod=0, grading_mod=0):
        """
        :param symbols: Dictionary mapping generator -> grading.
        :param differentials: Dictionary mapping generator -> polynomial.
        :param coeff_mod: m when using Z/mZ coeffs. Must be zero or prime.
        :param grading_mod: m when using Z/mZ grading.
        """
        self.gradings = gradings
        self.symbols = self.gradings.keys()
        self.differentials = differentials
        self.coeff_mod = coeff_mod
        self.grading_mod = grading_mod
        self._verify_init_args()
        self._correct_gradings()

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

    def __init__(self, gradings, differentials, coeff_mod=0, grading_mod=0):
        super(ChainComplex, self).__init__(
            gradings=gradings,
            differentials=differentials,
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
        self.linear_maps = dict()
        for i in self.gradings.values():
            j = i - 1
            if self.grading_mod != 0:
                j %= self.grading_mod
            range_symbols = [s for s in self.symbols if self.gradings[s] == j]
            domain_symbols = [s for s in self.symbols if self.gradings[s] == i]
            diffs = {s: self.differentials[s].expression for s in domain_symbols}
            self.linear_maps[i] = LinearMap(
                coeff_dict=diffs,
                range_symbols=range_symbols,
                coeff_mod=self.coeff_mod
            )

    def _set_poincare_poly(self):
        output_dict = dict()
        for i in range(len(self.linear_maps)):
            j = i - 1
            if self.grading_mod != 0:
                j %= self.grading_mod
            rank_im = self.linear_maps[i].matrix.rank_im
            rank_ker = self.linear_maps[i].matrix.rank_ker
            if i in output_dict.keys():
                output_dict[i] += rank_ker
            else:
                output_dict[i] = rank_ker
            if j in output_dict.keys():
                output_dict[j] -= rank_im
            else:
                output_dict[j] = -rank_im
        t = sympy.Symbol('t')
        output = 0
        for i in output_dict.keys():
            output += output_dict[i] * (t ** i)
        self.poincare_poly = output


class DGA(DGBase):

    def __init__(self, gradings, differentials, coeff_mod=0, grading_mod=0, lazy_bilin=False):
        super(DGA, self).__init__(
            gradings=gradings,
            differentials=differentials,
            coeff_mod=coeff_mod,
            grading_mod=grading_mod
        )
        self._set_augmentations()
        self.n_augs = len(self.augmentations)
        LOG.info(f"Found {self.n_augs} augmentations of DGA")
        self.bilin_polys = [[None for _ in range(self.n_augs)] for _ in range(self.n_augs)]
        self.lazy_bilin = lazy_bilin
        if not lazy_bilin:
            self.set_all_bilin()

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
        bilin_counter = 0
        bilin_diff_storage = list()
        for i in range(self.n_augs):
            for j in range(self.n_augs):
                aug_1 = self.augmentations[i]
                aug_2 = self.augmentations[j]
                bilin_diff = self.get_bilin_differential(aug_1, aug_2)
                storage_matches = [
                    i for i in range(len(bilin_diff_storage)) if bilin_diff_storage[i]["diff"] == bilin_diff
                ]
                if len(storage_matches) > 0:
                    i = storage_matches[0]
                    poincare_poly = bilin_diff_storage[i]["poincare_poly"]
                    bilin_diff_storage[i]["augs"].append([i, j])
                else:
                    cx = ChainComplex(
                        gradings=self.gradings,
                        differentials=bilin_diff,
                        grading_mod=self.grading_mod,
                        coeff_mod=self.coeff_mod
                    )
                    poincare_poly = cx.poincare_poly
                    bilin_diff_storage.append(
                        {
                            "diff": bilin_diff,
                            "poincare_poly": poincare_poly,
                            "augs": [[i, j]]
                        }
                    )
                self.bilin_polys[i][j] = poincare_poly
                bilin_counter += 1
                if bilin_counter % 20 == 0:
                    LOG.info(f"Computed {bilin_counter} bilinearized Poincare polys so far")
                    diff_frequency = dict(Counter([len(d["augs"]) for d in bilin_diff_storage]))
                    LOG.info(f"Frequency of repetition in bilin homologies so far: {diff_frequency}")
        diff_frequency = dict(Counter([len(d["augs"]) for d in bilin_diff_storage]))
        LOG.info(f"Frequency of repetition in bilin homologies: {diff_frequency}")

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
    def _set_augmentations(self):
        # this pattern is bad, because we allow 0 coeff_mod upon instance creation and then this method always runs
        if self.coeff_mod == 0:
            raise ValueError("We cannot search for augmentations over ZZ. It's too difficult :(")
        zero_graded_symbols = [g for g in self.symbols if self.gradings[g] == 0]
        LOG.info(f"DGA has {len(zero_graded_symbols)} generators with 0 grading")
        symbols_to_comm_symbols = {g: sympy.Symbol(str(g), commutative=True) for g in zero_graded_symbols}
        comm_symbols = list(symbols_to_comm_symbols.values())
        if len(comm_symbols) == 0:
            self.augmentations = []
            return
        d_expressions = [d.expression for d in self.differentials.values()]
        # Set all non-zero graded elements to zero.
        zero_substitutions = {g: 0 for g in self.symbols if self.gradings[g] != 0}
        # Make polynomials in commutative variables of the deg=0 generators.
        # This is required to apply zero_set which only works with commuting variables!
        polys = []
        for exp in d_expressions:
            polys.append(sympy.sympify(exp).subs(zero_substitutions))
        polys = [p.subs(symbols_to_comm_symbols) for p in polys]
        polys = utils.unique_elements(polys)
        if 0 in polys:
            polys.remove(0)
        comm_augmentations = polynomials.zero_set(polys=polys, symbols=comm_symbols, modulus=self.coeff_mod)
        # comm_augs will be a list of dicts whose keys are the commutative symbols.
        # We need to switch them back!
        comm_symbols_to_symbols = {v:k for k, v in symbols_to_comm_symbols.items()}
        augmentations = []
        for aug in comm_augmentations:
            augmentations.append({comm_symbols_to_symbols[k]: v for k, v in aug.items()})
        self.augmentations = augmentations
