from math import prod
import numpy as np
import sympy
import utils
import polynomials

LOG = utils.get_logger(__name__)


class Differential(object):
    """Differential of a DGA element.

    self.summands is a list of lists. Each list stores the factors in the order in which they appear. Eg.
    1 + x + z + x*y*z will be stored as [[1], [x], [z], [x,y,z]] whereas...
    1 + x + z + z*y*x will be stored as [[1], [x], [z], [z,y,x]]
    This non-commutative storage is needed for bilinearization.
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
        raise NotImplementedError()


class Matrix(object):
    """Attributes:

    lazy_ref: If False, compute row echelon form upon initialization
    values: np.array
    n_rows:
    n_cols:
    coeff_modulus: Consider as matrix with coeffs modulus this number
    ref: Row echelon form
    ref_q: Matrix q such that q*ref = self
    ref_q_inv: Inverse of ref_q
    rank_im: Dimension of the image
    rank_ker: Dimension of the kernel
    """

    def __init__(self, values, coeff_modulus=0, lazy_ref=False):
        self.lazy_ref = lazy_ref
        self.values = np.array(values)
        self.n_rows = self.values.shape[0]
        self.n_cols = self.values.shape[1]
        self.coeff_modulus = coeff_modulus
        if coeff_modulus != 0:
            self.values %= self.coeff_modulus
        self.ref = None
        self.ref_q = None
        self.ref_q_inv = None
        self.rank_im = None
        self.rank_ker = None
        if not lazy_ref:
            self.set_row_echelon()

    def __repr__(self):
        return str(self.values)

    def multiply(self, other):
        """Multiply self with another Matrix. Must have appropriate dimensions and have the same coeff modulus.

        :param other: Other instance of Matrix
        :return: Matrix
        """
        if not isinstance(other, Matrix):
            raise ValueError(f"Trying to multiply matrix with {other.__class__}")
        if not self.coeff_modulus == other.coeff_modulus:
            raise ValueError(f"Trying to multiply matrix having coeff_modulus={self.coeff_modulus} "
                             f"with other having coeff_modulus={other.coeff_modulus}")
        values = np.matmul(self.values, other.values)
        return Matrix(values=values, coeff_modulus=self.coeff_modulus, lazy_ref=self.lazy_ref)

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
        self.ref = Matrix(ref, coeff_modulus=self.coeff_modulus, lazy_ref=True)
        self.ref_q = Matrix(q, coeff_modulus=self.coeff_modulus, lazy_ref=True)
        self.ref_q_inv = Matrix(q_inv, coeff_modulus=self.coeff_modulus, lazy_ref=True)
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
            q = (ref[i, l] // ref[k, l])
            if self.coeff_modulus != 0:
                q %= self.coeff_modulus
            # row add i,k,-q
            ref[i] += (-q * ref[k])
            q_inv[i] += (-q * q_inv[k])  # row add
            q[:, k] += (q * q[:, i])  # column add (note i,k are switched)
            if self.coeff_modulus != 0:
                ref[i] %= self.coeff_modulus
                q_inv[i] %= self.coeff_modulus
                q[:, k] %= self.coeff_modulus
        return ref, q, q_inv


class LinearMap(object):
    """Interface to matrices for linear polynomial expressions"""

    def __init__(self, coeff_dict, coeff_mod=2):
        """
        :param coeff_dict: Dict of the form {symbol: sympy expression}
        :param coeff_mod: Coefficient modulus to use
        """
        self.coeff_dict = coeff_dict
        self.coeff_mod = coeff_mod
        self._set_var_mapping()
        self._set_matrix()

    def _set_var_mapping(self):
        """Encode input and output symbols as `onehot` vectors so that they could be turned into a matrix"""

    def _set_matrix(self):
        raise NotImplementedError()


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

    def _verify_linear_differentials(self):
        for k, v in self.differentials.items():
            if not v.is_linear():
                raise ValueError(f"Trying to instantiate ChainComplex with non-linear differential {k}: {v}")

    def _set_betti_numbers(self):
        raise NotImplementedError()


class DGA(DGBase):

    def __init__(self, gradings, differentials, coeff_mod=0, grading_mod=0):
        super(DGA, self).__init__(
            gradings=gradings,
            differentials=differentials,
            coeff_mod=coeff_mod,
            grading_mod=grading_mod
        )
        self._set_augmentations()

    def get_verbose_subs_from_aug(self, aug):
        """Extend aug to all the things it must be zero on

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

    def bilinearized_poincare_poly(self, aug_1, aug_2):
        """If we assume that the base field is a prime (so that the coeff ring is a field) then we don't need to
        understand the homology groups in terms of generators and relations. We can get the ith betti number by
        b_i = dim(ker del_i) - rank del_{i-1}. Both of these should be available in sympy.

        :param aug_1:
        :param aug_2:
        :return:
        """
        raise NotImplementedError()

    def linearized_poincare_poly(self, aug):
        return self.bilinearized_poincare_poly(aug, aug)

    def _set_augmentations(self):
        # this pattern is bad, because we allow 0 coeff_mod upon instance creation and then this method always runs
        if self.coeff_mod == 0:
            raise ValueError("We cannot search for augmentations over ZZ. It's too difficult :(")
        zero_graded_symbols = [g for g in self.symbols if self.gradings[g] == 0]
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
        if 0 in polys:
            polys.remove(0)
        polys = list(set(polys))
        comm_augmentations = polynomials.zero_set(polys=polys, symbols=comm_symbols, modulus=self.coeff_mod)
        # comm_augs will be a list of dicts whose keys are the commutative symbols.
        # We need to switch them back!
        comm_symbols_to_symbols = {v:k for k, v in symbols_to_comm_symbols.items()}
        augmentations = []
        for aug in comm_augmentations:
            augmentations.append({comm_symbols_to_symbols[k]: v for k, v in aug.items()})
        self.augmentations = augmentations
