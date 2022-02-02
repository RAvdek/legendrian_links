from math import prod
import sympy
import utils


LOG = utils.get_logger(__name__)


class Differential(object):
    """Differential of a DGA element.

    self.summands is a list of lists. Each list stores the factors in the order in which they appear. Eg.
    1 + x + z + x*y*z will be stored as [[1], [x], [z], [x,y,z]] whereas...
    1 + x + z + z*y*x will be stored as [[1], [x], [z], [z,y,x]]
    This non-commutative storage is needed for bilinearization.
    """

    def __init__(self, summands=[], coeff_mod=0, signs=[]):
        self.summands = summands
        self.coeff_mod = coeff_mod
        self.signs = signs

    def __repr__(self):
        return str(self.summands)

    def to_polynomial(self):
        # When char(field) != 2, we have to worry about how the orders of the elements are multiplied, as minus signs
        # may be created for swapping orders of odd-degree generators.
        if not self.coeff_mod == 2:
            raise NotImplementedError()
        output = 0
        for i in range(len(self.summands)):
            output += self.signs[i]*prod(self.summands[i])
        return output

    def linearize(self, subs):
        """Return linearized differential"""
        return self.bilinearize(subs, subs)

    def bilinearize(self, subs_1, subs_2):
        """Return linearized differential"""
        raise NotImplementedError()

class DGA(object):

    def __init__(self, gradings, differentials, coeff_mod=0, grading_mod=0):
        """
        :param symbols: Dictionary mapping generator -> grading.
        :param differentials: Dictionary mapping generator -> polynomial.
        :param coeff_mod: m when using Z/mZ coeffs. Must be zero or prime.
        :param grading_mod: m when using Z/mZ grading.
        """
        self.symbols = self.gradings.keys()
        self.gradings = gradings
        self.differentials = differentials
        self.coeff_mod = coeff_mod
        self.grading_mod = grading_mod
        self._verify_init_args()
        self._correct_gradings()
        self._set_augmentations()

    def reduced_dga(self, ceoff_mod, grading_mod):
        if self.coeff_mod % ceoff_mod != 0:
            raise ValueError(f"New coefficient modulus {ceoff_mod} doesn't divide old modulus {self.coeff_mod}")
        if self.grading_mod % grading_mod != 0:
            raise ValueError(f"New grading modulus {grading_mod} doesn't divide old modulus {self.grading_mod}")
        return DGA(
            gradings=self.gradings,
            differentials=self.differentials,
            coeff_mod=ceoff_mod,
            grading_mod=grading_mod
        )

    def are_homotopy_equivalent(self, aug_1, aug_2):
        """The augmentations are homotopy equivalent iff the induced map on bilinearized homology to the base field
        given by aug = (aug_1 - aug_2) is zero. How do we figure this out using only generators of the underlying vector
        spaces of generators? Since aug(del(x)) always vanishes, it is already zero on Image(del). Therefore aug == 0
        on homology iff it vanishes on ker del.

        This is actually annoying to do in sympy and easy to do in sage :(
        However, we can do a trick... Represent ker as the intersection of the varieties determined by the row vectors
        of a marix representation of del. Say these linear equations are p_i. Since aug != 0 iff 1 in Im(aug), we see
        that aug restricted to ker is non-zero iff there is a solution to [p_1, ..., p_k, aug + 1]. This can be worked
        out with the equation solver using existence_only=True.

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
        polys = [d.to_polynomial() for d in self.differentials.values()]
        polys = [sympy.Poly(p, *self.symbols, modulus=self.coeff_mod) for p in polys]
        # Set all non-zero graded elements to zero
        zero_substitutions = {g: 0 for g in self.symbols if self.gradings[g] != 0}
        polys = [p.subs(zero_substitutions) for p in polys]
        if 0 in polys:
            polys.remove(0)
        # remove duplicates
        polys = list(set(polys))
        symbols = [g for g in self.symbols if self.gradings[g] == 0]
        self.augmentations = zero_set(polys=polys, symbols=symbols, modulus=self.coeff_mod)

    def _verify_init_args(self):
        if self.coeff_mod != 0:
            if not sympy.is_prime(self.coeff_mod):
                raise ValueError(f"Coefficient modulus {self.coeff_mod} is not 0 or prime")
        if self.symbols.keys() != self.differentials.keys():
            raise ValueError("generators don't match in symbols and keys")
        for g, d in self.differentials:
            if not isinstance(d, Differential):
                raise ValueError(f"Differential for {g} is not instance of class Differential.")

    def _correct_gradings(self):
        corrected_gradings = {g: self.gradings[g] % self.grading_mod for g in self.gradings.keys()}
        self.gradings = corrected_gradings


def zero_set(polys, symbols, modulus=2, existence_only=False):
    """
    Enumerate points in affine space generated by `symbols` in the intersection of the varieties
    determined by the polys.

    Like most of our code, this is a knock-off of Sivek's lch.sage code.
    """
    if modulus == 0:
        raise ValueError("We can only solve for zero sets over Z/mZ with m!=0.")
    root = SolutionSearchNode(polys, symbols, modulus=modulus, subs_dict=dict())
    nodes = [root]
    spawned_nodes = []
    solution_nodes = []
    while True:
        unspawned_nodes = [node for node in nodes if node not in spawned_nodes]
        if len(unspawned_nodes) == 0:
            break
        for node in unspawned_nodes:
            if node.has_solution():
                solution_nodes.append(node)
                # Terminate as early as possible if only checking existence
                if existence_only:
                    return True
            nodes += node.get_spawn()
            spawned_nodes.append(node)
    if existence_only:
        return len(solution_nodes) > 0
    return [node.subs_dict for node in solution_nodes]


def ideal_basis(polys, symbols, modulus):
    """Get a simplified basis of the idea generated by polys in GF(modulus)[symbols]"""
    return sympy.GroebnerBasis(polys, *symbols, modulus=modulus).polys


def finite_field_elements(modulus):
    return [n for n in range(modulus)]


class SolutionSearchNode(object):

    def __init__(self, polys, symbols, modulus, subs_dict):
        self.polys = [sympy.Poly(p, *symbols, modulus=modulus) for p in polys]
        self.symbols = symbols
        self.modulus = modulus
        self.subs_dict = subs_dict
        self.TERMINAL = False
        self.UNSOLVEABLE = False
        self._update_subs()

    def get_unset_vars(self):
        return [g for g in self.symbols if g not in self.subs_dict.keys()]

    def has_solution(self):
        return self.TERMINAL and (not self.UNSOLVEABLE)

    def get_spawn(self):
        output = []
        if self.TERMINAL:
            return output
        # here we can be smarter about which g we pick
        g = self.get_unset_vars()[0]
        for c in finite_field_elements(modulus=self.modulus):
            subs = self.subs_dict.copy()
            subs[g] = c
            output.append(
                SolutionSearchNode(
                    polys=self.polys,
                    symbols=self.symbols,
                    modulus=self.modulus,
                    subs_dict=subs
                )
            )
        return output

    def _update_subs(self):
        if len(self.get_unset_vars()) == 0:
            self.TERMINAL = True
            return
        if len(self.polys) == 0:
            return
        self._apply_subs()
        self._simplify_polys()
        self._check_for_const_polys()
        modified = False
        if not self.UNSOLVEABLE:
            modified = self._check_for_linear_polys()
        if len(self.get_unset_vars()) == 0:
            self.TERMINAL = True
            return
        if modified:
            self._update_subs()

    def _apply_subs(self):
        self.polys = [p.subs(self.subs_dict) for p in self.polys]

    def _simplify_polys(self):
        self.polys = ideal_basis(polys=self.polys, symbols=self.symbols, modulus=self.modulus)

    def _check_for_const_polys(self):
        if 0 in self.polys:
            self.polys.remove(0)
        for c in finite_field_elements(modulus=self.modulus):
            if c in self.polys:
                if c != 0:
                    self.TERMINAL = True
                    self.UNSOLVEABLE = True

    def _check_for_linear_polys(self):
        modified = False
        for c in finite_field_elements(modulus=self.modulus):
            for g in self.get_unset_vars():
                if g + c in self.polys:
                    self.subs_dict[g] = c
                    modified = True
        return modified
