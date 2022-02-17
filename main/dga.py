from math import prod
import sympy
import utils
from zero_set import zero_set

LOG = utils.get_logger(__name__)


class Differential(object):
    """Differential of a DGA element.

    self.summands is a list of lists. Each list stores the factors in the order in which they appear. Eg.
    1 + x + z + x*y*z will be stored as [[1], [x], [z], [x,y,z]] whereas...
    1 + x + z + z*y*x will be stored as [[1], [x], [z], [z,y,x]]
    This non-commutative storage is needed for bilinearization.
    """

    def __init__(self, summands=[], coeff_mod=2, signs=[]):
        self.summands = summands
        self.coeff_mod = coeff_mod
        self.signs = signs

    def __repr__(self):
        if len(self.summands) == 0:
            return str(0)
        output_terms = []
        for i in range(len(self.summands)):
            output_terms.append('*'.join([str(f) for f in self.summands[i]]))
        return " + ".join(output_terms)

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
        """Return bilinearized differential"""
        # throw out contributions with 0 negative ends
        bilin_signs = []
        bilin_summands = []
        for i in range(len(self.summands)):
            if len(self.summands[i]) > 0:
                bilin_signs.append(self.signs[i])
                bilin_summands.append(self.summands[i])
        output_signs = []
        output_summands = []
        for j in range(len(bilin_summands)):
            sign = bilin_signs[i]
            bs = bilin_summands[i]
            for i in range(len(bs)):
                aug_1_factor = 1 if i == 0 else prod(bs[:i]).subs(subs_1)
                aug_2_factor = 1 if i == len(bs) - 1 else prod(bs[i+1:]).subs(subs_2)
                output_signs.append(sign*aug_1_factor*aug_2_factor)
                output_summands.append([bs[i]])
        return Differential(summands=output_summands, coeff_mod=self.coeff_mod, signs=output_signs)


class DGA(object):

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
        on homology iff it vanishes on ker del.

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
