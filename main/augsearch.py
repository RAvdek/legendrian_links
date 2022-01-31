from time import time
import sympy


def get_augs(polys, symbols, modulus=2):
    """This actually doesn't work yet for the intended application,
    because it assumes that all of the symbols are cycles.

    Applies a tree search to find augmentations of a polynomial ring generated by `symbols`
    quotiented out by the idea determined by the polynomials `polys`.
    We use coeffs in ZZ/(`modulus`)ZZ.

    IOW, we enumerate points in affine space in the intersection of the varieties determined by the polys.

    Like most of our code, this is a knock-off of Sivek's lch.sage code.
    """
    root = SolutionSearchNode(polys, symbols, modulus=modulus, subs_dict=dict())
    root.update()
    nodes = [root]
    spawned_nodes = []
    while True:
        unspawned_nodes = [node for node in nodes if node not in spawned_nodes]
        if len(unspawned_nodes) == 0:
            break
        for node in unspawned_nodes:
            node.update()
            nodes += node.get_spawn()
            spawned_nodes.append(node)
    solution_nodes = [node for node in nodes if node.TERMINAL and (not node.UNSOLVEABLE)]
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

    def get_unset_vars(self):
        return [g for g in self.symbols if g not in self.subs_dict.keys()]

    def update(self):
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
            self.update()

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
