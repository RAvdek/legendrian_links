import sympy
import utils


LOG = utils.get_logger(__name__)


def zero_set(polys, symbols, modulus=2, existence_only=False):
    """
    Enumerate points in affine space generated by `symbols` in the intersection of the varieties
    determined by the polys. Like most of our code, this is a knock-off of Sivek's lch.sage code

    :param polys: List of sympy polynomials.
    :param symbols: List of sympy symbols.
    :param modulus: Coefficient modulus.
    :param existence_only: Boolean. Are we just checking existence?
    :return: List of points in the intersection of varieties if existence_only==False else Boolean
        indicating if the zeroset is non-empy.
    """
    if modulus == 0:
        raise ValueError("We can only solve for zero sets over Z/mZ with m!=0.")
    if len(symbols) == 0:
        return []
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


def is_linear(poly):
    if isinstance(poly, int):
        if poly == 0:
            return True
        else:
            return False
    args = poly.args
    for monom in args:
        if isinstance(monom, int):
            return False
        if isinstance(monom, sympy.Symbol):
            return True
        LOG.info(args)
        LOG.info(monom)
        m_args = monom.args
        if isinstance(m_args[0], int):
            m_args = m_args[1:]
        if len(m_args) != 1:
            return False
    return True


def highest_frequency_symbol(polys, symbols):
    # Ensure that the order of symbols in input match those in polys by re-initializing
    polys = [sympy.Poly(p, *symbols) for p in polys]
    counter = {g: 0 for g in symbols}
    for p in polys:
        monoms = p.monoms()
        for m in monoms:
            for i in range(len(symbols)):
                if m[i] != 0:
                    # Here we are not weighting by powers! I think that this should not make a huge
                    # difference for LCH computations as there should be no monomials divisible
                    # by squares (eg. 1 + x + x*x*y) appearing in a plat.
                    counter[symbols[i]] += 1
    highest_freq = max(counter.values())
    highest_freq_symbols = [g for g in symbols if counter[g] == highest_freq]
    return highest_freq_symbols[0]


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
        # We choose the symbol which appears in the most monomials
        g = highest_frequency_symbol(polys=self.polys, symbols=self.get_unset_vars())
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
