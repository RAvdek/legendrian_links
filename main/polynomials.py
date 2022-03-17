from collections import Counter
import sympy
import utils


LOG = utils.LOG
# How long (in seconds) before we abandon computing Groebner bases?
GROEBNER_TIMEOUT_DEFAULT = .1
UNSET_VAR = None


def expand_zero_set_from_unset(zeros, modulus=2, fill_na=None):
    ff_elements = finite_field_elements(modulus)
    output = list()
    # Fill in the UNSET_VAR
    if fill_na is None:
        for z in zeros:
            to_process = [z]
            while len(to_process) > 0:
                p = to_process.pop()
                if any([v is UNSET_VAR for v in p.values()]):
                    new_ps = [p.copy() for _ in ff_elements]
                    to_modify = None
                    for k, v in p.items():
                        if v is UNSET_VAR:
                            to_modify = k
                    for i in ff_elements:
                        new_ps[i][to_modify] = i
                    to_process += new_ps
                else:
                    output.append(p)
    else:
        if fill_na not in finite_field_elements(modulus=modulus):
            raise ValueError(f"Attempting to fill_na with {fill_na} not in Z/{modulus}Z")
        for z in zeros:
            z_filled = {
                k: v if v is not UNSET_VAR else fill_na
                for k, v in z.items()
            }
            output.append(z_filled)
    # Fill in any affine equations
    zeros = output.copy()
    output = []
    for z in zeros:
        subs = {k: sympy.sympify(v) for k, v in z.items()}
        while not all([v.is_number for v in subs.values()]):
            subs_copy = subs.copy()
            subs_copy_clean = {k: v for k, v in subs_copy.items() if v.is_number}
            for k in [k for k in subs_copy.keys() if k not in subs_copy_clean.keys()]:
                subs[k] = subs_copy[k].subs(subs_copy_clean)
        output.append(subs)
    return output


@utils.log_start_stop
def zero_set(
        polys, symbols, modulus=2,
        subs_dicts=None, existence_only=False,
        groebner_timeout=GROEBNER_TIMEOUT_DEFAULT,
        allow_unset=True):
    """
    Enumerate points in affine space generated by `symbols` in the intersection of the varieties
    determined by the polys. Like most of our code, this is a knock-off of Sivek's lch.sage code.

    We use depth first search so that it if augmentations exist you can more quickly track their existence by watching
    logs from the command line.

    :param polys: List of sympy polynomials. Must have commutative variables.
    :param symbols: List of sympy symbols.
    :param modulus: Coefficient modulus.
    :param subs_dicts: List of dicts mapping symbols to values which we require to be satisfied.
    :param existence_only: Bool. Are we just checking existence?
    :param groebner_timeout: compute Groebner bases until timeout is reached. If failed try larger number of simpler problems.
    :param allow_unset: Bool. If a variable is arbitrary, we make it UNSET_VAR.
    :return: List of points in the intersection of varieties if existence_only==False else Boolean
        indicating if the zeroset is non-empy.
    """
    LOG.info(f"Searching for zero set of {len(polys)} polynomials in {len(symbols)} variables")
    if modulus == 0:
        raise ValueError("We can only solve for zero sets over Z/mZ with m!=0.")
    if subs_dicts is None:
        subs_dicts = [dict()]
    roots = []
    counter = Counter()
    for sd in subs_dicts:
        subbed_symbols = list(sd.keys())
        symbols_to_process = [sym for sym in symbols if sym not in subbed_symbols]
        subs_node = SubsNode(subs=sd, parent=None)
        roots.append(
            SolutionSearchNode(
                polys=polys,
                symbols=symbols_to_process,
                counter=counter,
                modulus=modulus,
                parent=None,
                initial_subs_node=subs_node,
                groebner_timeout=groebner_timeout,
                allow_unset=allow_unset
            )
        )
    living_nodes = roots
    solution_nodes = []
    loop_counter = 0
    n_analyzed_nodes = 0
    try:
        while True:
            LOG.info(f"{loop_counter}th loop of zero_set search: "
                     f"{len(living_nodes)} living nodes, "
                     f"{n_analyzed_nodes} analyzed nodes, "
                     f"{len(solution_nodes)} solution nodes")
            loop_counter += 1
            if len(living_nodes) == 0:
                break
            else:
                node = living_nodes.pop(0)
                node.run()
                n_analyzed_nodes += 1
                if node.has_solution():
                    solution_nodes.append(node)
                    # Terminate as early as possible if only checking existence
                    if existence_only:
                        return True
                else:
                    # Search depth first by spawning at the first unspawned node
                    living_nodes = node.spawn + living_nodes
        if existence_only:
            return len(solution_nodes) > 0
        output = [node.subs_node.full_subs() for node in solution_nodes]
        n_unique = len(utils.unique_elements(output))
        if n_unique != len(output):
            raise RuntimeError(f"zero_set augs are not unique! n_unique={n_unique}, n={len(output)}")
        LOG.info(f"zero_set search statistics: {dict(counter)}")
    except KeyboardInterrupt as e:
        LOG.info(f"zero_set search statistics: {dict(counter)}")
        raise e
    return output


def filtered_zero_set(polys, symbols, filtration_levels, modulus=2, groebner_timeout=GROEBNER_TIMEOUT_DEFAULT):
    if modulus == 0:
        raise ValueError("We can only solve for zero sets over Z/mZ with m!=0.")
    f_level_values = sorted(utils.unique_elements(list(filtration_levels.values())))
    f_to_sym = {
        f: {sym for sym in symbols if filtration_levels[sym] == f}
        for f in f_level_values
    }
    cum_sym = set()
    f_to_sym_cum = dict()
    for f in f_level_values:
        cum_sym.update(f_to_sym[f])
        f_to_sym_cum[f] = cum_sym.copy()
    poly_to_symbols = {p: poly_symbols(p) for p in polys}
    f_to_poly = {
        f: {k for k, v in poly_to_symbols.items() if v.issubset(f_to_sym_cum[f])}
        for f in f_level_values
    }
    subs_dicts = [dict()]
    running_symbols = set()
    for f in f_level_values:
        running_symbols.update(f_to_sym[f])
        f_polys = [
            {sympy.Poly(p, *running_symbols, modulus=modulus).subs(d) for p in f_to_poly[f]}
            for d in subs_dicts
        ]
        LOG.info(f"At filtration level={f}: {len(subs_dicts)} prior augmentations")
        extension_subs = [
            zero_set(
                polys=poly_set,
                symbols=f_to_sym[f],
                modulus=modulus,
                groebner_timeout=groebner_timeout
            )
            for poly_set in f_polys
        ]
        new_subs = []
        for i in range(len(subs_dicts)):
            d = subs_dicts[i]
            for ext in extension_subs[i]:
                new_subs.append({**d, **ext})
        subs_dicts = new_subs
    return subs_dicts


@utils.log_start_stop
def batch_zero_set(polys, symbols, modulus=2, groebner_timeout=GROEBNER_TIMEOUT_DEFAULT, batch_size=10):
    """
    Enumerate points in affine space generated by `symbols` in the intersection of the varieties determined by
    the polys. The batch implementation proceeds by iteratively adding new variables and new polynomials in batches.

    Should only be necessary when there are eg. > 100 polys.

    :param polys: List of sympy polynomials. Must have commutative variables.
    :param symbols: List of sympy symbols.
    :param modulus: Coefficient modulus.
    :param groebner_timeout: compute Groebner bases until timeout is reached.
        If failed try larger number of simpler problems.
    :param batch_size: Size of batches used.
    :return: List of points in the intersection of varieties if existence_only==False else Boolean
        indicating if the zeroset is non-empy.
    """
    LOG.info(f"Searching for (batch) zero set of {len(polys)} polynomials in {len(symbols)} variables")
    if modulus == 0:
        raise ValueError("We can only solve for zero sets over Z/mZ with m!=0.")
    if batch_size <= 0:
        raise ValueError(f"Need batch_size={batch_size} at least 1")
    poly_to_symbols = {p: poly_symbols(p) for p in polys}
    # Frequencies of symbols will help isolate those which occur most often
    ## Could just try to find zeros determined by polynomials having the most frequent symbols
    ## eg. take top 10 symbols appearing most frequently. then polys defined in terms of those symbols
    symbol_freq = Counter()
    for symbol_set in poly_to_symbols.values():
        symbol_freq.update(list(symbol_set))
    symbol_freq = dict(symbol_freq)
    for sym in symbols:
        if sym not in symbol_freq:
            symbol_freq[sym] = 0
    analyzed_polys = []
    analyzed_symbols = set()
    subs_dicts = [dict()]
    loop_counter = 0
    n_analyzed_symbols = 0
    while len(analyzed_polys) < len(polys) and len(subs_dicts) > 0:
        # We want to add new variables so that the ratio of new_polys / new_vars is the highest
        # This pattern gives a hueristic
        unanalyzed_symbols = [x for x in symbols if x not in analyzed_symbols]
        unanalyzed_symbols.sort(key=lambda k: symbol_freq[k])
        new_symbols = []
        new_polys = []
        while len(unanalyzed_symbols) > 0 and len(new_polys) <= batch_size:
            new_symbols.append(unanalyzed_symbols.pop())
            analyzed_symbols.update(set(new_symbols))
            new_polys = [
                k for k, v in poly_to_symbols.items()
                if v.issubset(analyzed_symbols) and k not in analyzed_polys
            ]
        LOG.info(f"Analysis on {loop_counter}th iteration of batch_zero_set:\n"
                 f"    n possible augmentations: {len(subs_dicts)}\n"
                 f"    n analyzed variables: {n_analyzed_symbols}\n"
                 f"    n analyzed polynomials: {len(analyzed_polys)}\n"
                 f"    n new variables: {len(analyzed_symbols) - n_analyzed_symbols}\n"
                 f"    n new polynomials: {len(new_polys)}\n"
                 f"    new polynomials: {new_polys}")
        analyzed_polys += new_polys
        n_analyzed_symbols = len(analyzed_symbols)
        loop_counter += 1
        subs_dicts = zero_set(
            polys=new_polys,
            symbols=list(analyzed_symbols),
            modulus=modulus,
            subs_dicts=subs_dicts,
            groebner_timeout=groebner_timeout,
            allow_unset=True
        )

    return subs_dicts



def ideal_basis(polys, symbols, modulus):
    """Get a simplified basis of the idea generated by polys in GF(modulus)[symbols]"""
    # grlex is much faster than lex! This unbroke some aug searches.
    # sympy has 'buchberger' and 'f5b' methods for Groebner bases. It seems f5b is more up-to-date from googling...
    # https://github.com/sympy/sympy/blob/master/sympy/polys/groebnertools.py
    return utils.unique_elements(
        sympy.GroebnerBasis(polys, *symbols, modulus=modulus, order='grlex', method='f5b').polys
    )


def finite_field_elements(modulus):
    return [n for n in range(modulus)]


def is_linear(poly):
    """Is the polynomial linear? This is not the same as being affine :)

    :param poly: polynomial
    :return: Bool... is linear or not
    """
    if poly == 0:
        return True
    poly = sympy.sympify(poly)
    coeff_dict = poly.as_coefficients_dict()
    for k in list(coeff_dict.keys()):
        if not k.is_symbol:
            return False
    return True


def is_affine(poly):
    """Is the polynomial affine?

    :param poly: polynomial
    :return: Bool... is affine or not
    """
    if poly == 0:
        return True
    poly = sympy.sympify(poly)
    coeff_dict = poly.as_coefficients_dict()
    for k in list(coeff_dict.keys()):
        if not (k.is_symbol or k.is_number):
            return False
    return True


def poly_symbols(poly):
    """Return a set of symbols which appears in a polynomial expression

    :param poly:
    :return:
    """
    # TODO: Requires testing
    output = set()
    poly = sympy.sympify(poly)
    coeff_dict = dict(poly.as_expr().as_coefficients_dict())
    for k in list(coeff_dict.keys()):
        if k.is_number:
            continue
        for x in k.as_terms()[-1]:
            output.add(x)
    return output

def many_poly_symbols(polys):
    """Return a set of symbols which appears in a tuple of polynomial expressions

    :param polys: list of polynomials
    :return: set of symbols
    """
    output = set()
    for p in polys:
        output.update(poly_symbols(p))
    return output


def highest_frequency_symbol(polys, symbols):
    """For a collection of polynomials, find the symbol (variable) which appears which greatest frequency

    :param polys: list of polynomials
    :param symbols: list of symbols
    :return: symbol
    """
    # Ensure that the order of symbols in input match those in polys by re-initializing
    polys_to_analyze = [sympy.Poly(p, *symbols) for p in polys]
    counter = {g: 0 for g in symbols}
    for p in polys_to_analyze:
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


class SubsNode(object):
    """
    To prevent replication of subs dicts (taking up tons of memory). These objects
    will patch together substitutions in a tree. Each node only stores the incremental
    changes to the subs dictionaries along with a pointer to the parent. From this data,
    we can construct the entire subs dict.
    """

    def __init__(self, subs, parent=None):
        self.subs = subs
        self.parent = parent

    def full_subs(self):
        output = self.subs
        current_node = self
        while True:
            current_node = current_node.parent
            if current_node is None:
                return output
            else:
                new_subs = current_node.subs
                for k, v in new_subs.items():
                    if k in output.keys():
                        old = output[k]
                        if old != v:
                            raise RuntimeError(f"Found mismatching duplicate subs, old={k}:{output[k]}, new={k}:{v}")
                    output[k] = v


class SolutionSearchNode(object):
    """
    Similar to SubsNode, instances of this class form a tree. Data is reconstructed by traversing to the root.
    This way we don't have to store the data of all polynomials and all of their substitutions.

    Usage pattern...
    node = SolutionSearchNode(...) # set root node
    node.run()                     # do a bunch of stuff
    children = nodes.spawn         # get child nodes

    """
    STATE_INIT = 0
    STATE_RUNNING = 1
    STATE_RUN_COMPLETE = 2
    STATE_TEARDOWN_COMPLETE = 3

    EVENT_ELIM_0 = 'elim_zero_poly'
    EVENT_NON0 = 'nonzero_poly'
    EVENT_AFFINE_SUB = 'affine_sub'
    EVENT_QUAD_SUB = 'quad_sub'

    def __init__(self, initial_subs_node, polys, counter, symbols=None, modulus=None, parent=None,
                 groebner_timeout=GROEBNER_TIMEOUT_DEFAULT, allow_unset=False):
        self.initial_subs_node = initial_subs_node
        self.parent = parent
        if self.parent is not None:
            self._set_from_parents()
        else:
            if symbols is None or modulus is None or allow_unset is None:
                raise ValueError("Node not properly initialized")
            self.symbols = symbols
            self.modulus = modulus
            self.allow_unset = allow_unset
            self.depth = 0
            self.groebner_timeout = groebner_timeout
            self.ff_elements = finite_field_elements(self.modulus)
        self.polys = [sympy.Poly(p, *self.symbols, modulus=self.modulus) for p in polys]
        self.counter = counter
        self._set_id()
        self.TERMINAL = False
        self.UNSOLVEABLE = False
        self.state = self.STATE_INIT

    def run(self):
        self.state = self.STATE_RUNNING
        self._set_subs_dict_from_subs_node()
        self._cleanup_subs_and_symbols()
        self._apply_subs()
        if len(self.polys) > 0:
            self._unique_polys()
            self._setup_quads()
            self._update_subs_and_polys()
        self._cleanup_subs_and_symbols()
        self._set_subs_node()
        self._check_unset_vars()
        self._set_spawn()
        self.state = self.STATE_RUN_COMPLETE
        self._teardown()
        self.state = self.STATE_TEARDOWN_COMPLETE

    def get_unset_vars(self):
        """
        :return: list of symbols not already set by self.subs_dict
        """
        return [g for g in self.symbols if g not in self.subs_dict.keys()]

    def has_solution(self):
        """
        :return: Bool for if the node has an already set solution
        """
        return self.TERMINAL and (not self.UNSOLVEABLE)

    def get_symbols_in_polys(self):
        """Get a set of symbols which appear in self.polys

        :return: set
        """
        return many_poly_symbols(self.polys)

    def _teardown(self):
        if self.parent is not None:
            del self.symbols
            del self.polys
            del self.modulus
            del self.groebner_timeout
            del self.allow_unset
            del self.ff_elements
            del self.depth
            del self.subs_dict

    def _set_from_parents(self):
        node = self
        depth = 0
        while node.parent is not None:
            node = node.parent
            depth += 1
        self.symbols = node.symbols
        self.modulus = node.modulus
        self.allow_unset = node.allow_unset
        self.groebner_timeout = node.groebner_timeout
        self.ff_elements = node.ff_elements
        self.depth = depth

    def _set_subs_dict_from_subs_node(self):
        self.subs_dict = self.initial_subs_node.full_subs()

    def _set_subs_node(self):
        subs = {
            k: v for k, v in self.subs_dict.items()
            if k not in self.initial_subs_node.full_subs().keys()
        }
        self.subs_node = SubsNode(subs=subs, parent=self.initial_subs_node)

    def _set_spawn(self):
        output = []
        if not self.TERMINAL:
            # We choose the symbol which appears in the most monomials
            g = highest_frequency_symbol(polys=self.polys, symbols=self.get_unset_vars())
            for c in self.ff_elements:
                subs_copy = self.subs_dict.copy()
                subs_copy[g] = c
                old_subs = self.initial_subs_node.full_subs()
                new_subs = dict()
                for k, v in subs_copy.items():
                    if k not in old_subs.keys():
                        new_subs[k] = v
                next_subs_node = SubsNode(subs=new_subs, parent=self.subs_node)
                output.append(SolutionSearchNode(
                    initial_subs_node=next_subs_node, polys=self.polys, counter=self.counter, parent=self))
        self.spawn = output

    def _set_id(self):
        self.id = "{" + f"D={self.depth}:id=" + utils.tiny_id() + "}"

    def _setup_quads(self):
        if self.modulus == 2:
            self._quads = set()
            unset_symbols = self.get_unset_vars()
            for s_1 in unset_symbols:
                for s_2 in unset_symbols:
                    self._quads.add(s_1 * s_2)

    def _cleanup_subs_and_symbols(self):
        if self.allow_unset:
            used_symbols = self.get_symbols_in_polys()
            output_symbols = list()
            for g in self.symbols:
                if g in used_symbols:
                    output_symbols.append(g)
                else:
                    if g not in self.subs_dict:
                        self.subs_dict[g] = UNSET_VAR
            self.symbols = output_symbols


    def _check_unset_vars(self):
        if len(self.get_unset_vars()) == 0:
            self.TERMINAL = True
        if self.allow_unset:
            if len(self.polys) == 0:
                self.TERMINAL = True
                for g in self.symbols:
                    if g not in self.subs_dict.keys():
                        self.subs_dict[g] = UNSET_VAR

    def _update_subs_and_polys(self):
        """Update self.subs_dict and self.polys

        :return: None
        """
        n_unset_vars = len(self.get_unset_vars())
        n_polys = len(self.polys)
        if n_polys == 0:
            return
        LOG.info(f"SolutionSearchNode={self.id}: "
                 f"{n_unset_vars} unset vars & {n_polys} polys")
        while not self.UNSOLVEABLE:
            self._manually_reduce_polys()
            n_old_unset_vars = n_unset_vars
            n_unset_vars = len(self.get_unset_vars())
            # Break if we didn't eliminate any variables
            if n_unset_vars == n_old_unset_vars:
                break
        if self.UNSOLVEABLE:
            self.TERMINAL = True
            return
        if n_unset_vars == 0:
            self.TERMINAL = True
            return
        modified = self._simplify_polys()
        self._cleanup_subs_and_symbols()
        if modified:
            self._update_subs_and_polys()

    def _apply_subs(self, specific_subs=None):
        """Update self.polys by applying substitutions in self.subs_dict or specific_subs.

        Contract is that we do not remove items from or add items to self.polys.

        :param specific_subs: dict or None... substitutions to apply to self.polys
        :return: None
        """
        # The following can lead to hitting a recursion limit!
        # self.polys = [p.subs(specific_subs) for p in self.polys]
        if specific_subs is not None:
            subs_to_apply = specific_subs
        else:
            subs_to_apply = self.subs_dict
        clean_subs_to_apply = {k: v for k, v in subs_to_apply.items() if v is not UNSET_VAR}
        unset_vars = self.get_unset_vars()
        output = []
        for p in self.polys:
            expr = p.as_expr().subs(clean_subs_to_apply)
            if len(unset_vars) > 0:
                output.append(sympy.Poly(expr, *unset_vars, modulus=self.modulus))
            else:
                # sympy will complain if you try to make a poly with no variables...
                output.append(expr % self.modulus)
        self.polys = output

    def _unique_polys(self):
        self.polys = utils.unique_elements(self.polys)

    def _simplify_polys(self):
        """Attempt to simplify ideal for polynomials using Groebner bases until timeout.

        Returns True/False for if self.polys has been modified.

        :return: Bool
        """
        modified = False
        in_main_thread = utils.in_main_thread()
        thread_timeout = False
        # If not running from the main thread, skip usage of Groebner.
        # For app usage (when we'll be in a worker thread), we have to assume that zero_set args are pretty simple...
        # We need to be in the main thread in order to use utils.timeout_ctx
        if not in_main_thread:
            return False
        def new_ideal_basis():
            with utils.timeout_ctx(self.groebner_timeout):
                return ideal_basis(polys=self.polys, symbols=self.get_unset_vars(), modulus=self.modulus)
        try:
            grob_polys = new_ideal_basis()
        except sympy.polys.polyerrors.CoercionFailed as e:
            LOG.info("BIG PROBS")
            LOG.info('\n'.join([str(p.as_expr()) for p in self.polys]))
            LOG.info(self.get_unset_vars())
            LOG.info(self.symbols)
            LOG.info(self.subs_dict)
            LOG.info(self.modulus)
            raise e
        if grob_polys is None:
            thread_timeout = True
        if thread_timeout:
            LOG.info(f"SolutionSearchNode={self.id}: Failed to compute Groebner basis")
        else:
            n_polys = len(self.polys)
            n_grob_polys = len(grob_polys)
            if set(self.polys) != set(grob_polys):
                LOG.info(f"SolutionSearchNode={self.id}: "
                         f"Groebner simplified {n_polys} polys to {n_grob_polys}")
                # Groebner can increase sizes of generating sets
                if n_grob_polys <= n_polys:
                    self.polys = grob_polys
                    modified = True
        return modified

    def _manually_reduce_polys(self):
        """Make a single pass over self.polys, checking for affine terms (ax + b)
        Currently this works over Z/2Z. Need to divide by a in general.

        :return: Bool. Are we modifying the set of polynomials?
        """
        n_polys = len(self.polys)
        del_indicies = set()
        modified = False
        for i in range(n_polys):
            p = self.polys[i]
            p_exp = p.as_expr()
            # remove zero polynomials
            if p_exp == 0:
                del_indicies.add(i)
                modified = True
                self.counter.update([self.EVENT_ELIM_0])
                continue
            # if a poly is non-zero constant, we can never find a zero
            if p_exp.is_number:
                # Must be non-zero due to previous `if`
                self.TERMINAL = True
                self.UNSOLVEABLE = True
                modified = True
                self.counter.update([self.EVENT_NON0])
                break
            # check for affine terms... linear_combination + const
            if is_affine(p_exp):
                coeff_dict = p_exp.as_coefficients_dict()
                sym_keys = [k for k in coeff_dict.keys() if k in self.get_unset_vars()]
                if len(sym_keys) > 0:
                    g = sym_keys[0]
                    # We'd have to clean this up for non-ZZ2. This is always true in this case...
                    if coeff_dict[g] == 1:
                        p_min_g = p_exp - g
                        del_indicies.add(i)
                        self.subs_dict[g] = p_min_g
                        self._apply_subs(specific_subs={g: p_min_g})
                        modified = True
                        self.counter.update([self.EVENT_AFFINE_SUB])
            # check for terms of the form p = g - c for g a symbol and c a constant
            # for c in self.ff_elements:
            #    p_min_c = p_exp - c
            #    if p_min_c.is_symbol:
            #        del_indicies.add(i)
            #        self.subs_dict[p_min_c] = c
            #        self._apply_subs(specific_subs={p_min_c: c})
            #        modified = True
            #        break
            # over Z/2Z we can also manually remove quadratic terms.
            # I've seen these come up frequently.
            # p = sym_1 * sym_2 + 1 = 0 => both are 1.
            # we should be able to generalize this to things of the form p = 1 + x_1 * ... * x_k
            # though I don't know how useful it would be
            if self.modulus == 2:
                # this is poorly designed: p_exp + 1 won't be in self._quads,
                # even though elements may agree over Z/2Z
                if p_exp - 1 in self._quads:
                    del_indicies.add(i)
                    # this will extract the individual factors
                    sym_1, sym_2 = (p_exp - 1).as_ordered_factors()
                    self.subs_dict[sym_1] = 1
                    self.subs_dict[sym_2] = 1
                    self._apply_subs(specific_subs={sym_1: 1, sym_2: 1})
                    modified = True
                    self.counter.update([self.EVENT_QUAD_SUB])
        if modified:
            polys_copy = [v for i, v in enumerate(self.polys) if i not in del_indicies]
            self.polys = polys_copy
            self._unique_polys()
        return modified
