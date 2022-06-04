import json
import os
import tempfile
import unittest
import numpy as np
import sympy
import polynomials
import utils
import algebra
import legendrian_links as ll


LOG = utils.get_logger(__name__)


def comparable_list_of_dicts(l):
    # Have to convert to sets of json strings as we cannot make sets of dicts in TestCase.AssertEqual
    return sorted([json.dumps({str(k): str(v) for k, v in d.items()}, sort_keys=True) for d in l])


class TestUtils(unittest.TestCase):

    def test_num_inverse(self):
        self.assertEqual(sympy.Rational(1, 2), utils.num_inverse(2, 0))
        self.assertEqual(sympy.Rational(3, 2), utils.num_inverse(sympy.Rational(2, 3), 0))
        self.assertEqual(1, utils.num_inverse(1, 2))
        self.assertEqual(2, utils.num_inverse(3, 5))

    def test_one_hot(self):
        self.assertTrue(np.array_equal(utils.one_hot_array(1, 2), [0, 1]))
        self.assertTrue(np.array_equal(utils.one_hot_array(3, 5), [0, 0, 0, 1, 0]))


class TestLinks(unittest.TestCase):

    def _test_link(self, n_strands, front_crossings, n_knots, n_disks, tb_rot=None, lch_gradings=None):

        pd = ll.PlatDiagram(n_strands=n_strands, front_crossings=front_crossings, lazy_lch=False, lazy_rsft=False)

        # There is exactly one right close and one left close
        self.assertEqual(len([ps for ps in pd.plat_segments if ps.left_close]), 1)
        self.assertEqual(len([ps for ps in pd.plat_segments if ps.right_close]), 1)

        # Number of knots, chords, and disks are expected
        self.assertEqual(len(pd.knots), n_knots)
        self.assertEqual(len(pd.chords), len(front_crossings) + (n_strands // 2))
        self.assertEqual(len(pd.disks), n_disks)

        # There is exactly one t_label for each component of the link
        for knot in pd.knots:
            t_labeled_segments = [ls for ls in knot['line_segments'] if ls.t_label]
            self.assertTrue(len(t_labeled_segments), 1)

        # Each capping path should have rotation a multiple of 1/2
        for cp in pd.capping_paths:
            self.assertEqual(int(2 * cp.rotation_number) % 2, 1)

        if lch_gradings is not None:
            self.assertEqual({g.grading for g in pd.lch_generators}, lch_gradings)

        # The x-values of disk segments in each disk are all distinct
        for d in pd.disks:
            ds_x_values = [ds.x for ds in d.disk_segments]
            self.assertEqual(len(ds_x_values), len(set(ds_x_values)))

        # tb and rot
        if tb_rot is not None:
            tb, rot = tb_rot
            if len(pd.knots) == 1:
                self.assertEqual(pd.knots[0]['tb'], tb)
                self.assertTrue(pd.knots[0]['rot'] in rot)

    def test_unknot(self):
        self._test_link(2, [], 1, 2, tb_rot=[-1, [0]], lch_gradings={1})

    def test_stabilized_unknot(self):
        self._test_link(2, [0], 1, 3, tb_rot=[-2, [1, -1]])

    def test_twisted_unknot(self):
        self._test_link(4, [1], 1, 4, lch_gradings={0, 1})

    def test_very_twisted_unknot(self):
        self._test_link(2, [0]*5, 1, 7, tb_rot=[-6, [1, -1]])

    def test_double_unlink(self):
        self._test_link(4, [], 2, 4)

    def test_hopf_link(self):
        self._test_link(4, [1]*2, 2, 7)

    def test_trefoil(self):
        self._test_link(4, [1]*3, 1, 10, tb_rot=[1, [0]], lch_gradings={0, 1})

    def test_composable(self):
        trefoil = ll.PlatDiagram(n_strands=4, front_crossings=[1, 1, 1])
        chord_1 = trefoil.get_chord_from_x(1)
        chord_2 = trefoil.get_chord_from_x(2)
        self.assertTrue(chord_1.is_composable_with(chord_1))
        self.assertTrue(chord_1.is_composable_with(chord_2))

        hopf = ll.PlatDiagram(n_strands=4, front_crossings=[1, 1])
        chord_1 = hopf.get_chord_from_x(1)
        chord_2 = hopf.get_chord_from_x(2)
        self.assertFalse(chord_1.is_composable_with(chord_1))
        self.assertTrue(chord_1.is_composable_with(chord_2))

    def test_word_is_admissible(self):
        trefoil = ll.PlatDiagram(n_strands=4, front_crossings=[1, 1, 1])
        chord_1 = trefoil.get_chord_from_x(1)
        chord_2 = trefoil.get_chord_from_x(2)
        self.assertTrue(trefoil.word_is_admissible([chord_1]))
        self.assertFalse(trefoil.word_is_admissible([chord_1, chord_1]))
        self.assertFalse(trefoil.word_is_admissible([chord_1, chord_2]))

        hopf = ll.PlatDiagram(n_strands=4, front_crossings=[1, 1])
        chord_1 = hopf.get_chord_from_x(1)
        chord_2 = hopf.get_chord_from_x(2)
        self.assertTrue(hopf.word_is_admissible([chord_1], cyclic=False))
        self.assertFalse(hopf.word_is_admissible([chord_1], cyclic=True))
        self.assertTrue(hopf.word_is_admissible([chord_1, chord_2], cyclic=True))
        self.assertTrue(hopf.word_is_admissible([chord_1, chord_2], cyclic=False))
        self.assertFalse(hopf.word_is_admissible([chord_1, chord_1], cyclic=True))
        self.assertFalse(hopf.word_is_admissible([chord_1, chord_1], cyclic=False))
        self.assertFalse(hopf.word_is_admissible([chord_1, chord_2, chord_1], cyclic=False))
        self.assertFalse(hopf.word_is_admissible([chord_1, chord_2, chord_1], cyclic=True))
        self.assertFalse(hopf.word_is_admissible([chord_1, chord_2, chord_1, chord_2], cyclic=True))


class TestZeroSet(unittest.TestCase):

    def setUp(self):
        self.a, self.b, self.c = sympy.symbols('a,b,c')
        self.symbols = [self.a, self.b, self.c]

    def test_is_linear(self):
        self.assertTrue(polynomials.is_linear(0))
        self.assertFalse(polynomials.is_linear(1 + self.a))
        self.assertTrue(polynomials.is_linear(self.a + self.b))
        self.assertFalse(polynomials.is_linear(1 + self.a + self.b))
        self.assertFalse(polynomials.is_linear(self.a * self.b))

    def test_is_affine(self):
        self.assertTrue(polynomials.is_affine(0))
        self.assertTrue(polynomials.is_affine(1 + self.a))
        self.assertTrue(polynomials.is_affine(self.a + self.b))
        self.assertTrue(polynomials.is_affine(1 + self.a + self.b))
        self.assertFalse(polynomials.is_affine(self.a * self.b))

    def test_poly_symbols(self):
        z = sympy.Poly(1 + self.a * self.b * self.c, *self.symbols)
        self.assertEqual(polynomials.poly_symbols(z), set(self.symbols))

        z = sympy.Poly(1, *self.symbols)
        self.assertEqual(polynomials.poly_symbols(z), set())

        z = sympy.Poly(1 + self.a * self.b, *self.symbols)
        self.assertEqual(polynomials.poly_symbols(z), {self.a, self.b})

    def test_expand_zero_set_from_unset(self):
        z = []
        output = polynomials.expand_zero_set_from_unset(z, 2)
        self.assertEqual(output, [])

        z = [{'a': polynomials.UNSET_VAR}]
        output = polynomials.expand_zero_set_from_unset(z, 2)
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts([{'a': 0}, {'a': 1}]))

        z = [{'a': 1}]
        output = polynomials.expand_zero_set_from_unset(z, 2)
        self.assertEqual(output, z)

        z = [{'a': 1, 'b': polynomials.UNSET_VAR}]
        output = polynomials.expand_zero_set_from_unset(z, 2)
        expected = [{'a': 1, 'b': 0}, {'a': 1, 'b': 1}]
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts(expected))

        z = [{self.a: 1, self.b: self.a}]
        output = polynomials.expand_zero_set_from_unset(z, 2)
        expected = [{self.a: 1, self.b: 1}]
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts(expected))

        z = [{self.a: 1, self.b: self.a, self.c: self.b}]
        output = polynomials.expand_zero_set_from_unset(z, 2)
        expected = [{self.a: 1, self.b: 1, self.c: 1}]
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts(expected))

        z = [{'a': 1, 'b': polynomials.UNSET_VAR}]
        output = polynomials.expand_zero_set_from_unset(z, modulus=2, fill_na=0)
        expected = [{'a': 1, 'b': 0}]
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts(expected))

        z = [{'a': 1, 'b': polynomials.UNSET_VAR}]
        output = polynomials.expand_zero_set_from_unset(z, modulus=2, fill_na=1)
        expected = [{'a': 1, 'b': 1}]
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts(expected))

        z = [{self.a: polynomials.UNSET_VAR, self.b: self.a, self.c: self.b}]
        output = polynomials.expand_zero_set_from_unset(z, modulus=2, fill_na=0)
        expected = [{self.a: 0, self.b: 0, self.c: 0}]
        self.assertEqual(comparable_list_of_dicts(output), comparable_list_of_dicts(expected))

    def test_non_zero_const(self):
        polys = [1]
        output = polynomials.zero_set(polys=polys, symbols=self.symbols)
        self.assertEqual(len(output), 0)

        output = polynomials.zero_set(polys=polys, symbols=[])
        self.assertEqual(len(output), 0)

    def test_zero_const(self):
        polys = [0]
        output = polynomials.zero_set(polys=polys, symbols=self.symbols, allow_unset=False)
        self.assertEqual(len(output), 8)
        output = polynomials.zero_set(polys=polys, symbols=self.symbols, allow_unset=True)
        self.assertEqual(len(output), 1)

        output = polynomials.zero_set(polys=polys, symbols=[])
        self.assertEqual(output, [dict()])

    def test_quad(self):
        polys = [1 + self.a*self.b]
        expected_output = [{self.a: 1, self.b: 1}]
        results = polynomials.zero_set(polys=polys, symbols=[self.a, self.b], allow_unset=False)
        self.assertEqual(expected_output, results)

        expected_output = [{self.a: 1, self.b: 1, self.c: polynomials.UNSET_VAR}]
        results = polynomials.zero_set(polys=polys, symbols=self.symbols, allow_unset=True)
        self.assertEqual(comparable_list_of_dicts(expected_output), comparable_list_of_dicts(results))

    def test_deg_three(self):
        polys = [1 + self.a * self.b * self.c]
        expected_output = [{self.a: 1, self.b: 1, self.c: 1}]
        output = polynomials.zero_set(polys=polys, symbols=self.symbols, allow_unset=False)
        self.assertEqual(expected_output, output)

    def test_trefoil(self):
        polys = [1 + self.a + self.c + self.a * self.b * self.c]
        results = polynomials.zero_set(polys=polys, symbols=self.symbols, modulus=2)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        expected_results = [
            {self.a: 0, self.c: 1, self.b: 0},
            {self.a: 0, self.c: 1, self.b: 1},
            {self.a: 1, self.b: 0, self.c: 0},
            {self.a: 1, self.b: 1, self.c: 0},
            {self.a: 1, self.b: 1, self.c: 1}
        ]
        self.assertEqual(comparable_list_of_dicts(results), comparable_list_of_dicts(expected_results))

    def test_fivefold(self):
        d, e = sympy.symbols('d, e')
        polys = [1 + self.a + e + (self.a * self.b * self.c * d * e)]
        symbols = self.symbols + [d, e]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2, allow_unset=False)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        self.assertEqual(len(results), 17)

    def test_fivefold_plus_linear(self):
        d, e = sympy.symbols('d, e')
        polys = [1 + self.a + e + (self.a * self.b * self.c * d * e), self.a]
        symbols = self.symbols + [d, e]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2, allow_unset=False)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        self.assertEqual(len(results), 8)

        polys = [1 + self.a + e + (self.a * self.b * self.c * d * e), self.a + 1]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2, allow_unset=False)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        self.assertEqual(len(results), 9)

        polys = [1 + self.a + e + (self.a * self.b * self.c * d * e), self.a + 1, e]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2, allow_unset=False)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        self.assertEqual(len(results), 8)

        polys = [1 + self.a + e + (self.a * self.b * self.c * d * e), self.a + 1, e + 1]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2, allow_unset=False)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        self.assertEqual(len(results), 1)

    def test_hopf_link(self):
        polys = [self.a * self.b]
        results = polynomials.zero_set(polys=polys, symbols=[self.a, self.b], modulus=2)
        results = polynomials.expand_zero_set_from_unset(results, modulus=2)
        expected_results = [
            {self.a: 1, self.b: 0},
            {self.a: 0, self.b: 0},
            {self.a: 0, self.b: 1}
        ]
        self.assertEqual(comparable_list_of_dicts(results), comparable_list_of_dicts(expected_results))


class TestMatrix(unittest.TestCase):

    def test_modulus(self):
        self.assertFalse(
            np.any(
                algebra.Matrix([[4]], coeff_mod=3).values
                - algebra.Matrix([[1]], coeff_mod=3).values
            )
        )
        self.assertFalse(
            np.any(
                algebra.Matrix([[1, 2], [2, 2]], coeff_mod=2).values
                - algebra.Matrix([[1, 0], [0, 0]], coeff_mod=2).values
            )
        )

    def test_ref_methods(self):
        ref_q = algebra.Matrix([[1, 0], [1, 1]], coeff_mod=2)
        ref_q_inv = ref_q
        ref = algebra.Matrix([[1, 1, 1, 0], [0, 1, 0, 0]], coeff_mod=2)
        mat = ref_q_inv.multiply(ref)
        mat.set_row_echelon()
        self.assertTrue(np.array_equal(ref.values, mat.ref.values))
        self.assertTrue(np.array_equal(ref_q_inv.values, mat.ref_q_inv.values))
        self.assertTrue(np.array_equal(ref_q.values, mat.ref_q.values))

        ref_q = algebra.Matrix([[1, 0], [1, 1]], coeff_mod=2)
        ref_q_inv = ref_q
        ref = algebra.Matrix([[1, 0, 1, 0], [0, 1, 0, 0]], coeff_mod=2)
        mat = ref_q_inv.multiply(ref)
        mat.set_red_row_echelon()
        self.assertTrue(np.array_equal(ref.values, mat.ref.values))
        self.assertTrue(np.array_equal(ref_q_inv.values, mat.ref_q_inv.values))
        self.assertTrue(np.array_equal(ref_q.values, mat.ref_q.values))

    def test_kernel(self):
        mat = algebra.Matrix([[0, 1], [0, 0]])
        ker = mat.kernel()
        self.assertEqual(len(ker), 1)
        self.assertTrue(np.array_equal([1, 0], ker[0]))

        mat = algebra.Matrix([[1, 5], [0, 0]])
        ker = mat.kernel()
        self.assertEqual(len(ker), 1)
        self.assertTrue(np.array_equal([-5, 1], ker[0]))

        mat = algebra.Matrix(np.zeros([2, 2]))
        ker = mat.kernel()
        self.assertEqual(len(ker), 2)


class TestMatrixChainComplex(unittest.TestCase):

    def test_trivial(self):
        ranks = {0: 1}
        differentials = dict()
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=2)
        self.assertEqual(cx.poincare_poly(), 1)

    def test_stabilization(self):
        ranks = {0: 1, 1: 1}
        differentials = {1: algebra.Matrix([[1]], coeff_mod=0)}
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=0)
        self.assertEqual(cx.poincare_poly(), 0)

    def test_circle(self):
        ranks = {0: 1, 1: 1}
        differentials = {1: algebra.Matrix([[0]], coeff_mod=0)}
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=0)
        self.assertEqual(cx.poincare_poly(), 1 + algebra.DEGREE_VAR)

    def test_rp2(self):
        ranks = {0: 1, 1: 1, 2: 1}
        differentials = {
            2: algebra.Matrix([[2]], coeff_mod=0),
            1: algebra.Matrix([[0]], coeff_mod=0)
        }
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=0)
        self.assertEqual(cx.poincare_poly(), 1)

        ranks = {0: 1, 1: 1, 2: 1}
        differentials = {
            2: algebra.Matrix([[2]], coeff_mod=2),
            1: algebra.Matrix([[0]], coeff_mod=2)
        }
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=2)
        self.assertEqual(cx.poincare_poly(), 1 + algebra.DEGREE_VAR + algebra.DEGREE_VAR ** 2)

    def test_non_diag(self):
        ranks = {0: 1, 1: 2}
        differentials = {
            1: algebra.Matrix([[1, 0]], coeff_mod=0)
        }
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=0)
        self.assertEqual(cx.poincare_poly(), algebra.DEGREE_VAR)

        ranks = {0: 2, 1: 2}
        differentials = {
            1: algebra.Matrix([[1, 0], [0, 0]], coeff_mod=0)
        }
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=0)
        self.assertEqual(cx.poincare_poly(), 1 + algebra.DEGREE_VAR)

        ranks = {0: 2, 1: 2}
        differentials = {
            1: algebra.Matrix([[0, 1], [0, 0]], coeff_mod=0)
        }
        cx = algebra.MatrixChainComplex(ranks=ranks, differentials=differentials, coeff_mod=0)
        cx.set_qrs_decomposition()
        self.assertEqual(cx.poincare_poly(), 1 + algebra.DEGREE_VAR)


class TestSpectralSequence(unittest.TestCase):

    def test_trivial(self):
        x = sympy.symbols('x')
        gradings = {x: 0}
        differentials = {x: algebra.Differential(0, coeff_mod=0)}
        filtration_levels = {x: 1}
        specseq = algebra.SpectralSequence(
            gradings=gradings, differentials=differentials, filtration_levels=filtration_levels, coeff_mod=0)
        p_poly = specseq.poincare_poly()
        self.assertEqual(p_poly, algebra.FILTRATION_VAR)

    def test_circle(self):
        # Morse function on a circle
        n, s = sympy.symbols('n, s')
        gradings = {n: 1, s: 0}
        differentials = {
            n: algebra.Differential(0, coeff_mod=0),
            s: algebra.Differential(0, coeff_mod=0)}
        filtration_levels = {n: 1, s: 1}
        specseq = algebra.SpectralSequence(
            gradings=gradings, differentials=differentials, filtration_levels=filtration_levels, coeff_mod=0)
        p_poly = specseq.poincare_poly()
        self.assertTrue(utils.poly_equal(p_poly, algebra.FILTRATION_VAR * (1 + algebra.DEGREE_VAR)))

    def test_sphere(self):
        # Morse function with a pair of critical points on S2 for north and south poles
        # This converges at the zeroth page
        n, s = sympy.symbols('n, s')
        gradings = {n: 2, s: 0}
        differentials = {
            n: algebra.Differential(0, coeff_mod=0),
            s: algebra.Differential(0, coeff_mod=0)}
        filtration_levels = {n: 2, s: 1}
        specseq = algebra.SpectralSequence(
            gradings=gradings, differentials=differentials, filtration_levels=filtration_levels, coeff_mod=0)
        p_poly = specseq.poincare_poly()
        expected = (
                (algebra.PAGE_VAR + 1) * algebra.FILTRATION_VAR
                * (1 + algebra.FILTRATION_VAR * algebra.DEGREE_VAR ** 2)
        )
        LOG.info('sphere')
        LOG.info(p_poly)
        LOG.info(expected)
        self.assertTrue(utils.poly_equal(p_poly, expected))

    def test_heart_sphere(self):
        # Try the heart shaped sphere, wth one of the humps higher than the other,
        # each critical point in its own filtration level
        w, x, y, z = sympy.symbols('w,x,y,z')
        gradings = {w: 2, x: 2, y: 1, z: 0}
        differentials = {
            w: algebra.Differential(y, coeff_mod=2),
            x: algebra.Differential(y, coeff_mod=2),
            y: algebra.Differential(0, coeff_mod=2),
            z: algebra.Differential(0, coeff_mod=2)}
        filtration_levels = {w: 4, x: 3, y: 2, z: 1}
        specseq = algebra.SpectralSequence(
            gradings=gradings, differentials=differentials, filtration_levels=filtration_levels, coeff_mod=2)
        p_poly = specseq.poincare_poly()
        # Checked by hand
        p = algebra.FILTRATION_VAR
        t = algebra.DEGREE_VAR
        r = algebra.PAGE_VAR
        expected = p + (t * p**2) + (p**3 + p**4)*t**2 + (r + r**2 + r**3)*(p + (p**4)*(t**2))
        LOG.info('heart sphere')
        LOG.info(p_poly)
        LOG.info(expected)
        self.assertTrue(utils.poly_equal(p_poly, expected))


class TestDGA(unittest.TestCase):

    def test_pickle(self):
        pd = ll.PlatDiagram(n_strands=4, front_crossings=[1, 1, 1], lazy_disks=False, lazy_lch=False, lazy_rsft=True)
        dga = pd.lch_dga
        with tempfile.TemporaryDirectory() as d:
            temp_name = os.path.join(d, 'blah.pk')
            dga.pickle(temp_name, compress=False)
            new_dga = algebra.DGA.from_pickle(temp_name, decompress=True)
        self.assertEqual(dga.augmentations, new_dga.augmentations)


if __name__ == '__main__':
    unittest.main()
