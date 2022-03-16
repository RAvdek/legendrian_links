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


LOG = utils.LOG


def comparable_list_of_dicts(l):
    # Have to convert to sets of json strings as we cannot make sets of dicts in TestCase.AssertEqual
    return sorted([json.dumps({str(k): str(v) for k, v in d.items()}, sort_keys=True) for d in l])


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

    def test_non_zero_const(self):
        polys = [1]
        output = polynomials.zero_set(polys=polys, symbols=self.symbols)
        self.assertEqual(len(output), 0)

    def test_zero_const(self):
        polys = [0]
        output = polynomials.zero_set(polys=polys, symbols=self.symbols, allow_unset=False)
        self.assertEqual(len(output), 8)
        output = polynomials.zero_set(polys=polys, symbols=self.symbols, allow_unset=True)
        self.assertEqual(len(output), 1)

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
