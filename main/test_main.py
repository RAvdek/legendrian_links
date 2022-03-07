import json
import unittest
import numpy as np
import sympy
import polynomials
import legendrian_links as ll
from algebra import Matrix


def comparable_list_of_dicts(l):
    # Have to convert to sets of json strings as we cannot make sets of dicts in TestCase.AssertEqual
    return sorted([json.dumps({str(k): v for k, v in d.items()}, sort_keys=True) for d in l])


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


class TestAugsearch(unittest.TestCase):

    def test_trefoil(self):
        x, y, z = sympy.symbols('x,y,z')
        symbols = [x, y, z]
        polys = [1 + x + z + x * y * z]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2)
        expected_results = [
            {x: 0, z: 1, y: 0},
            {x: 0, z: 1, y: 1},
            {x: 1, y: 0, z: 0},
            {x: 1, y: 1, z: 0},
            {x: 1, y: 1, z: 1}
        ]
        self.assertEqual(comparable_list_of_dicts(results), comparable_list_of_dicts(expected_results))

    def test_hopf_link(self):
        x, y = sympy.symbols('x,y')
        symbols = [x, y]
        polys = [y * x]
        results = polynomials.zero_set(polys=polys, symbols=symbols, modulus=2)
        expected_results = [
            {x: 1, y: 0},
            {x: 0, y: 0},
            {x: 0, y: 1}
        ]
        self.assertEqual(comparable_list_of_dicts(results), comparable_list_of_dicts(expected_results))


class TestMatrix(unittest.TestCase):

    def test_modulus(self):
        self.assertFalse(
            np.any(Matrix([[4]], coeff_mod=3).values - Matrix([[1]], coeff_mod=3).values)
        )
        self.assertFalse(
            np.any(Matrix([[1, 2], [2, 2]], coeff_mod=2).values - Matrix([[1, 0], [0, 0]], coeff_mod=2).values)
        )


if __name__ == '__main__':
    unittest.main()
