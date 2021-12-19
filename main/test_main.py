import unittest
import legendrian_links as ll


class TestLinks(unittest.TestCase):

    def _test_link(self, n_strands, front_crossings, n_knots, n_disks, tb_rot=None):
        # TODO: Should be able to test for rot and tb.

        pd = ll.PlatDiagram(n_strands=n_strands, front_crossings=front_crossings)

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
        self._test_link(2, [], 1, 2, tb_rot=[-1, [0]])

    def test_stabilized_unknot(self):
        self._test_link(2, [0], 1, 3, tb_rot=[-2, [1, -1]])

    def test_twisted_unknot(self):
        self._test_link(4, [1], 1, 4)

    def test_very_twisted_unknot(self):
        self._test_link(2, [0]*5, 1, 7, tb_rot=[-6, [1, -1]])

    def test_double_unlink(self):
        self._test_link(4, [], 2, 4)

    def test_hopf_link(self):
        self._test_link(4, [1]*2, 2, 7)

    def test_trefoil(self):
        self._test_link(4, [1]*3, 1, 10, tb_rot=[1, [0]])


if __name__ == '__main__':
    unittest.main()
