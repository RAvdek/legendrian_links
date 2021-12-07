import unittest
import legendrian_links as ll


class TestLinks(unittest.TestCase):

    def _test_link(self, n_strands, crossings, n_crossings, n_disks):

        pd = ll.PlatDiagram(n_strands=n_strands, front_crossings=crossings)

        # There is exactly one right close and one left close
        self.assertEqual(len([ps for ps in pd.plat_segments if ps.left_close]), 1)
        self.assertEqual(len([ps for ps in pd.plat_segments if ps.right_close]), 1)

        # The x-values of crossings are all distinct
        crossing_x_values = [c.x for c in pd.crossings]
        self.assertEqual(len(crossing_x_values), len(crossing_x_values))

        # Number of disks and crossings are expected
        self.assertEqual(len(pd.crossings), n_crossings)
        self.assertEqual(len(pd.disks), n_disks)

        # The x-values of disk segments in each disk are all distinct
        for d in pd.disks:
            ds_x_values = [ds.x for ds in d]
            self.assertEqual(len(ds_x_values), len(set(ds_x_values)))

    def test_unknot(self):
        self._test_link(2, None, 1, 2)

    def test_stabilized_unknot(self):
        self._test_link(2, [0], 2, 3)

    def test_very_twisted_unknot(self):
        self._test_link(2, [0]*5, 6, 7)

    def test_twisted_unknot(self):
        self._test_link(4, [1], 3, 4)

    def test_double_unlink(self):
        self._test_link(4, None, 2, 4)

    def test_hopf_link(self):
        self._test_link(4, [1]*2, 4, 7)

    def test_trefoil(self):
        self._test_link(4, [1]*3, 5, 10)


if __name__ == '__main__':
    unittest.main()
