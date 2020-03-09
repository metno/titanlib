from __future__ import print_function
import unittest
import titanlib
import numpy as np


lats = [60, 60, 60, 60, 60, 70]
lons = [10,10.1,10.2,10.3,10.4, 10]


class KDTreeTest(unittest.TestCase):
    def test_get_neighbours(self):
        """Check that the test doesn't fail"""
        tree = titanlib.KDTree(lats, lons)
        neighbours = list(tree.get_neighbours(60, 10.101, 10000, 0, False))
        neighbours.sort()
        self.assertListEqual(neighbours, [0,1,2])
        neighbours = list(tree.get_neighbours(60, 10.101, 10000, 0, True))
        neighbours.sort()
        self.assertListEqual(neighbours, [0,1,2])

    def test_get_neighbours_with_num(self):
        tree = titanlib.KDTree(lats, lons)
        neighbours = list(tree.get_neighbours(60, 10.101, 10000, 1, False))
        neighbours.sort()
        self.assertListEqual(neighbours, [1])

    def test_get_neighbours_with_match(self):
        tree = titanlib.KDTree(lats, lons)
        neighbours = list(tree.get_neighbours(60, 10.1, 10000, 0, False))
        neighbours.sort()
        self.assertListEqual(neighbours, [0, 2])
        neighbours = list(tree.get_neighbours(60, 10.1, 10000, 0, True))
        neighbours.sort()
        self.assertListEqual(neighbours, [0, 1, 2])

    def test_get_neighbours_with_distance(self):
        tree = titanlib.KDTree(lats, lons)
        neighbours, distances = tree.get_neighbours_with_distance(60, 10.101, 10000, 0, False)
        self.assertEqual(len(distances), 3)
        self.assertAlmostEqual(np.min(distances), 55.506, 2)
        self.assertAlmostEqual(np.max(distances), 5614, 0)

    def test_get_num_neighbours(self):
        tree = titanlib.KDTree(lats, lons)
        self.assertEqual(tree.get_num_neighbours(60, 10.101, 10000, 0, False), 3)
        self.assertEqual(tree.get_num_neighbours(60, 10.101, 100, 0, False), 1)
        self.assertEqual(tree.get_num_neighbours(60, 10.101, 2000, 0, False), 1)

    def test_get_itself(self):
        tree = titanlib.KDTree(lats, lons)
        self.assertEqual(tree.get_nearest_neighbour(60, 10.1, True), 1)
        self.assertEqual(tree.get_nearest_neighbour(60, 10.101, True), 1)

    def test_get_neearest_neighbour(self):
        tree = titanlib.KDTree(lats, lons)
        index = tree.get_nearest_neighbour(60, 10.1, True)
        self.assertEqual(index, 1)
        indices = tree.get_nearest_neighbour([60], [10.1], True)
        self.assertListEqual(list(indices), [1])

    def test_circle(self):
        """If a circle is used, only the first and last point should match"""
        tree = titanlib.KDTree([0, 1, 1] , [1, 0, 1])
        self.assertEqual(tree.get_num_neighbours(0, 0, 120000, 0, False), 2)

    def test_matching(self):
        tree = titanlib.KDTree([0, 0, 1, 1] , [0, 1.0001, 0, 1])
        self.assertEqual(tree.get_num_neighbours(0, 0, 120000, 0, False), 2)
        self.assertEqual(tree.get_num_neighbours(0, 0, 120000, 0, True), 3)
        """Check that match work when limiting the number of stations"""
        self.assertEqual(tree.get_num_neighbours(0, 0, 120000, 1, False), 1)
        self.assertEqual(tree.get_num_neighbours(0, 0, 120000, 1, True), 1)
        self.assertEqual(tree.get_neighbours(0, 0, 120000, 1, False)[0], 2)
        self.assertEqual(tree.get_neighbours(0, 0, 120000, 1, True)[0], 0)
        """Check that match works When limiting the radius"""
        self.assertEqual(tree.get_num_neighbours(0, 0, 0, 1, False), 1)
        self.assertEqual(tree.get_num_neighbours(0, 0, 0, 1, True), 1)
        self.assertEqual(tree.get_neighbours(0, 0, 0, 1, False)[0], 2)
        self.assertEqual(tree.get_neighbours(0, 0, 0, 1, True)[0], 0)


if __name__ == '__main__':
    unittest.main()
