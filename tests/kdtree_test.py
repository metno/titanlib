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
        neighbours = list(tree.get_neighbours(60, 10.101, 10000))
        neighbours.sort()
        self.assertListEqual(neighbours, [0,1,2])

    def test_get_neighbours_with_num(self):
        tree = titanlib.KDTree(lats, lons)
        neighbours = list(tree.get_neighbours(60, 10.101, 10000, 1))
        neighbours.sort()
        self.assertListEqual(neighbours, [1])

    def test_get_neighbours_with_distance(self):
        tree = titanlib.KDTree(lats, lons)
        neighbours, distances = tree.get_neighbours_with_distance(60, 10.101, 10000)
        self.assertEqual(len(distances), 3)
        self.assertAlmostEqual(np.min(distances), 55.506, 2)
        self.assertAlmostEqual(np.max(distances), 5614, 0)

    def test_get_num_neighbours(self):
        tree = titanlib.KDTree(lats, lons)
        self.assertEqual(tree.get_num_neighbours(60, 10.101, 10000), 3)
        self.assertEqual(tree.get_num_neighbours(60, 10.101, 100), 1)
        self.assertEqual(tree.get_num_neighbours(60, 10.101, 2000), 1)

    def test_get_closest_neighbours(self):
        tree = titanlib.KDTree(lats, lons)
        neighbours = list(tree.get_closest_neighbours(60, 10.101, 1))
        neighbours.sort()
        self.assertListEqual(neighbours, [1])
        neighbours = list(tree.get_closest_neighbours(60, 10.09, 2))
        neighbours.sort()
        self.assertListEqual(neighbours, [0,1])
        neighbours = list(tree.get_closest_neighbours(60, 10.1, 6))
        neighbours.sort()
        self.assertListEqual(neighbours, [0, 2, 3, 4, 5])
        neighbours = list(tree.get_closest_neighbours(60, 10.101, 6))
        neighbours.sort()
        self.assertListEqual(neighbours, [0, 1, 2, 3, 4, 5])
        neighbours = list(tree.get_closest_neighbours(60, 10.1, 12))
        neighbours.sort()
        self.assertListEqual(neighbours, [0, 2, 3, 4, 5])
        neighbours = list(tree.get_closest_neighbours(60, 10.101, 12))
        neighbours.sort()
        self.assertListEqual(neighbours, [0, 1, 2, 3, 4, 5])

    def test_get_itself(self):
        tree = titanlib.KDTree(lats, lons)
        self.assertEqual(tree.get_nearest_neighbour(60, 10.1), 2)
        self.assertEqual(tree.get_nearest_neighbour(60, 10.101), 1)


if __name__ == '__main__':
    unittest.main()
