from __future__ import print_function
import unittest
import titanlib
import numpy as np
import scipy.spatial
import time
import util


"""Convenient vectors used as inputs"""


class IsolationTest(unittest.TestCase):
    def test_no_elevation(self):
        # These points are roughly 5 km apart, except for the last
        lats = [60, 60, 60, 60, 60, 70]
        lons = [10,10.1,10.2,10.3,10.4, 10]
        elevs = [0, 50, 50, 50, 400, 500]
        points = titanlib.Points(lats, lons, elevs)

        radius = 15000

        # All except the last one have at least 2 neighbours
        flags = titanlib.isolation_check(points, 2, radius)
        self.assertListEqual(list(flags), [0, 0, 0, 0, 0, 1])

        # Test restricting by vertical radius
        vertical_radius = 100
        flags = titanlib.isolation_check(points, 2, radius, vertical_radius)
        self.assertListEqual(list(flags), [0, 0, 0, 0, 1, 1])
        vertical_radius = 0
        flags = titanlib.isolation_check(points, 2, radius, vertical_radius)
        self.assertListEqual(list(flags), [1, 0, 0, 0, 1, 1])

        # Only the middle one has 4 neighbours
        flags = titanlib.isolation_check(points, 4, radius)
        self.assertListEqual(list(flags), [1, 1, 0, 1, 1, 1])

        # None have neighbours within 1 km
        flags = titanlib.isolation_check(points, 1, 1000)
        self.assertListEqual(list(flags), [1, 1, 1, 1, 1, 1])

    def test_vector_version(self):
        lats = [60, 60, 60, 60, 60, 70]
        lons = [10,10.1,10.2,10.3,10.4, 10]
        elevs = [0, 100, 200, 100, 0, 100]
        N = len(lats)
        points = titanlib.Points(lats, lons, elevs)

        radius = [12000, 8000, 12000, 8000, 12000, 10000000]
        vertical_radius = [110, 110, 110, 200, 100, 50]
        # Expected num neighbours:
        expected_num = np.array([1, 2, 2, 2, 1, 2])

        num_min = np.array([3, 4, 4, 4, 3, 4])
        flags = titanlib.isolation_check(points, num_min, radius, vertical_radius)
        np.testing.assert_array_almost_equal(list(flags), expected_num < num_min)

        num_min = np.array([1, 2, 2, 2, 1, 2])
        flags = titanlib.isolation_check(points, num_min, radius, vertical_radius)
        np.testing.assert_array_almost_equal(list(flags), expected_num < num_min)

        num_min = expected_num
        flags = titanlib.isolation_check(points, num_min, radius, vertical_radius)
        np.testing.assert_array_almost_equal(list(flags), expected_num < num_min)

        num_min = np.array([1, 1, 1, 1, 1, 1])
        flags = titanlib.isolation_check(points, num_min, radius, vertical_radius)
        np.testing.assert_array_almost_equal(list(flags), expected_num < num_min)

    def test_empty_vertical_radius(self):
        lats = [60, 60, 60, 60, 60, 70]
        lons = [10,10.1,10.2,10.3,10.4, 10]
        elevs = [0, 100, 200, 100, 0, 100]
        N = len(lats)
        points = titanlib.Points(lats, lons, elevs)

        radius = [12000, 8000, 12000, 8000, 12000, 10000000]
        # Expected num neighbours:
        expected_num = np.array([2, 2, 4, 2, 2, 5])

        num_min = np.array([3, 4, 4, 4, 3, 4])
        flags = titanlib.isolation_check(points, num_min, radius)
        np.testing.assert_array_almost_equal(list(flags), expected_num < num_min)

    def test_vector_invalid_arguments(self):
        lats = [60, 60, 60, 60, 60, 70]
        lons = [10,10.1,10.2,10.3,10.4, 10]
        elevs = [0, 100, 200, 100, 0, 100]
        N = len(lats)
        points = titanlib.Points(lats, lons, elevs)

        num_min = np.full(N, 2)
        radius = np.full(N, 15000)
        vertical_radius = np.full(N, 100)

        with self.assertRaises(ValueError):
            flags = titanlib.isolation_check(points, num_min, radius, [1])

        with self.assertRaises(ValueError):
            flags = titanlib.isolation_check(points, [1], radius, vertical_radius)

        with self.assertRaises(ValueError):
            flags = titanlib.isolation_check(points, num_min, [1], vertical_radius)


if __name__ == '__main__':
    unittest.main()
