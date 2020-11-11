from __future__ import print_function
import unittest
import titanlib
import numpy as np
import util
import time


"""Convenient vectors used as inputs"""
N = 10
lats = [60]*N
lons = np.linspace(60, 60.001, N)
elevs = [0]*10
points = titanlib.Points(lats, lons, elevs)
values = [0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1]
N = len(lats)
radius = [10000]
num_min = [1]
threshold = 1
elev_gradient = -0.0065
max_elev_diff = 200
min_std = 0.01
num_iterations = 2


class BuddyCheckTest(unittest.TestCase):
    def test_1(self):
        print(values)
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        print(flags)
        np.testing.assert_array_equal(flags, [0]*8 + [1]*2)

    def test_min_std(self):
        """Check that when min_std is high enough, nothing is flagged"""
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, 0.3, num_iterations)
        np.testing.assert_array_equal(flags, [0]*9+[1])
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, 1, num_iterations)
        np.testing.assert_array_equal(flags, [0]*N)

    def test_min_num(self):
        """Check that when min_num is high enough, nothing is flagged"""
        flags = titanlib.buddy_check(points, values, radius, [20], threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*N)

    def test_num_iterations(self):
        """Check that min_iterations affects the flags"""
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, 1)
        np.testing.assert_array_equal(flags, [0]*9 + [1])

    def test_elev_gradient(self):
        elevs0 = [0]*9 + [-153.8]
        points0 = titanlib.Points(lats, lons, elevs0)
        flags = titanlib.buddy_check(points0, values, radius, num_min, threshold,
                max_elev_diff, 0, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*8 + [1, 1])

        flags = titanlib.buddy_check(points0, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*8 + [1, 0])

    def test_max_elev_diff(self):
        """Check that test is not run on a point which has no other points within the elevation range"""
        elevs0 = [0]*9 + [100]
        points0 = titanlib.Points(lats, lons, elevs0)
        flags = titanlib.buddy_check(points0, values, radius, num_min, threshold,
                1, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*8 + [1, 0])

if __name__ == '__main__':
    unittest.main()
