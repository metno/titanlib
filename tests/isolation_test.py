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
        points = titanlib.Points(lats, lons)

        radius = 15000

        # All except the last one have at least 2 neighbours
        flags = titanlib.isolation_check(points, 2, radius)
        self.assertListEqual(list(flags), [0, 0, 0, 0, 0, 1])

        # Only the middle one has 4 neighbours
        flags = titanlib.isolation_check(points, 4, radius)
        self.assertListEqual(list(flags), [1, 1, 0, 1, 1, 1])

        # None have neighbours within 1 km
        flags = titanlib.isolation_check(points, 1, 1000)
        self.assertListEqual(list(flags), [1, 1, 1, 1, 1, 1])

    def test_summer_case(self):
        """Check that the test doesn't fail"""
        nmin = 5
        radius = 15000
        lats, lons, elevs, values = util.summer_temperature_example()
        points = titanlib.Points(lats, lons, elevs)
        dz = 100
        s_time = time.time()
        flags = titanlib.isolation_check(points, nmin, radius, dz)
        e_time = time.time()
        print("%.2f %d" %(e_time - s_time, np.sum(flags)))
        s_time = time.time()
        flags = titanlib.isolation_check(points, nmin, radius)
        e_time = time.time()
        print("%.2f %d" %(e_time - s_time, np.sum(flags)))
        print(np.sum(flags))


if __name__ == '__main__':
    unittest.main()
