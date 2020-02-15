from __future__ import print_function
import unittest
import titanlib
import numpy as np
import scipy.spatial
import time
import util


"""Convenient vectors used as inputs"""


class IsolationTest(unittest.TestCase):
    def test_1(self):
        """Check that the test doesn't fail"""
        nmin = 5
        radius = 15000
        lats, lons, elevs, values = util.summer_temperature_example()
        dz = 100
        s_time = time.time()
        status, flags = titanlib.isolation_check(lats, lons, elevs, nmin, radius, dz)
        e_time = time.time()
        print("%.2f %d" %(e_time - s_time, np.sum(flags)))
        s_time = time.time()
        status, flags = titanlib.isolation_check(lats, lons, nmin, radius)
        e_time = time.time()
        print("%.2f %d" %(e_time - s_time, np.sum(flags)))
        print(np.sum(flags))


if __name__ == '__main__':
    unittest.main()
