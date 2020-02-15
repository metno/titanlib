from __future__ import print_function
import unittest
import titanlib
import numpy as np
import scipy.spatial


"""Convenient vectors used as inputs"""


class UtilTest(unittest.TestCase):
    def test_1(self):
        """Check that the test doesn't fail"""
        N = 5
        np.random.seed(0)
        lats = [0, 1];
        lons = [0, 0];
        status, x, y, z = titanlib.convert_coordinates(lats, lons)
        print(x, y, z)
        print(np.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2))

if __name__ == '__main__':
    unittest.main()
