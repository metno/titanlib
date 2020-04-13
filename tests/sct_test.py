from __future__ import print_function
import unittest
import titanlib
import numpy as np


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
len100 = np.random.randn(100)
z100 = np.zeros(100)
o100 = np.ones(100)


class SctTest(unittest.TestCase):
    def test_valid_input(self):
        """Check that the test doesn't fail"""
        N = 100
        lats = np.random.randn(N) * 10000;
        lons = np.random.randn(N) * 10000;
        elevs = np.random.rand(N) * 100;
        values = np.random.randn(N) * 6;
        t2pos = np.ones(N) * 4;
        t2neg = np.ones(N) * 4;
        eps2 = np.ones(N) * 0.5;
        sct, rep, flags = titanlib.sct(lats, lons, elevs, values, 100, 1000, 100, 100, 10000,
                100, t2pos, t2neg, eps2)

    def test_nminprof(self):
        nminprof = 20
        sct, rep, flags = titanlib.sct(len100, len100, len100, z100, 100, 1000,
                nminprof, 100, 10000, 100, 2*o100, 2*o100, 0.5*o100)

if __name__ == '__main__':
    unittest.main()
