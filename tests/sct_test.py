from __future__ import print_function
import unittest
import titanlib
import numpy as np
import time


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
len100 = np.random.randn(100)
z100 = np.zeros(100)
o100 = np.ones(100)


class SctTest(unittest.TestCase):
    def test_time(self):
        s_time = time.time()

        N = 30000
        lats = np.random.randn(N) * 1;
        lons = np.random.randn(N) * 1;
        elevs = np.random.rand(N) * 100;
        values = np.random.randn(N) * 6;
        t2pos = np.ones(N) * 4;
        t2neg = np.ones(N) * 4;
        eps2 = np.ones(N) * 0.5;
        flags, sct, rep = titanlib.sct(lats, lons, elevs, values, 5, 100, 4000, 10000, 2, 30, 100, 1000, 100, t2pos, t2neg, eps2)

        print("%.1fs" % (time.time() - s_time))

if __name__ == '__main__':
    unittest.main()
