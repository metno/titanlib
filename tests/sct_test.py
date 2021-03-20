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

lats = [60, 60, 60]
lons = [10, 10.01, 10.02]
elevs = [0, 0, 0]
values = [0, 0, -100]
N = len(lats)
num_min = 3
num_max = 10
inner_radius = 10000
outer_radius = 10000
num_iterations = 1
num_min_prof = 0
min_elev_diff = 100
min_horizontal_scale = 10000
vertical_scale = 200
pos = np.ones(N) * 2;
neg = np.ones(N) * 2;
eps2 = np.ones(N) * 0.5;


class SctTest(unittest.TestCase):
    def test_simple(self):
        elevs0 = [0, 0, 0]
        values0 = [0, 1, 100]
        points0 = titanlib.Points(lats, lons, elevs0)

        flags, sct, rep = titanlib.sct(points0, values0, num_min, num_max, inner_radius,
                outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale,
                vertical_scale, pos, neg, eps2)
        np.testing.assert_array_equal(flags, [0, 0, 1])

    def test_time(self):
        s_time = time.time()

        N = 3000
        lats = np.random.randn(N) * 1;
        lons = np.random.randn(N) * 1;
        elevs = np.random.rand(N) * 100;
        points = titanlib.Points(lats, lons, elevs)
        values = np.random.randn(N) * 6;
        pos = np.ones(N) * 4;
        neg = np.ones(N) * 4;
        eps2 = np.ones(N) * 0.5;
        flags, sct, rep = titanlib.sct(points, values, 5, 100, 4000, 10000, 2, 30, 100, 1000, 100, pos, neg, eps2)

        print("%.1fs" % (time.time() - s_time))

if __name__ == '__main__':
    unittest.main()
