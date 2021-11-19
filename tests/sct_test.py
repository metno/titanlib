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
    def test_invalid_arguments(self):
        points = titanlib.Points(lats, lons, elevs)
        N = 3
        args0 = [points, values, num_min, num_max, inner_radius,
                outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale,
                vertical_scale, pos, neg, eps2]
        invalid_arguments = [(2, -1), (3, -1), (4, -1000),
                             (5, -1000), (6, -1), (7, -1), (8, -1), (9, -1),
                             (10, -1), (11, [-1] * N), (12, [-1] * N), (13, [-1] * N), (13, [0] * N)]
        for invalid_argument in invalid_arguments:
            with self.subTest(argument=invalid_argument[0]):
                with self.assertRaises(ValueError):
                    args = args0
                    args[invalid_argument[0]] = invalid_argument[1]
                    titanlib.sct(*args)

    def test_simple(self):
        elevs0 = [0, 0, 0]
        values0 = [0, 1, 100]
        points0 = titanlib.Points(lats, lons, elevs0)

        flags, sct, rep = titanlib.sct(points0, values0, num_min, num_max, inner_radius,
                outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale,
                vertical_scale, pos, neg, eps2)
        np.testing.assert_array_equal(flags, [0, 0, 1])

    def test_simple_dataset(self):
        elevs0 = [0, 0, 0]
        values0 = [0, 1, 100]
        points0 = titanlib.Points(lats, lons, elevs0)

        dataset = titanlib.Dataset(points0, values0)

        sct, rep = dataset.sct(num_min, num_max, inner_radius,
                outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale,
                vertical_scale, pos, neg, eps2)
        np.testing.assert_array_equal(dataset.flags, [0, 0, 1])


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

    def test_obs_to_check(self):
        elevs0 = [0, 0, 0, 0, 0, 0]
        values0 = [0, 1, 1, 1, 100, 100]
        points0 = titanlib.Points(lats*2, lons*2, elevs0)
        obs_to_check0 = [0, 1, 1, 1, 1, 0]

        flags, sct, rep = titanlib.sct(points0, values0, num_min, num_max, inner_radius,
                outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale,
                vertical_scale, np.ones(6), np.ones(6), 0.5 * np.ones(6), obs_to_check0)
        np.testing.assert_array_equal(flags, [0, 0, 0, 0, 1, 0])

    def test_sct_identical_elevations(self):
        lats = [0, 0]
        lons = [0, 0]
        elevs = [0, 0]
        values = [0, 0]
        num_min_prof = 0
        min_elev_diff = 0
        points = titanlib.Points(lats, lons, elevs)
        t2pos = [1, 1]
        t2neg = [1, 1]
        eps2 = [1, 1]


if __name__ == '__main__':
    unittest.main()
