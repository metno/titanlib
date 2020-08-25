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
obs_to_check = [1, 1, 1]
background_values = [0, 0, 0]
background_elab_type = "vertical_profile"
N = len(lats)
num_min = 3
num_max = 10
inner_radius = 10000
outer_radius = 10000
num_iterations = 1
num_min_prof = 0
min_elev_diff = 100
min_horizontal_scale = 10000
max_horizontal_scale = 100000
kth_closest_obs_horizontal_scale = 2
vertical_scale = 200
tpos_score = np.ones(N) * 2;
tneg_score = np.ones(N) * 2;
t_sod = np.ones(N) * 4;
eps2 = np.ones(N) * 0.5;


class SctTest(unittest.TestCase):
    def test_simple(self):
        elevs0 = [0, 0, 0]
        values0 = [0, 1, 100]

        flags, sct, rep, sod, num_inner, horizontal_scale, an_inc, an_res, cv_res, innov, idi, idiv, sig2o = titanlib.sct(lats, lons, elevs0, values0, obs_to_check, background_values, background_elab_type, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, eps2, tpos_score, tneg_score, t_sod)
        np.testing.assert_array_equal(flags, [0, 0, 1])

    def test_time(self):
        s_time = time.time()

        N = 3000
        lats = np.random.randn(N) * 1;
        lons = np.random.randn(N) * 1;
        elevs = np.random.rand(N) * 100;
        values = np.random.randn(N) * 6;
        obs_to_check = np.ones(N);
        background_values = np.random.randn(N) * 6;
        background_elab_type = "vertical_profile";
        tpos_score = np.ones(N) * 4;
        tneg_score = np.ones(N) * 4;
        t_sod = np.ones(N) * 2;
        eps2 = np.ones(N) * 0.5;
        
        flags, sct, rep, sod, num_inner, horizontal_scale, an_inc, an_res, cv_res, innov, idi, idiv, sig2o = titanlib.sct(lats, lons, elevs, values, obs_to_check, background_values, background_elab_type, 5, 100, 4000, 10000, 2, 30, 100, 1000, 10000, 3, 100, eps2, tpos_score, tneg_score, t_sod)

        print("%.1fs" % (time.time() - s_time))

if __name__ == '__main__':
    unittest.main()
