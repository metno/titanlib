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
value_min = -50;
value_max = 50;
sig2o_min = 0.01;
sig2o_max = 20;
debug = True;

class SctTest(unittest.TestCase):
    def test_simple(self):
        elevs0 = [0, 0, 0]
        values0 = [0, 1, 100]

        flags, sct, rep, sod, num_inner, horizontal_scale, an_inc, an_res, cv_res, innov, idi, idiv, sig2o = titanlib.sct(lats, lons, elevs0, values0, obs_to_check, background_values, background_elab_type, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, value_min, value_max, sig2o_min, sig2o_max, eps2, tpos_score, tneg_score, t_sod, debug)
        np.testing.assert_array_equal(flags, [0, 0, 1])

    def test_time(self):
        s_time = time.time()

        P = 3000
        lats1 = np.random.randn(P) * 1;
        lons1 = np.random.randn(P) * 1;
        elevs1 = np.random.rand(P) * 100;
        values1 = np.random.randn(P) * 6;
        obs_to_check1 = np.ones(P) * 1;
        background_values1 = np.random.randn(P) * 6;
        background_elab_type1 = "vertical_profile";
        tpos_score1 = np.ones(P) * 4;
        tneg_score1 = np.ones(P) * 4;
        t_sod1 = np.ones(P) * 2;
        eps21 = np.ones(P) * 0.5;
        value_min1 = -50;
        value_max1 = 50;
        sig2o_min1 = 0.01;
        sig2o_max1 = 20;
        debug1 = False;

        
        flags, sct, rep, sod, num_inner, horizontal_scale, an_inc, an_res, cv_res, innov, idi, idiv, sig2o = titanlib.sct(lats1, lons1, elevs1, values1, obs_to_check1, background_values1, background_elab_type1, 5, 100, 4000, 50000, 2, 30, 100, 1000, 10000, 3, 100, value_min1, value_max1, sig2o_min1, eps21, sig2o_max1, tpos_score1, tneg_score1, t_sod1, debug1)

        print("%.1fs" % (time.time() - s_time))

if __name__ == '__main__':
    unittest.main()
