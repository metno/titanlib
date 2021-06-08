from __future__ import print_function
import unittest
import titanlib
import numpy as np
import time



class SctResistantTest(unittest.TestCase):
#    def test_simple(self):

#        """Convenient vectors used as inputs"""
#        lats = [60, 60, 60, 60, 60]
#        lons = [10, 10.005, 10.01, 10.015, 10.02]
#        elevs = [0, 0, 0, 0, 0]
#        values = [0, 1, 10, -10, -100]
#        obs_to_check = [1, 1, 1, 1, 1]
#        background_values = [0, 0, 0, 0, 0]
#        background_elab_type = "vertical_profile"
#        N = len(lats)
#        print(" -- Test over", N, "observations --") 
#        num_min_outer = 3
#        num_max_outer = 10
#        inner_radius = 50000
#        outer_radius = 50000
#        num_iterations = 10
#        num_min_prof = 0
#        min_elev_diff = 100
#        min_horizontal_scale = 10000
#        max_horizontal_scale = 100000
#        kth_closest_obs_horizontal_scale = 2
#        vertical_scale = 200
#        tpos = np.ones(N) * 16
#        tneg = np.ones(N) * 16
#        eps2 = np.ones(N) * 0.5
#        values_mina = values - 20 * np.ones(N)
#        values_maxa = values + 20 * np.ones(N)
#        values_minv = values - 1 * np.ones(N)
#        values_maxv = values + 1 * np.ones(N)
#        debug = False
#
#        points = titanlib.Points(lats, lons, elevs)
#        flags, score = titanlib.sct_resistant(points, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, values_mina, values_maxa, values_minv, values_maxv, eps2, tpos, tneg, debug)
#        print(flags)
#
#        np.testing.assert_array_equal(flags, [0, 0, 0, 0, 1])
#
#        print("===============================================================")

    def test_time(self):
        s_time = time.time()
        N = 5000
        print(" -- Test over", N, "observations --") 
        lats = np.ones(N) * 55 + np.random.random(N) * ( 70 - 55)
        lons = np.ones(N) *  5 + np.random.random(N) * ( 30 -  5)
        elevs = np.random.random(N) * 2500
        values =  np.ones(N) * 30 -0.0065 * elevs
        pGE = 0.3
        nGE = int( np.ceil( N * pGE))
        idx = np.random.randint(0, high=(N-1), size=nGE)
        obs_to_check = [1] * N
        true_GE = [0] * N
        for i in idx: 
            values[i] = np.random.random(1) * 100 - 50
            true_GE[i] = 1 
        for i in range(1, N):
            obs_to_check[i] = int(1)
        background_values = np.random.randn(N) * 6
        background_elab_type = titanlib.VerticalProfileTheilSen
        num_min_outer = 3
        num_max_outer = 50
        inner_radius = 30000
        outer_radius = 50000
        num_iterations = 100
        num_min_prof = 10
        min_elev_diff = 500
        min_horizontal_scale = 500
        max_horizontal_scale = 10000
        kth_closest_obs_horizontal_scale = 3
        vertical_scale = 600
        tpos = np.ones(N) * 3
        tneg = np.ones(N) * 3
        eps2 = np.ones(N) * 0.5
        values_mina = values - 20 * np.ones(N)
        values_maxa = values + 20 * np.ones(N)
        values_minv = values - 1 * np.ones(N)
        values_maxv = values + 1 * np.ones(N)
        debug = False
        basic = False
        flags = np.ones(N) * (-999.)

        points = titanlib.Points(lats, lons, elevs)
        flags, score = titanlib.sct_resistant(points, values, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, values_mina, values_maxa, values_minv, values_maxv, eps2, tpos, tneg, debug, basic)

        print("%.1fs" % (time.time() - s_time))

        a = 0
        b = 0
        c = 0
        d = 0 
        for i in range(1, N):
           if ( flags[i] == 1 and true_GE[i] == 1):
               a = a + 1
           if ( flags[i] == 1 and true_GE[i] == 0):
               b = b + 1
           if ( flags[i] == 0 and true_GE[i] == 1): 
               c = c + 1
           if ( flags[i] == 0 and true_GE[i] == 0): 
               d = d + 1
#           print( "flags[i] true_GE[i] a(bad) b c d:", flags[i], true_GE[i], a,"(",nGE,")", b, c, d, a+b+c+d)
        rand = (a+c) * (a+b) / (a+b+c+d)
        ets = (a-rand) / (a+b+c-rand)
        acc = (a+d)/(a+b+c+d)
        pod = a/(a+c)
        pofa = b/(b+d)
        print( "a(bad) b c d:", a,"(",nGE,")", b, c, d, a+b+c+d)
#        print( "acc pod pofa ets", acc, round(pod,2), round(pofa,2), round(ets,2))
        print( "acc pod pofa ets", round(acc,2), pod, pofa, ets)
 
        print("===============================================================")

if __name__ == '__main__':
    unittest.main()
