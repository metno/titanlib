from __future__ import print_function
import unittest
import titanlib
import numpy as np
import time



class SctDualTest(unittest.TestCase):

    def test_time(self):
        s_time = time.time()
        N = 100000
        print(" -- Test over", N, "observations --") 
        lats = np.ones(N) * 55 + np.random.random(N) * ( 70 - 55)
        lons = np.ones(N) *  5 + np.random.random(N) * ( 30 -  5)
        elevs = np.random.random(N) * 2500
        values = np.ones(N) * 10
        for i in range(1, N):
            if ( lons[i] >= 17.5):
                values[i] = 0
        pGE = 0.01
        nGE = int( np.ceil( N * pGE))
        idx = np.random.randint(0, high=(N-1), size=nGE)
        obs_to_check = [1] * N
        true_GE = [0] * N
        for i in idx:
            if ( values[i] == 10):
                values[i] = 0
            if ( values[i] == 0):
                values[i] = 10
            true_GE[i] = 1 
        for i in range(1, N):
            obs_to_check[i] = int(1)
        event_thresholds = np.ones(N) * 0.1
        condition = titanlib.Gt
        test_thresholds = np.ones(N) * 0.5
        num_min_outer = 3
        num_max_outer = 50
        inner_radius = 30000
        outer_radius = 50000
        num_iterations = 100
        min_horizontal_scale = 500
        max_horizontal_scale = 10000
        kth_closest_obs_horizontal_scale = 3
        vertical_scale = 600
        debug = False
        flags = np.ones(N) * (-999.)

        points = titanlib.Points(lats, lons, elevs)
        flags = titanlib.sct_dual(points, values, obs_to_check, event_thresholds, condition, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, test_thresholds, debug)

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
