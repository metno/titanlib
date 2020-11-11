from __future__ import print_function
import unittest
import titanlib
import numpy as np
import util
import time


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
ut = 1546300800


class BuddyEventCheckTest(unittest.TestCase):
    def test_1(self):
        """Check that the test doesn't fail"""
        lats, lons, elevs, values = util.summer_temperature_example()
        points = titanlib.Points(lats, lons, elevs)
        I = slice(0, len(lats))
        lats = lats[I]
        lons = lons[I]
        elevs = elevs[I]
        values = values[I]
        s_time = time.time()
        event_thr = 0.2
        thr = 0.25
        elev_gradient = -0.0065
        num_iterations = 1
        flags = titanlib.buddy_event_check(points, values, [10000], [5], event_thr, thr, 5, elev_gradient, num_iterations)
        e_time = time.time()
        print(e_time - s_time)
        print("Fraction of stations removed: %.1f %%" % (np.mean(flags) * 100))


if __name__ == '__main__':
    unittest.main()
