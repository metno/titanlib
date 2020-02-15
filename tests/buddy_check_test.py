from __future__ import print_function
import unittest
import titanlib
import numpy as np
import util


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
ut = 1546300800


class BuddyCheckTest(unittest.TestCase):
    def test_1(self):
        """Check that the test doesn't fail"""
        #lats, lons, elevs, values = util.summer_temperature_example()
        #status, flags = titanlib.buddy_check(lats, lons, elevs, values, [100000], [0], [5], [5], 5, False)
        #print("Number of stations: %d" % np.sum(flags))
        pass


if __name__ == '__main__':
    unittest.main()
