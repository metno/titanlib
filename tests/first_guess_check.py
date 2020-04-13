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


class FirstGuessCheckTest(unittest.TestCase):
    def test_1(self):
        """Check that the test doesn't fail"""
        values = [10, 13]
        fg = [10, 10]
        fg_min = [3]
        fg_max = [3, 2]

        s_time = time.time()
        flags = titanlib.first_guess_check(values, fg, fg_min, fg_max)
        e_time = time.time()
        print(e_time - s_time)
        print("Fraction of stations removed: %.1f %%" % (np.mean(flags) * 100))


if __name__ == '__main__':
    unittest.main()
