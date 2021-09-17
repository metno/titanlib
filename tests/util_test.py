from __future__ import print_function
import unittest
import titanlib
import numpy as np
import scipy.spatial


"""Convenient vectors used as inputs"""


class UtilTest(unittest.TestCase):
    def test_interpolate_to_points(self):
        olats = [70.1, 70]
        olons = [10.7, 10]
        ilats, ilons = np.meshgrid(np.linspace(50,80,301), np.linspace(10,11, 101))
        ivalues = ilats + ilons
        ovalues = titanlib.interpolate_to_points(ilats, ilons, ivalues, olats, olons)
        np.testing.assert_almost_equal(list(ovalues), [80.8, 80], 4)

    def test_is_valid(self):
        self.assertFalse(titanlib.is_valid(np.nan))
        self.assertFalse(titanlib.is_valid(np.inf))
        self.assertFalse(titanlib.is_valid(float('nan')))
        self.assertTrue(titanlib.is_valid(0))
        self.assertTrue(titanlib.is_valid(1))
        self.assertTrue(titanlib.is_valid(-1))

    def test_subset(self):
        array = [0, 1, 2]
        indices = [1]
        np.testing.assert_almost_equal(titanlib.subset(array, indices), [1])
        np.testing.assert_almost_equal(titanlib.subset([], indices), [])

if __name__ == '__main__':
    unittest.main()
