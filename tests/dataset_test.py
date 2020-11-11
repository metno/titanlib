from __future__ import print_function
import unittest
import titanlib
import numpy as np


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]


class DatasetTest(unittest.TestCase):
    def test_1(self):
        """Check that the test doesn't fail"""
        N = 5
        np.random.seed(0)
        lats = np.random.randn(N) * 10000;
        lons = np.random.randn(N) * 10000;
        elevs = np.random.rand(N) * 100;
        points = titanlib.Points(lats, lons, elevs)
        values = np.random.randn(N) * 6;
        dataset = titanlib.Dataset(points, values);
        dataset.range_check([-100], [-100], [0, 1, 2])
        dataset.range_check([-100], [-100], [1, 2, 4])
        self.assertListEqual([i for i in dataset.flags], [1, 1, 1, 0, 1])


if __name__ == '__main__':
    unittest.main()
