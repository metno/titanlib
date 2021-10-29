from __future__ import print_function
import unittest
import titanlib
import numpy as np


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]


class DatasetTest(unittest.TestCase):
    def get_dataset(self):
        N = 6
        lats = np.linspace(0, 0.1, N);
        lons = np.zeros(N);
        elevs = np.ones(N) * 100
        values = np.zeros(N);
        points = titanlib.Points(lats, lons, elevs, elevs)
        dataset = titanlib.Dataset(points, values);
        return dataset

    def test_1(self):
        """Check that the test doesn't fail"""
        N = 5
        np.random.seed(0)
        lats = np.random.randn(N) * 10000;
        lons = np.random.randn(N) * 10000;
        elevs = np.random.rand(N) * 100;
        points = titanlib.Points(lats, lons, elevs)
        values = np.zeros(N);
        dataset = titanlib.Dataset(points, values);
        dataset.range_check([-100], [-100], [0, 1, 2])
        dataset.range_check([-100], [-100], [1, 2, 4])
        self.assertListEqual([i for i in dataset.flags], [1, 1, 1, 0, 1])

    def test_empty_indices(self):
        dataset = self.get_dataset()
        indices = []
        dataset.isolation_check(2, 2500, 0)
        # self.assertListEqual([i for i in dataset.flags], [1, 0, 0, 0, 0, 1])
        dataset.isolation_check(2, 2500, 100, indices)


if __name__ == '__main__':
    unittest.main()
