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
        N = 5
        np.random.seed(0)
        lats = np.random.randn(N) * 10000;
        lons = np.random.randn(N) * 10000;
        elevs = np.random.rand(N) * 100;
        points = titanlib.Points(lats, lons, elevs)
        values = np.zeros(N);
        dataset = titanlib.Dataset(points, values);
        return dataset

    def test_indices(self):
        """Check that the test doesn't fail"""
        dataset = self.get_dataset()
        dataset.range_check([-100], [-100], [0, 1, 2])
        self.assertListEqual([i for i in dataset.flags], [1, 1, 1, 0, 0])
        dataset.range_check([-100], [-100], [1, 2, 4])
        self.assertListEqual([i for i in dataset.flags], [1, 1, 1, 0, 1])

    def test_merge(self):
        dataset = self.get_dataset()
        dataset.isolation_check(3, 1000000000)
        self.assertListEqual([i for i in dataset.flags], [0, 0, 0, 0, 0])

        dataset = self.get_dataset()
        N = dataset.points.size()
        dataset.isolation_check(3 * np.ones(N), 1000000000 * np.ones(N))
        self.assertListEqual([i for i in dataset.flags], [0, 0, 0, 0, 0])

    def test_empty_indices(self):
        dataset = self.get_dataset()
        dataset.range_check([-100], [-100], [])
        self.assertListEqual([i for i in dataset.flags], [0, 0, 0, 0, 0])

    def test_subset(self):
        dataset = self.get_dataset()
        N = len(dataset.values)
        dataset.range_check([-100] * N, [-100] * N, [])

    def test_buddy_check(self):
        N = 10
        lats = [60]*N
        lons = np.linspace(60, 60.001, N)
        elevs = [0]*10
        points = titanlib.Points(lats, lons, elevs)
        values = [0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1]
        N = len(lats)
        radius = [10000]
        num_min = [1]
        threshold = 1
        elev_gradient = -0.0065
        max_elev_diff = 200
        min_std = 0.01
        num_iterations = 2
        dataset = titanlib.Dataset(points, values)
        dataset.buddy_check(radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(dataset.flags, [0]*8 + [1]*2)

    def test_get_values(self):
        N = 10
        lats = [60]*N
        lons = np.linspace(60, 60.001, N)
        elevs = [0]*10
        points = titanlib.Points(lats, lons, elevs)
        values = [0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1]
        N = len(lats)
        radius = [10000]
        num_min = [1]
        threshold = 1
        elev_gradient = -0.0065
        max_elev_diff = 200
        min_std = 0.01
        num_iterations = 2
        dataset = titanlib.Dataset(points, values)
        dataset.values[:] = [1]
        print(dataset.get_values())
        dataset.points = points



if __name__ == '__main__':
    unittest.main()
