from __future__ import print_function
import unittest
import titanlib
import numpy as np
import time


class Test(unittest.TestCase):
    def test_compute_vertical_profile(self):
        lats = [0, 0]
        lons = [0, 0]
        elevs = [0, 10]
        oelevs = [20]
        values = [0, 10]
        q = titanlib.compute_vertical_profile(elevs, oelevs, values, 30, 11)
        # Use a mean temperature
        self.assertEqual(q, 5)

        q = titanlib.compute_vertical_profile(elevs, oelevs, values, 30, 9)
        # Use a linear temperature gradient
        self.assertTrue(np.abs(q -  20 < 1))


if __name__ == '__main__':
    unittest.main()
