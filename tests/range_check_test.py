from __future__ import print_function
import unittest
import titanlib


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
ut = 1546300800


class RangeCheckTest(unittest.TestCase):
    def test_valid_input(self):
        """Check that the test doesn't fail"""
        self.assertTrue(titanlib.range_check(len1, len1, len1)[0])
        self.assertTrue(titanlib.range_check(len3, len1, len1)[0])
        self.assertTrue(titanlib.range_check(len3, len1, len3)[0])
        self.assertTrue(titanlib.range_check(len3, len1, len1)[0])

    def test_invalid_input(self):
        """Check that the test fails"""
        # Min/max inconsistent with inputs
        status, flags  = titanlib.range_check(len3, len2, len2)
        self.assertFalse(status)

    def test_simple(self):
        """Check that the test returns the correct flags"""
        # Min/max inconsistent with inputs
        values = [4, 2, 2]
        min = [5, 1, 0]
        max = [6, 3, 1]

        status, flags  = titanlib.range_check(values, min, max)
        self.assertTrue(status)
        self.assertEqual(flags, (1, 0, 1))


if __name__ == '__main__':
    unittest.main()
