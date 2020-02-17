from __future__ import print_function
import unittest
import titanlib


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
ut = 1546300800


class ClimatologyCheckTest(unittest.TestCase):
    def test_valid_input(self):
        """Check that the test doesn't fail"""
        self.assertTrue(titanlib.range_check_climatology(len1, len1, len1, len1, ut, len1, len1)[0])
        self.assertTrue(titanlib.range_check_climatology(len3, len3, len3, len3, ut, len1, len1)[0])
        self.assertTrue(titanlib.range_check_climatology(len3, len3, len3, len3, ut, len1, len3)[0])
        self.assertTrue(titanlib.range_check_climatology(len3, len3, len3, len3, ut, len1, len1)[0])

    def test_invalid_input(self):
        """Check that the test fails"""
        # Min/max inconsistent with inputs
        status, flags  = titanlib.range_check_climatology(len3, len3, len3, len3, ut, len2, len2)
        self.assertFalse(status)

        # Lat / lon inconsistent
        status, flags  = titanlib.range_check_climatology(len3, len1, len3, len3, ut, len1, len1)
        self.assertFalse(status)


if __name__ == '__main__':
    unittest.main()
