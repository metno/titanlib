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
        # TODO: Reenable these
        #self.assertEqual(titanlib.range_check_climatology(len1, len1, len1, len1, ut, len1, len1)[0],0)
        #self.assertEqual(titanlib.range_check_climatology(len3, len3, len3, len3, ut, len1, len1)[0],0)
        #self.assertEqual(titanlib.range_check_climatology(len3, len3, len3, len3, ut, len1, len3)[0],0)
        #self.assertEqual(titanlib.range_check_climatology(len3, len3, len3, len3, ut, len1, len1)[0],0)
        pass

    def test_invalid_input(self):
        """Check that the test fails"""
        # TODO: Reenable these
        # Min/max inconsistent with inputs
        #flags  = titanlib.range_check_climatology(len3, len3, len3, len3, ut, len2, len2)
        #self.assertEqual(status,1)

        # Lat / lon inconsistent
        #flags  = titanlib.range_check_climatology(len3, len1, len3, len3, ut, len1, len1)
        #self.assertEqual(status,1)


if __name__ == '__main__':
    unittest.main()
