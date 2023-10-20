from __future__ import print_function
import unittest
import titanlib
import numpy as np
import util
import time


"""Convenient vectors used as inputs"""
values = [0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1]
N = len(values)

def get_points(elevs):
    N = len(elevs)
    y = [0]*N
    x = np.linspace(0, 100, N)
    lafs = [0]*N
    return titanlib.Points(y, x, elevs, lafs, titanlib.Cartesian)

points  = get_points([0]*N)
radius = [10000]
num_min = [1]
threshold = 1
elev_gradient = -0.0065
max_elev_diff = 200
min_std = 0.01
num_iterations = 2

class BuddyCheckTest2iter(unittest.TestCase):
    def test_bad_data_iterative_effect(self):
        ''' Check that after 1 iterations, the last bad data is removed and
        that after 2 iterations the last bad data is not anymore taken into account
        for the variance which makes the variance being lower at the 2nd iteration,
        and thus make the the (N-2)th element being also flagged as bad data'''
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, 1)
        np.testing.assert_array_equal(flags, [0]*(N-1) + [1])
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, 2)
        np.testing.assert_array_equal(flags, [0]*(N-2) + [1]*2)
    
    def test_empty_vector_inputs(self):
        '''Check that empty vectors as inputs are correctly processed'''
        flags = titanlib.buddy_check(titanlib.Points([], [], []), [], [], [],
		                        threshold, max_elev_diff, elev_gradient,
                                        min_std, num_iterations)
        np.testing.assert_array_equal(flags, [])

    def test_empty_vector_mixed_size1_input(self):
        '''Check that empty points and values used together
           with 1-size radius and num_min vectors are correctly processed'''
        flags = titanlib.buddy_check(titanlib.Points([], [], []), [],
                radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [])

    def test_inconsistent_size_input(self):
        '''Check that error is raised if points/values points/num_min lengths are not the same'''
        with self.assertRaises(ValueError):
            titanlib.buddy_check(titanlib.Points([], [], []), values, 
                                 radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        with self.assertRaises(ValueError):
            titanlib.buddy_check(points, values, radius, [], threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
            
    def test_use_all_elements_radius(self):
        '''When using N dimension radius, check that all the elements for the radius are 
           used and that the small radius makes more stations being accepeted'''
        flags = titanlib.buddy_check(points, values, [10000]*(N-1) + [1], num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*N)

    def test_equivalence_size(self):
        '''Check that we get the same result when passing 1 dimension or N dimensions radius/num_min vectors'''
        flags1 = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        flags2 = titanlib.buddy_check(points, values, [radius[0]]*N, [num_min[0]]*N, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags1, flags2)

    def test_min_std(self):
        """Check that when min_std is high enough, nothing is flagged"""
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, 0.3, num_iterations)
        np.testing.assert_array_equal(flags, [0]*(N-1)+[1])
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, 1, num_iterations)
        np.testing.assert_array_equal(flags, [0]*N)

    def test_min_num(self):
        """Check that when min_num is high enough, nothing is flagged"""
        flags = titanlib.buddy_check(points, values, radius, [20], threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*N)

    def test_num_iterations(self):
        """Check that num_iterations affects the flags"""
        flags = titanlib.buddy_check(points, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, 1)
        np.testing.assert_array_equal(flags, [0]*(N-1) + [1])

    def test_elev_gradient(self):
        '''Check flags takes into account the altitude gradient correction'''
        elevs0 = [0]*(N-1) + [-153.8]
        points0 = get_points(elevs0)
        flags = titanlib.buddy_check(points0, values, radius, num_min, threshold,
                max_elev_diff, 0, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*(N-2) + [1, 1])

        flags = titanlib.buddy_check(points0, values, radius, num_min, threshold,
                max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*(N-2) + [1, 0])

    def test_max_elev_diff(self):
        """Check that test is not run on a point which has no other points within the elevation range,
        done by taking a small value for the max_elev_diff """
        elevs0 = [0]*(N-1) + [100]
        points0 = get_points(elevs0)
        small_max_elev_diff = 1
        flags = titanlib.buddy_check(points0, values, radius, num_min, threshold,
                small_max_elev_diff, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [0]*(N-2) + [1, 0])

    def test_missing(self):
        '''Check that NaNs are still flagged, eventhough they are not checked'''
        values0 = np.copy(values)
        values0[0] = np.nan
        flags = titanlib.buddy_check(points, values0, radius, num_min, 0.0001,
                1, elev_gradient, min_std, num_iterations)
        np.testing.assert_array_equal(flags, [1]*N)

if __name__ == '__main__':
    unittest.main()
