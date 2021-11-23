from __future__ import print_function
import unittest
import titanlib
import numpy as np
import csv
import glob
import yaml
import os


"""Checks that titanlib works on bigger files"""

class BulkTest(unittest.TestCase):
    def load(self, filename):
        """Reads and parses a test file
            Returns:
                parameters: dictionary of parameter->value corresponding to a test's options
                data: dictionary of input->values 
        """
        data = dict()
        parameters = dict()
        with open(filename, 'r') as file:
            data = yaml.load(file, Loader=yaml.SafeLoader)
            return data

    def get_points(self, dataset):
        args = [dataset["lats"], dataset["lons"]]
        if "elevs" in dataset:
            args += [dataset["elevs"]]
        if "lafs" in dataset:
            if "elevs" not in dataset:
                args += [np.nan * np.zeros(len(dataset["lafs"]))]
            args += [dataset["lafs"]]
        assert(len(args) <= 4)

        return titanlib.Points(*args)

    def get_dataset(self, dataset):
        points = self.get_points(dataset)
        if "values" in dataset:
            values = dataset["values"]
        else:
            values = np.nan * np.zeros(points.size())
        dataset = titanlib.Dataset(points, values)
        return dataset

    def run_checks(self, filename):
        """Check that the test doesn't fail"""
        alldata = self.load(filename)
        for i, data in enumerate(alldata):
            base_args = []
            base_test = None
            if "base_args" in data:
                base_args = data["base_args"]
            if "base_test" in data:
                base_test = data["base_test"]

            for j, test in enumerate(data["tests"]):
                description = data["description"] if "description" in data else ""
                test_type = test['type'] if base_test is None else base_test
                with self.subTest(filename=filename, dataset=i, test=test_type, test_index=j, description=description):
                    dataset = self.get_dataset(data["dataset"])
                    # Get args
                    needs_points = True
                    needs_values = True
                    func = getattr(titanlib, test_type)
                    func_dataset = getattr(dataset, test_type)
                    optional_args = []
                    required_args = []
                    if test_type == 'range_check':
                        required_args = ['min', 'max']
                        needs_points = False
                    elif test_type == 'range_check_climatology':
                        required_args = ['unixtime', 'pos', 'neg']
                    elif test_type == "isolation_check":
                        required_args = ['num_min', 'radius']
                        optional_args = ['vertical_radius']
                        needs_values = False
                    elif test_type == 'metadata_check':
                        optional_args = ['check_lat', 'check_lon', 'check_elev', 'check_laf']
                        needs_values = False
                    elif test_type == 'duplicate_check':
                        required_args = ['radius']
                        optional_args = ['vertical_range']
                        needs_values = False
                    elif test_type == 'buddy_check':
                        required_args = ['radius', 'num_min', 'threshold', 'max_elev_diff', 'elev_gradient', 'min_std', 'num_iterations']
                        optional_args = ['obs_to_check']
                    elif test_type == 'buddy_event_check':
                        required_args = ['radius', 'num_min', 'event_threshold', 'threshold', 'max_elev_diff', 'elev_gradient', 'num_iterations']
                        optional_args = ['obs_to_check']
                    else:
                        raise NotImplementedError
                    dataset_args = list()

                    for arg in required_args:
                        if 'args' in test and arg in test['args']:
                            dataset_args += [test['args'][arg]]
                        elif arg in base_args:
                            dataset_args += [base_args[arg]]
                        else:
                            raise Exception("Test '%s' requires argument '%s'" % (test_type, arg))

                    for arg in optional_args:
                        if 'args' in test and arg in test['args']:
                            dataset_args += [test['args'][arg]]
                        elif arg in base_args:
                            dataset_args += [base_args[arg]]

                    # Check expected using the dataset interface
                    func_dataset(*dataset_args)
                    flags = [int(i) for i in dataset.flags]
                    np.testing.assert_array_almost_equal(flags, test['expected'])

                    # Check expected using the function interface
                    args = list()
                    if needs_points:
                        args += [dataset.points]
                    if needs_values:
                        args += [[float(i) for i in dataset.values]]
                    args += dataset_args
                    func(*args)
                    np.testing.assert_array_almost_equal(flags, test['expected'])

    def test_run_all(self):
        directory = os.path.join(os.path.dirname(__file__))
        filenames = glob.glob('%sfiles/*.yml' % directory)
        for filename in filenames:
            with self.subTest(filename=filename):
                self.run_checks(filename)


def is_numeric(value):
    try:
        ret = float(value)
        return True
    except Exception as e:
        return False


if __name__ == '__main__':
    unittest.main()
