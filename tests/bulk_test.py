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
        return titanlib.Points(*args)

    def get_dataset(self, dataset):
        points = self.get_points(dataset)
        dataset = titanlib.Dataset(points, dataset["values"])
        return dataset

    def run_checks(self, filename):
        """Check that the test doesn't fail"""
        alldata = self.load(filename)
        for i, data in enumerate(alldata):
            for j, test in enumerate(data["tests"]):
                description = data["description"] if "description" in data else ""
                with self.subTest(filename=filename, dataset=i, test=test['type'], test_index=j, description=description):
                    dataset = self.get_dataset(data["dataset"])
                    # Get args
                    needs_points = True
                    needs_values = True
                    test_type = test['type']
                    func = getattr(titanlib, test_type)
                    func_dataset = getattr(dataset, test_type)
                    optional_args = []
                    if test['type'] == 'range_check':
                        required_args = ['min', 'max']
                        needs_points = False
                    elif test['type'] == "isolation_check":
                        required_args = ['num_min', 'radius']
                        optional_args = ['vertical_radius']
                        needs_values = False
                    elif test['type'] == 'metadata_check':
                        optional_args = ['check_lat', 'check_lon', 'check_elev', 'check_laf']
                        needs_values = False
                    elif test['type'] == 'duplicate_check':
                        required_args = ['radius']
                        optional_args = ['vertical_range']
                        needs_values = False
                    else:
                        raise NotImplementedError
                    dataset_args = list()
                    for arg in required_args:
                        if arg not in test['args']:
                            raise Exception("Test '%s' requires argument '%s'" % (test_type, arg))
                        dataset_args += [test['args'][arg]]

                    for arg in optional_args:
                        if arg in test['args']:
                            dataset_args += [test['args'][arg]]

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
