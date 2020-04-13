from __future__ import print_function
import unittest
import titanlib
import numpy as np
import csv
import glob


"""Checks that titanlib works on bigger files"""

def load(filename):
    """Reads and parses a test file
        Returns:
            parameters: dictionary of parameter->value corresponding to a test's options
            data: dictionary of input->values 
    """
    data = dict()
    parameters = dict()
    with open(filename, 'r') as file:
        header = None
        for line in file:
            line = line.strip()
            if ":" in line:
                line = ''.join([i for i in line if i != ' '])
                words = line.split(':')
                key = words[0]
                value = words[1]
                if is_numeric(value):
                    value = float(value)
                parameters[key] = value
            elif line[0] == '#':
                continue
            elif header is None:
                header = line.split(' ')
                for h in header:
                    data[h] = list()
            else:
                words = line.split(" ")
                for i in range(len(words)):
                    h = header[i]
                    data[h] += [float(words[i])]
    for key in data:
        if key == 'flag':
            data[key] = np.array(data[key], 'int')
        else:
            data[key] = np.array(data[key], 'float')
    assert('test' in parameters)
    assert('flag' in data)
    return parameters, data


class BulkTest(unittest.TestCase):
    def run_check(self, filename):
        """Check that the test doesn't fail"""
        parameters, data = load(filename)
        if parameters['test'] == 'range':
            flags  = titanlib.range_check(data['value'], [parameters['min']], [parameters['max']])
            self.assertListEqual(list(flags), data['flag'].tolist())
        else:
            raise NotImplementedError

    def test_run_all(self):
        directory = '/'.join(__file__.split('/')[0:-1])
        filenames = glob.glob('%s/files/*.txt' % directory)
        for filename in filenames:
            self.run_check(filename)


def is_numeric(value):
    try:
        ret = float(value)
        return True
    except Exception as e:
        return False


if __name__ == '__main__':
    unittest.main()
