# How to set up tests in YAML files

To add tests, create a yaml file (must have extension .yml) in files/ with the following format. The file is a
list of dataset and test combinations. Each item in the list has the keys:

```
description: Description of dataset (optional)
dataset:                 # Dataset of observations
   lats: [60]            # Latitude in degrees
   lons: [10]            # Longitudes in degrees
   elevs: [10]           # Elevations in meters
   values: [25]          # Observational value
tests:                   # A list of checks to run on the above dataset
    - type: range_check  # Which QC test to run (must be a function name in titanlib)
      args:              # The list of arguments to the function. Must be the same names as in the titanlib.h file
         min: [0]
         max: [1]
      expected: [1]      # The expected flags for each value (1 for flagged, 0 for unflagged)
```

# More advanced options
Sometimes its useful to run a series of tests on almost identical input, where one parameter is changed at at
time. This can become verbose with the above approach. To facilitate this use case, add the `base_test` and
`base_args` keywords like this:

```
dataset:
   ...
base_test: range_check
base_args:
   min: [0]
   max: [1]
tests:
    - args:
         max: [25]
         expected: [0]
    - args:
        min: [-1]
        expected: [1]
```

In this case, the base_args are used when a test does not provide one. Also base_test is used as the test for
all tests.
