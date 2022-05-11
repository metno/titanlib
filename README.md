# Titanlib 
[!["Latest release"](https://img.shields.io/github/v/release/metno/titanlib.svg)](https://github.com/metno/titanlib/releases)
[![C/C++ CI](https://github.com/metno/titanlib/workflows/C/C++%20CI/badge.svg)](https://github.com/metno/titanlib/actions)

Titanlib is a library of **automatic quality control** routines for weather observations. It emphases **spatial
checks** and is suitable for use with dense observation networks, such as citizen weather observations. It is
written in C++ and has bindings for python and R. The library consists of a set of functions that perform
tests on data.

Titanlib is currently under active development and the current version is a prototype for testing. Feedback
is welcome, either by using the issue tracker in Github, or by contacting Thomas Nipen (thomasn@met.no).

![Example of titanlib](docs/image.jpg)

## Documentation

For more information on how to use Titanlib, check out the wiki at https://github.com/metno/titanlib/wiki.

## Features

- A wide variety of spatial checks, such as **spatial consistency test**, **buddy check**, **isolation check**.
- Plausability tests such as **range check** and **climatology check**.
- Fast C++ implementation for efficient processing of large observation datasets 

## Quick-start in python

The easiest is to install the latest release of the package using pip. Provided you have installed the dependencies listed above, you can install the most recent release of the python package as follows:
```bash
pip install titanlib
```

## Installation

For full installation instructions, see the following [wiki page](https://github.com/metno/titanlib/wiki/Installation).

## Python example

Here is an example using the buddy check, which has the following function signature:
```python
buddy_check(points, values, radius, num_min, threshold, max_elev_diff, elev_gradient, min_std, num_iterations)
```

The test reveals that the last observation (-111) is likely faulty:

```python
import titanlib
lats = [60,60.1,60.2]
lons = [10,10,10]
elevs = [0,0,0]
obs = [0, 1, -111]
points = titanlib.Points(lats, lons, elevs)
radius = 50000
num_min = 2

flags = titanlib.buddy_check(points,
   obs,
   [radius],
   [num_min],
   2,  # threshold
   200,  # max_elev_diff
   0,  # elev_gradient
   1,  # min_std
   2,  # num_iterations
)

print(flags)

>>> [0 0 1]
```

## R example

Run the following code in R from the build directory, or if you want to run from any other directory, just
put in the proper paths for rtitanlib and titanlib.R

```
dyn.load(paste("extras/SWIG/R/titanlib", .Platform$dynlib.ext, sep=""))
source("extras/SWIG/R/titanlib.R")
cacheMetaData(1)

sct(c(60,60.1,60.2), c(10,10,10), c(0,0,0), c(0,1,-111),50000,2,2,100,0,1,2)
```

See also the _Installation tips and tricks_ on the [wiki](https://github.com/metno/titanlib/wiki/R-interface) .

## Copyright and license

Copyright © 2019-2022 Norwegian Meteorological Institute. Titanlib is licensed under The GNU Lesser General
Public License (LGPL). See LICENSE file.

## Contact

E-mail: Thomas Nipen (thomasn@met.no)
