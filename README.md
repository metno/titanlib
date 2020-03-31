# Titanlib 

![C/C++ CI](https://github.com/metno/titanlib/workflows/C/C++%20CI/badge.svg)

Titanlib is a library of automatic quality control routines for weather observations. It emphases spatial
checks and is suitable for use with dense observation networks, such as citizen weather observations. It is
written in C++ and has bindings for python and R. The library consists of a set of functions that perform
tests on data.

Titanlib is currently under development.

![Example of titanlib](extras/image.jpg)

## Features

- A wide variety of spatial checks, such as **spatial consistency test**, **buddy check**, **isolation check**.
- Plausability tests such as **range check** and **climatology check**.
- A graphical interface for tuning the parameters of the checks
- Fast C++ implementation for efficient processing of large observation datasets 

## Documentation

For more information on how to use Titanlib, check out the wiki at https://github.com/metno/titanlib/wiki.

## Installation

Here are installation instructions for Ubuntu Bionic. First, install the pre-requisites:

```
sudo apt install swig doxygen r-base-core libboost-dev libproj-dev cmake libgsl-dev python3-setuptools python3-nose python3-numpy python3-scipy
```

Next, configure the installation cmake. Create a build directory and perform the
installation there:

```
mkdir build
cd build
cmake ../
```

To compile the library and bindings, run

```bash
make all
sudo make install
```

To install the bindings, run
```bash
make install-python-user
make build-r
```

This will install the library in `/usr/local/lib/libtitanlib.so`, the python bindings in
`~/local/lib/python3.6/site-packages/titanlib.py`. To install the python bindings
system-wide, use `sudo make install-python`.  Currently, the R package is not installed
centrally, but instead is placed in `SWIG/R/titanlib.R`.

## Python example

```python
import titanlib
status, flags = titanlib.sct([60,61],[10,11],[0,100],[-4,2],100,300,100,100,100,100,[2,2],[2,2],[2,2])
print(flags)
```

## R example

Run the following code in R from the build directory, or if you want to run from any other directory, just
put in the proper paths for rtitanlib and titanlib.R

```
dyn.load(paste("SWIG/R/titanlib", .Platform$dynlib.ext, sep=""))
source("SWIG/R/titanlib.R")
cacheMetaData(1)

sct(c(60,61), c(10,11), c(0,100), c(-4,2),100,300,100,100,100,100,c(2,2),c(2,2),c(2,2))
```

## Copyright and license

Copyright © 2019-2020 Norwegian Meteorological Institute. Titanlib is licensed under The GNU Lesser General
Public License (LGPL). See LICENSE file.
