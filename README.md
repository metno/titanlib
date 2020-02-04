# titanlib

Automatic quality control of in-situ observations with an emphasis on spatial controls.

## Required dependencies

You need to install swig:

```
sudo apt install swig
```

## Installation of library

Install the library and packages using cmake. Create a build directory and perform the
installation there:

```
mkdir build
cd build
cmake ../
```

To compile just the titan library, do this:
```
make titanlib
```

This creates the file libtitanlib.so

## Installation and use of python module

To **install** the python module, you only need the following (i.e. you do not need to also
install the titan library)

```
sudo make install-python
```

The installation above should install the python package in a central place, which should
automatically be accessible in python3.

If you instead just want to **build** the package, do this:

```
make build-python
```

The package is then available in SWIG/python/titanlib.py. Here is an example of how to use the package:

```python
import titanlib
titanlib.calc_gamma(1, 2)
status, flags = titanlib.sct([60,61],[10,11],[0,100],[-4,2],0,0,0,0,0,0,[2,2],[2,2],[2,2])
print(flags)
```

## Installation and use of R module

To build the R module, do this:

```
make build-r
```

Run the following code in R from the build directory:

```
dyn.load(paste("SWIG/R/titanlib", .Platform$dynlib.ext, sep=""))
source("SWIG/R/titanlib.R")
cacheMetaData(1)

calc_gamma(1, 2)
sct(c(60,61), c(10,11), c(0,100), c(-4,2),0,0,0,0,0,0,c(2,2),c(2,2),c(2,2))
```

or if you want to run from any other directory, just put in the proper paths for rtitanlib and titanlib.R

## Run the test suite

Tests are written in src/tests. To run these, first install the python library, then run the following:

```bash
pip3 install nose

cd src
nosetests3
```
