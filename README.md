# Titanlib 

Titanlib is a library of automatic quality control routines for in-situ observations with an emphasis on spatial checks. It is written in C++ and has bindings for python and R.

![Example of titanlib](extras/image.jpg)

![C/C++ CI](https://github.com/metno/titanlib/workflows/C/C++%20CI/badge.svg)


## Installation

Here are installation instructions for Ubuntu Bionic. First, install the pre-requisites:

```
sudo apt install swig doxygen R-base-core libboost-dev libproj-dev cmake libgsl-dev
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
sudo make install-python
make build-r
```

This will install the library in `/usr/local/lib/libtitanlib.so`, the python bindings in
`/usr/lib/python3/dist-packages/titanlib.py`. Currently, the R package is not installed centrally, but
instead is placed in `SWIG/R/titanlib.R`.

## Python example

```python
import titanlib
titanlib.calc_gamma(1, 2)
status, flags = titanlib.sct([60,61],[10,11],[0,100],[-4,2],0,0,0,0,0,0,[2,2],[2,2],[2,2])
print(flags)
```

## R example

Run the following code in R from the build directory, or if you want to run from any other directory, just
put in the proper paths for rtitanlib and titanlib.R

```
dyn.load(paste("SWIG/R/titanlib", .Platform$dynlib.ext, sep=""))
source("SWIG/R/titanlib.R")
cacheMetaData(1)

calc_gamma(1, 2)
sct(c(60,61), c(10,11), c(0,100), c(-4,2),0,0,0,0,0,0,c(2,2),c(2,2),c(2,2))
```

## Development instructions

The C++ library and python bindings can be installed separately like this:

```
sudo make install
sudo make install-python
```

To build the library and packages (and not install them), run:

```
make titanlib
make build-python
make build-r
```

*Unit tests* are written in the tests/ directory. The tests rely on the python interface to titanlib being installed (see above). Use nosetests3 to run all tests. First install the python library  and the nose package (pip3 install nose) then run the following:

```bash
nosetests3
```

To *build the documentation*, run the following from the build directory:

```bash
make doc_doxygen
```

The html documentation will appear in build/doxygen/html.

To build a *debug version* of the library, use the following options when running cmake:

```bash
cmake .. -DCMAKE_BUILD_TYPE=DEBUG
```
