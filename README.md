# titanlib

Automatic quality control of in-situ observations with an emphasis on spatial controls.

## Installation

Install the library and packages using cmake. Create a build directory and perform the installation
there:

```
mkdir build
cd build
cmake ../
```

To install just the titan library, do this:
```
sudo make install
```

This will install it in /usr/local/lib/libtitanlib.so

To install the python module, you only need the following (i.e. you do not need to also
install the titan library)

```
sudo make install-python
```

## Use in python

The installation above should install the python package in a central place, which should
automatically be accessible in python3:

```python
import titanlib
titanlib.calc_gamma(1, 2)
```

## Use in R

```
dyn.load(paste("titanlib", .Platform$dynlib.ext, sep=""))
source("titanlib.R")
cacheMetaData(1)

calc_gamma(1, 2)
```
