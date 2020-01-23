Installation
------------

Install the library and packages using cmake. Create a build directory and perform the installation
there:

.. code-block:: bash

  mkdir build
  cd build
  cmake ../
  make

  make install-python


Use in python
-------------

The installation above should install the python package in a central place, which should
automatically be accessible in python3:

.. code-block:: python

  import titanlib
  titanlib.calc_gamma(1, 2)


Use in R
--------

.. code-block:: bash

  dyn.load(paste("titanlib", .Platform$dynlib.ext, sep=""))
  source("titanlib.R")
  cacheMetaData(1)

  calc_gamma(1, 2)
