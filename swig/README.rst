Installation
------------

```
make python
make r
```

Use in python
-------------

```
python

import titanlib
titanlib.calc_gamma(1, 2)
```

Use in R
--------

```
R

dyn.load(paste("titanlib", .Platform$dynlib.ext, sep=""))
source("titanlib.R")
cacheMetaData(1)

calc_gamma(1, 2)

```
