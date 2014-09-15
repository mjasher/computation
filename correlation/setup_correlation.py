#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# inplace extension module
mod_correlation = Extension("_correlation",
                   ["correlation_wrap.cxx","correlation.cpp"],
                   include_dirs = [numpy_include],
                   extra_compile_args = ['-fopenmp','-fpic'],
                   extra_link_args = ['-lgomp']
                   )

# NumyTypemapTests setup
setup(  name        = "correlation",
        description = "correlation matrix builder",
        author      = "mjasher",
        version     = "1.0",
        ext_modules = [mod_correlation]
        )
