
###  Run this in Terminal to compile code:
###  python setup.py build_ext --inplace

import numpy
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("gen_chi2_cube.pyx"),
    include_dirs=[numpy.get_include()]
)
