from setuptools import setup
from Cython.Build import cythonize
import numpy
setup(
    name='orthonormal_test',
    ext_modules=cythonize("orthonormal_test.pyx"),include_dirs=[numpy.get_include()]
)
