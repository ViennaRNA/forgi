from setuptools import setup
from Cython.Build import cythonize
import numpy
setup(
    name='orthonormal_test',
    ext_modules=cythonize(["orthonormal_test.pyx","brokenml1.pyx", "brokenml2.pyx", "cppvect.pyx", "broken_ml_core.cpp"], language="C++"),
    include_dirs=[numpy.get_include()]
)
