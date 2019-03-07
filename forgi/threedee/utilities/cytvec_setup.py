from setuptools import setup
from Cython.Build import cythonize
import numpy
setup(
    name='cytvec',
    ext_modules=cythonize("cytvec.pyx"),include_dirs=[numpy.get_include()]
)
