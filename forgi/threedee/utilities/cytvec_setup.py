from setuptools import setup
from Cython.Build import cythonize

setup(
  name = 'cytvec',
  ext_modules = cythonize("cytvec.pyx"),
)
