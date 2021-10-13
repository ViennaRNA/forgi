from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

setup(
  name = 'cppvect',
  ext_modules=
    cythonize(Extension('cpp_vect',
              sources=['cppvect.pyx'],
              extra_compile_args=['-O3', "-std=c++11"],
              language='c++')),
    include_dirs=[numpy.get_include()]

)

