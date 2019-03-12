from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
setup(
    name='cytvec',
    ext_modules=cythonize(Extension("cytvec",
                                    sources=["cytvec.pyx"],
                                    extra_compile_args=['-O3', "-std=c++11"],
                                    language='c++')),
    include_dirs=[numpy.get_include()]
)
