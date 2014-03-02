from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(name="cytvec", 
                sources = ["cytvec.pyx"])]

setup(
  name = 'Cython Vector Functions',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
