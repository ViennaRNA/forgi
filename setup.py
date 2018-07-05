from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
import subprocess
import os
import warnings

def try_cythonize(arg):
  try:
    from Cython.Build import cythonize
    cythonize(arg)
  except Exception as e:
    warnings.warn("Could not use cython. Exception of type {} occurred: {}".format(type(e), e))

try: #If we are in a git-repo, get git-describe version.
    path = os.path.abspath(os.path.dirname(__file__))
    forgi_version = subprocess.check_output(["git", "describe", "--always"], universal_newlines=True).strip()
    try:
        subprocess.check_call(['git', 'diff-index', '--quiet', 'HEAD', '--'], universal_newlines=True)
    except subprocess.CalledProcessError:
        forgi_version+="+uncommited_changes"
    #Use a subclass of build_py from distutils to costumize the build.
    class build_py(_build_py):
        def run(self):
            """
            During building, adds a variable with the complete version (from git describe)
            to forgi/__init__.py.
            """
            outfile = self.get_module_outfile(self.build_lib, ["forgi"], "__init__")
            try:
                os.remove(outfile) #If we have an old version, delete it, so _build_py will copy the original version into the build directory.
            except:
                pass
            # Superclass build
            _build_py.run(self)
            # Apped the version number to init.py
            with open(outfile, "a") as of:
                of.write('\n__complete_version__ = "{}"'.format(forgi_version))
except: #Outside of a git repo, do nothing.
    build_py = _build_py


setup(
      cmdclass={'build_py': build_py},
      name='forgi',
      version='2.0-alpha',
      description='RNA Graph Library',
      author='Peter Kerpedjiev, Bernhard Thiel',
      author_email='pkerp@tbi.univie.ac.at, thiel@tbi.univie.ac.at',
      license='GNU Affero GPL 3.0',
      url='http://www.tbi.univie.ac.at/~pkerp/forgi/',
      ext_modules = try_cythonize("forgi/threedee/utilities/cytvec.pyx"),
      packages=['forgi', 'forgi.graph', 'forgi.threedee',
                'forgi.threedee.model', 'forgi.utilities',
                'forgi.threedee.utilities',
                'forgi.threedee.classification',
                'forgi._k2n_standalone', 'forgi.threedee.visual',
                'forgi.visual', 'forgi.projection'],
      package_data={'forgi.threedee': ['data/*.pdb', 'data/stats/temp.stats', 'data/average_atom_positions.json', 'data/aminor_geometries.csv', 'data/aminor_params.json']},
      scripts=['examples/rnaConvert.py',
               'examples/describe_cg.py',
               'examples/compare_RNA.py',
               'examples/visualize_cg.py',
               'examples/visualize_pdb.py'],
      install_requires=[
		'numpy>=1.10.0',
                'scipy>=0.19.1',
                'networkx>=2.0',
                'future',
                'biopython',
                'pandas>=0.20',
                'appdirs>=1.4',
                'logging_exceptions>=0.1.8',
                'beautifulsoup4>=4.6'
	],
      extras_require={
        'forgi.visual': ["matplotlib>=2.0"],
        'faster vector': ["cython"]
      },
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU Affero General Public License v3',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
         ],

     )
