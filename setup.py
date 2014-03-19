from distutils.core import setup

setup(name='forgi',
      version='0.1',
      description='RNA Graph Library',
      author='Peter Kerpedjiev',
      author_email='pkerp@tbi.univie.ac.at',
      url='http://www.tbi.univie.ac.at/~pkerp/forgi/',
      packages=['forgi', 'forgi.graph', 'forgi.utilities'],
      package_data={'forgi': ['data/*.pdb']}
     )
