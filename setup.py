from distutils.core import setup

from distutils.command.build_py import build_py as _build_py
import subprocess
import os

try: #If we are in a git-repo, get git-describe version.
    path = os.path.abspath(os.path.dirname(__file__))
    forgi_version = subprocess.check_output([os.path.join(path, "git"), "describe", "--always"]).strip()
    try:
        subprocess.check_call([os.path.join(path, 'git'), 'diff-index', '--quiet', 'HEAD', '--'])
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
except OSError: #Outside of a git repo, do nothing.
    build_py = _build_py


setup(
      cmdclass={'build_py': build_py},
      name='forgi',
      version='0.40',
      description='RNA Graph Library',
      author='Peter Kerpedjiev, Bernhard Thiel',
      author_email='pkerp@tbi.univie.ac.at, thiel@tbi.univie.ac.at',
      license='GNU Affero GPL 3.0',
      url='http://www.tbi.univie.ac.at/~pkerp/forgi/',
      packages=['forgi', 'forgi.graph', 'forgi.threedee', 
                'forgi.threedee.model', 'forgi.utilities', 
                'forgi.threedee.utilities', 'forgi.aux', 
                'forgi.aux.k2n_standalone', 'forgi.threedee.visual', 
                'forgi.visual', 'forgi.projection'],
      package_data={'forgi.threedee': ['data/*.pdb', 'data/stats/temp.stats']},
      scripts=['examples/visualize_cg.py', 
               'examples/visualize_pdb.py', 
               'examples/pdb_rmsd.py',
               'examples/cg_rmsd.py',
               'examples/cg_rog.py',
               'examples/cg_to_bpseq_string.py',
               'examples/cg_to_fornac.py',
               'examples/cg_to_ss_fasta.py',
               'examples/pdb_to_cg.py', 
               'examples/pdb_to_ss_fasta.py',
               'examples/bpseq_to_bulge_graph.py', 
               'examples/ss_to_bulge_graph.py', 
               'examples/cg_to_ss_fasta.py',
               'examples/cg_to_bpseq_string.py', 
               'examples/average_atom_positions.py',
               'examples/dotbracket_to_element_string.py', 
               'examples/dotbracket_to_bulge_graph.py',
               'examples/compare_projections.py',
               'examples/graph_to_neato.py',
               'examples/investigate_projection.py',
               'examples/multipleProjectionRMSD.py',
               'examples/plot_projection.py'],

     )
