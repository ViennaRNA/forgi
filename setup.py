from distutils.core import setup

setup(name='forgi',
      version='0.20',
      description='RNA Graph Library',
      author='Peter Kerpedjiev',
      author_email='pkerp@tbi.univie.ac.at',
      url='http://www.tbi.univie.ac.at/~pkerp/forgi/',
      packages=['forgi', 'forgi.graph', 'forgi.threedee', 'forgi.threedee.model', 
                'forgi.utilities', 'forgi.threedee.utilities', 'forgi.aux', 
                'forgi.aux.k2n_standalone', 'forgi.threedee.visual', 'forgi.visual'],
      package_data={'forgi.threedee': ['data/*.pdb', 'data/stats/temp.stats']},
      scripts=['examples/visualize_cg.py', 'examples/visualize_pdb.py', 
               'examples/pdb_to_cg.py', 'examples/pdb_to_ss_fasta.py',
               'examples/bpseq_to_bulge_graph.py', 
               'examples/ss_to_bulge_graph.py', 'examples/cg_to_ss_fasta.py',
               'examples/cg_to_bpseq_string.py', 'examples/average_atom_positions.py',
              'examples/dotbracket_to_element_string.py'],

     )
