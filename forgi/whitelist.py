from .graph import bulgegraph as fgb
from .threedee.model import stats as ftms
from .threedee.utilities import graph_pdb as ftug
from .threedee.utilities import pdb as ftup
from .threedee.utilities import vector as ftuv
"""
This is used for vulture (dead code detection tool) 
to mark library-functions as useful, even if they
are not used within the forgi project.
"""

fgb.print_brackets
bg = fgb.BulgeGraph()
bg.interior_loop_iterator()
bg.is_single_stranded()

s = ftms.ContinuousAngleStats()
ftug.spos_to_pos
ftup.output_chain.HSelect.accept_atom

ftuv.null_array
ftuv.tau
ftuv.is_almost_colinear([1], [2])
ftuv.get_standard_basis(3)
ftuv.vec_dot([0, 2, 2], [3, 4, 5])
