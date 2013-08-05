import corgy.graph.bulge_graph as cgb
import subprocess as sp
import corgy.utilities.mcannotate as cum
import corgy.utilities.debug as cud

import tempfile as tf
import time
import sys
import os.path as op

def from_pdb(pdb_filename):
    '''
    Load a coarse grain model from a PDB file, by extracing
    the bulge graph.
    '''
    pdb_base = op.splitext(op.basename(pdb_filename))[0]
    # first we annotate the 3D structure
    p = sp.Popen(['MC-Annotate', pdb_filename], stdout=sp.PIPE)
    out, err = p.communicate()
    lines = out.strip().split('\n')
    # convert the mcannotate output into bpseq format
    dotplot = cum.get_dotplot(lines)

    f = tf.NamedTemporaryFile()
    f.write(dotplot)
    f.flush()

    # remove pseudoknots
    p = sp.Popen(['aux/k2n_standalone/knotted2nested.py', '-f', 'bpseq', 
                  '-F', 'vienna', f.name], stdout = sp.PIPE)
    out, err = p.communicate()
    out = out.replace(' Nested structure', pdb_base)

    bg = cgb.BulgeGraph()
    bg.from_fasta(out)

    print bg.to_bg_string()

    f.close()
    

class CoarseGrainRNA:
    '''
    A coarse grain model of RNA structure based on the
    bulge graph representation.

    Each stem is represented by four parameters (two endpoints)
    and two twist vetors pointing towards the centers of the base
    pairs at each end of the helix.
    '''
    def __init__(self):
        pass

