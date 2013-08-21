import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg

import corgy.utilities.debug as cud
import corgy.utilities.mcannotate as cum
import corgy.utilities.pdb as cup

import Bio.PDB as bpdb
import collections as c
import os.path as op
import subprocess as sp
import sys
import tempfile as tf
import time
import warnings

def load_cg_from_pdb(pdb_filename):
    '''
    Load a coarse grain model from a PDB file, by extracing
    the bulge graph.
    '''
    # create a temporary file to hold the biggest chain
    with tf.NamedTemporaryFile() as f:
        cup.get_biggest_chain(pdb_filename, f) 
        f.flush()

        # f1 will store the rosetta-fied pdb file
        with tf.NamedTemporaryFile() as f1:
            # Use rosetta's rna preparation script to rename some of the atoms
            p = sp.Popen(['aux/make_rna_rosetta_ready.py', f.name], 
                         stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = p.communicate()

            # Change rosetta's residue naming a little bit
            line = out
            line = line.replace('rA', 'A')
            line = line.replace('rC', 'C')
            line = line.replace('rG', 'G')
            line = line.replace('rU', 'U')
            f1.write(line)

            f1.flush()

            pdb_base = op.splitext(op.basename(pdb_filename))[0]
            # first we annotate the 3D structure
            p = sp.Popen(['MC-Annotate', f1.name], stdout=sp.PIPE)
            out, err = p.communicate()
            lines = out.strip().split('\n')
            # convert the mcannotate output into bpseq format
            dotplot = cum.get_dotplot(lines)

            # f2 will store the dotbracket notation
            with tf.NamedTemporaryFile() as f2:
                f2.write(dotplot)
                f2.flush()

                # remove pseudoknots
                p = sp.Popen(['aux/k2n_standalone/knotted2nested.py', '-f', 'bpseq', 
                              '-F', 'vienna', f2.name], stdout = sp.PIPE)

                out, err = p.communicate()
                out = out.replace(' Nested structure', pdb_base)

                bg = cgb.BulgeGraph()
                bg.from_fasta(out, dissolve_length_one_stems=1)

            # Add the 3D information about the starts and ends of the stems
            # and loops
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = bpdb.PDBParser().get_structure('temp', f1.name)
                chain = list(s.get_chains())[0]

            cg = CoarseGrainRNA(bg)
            cgg.add_stem_information_from_pdb_chain(cg, chain)
            cgg.add_bulge_information_from_pdb_chain(cg, chain)
            cgg.add_loop_information_from_pdb_chain(cg, chain)

    return cg

    
def from_pdb(pdb_filename):
    cg = load_cg_from_pdb(pdb_filename)

    return cg

class CoarseGrainRNA():
    '''
    A coarse grain model of RNA structure based on the
    bulge graph representation.

    Each stem is represented by four parameters (two endpoints)
    and two twist vetors pointing towards the centers of the base
    pairs at each end of the helix.
    '''
    def __init__(self, bg):
        '''
        Initialize the new structure using an existing BulgeGraph.
        Without a BulgeGraph to describe the topology of the RNA
        structure, this class is useless.
        '''
        self.bg = bg
        self.coords = dict()
        self.twists = dict()
        self.sampled = dict()

        self.vposs = c.defaultdict( dict )
        self.vvecs = c.defaultdict( dict )
        self.v3dposs = c.defaultdict( dict )
        self.vbases = c.defaultdict( dict )
        self.vinvs = c.defaultdict( dict )

        self.longrange = c.defaultdict( set )

        pass

    def get_coord_str(self):
        '''
        Place the start and end coordinates of each stem into a string.

        The format is:

            coord s1 x1 y1 z1 x2 y2 z2

        Where s1 is the name of the stem, (x1,y1,z1) and (x2,y2,z2) are
        the two endpoints of the stem.

        @return: A string containing the coordinates for all of the stems.
        '''
            
        out_str = ''
        for key in self.coords.keys():
            [p, n] = self.coords[key]
            out_str += "coord %s %s %s" % (key, " ".join([str(pt) for pt in p]),
                                           " ".join([str(pt) for pt in n]))
            out_str += '\n'
        return out_str

    def get_twist_str(self):
        '''
        Place the twist vectors into a string. 

        The format is:

            twist s1 x1 y1 z1 x2 y2 z2

        Where s1 is the name of the stem and (x1,y1,z1) and (x2,y2,z2) are
        the unit vectors from the start of each end of the stem to the middle
        of the base pairs.
        '''
        out_str = ''
        for key in self.twists.keys():
            [p, n] = self.twists[key]
            out_str += "twist %s %s %s" % (key, " ".join([str(pt) for pt in p]),
                                           " ".join([str(pt) for pt in n]))
            out_str += '\n'
        return out_str

    def get_long_range_str(self):
        out_str = ''
        for key1 in self.longrange.keys():
            for key2 in self.longrange[key1]:
                out_str += "longrange %s %s\n" % (key1, key2)

        return out_str

    def to_cg_string(self):
        '''
        Output this structure in string form.
        '''
        curr_str = self.bg.to_bg_string()
        curr_str += self.get_coord_str()
        curr_str += self.get_twist_str()
        curr_str += self.get_long_range_str()

        return curr_str

    def from_cg_string(self, cg_string):
        '''
        Populate this structure from the string
        representation of a graph.
        '''
        lines = cg_string.split('\n')
        for line in lines:
            line = line.strip()
            parts = line.split()
            if len(parts) == 0:
                continue
            if parts[0] == 'coord':
                name = parts[1]
                self.coords[name] = [map(float, parts[2:5]), map(float, parts[5:8])]
            if parts[0] == 'twist':
                name = parts[1]
                self.twists[name] = [map(float, parts[2:5]), map(float, parts[5:8])]
