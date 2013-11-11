import forgi.graph.bulge_graph as cgb
import forgi.threedee.utilities.graph_pdb as cgg

import forgi.aux.k2n_standalone.knotted2nested as cak
import forgi.utilities.debug as cud
import forgi.threedee.utilities.mcannotate as cum
import forgi.threedee.utilities.pdb as cup

import Bio.PDB as bpdb
import collections as c
import contextlib
import numpy as np
import os
import os.path as op
import shutil
import subprocess as sp
import sys
import tempfile as tf
import time
import warnings

def remove_hetatm(lines):
    '''
    Go through the lines of a pdb file and remove any which refer to a
    HETATM.

    @param lines: A an array of lines of text from a pdb file.
    '''
    new_lines = []

    for line in lines:
        if line.find('HETATM') == 0:
            if line.find('5MU') > 0:
                line = line.replace('5MU', '  U')
            elif line.find('PSU') > 0:
                line = line.replace('PSU', '  U')
            elif line.find('5MC') > 0:
                line = line.replace('5MC', '  C')
            elif line.find('1MG') > 0:
                line = line.replace('1MG', '  G')
            elif line.find('H2U') > 0:
                line = line.replace('H2U', '  G')
            else:
                continue
        
        line = line.replace('HETATM', 'ATOM  ')
        new_lines += [line]

    return new_lines

@contextlib.contextmanager
def make_temp_directory():
    '''
    Yanked from:

    http://stackoverflow.com/questions/13379742/right-way-to-clean-up-a-temporary-folder-in-python-class
    '''
    temp_dir = tf.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

def load_cg_from_pdb_in_dir(pdb_filename, output_dir, secondary_structure=''):
    '''
    Create the coarse grain model from a pdb file and store all
    of the intermediate files in the given directory.

    @param pdb_filename: The name of the pdb file to be coarseified
    @param output_dir: The name of the output directory
    @param secondary_structure: Specify a particular secondary structure
                                for this coarsification.
    '''
    chain = cup.load_structure(pdb_filename)
    pdb_base = op.splitext(op.basename(pdb_filename))[0]
    output_dir = op.join(output_dir, pdb_base)

    if not op.exists(output_dir):
        os.makedirs(output_dir)

    with open(op.join(output_dir, 'temp.pdb'), 'w') as f:
        # output the biggest RNA chain
        cup.output_chain(chain, f.name)
        f.flush()

        pdb_base = op.splitext(op.basename(pdb_filename))[0]
        # first we annotate the 3D structure
        p = sp.Popen(['MC-Annotate', f.name], stdout=sp.PIPE)
        out, err = p.communicate()
        lines = out.strip().split('\n')
        # convert the mcannotate output into bpseq format
        dotplot = cum.get_dotplot(lines)

        # f2 will store the dotbracket notation
        with open(op.join(output_dir, 'temp.dotplot'), 'w') as f2:
            f2.write(dotplot)
            f2.flush()

            # remove pseudoknots
            '''
            p = sp.Popen(['aux/k2n_standalone/knotted2nested.py', '-f', 'bpseq', 
                          '-F', 'vienna', f2.name], stdout = sp.PIPE)

            out, err = p.communicate()
            '''
            out = cak.k2n_main(f2.name, input_format='bpseq',
                               output_format = 'vienna',
                               method = cak.DEFAULT_METHOD,
                               opt_method = cak.DEFAULT_OPT_METHOD,
                               verbose = cak.DEFAULT_VERBOSE,
                               removed= cak.DEFAULT_REMOVED)

            out = out.replace(' Nested structure', pdb_base)

            if secondary_structure != '':
                lines = out.split('\n')

                if len(secondary_structure) != len(lines[1].strip()):
                    print >>sys.stderr, "The provided secondary structure \
                            does not match the length of the 3D structure"
                    print >>sys.stderr, "Sequence:", lines[1]
                    print >>sys.stderr, "ss_struc:", secondary_structure
                    sys.exit(1)

                lines[-1] = secondary_structure
                out = "\n".join(lines)

            # Add the 3D information about the starts and ends of the stems
            # and loops
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = bpdb.PDBParser().get_structure('temp', f.name)
                chain = list(s.get_chains())[0]

            cg = CoarseGrainRNA()
            cg.from_fasta(out, dissolve_length_one_stems=1)
            cgg.add_stem_information_from_pdb_chain(cg, chain)
            cgg.add_bulge_information_from_pdb_chain(cg, chain)
            cgg.add_loop_information_from_pdb_chain(cg, chain)

            with open(op.join(output_dir, 'temp.cg'), 'w') as f3:
                f3.write(cg.to_cg_string())
                f3.flush()

            return cg
    print >>sys.stderr, "Uh oh... couldn't generate the coarse-grain structure."
    print >>sys.stderr, "Prepare for an incoming exception."

def load_cg_from_pdb(pdb_filename, secondary_structure='', 
                     intermediate_file_dir=''):
    '''
    Load a coarse grain model from a PDB file, by extracing
    the bulge graph.

    @param pdb_filename: The filename of the 3D model
    @param secondary_structure: A dot-bracket string encoding the secondary
                                structure of this molecule
    '''

    if intermediate_file_dir != '':
        output_dir = intermediate_file_dir

        cg = load_cg_from_pdb_in_dir(pdb_filename, output_dir, secondary_structure)
    else:
        with make_temp_directory() as output_dir:
            print "output_dir"
            print op.exists(output_dir)
            cg = load_cg_from_pdb_in_dir(pdb_filename, output_dir, secondary_structure)

    return cg

def from_file(cg_filename):
    '''
    Read a coarse-grain structure file.

    @param cg_filename: The filename.
    @return: A CoarseGrainRNA from a file.
    '''
    with open(cg_filename, 'r') as f:
        lines = "".join(f.readlines())

        cg = CoarseGrainRNA()
        cg.from_cg_string(lines)

        return cg
    
def from_pdb(pdb_filename, secondary_structure='', intermediate_file_dir=''):
    cg = load_cg_from_pdb(pdb_filename, secondary_structure, intermediate_file_dir)

    return cg

class CoarseGrainRNA(cgb.BulgeGraph):
    '''
    A coarse grain model of RNA structure based on the
    bulge graph representation.

    Each stem is represented by four parameters (two endpoints)
    and two twist vetors pointing towards the centers of the base
    pairs at each end of the helix.
    '''
    def __init__(self):
        '''
        Initialize the new structure.
        '''
        super(CoarseGrainRNA, self).__init__()

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
        curr_str = self.to_bg_string()
        curr_str += self.get_coord_str()
        curr_str += self.get_twist_str()
        curr_str += self.get_long_range_str()

        return curr_str


    def from_cg_string(self, cg_string):
        '''
        Populate this structure from the string
        representation of a graph.
        '''
        self.from_bg_string(cg_string)

        lines = cg_string.split('\n')
        for line in lines:
            line = line.strip()
            parts = line.split()
            if len(parts) == 0:
                continue
            if parts[0] == 'coord':
                name = parts[1]
                self.coords[name] = np.array([map(float, parts[2:5]), map(float, parts[5:8])])
            if parts[0] == 'twist':
                name = parts[1]
                self.twists[name] = np.array([map(float, parts[2:5]), map(float, parts[5:8])])
