from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

from ...graph import bulge_graph as fgb
from ..utilities import graph_pdb as ftug
from ..model import stats as ftms

from ...aux.k2n_standalone import knotted2nested as cak
from ..utilities import mcannotate as ftum
from ..utilities import pdb as ftup
from ..utilities import rmsd as ftur
from ..utilities import vector as ftuv
from ...utilities import debug as fud
from ...utilities import stuff as fus
from ...utilities.observedDict import observedDict

import Bio.PDB as bpdb
import collections as c
import contextlib
import numpy as np
import scipy.spatial
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

    :param lines: A an array of lines of text from a pdb file.
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


def add_longrange_interactions(cg, lines):
    '''
    Iterate over the lines in an MC-Annotate file and add information
    about interactions between non-adjacent elements.

    :param cg: A CoarseGrainRNA structure
    :param lines: All the lines in an MC-Annotate file
    '''
    for line in ftum.iterate_over_interactions(lines):
        (from_chain, from_base, to_chain, to_base) =  ftum.get_interacting_base_pairs(line)

        seq_id1 = cg.seq_ids.index(ftum.parse_resid(from_base)) + 1
        seq_id2 = cg.seq_ids.index(ftum.parse_resid(to_base)) + 1

        node1 = cg.get_node_from_residue_num(seq_id1)
        node2 = cg.get_node_from_residue_num(seq_id2)

        if abs(seq_id2 - seq_id1) > 1 and node1 != node2 and not cg.has_connection(node1, node2):
            cg.longrange[node1].add(node2)
            cg.longrange[node2].add(node1)

def load_cg_from_pdb_in_dir(pdb_filename, output_dir, secondary_structure='', 
                            chain_id=None, remove_pseudoknots=True, parser=None):
    '''
    Create the coarse grain model from a pdb file and store all
    of the intermediate files in the given directory.

    :param pdb_filename: The name of the pdb file to be coarseified
    :param output_dir: The name of the output directory
    :param secondary_structure: Specify a particular secondary structure
                                for this coarsification.
    :param chain_id: The id of the chain to create the CG model from
    '''
    if parser == None:
        parser = bpdb.PDBParser()

    #chain = ftup.load_structure(pdb_filename)
    if chain_id == None:
        chain = ftup.get_biggest_chain(pdb_filename, parser=parser)
    else:
        chain = ftup.get_particular_chain(pdb_filename, chain_id, parser=parser)

    chain = ftup.rename_modified_ress(chain)
    chain = ftup.rename_rosetta_atoms(chain)
    chain = ftup.remove_hetatm(chain)

    # output the biggest RNA chain
    pdb_base = op.splitext(op.basename(pdb_filename))[0]
    output_dir = op.join(output_dir, pdb_base + "_" + chain.id)

    if not op.exists(output_dir):
        os.makedirs(output_dir)

    with open(op.join(output_dir, 'temp.pdb'), 'w') as f:
        # TODO: the following should be changed to take the input parser
        # and use that to output the chain
        ftup.output_chain(chain, f.name)
        f.flush()

        pdb_base = op.splitext(op.basename(pdb_filename))[0]
        pdb_base += "_" + chain.id

        cg = CoarseGrainRNA()
        cg.name = pdb_base

        if len(chain.get_list()) == 0:
            return cg

        # first we annotate the 3D structure
        p = sp.Popen(['MC-Annotate', f.name], stdout=sp.PIPE)
        out, err = p.communicate()

        with open(op.join(output_dir, 'temp.mcannotate'), 'w') as f3:
            f3.write(out)

        lines = out.strip().split('\n')
        # convert the mcannotate output into bpseq format
        try:
            (dotplot, residue_map) = ftum.get_dotplot(lines)
        except Exception as e:
            print (e, file=sys.stderr)
            return cg

        # f2 will store the dotbracket notation
        with open(op.join(output_dir, 'temp.bpseq'), 'w') as f2:
            f2.write(dotplot)
            f2.flush()

            # remove pseudoknots
            '''
            p = sp.Popen(['aux/k2n_standalone/knotted2nested.py', '-f', 'bpseq', 
                          '-F', 'vienna', f2.name], stdout = sp.PIPE)

            out, err = p.communicate()
            '''
            if remove_pseudoknots:
                out = cak.k2n_main(f2.name, input_format='bpseq',
                                   #output_format = 'vienna',
                                   output_format = 'bpseq',
                                   method = cak.DEFAULT_METHOD,
                                   opt_method = cak.DEFAULT_OPT_METHOD,
                                   verbose = cak.DEFAULT_VERBOSE,
                                   removed= cak.DEFAULT_REMOVED)

                out = out.replace(' Nested structure', pdb_base)
            else:
                out = dotplot
            #(out, residue_map) = add_missing_nucleotides(out, residue_map)

            '''
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
            '''

            # Add the 3D information about the starts and ends of the stems
            # and loops
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = bpdb.PDBParser().get_structure('temp', f.name)
                chains = list(s.get_chains())
                if len(chains) < 1:
                    raise Exception("No chains in the PDB file")

                chain = chains[0]

            #cg.from_fasta(out, dissolve_length_one_stems=1)
            cg.from_bpseq_str(out, dissolve_length_one_stems=False)
            cg.name = pdb_base
            cg.seqids_from_residue_map(residue_map)
            ftug.add_stem_information_from_pdb_chain(cg, chain)
            ftug.add_bulge_information_from_pdb_chain(cg, chain)
            ftug.add_loop_information_from_pdb_chain(cg, chain)

            cg.chain = chain

            add_longrange_interactions(cg, lines)

            with open(op.join(output_dir, 'temp.cg'), 'w') as f3:
                f3.write(cg.to_cg_string())
                f3.flush()

            return cg
    print ("Uh oh... couldn't generate the coarse-grain structure.", file=sys.stderr)
    print ("Prepare for an incoming exception.", file=sys.stderr)

def load_cg_from_pdb(pdb_filename, secondary_structure='', 
                     intermediate_file_dir=None, chain_id=None,
                    remove_pseudoknots=True, parser=None):
    '''
    Load a coarse grain model from a PDB file, by extracing
    the bulge graph.

    :param pdb_filename: The filename of the 3D model
    :param secondary_structure: A dot-bracket string encoding the secondary
                                structure of this molecule
    '''

    if intermediate_file_dir is not None:
        output_dir = intermediate_file_dir

        cg = load_cg_from_pdb_in_dir(pdb_filename, output_dir, 
                                     secondary_structure, chain_id=chain_id,
                                    remove_pseudoknots=remove_pseudoknots, parser=parser)
    else:
        with fus.make_temp_directory() as output_dir:
            cg = load_cg_from_pdb_in_dir(pdb_filename, output_dir, 
                                         secondary_structure, chain_id = chain_id,
                                        remove_pseudoknots=remove_pseudoknots, parser=parser)

    return cg

def from_file(cg_filename):
    '''
    Read a coarse-grain structure file.

    :param cg_filename: The filename.
    :return: A CoarseGrainRNA from a file.
    '''
    with open(cg_filename, 'r') as f:
        lines = "".join(f.readlines())

        cg = CoarseGrainRNA()
        cg.from_cg_string(lines)

        return cg
    
def from_pdb(pdb_filename, secondary_structure='', intermediate_file_dir=None, 
             chain_id=None, remove_pseudoknots=True, parser=None):
    cg = load_cg_from_pdb(pdb_filename, secondary_structure, 
                          intermediate_file_dir, chain_id=chain_id,
                         remove_pseudoknots=remove_pseudoknots, parser=parser)

    return cg

class RnaMissing3dError(LookupError):
    pass

class CoarseGrainRNA(fgb.BulgeGraph):
    '''
    A coarse grain model of RNA structure based on the
    bulge graph representation.

    Each stem is represented by four parameters (two endpoints)
    and two twist vetors pointing towards the centers of the base
    pairs at each end of the helix.
    '''
    def __init__(self, cg_file=None, dotbracket_str='', seq=''):
        '''
        Initialize the new structure.
        '''
        super(CoarseGrainRNA, self).__init__(cg_file, dotbracket_str, seq)

        self._virtual_atom_cache={}
        #: Keys are element identifiers (e.g.: "s1" or "i3"), values are 2-tuples of vectors
        #: The first value corresponds to the start of the stem 
        #: (the one with the lowest nucleotide number),
        #: The second value to the end of the stem.
        #: If the coordinates for an element change, the virtual atom and virtual residue 
        #: coordinates are automatically invalidated.
        self.coords = observedDict(on_change=self.reset_vatom_cache)
        self.twists = observedDict(on_change=self.reset_vatom_cache)
        self.sampled = dict()
        
        #: Global (carthesian) position of the virtual residue 
        #: (=offset of the residue's coordinate-system)
        #: generated by self.add_all_virtual_residues()
        self.vposs = c.defaultdict( dict )    
        #: The coordinate system specific to each virtual residue 
        #: (3x3 matrix, carthesian coordiantes; each row is one unit vector)
        #: generated by self.add_all_virtual_residues()
        self.vbases = c.defaultdict( dict )
        #: generated by self.add_all_virtual_residues()
        self.vvecs = c.defaultdict( dict )
        #: generated by self.add_all_virtual_residues()
        self.v3dposs = c.defaultdict( dict )
        #: generated by self.add_all_virtual_residues()
        self.vinvs = c.defaultdict( dict )

        #: A 3D vector. Used as hint, from what direction the Projection2D object 
        #: should be generated in the default case.
        self.project_from = None

        self.longrange = c.defaultdict( set )
        self.chain = None #the PDB chain if loaded from a PDB file

        if cg_file is not None:
            self.from_file(cg_file)
        pass

    def get_coord_str(self):
        '''
        Place the start and end coordinates of each stem into a string.

        The format is:

            coord s1 x1 y1 z1 x2 y2 z2

        Where s1 is the name of the stem, (x1,y1,z1) and (x2,y2,z2) are
        the two endpoints of the stem.

        :return: A string containing the coordinates for all of the stems.
        '''
            
        out_str = ''
        for key in self.coords.keys():
            [p, n] = self.coords[key]
            out_str += "coord %s %s %s" % (key, " ".join([str(pt) for pt in p]),
                                           " ".join([str(pt) for pt in n]))
            out_str += '\n'
        return out_str

    def add_all_virtual_residues(self):
        """
        Calls ftug.add_virtual_residues() for all stems of this RNA.
    
        .. note::
           Don't forget to call this again if you changed the structure of the RNA,
           to avoid leaving it in an inconsistent state.

        .. warning::
           Virtual residues are only added to stems, not to loop regions.
           The position of residues in loops is much more flexible, which is why virtual 
           residue positions for loops usually do not make sense.
        """
        for stem in self.stem_iterator():
            try:
                ftug.add_virtual_residues(self, stem)
            except KeyError:
                if stem not in self.coords:
                    raise RnaMissing3dError("No 3D coordinates available for stem {}".format(stem))
                elif stem not in self.twists:
                    raise RnaMissing3dError("No twists available for stem {}".format(stem))
                else: 
                    raise

    def get_ordered_stem_poss(self):
        points = []
        for s in self.sorted_stem_iterator():
            points += [self.coords[s][0]]
            points += [self.coords[s][1]]
        return np.array(points)

    def get_ordered_virtual_residue_poss(self):
        """
        Get the coordinates of the virtual residues in a consistent order.
        
        This is used for RMSD calculation and is ment to replace ftug.bg_virtual_residues
        If no virtual_residue_positions are known, self.add_all_virtual_residues() is called 
        automatically.
    
        :returns: A numpy array.
        """
        if not self.v3dposs:
            self.add_all_virtual_residues()
        vress = []
        for s in self.sorted_stem_iterator():
            for i in range(self.stem_length(s)):
                vres=self.v3dposs[s][i]
                vress += [vres[0] + vres[2], vres[0] + vres[3]]
        return np.array(vress)

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
        printed = set()

        for key1 in self.longrange.keys():
            for key2 in self.longrange[key1]:
                k = [key1, key2]
                k.sort()
                k = tuple(k)
                if k in printed:
                    continue
                printed.add(k)

                out_str += "longrange %s %s\n" % (key1, key2)

        return out_str

    def get_sampled_stems_str(self):
        out_str = ''
        for key in self.sampled.keys():
            out_str += 'sampled %s %s\n' % (key, " ".join(map(str, self.sampled[key])))
        return out_str

    def to_cg_string(self):
        '''
        Output this structure in string form.
        '''
        curr_str = self.to_bg_string()
        curr_str += self.get_coord_str()
        curr_str += self.get_twist_str()
        curr_str += self.get_sampled_stems_str()
        curr_str += self.get_long_range_str()
        if self.project_from is not None:
            curr_str+="project {} {} {}\n".format(*self.project_from)
        return curr_str

    def to_file(self, filename):
        with open(filename, 'w') as f:
            cg_str = self.to_cg_string()
            f.write(cg_str)

    def get_bulge_angle_stats_core(self, define, connections):
        '''
        Return the angle stats for a particular bulge. These stats describe the
        relative orientation of the two stems that it connects.

        :param define: The name of the bulge.
        :param connections: The two stems that are connected by it.
        :return: ftms.AngleStat object
        '''
        (stem1, twist1, stem2, twist2, bulge) = ftug.get_stem_twist_and_bulge_vecs(self, define, connections)

        # Get the orientations for orienting these two stems
        (r, u, v, t) = ftug.get_stem_orientation_parameters(stem1, twist1, 
                                                            stem2, twist2)
        (r1, u1, v1) = ftug.get_stem_separation_parameters(stem1, twist1, bulge)

        dims =self.get_bulge_dimensions(define)
        ang_type = self.connection_type(define, connections)
        seqs = self.get_define_seq_str(define, adjacent=True)

        angle_stat = ftms.AngleStat(self.name, dims[0], dims[1], u, v, t, r1, 
                                    u1, v1, ang_type, self.defines[define], 
                                    seqs)

        return angle_stat

    def get_loop_stat(self, d):
        '''
        Return the statistics for this loop.

        These stats describe the relative orientation of the loop to the stem
        to which it is attached.

        :param d: The name of the loop
        '''
        loop_stat = ftms.LoopStat()
        loop_stat.pdb_name = self.name

        loop_stat.bp_length = self.get_length(d)
        loop_stat.phys_length = ftuv.magnitude(self.coords[d][1] - self.coords[d][0])

        stem1 = list(self.edges[d])[0]
        (s1b, s1e) = self.get_sides(stem1, d)

        stem1_vec = self.coords[stem1][s1b] - self.coords[stem1][s1e]
        twist1_vec = self.twists[stem1][s1b]
        bulge_vec = self.coords[d][1] - self.coords[d][0] + 0.1 * (stem1_vec / ftuv.magnitude(stem1_vec))

        (r,u,v) = ftug.get_stem_separation_parameters(stem1_vec, twist1_vec, 
                                                      bulge_vec)
        (loop_stat.r, loop_stat.u, loop_stat.v) = (r,u,v)

        return loop_stat

    def get_bulge_angle_stats(self, bulge):                                                                           
        '''
        Return the angle stats for a particular bulge. These stats describe the                                       
        relative orientation of the two stems that it connects.                                                       

        :param bulge: The name of the bulge.
        :param connections: The two stems that are connected by it.
        :return: The angle statistics in one direction and angle statistics in
                 the other direction                    
        '''  
        if bulge == 'start':
            return (ftms.AngleStat(), cbs.AngleStat())

        connections = self.connections(bulge)

        angle_stat1 = self.get_bulge_angle_stats_core(bulge, connections)
        angle_stat2 = self.get_bulge_angle_stats_core(bulge, list(reversed(connections)))

        return (angle_stat1, angle_stat2)
    
    def get_stem_stats(self, stem):                                                                                   
        '''
        Calculate the statistics for a stem and return them. These statistics will describe the                       
        length of the stem as well as how much it twists.                                                             

        :param stem: The name of the stem.                                                                            

        :return: A StemStat structure containing the above information.                                               
        '''                                                                                                           
        ss = ftms.StemStat()
        ss.pdb_name = self.name
        #ss.bp_length = abs(self.defines[stem][0] - self.defines[stem][1])                                            
        ss.bp_length = self.stem_length(stem)
        ss.phys_length = ftuv.magnitude(self.coords[stem][0] - self.coords[stem][1])
        ss.twist_angle = ftug.get_twist_angle(self.coords[stem], self.twists[stem])
        ss.define = self.defines[stem]

        return ss  

    #def get_loop_from_residue(self, residue) ->  use BulgeGraph.get_node_from_residue_num()!

    def from_file(self, cg_filename):
        '''
        Load this data structure from a file.
        '''
        with open(cg_filename, 'r') as f:
            lines = "".join(f.readlines())

            self.from_cg_string(lines)

    def from_cg_string(self, cg_string):
        '''
        Populate this structure from the string
        representation of a graph.
        '''
        # Reading the bulge_graph-part of the file
        self.from_bg_string(cg_string)

        #Reading the part of the file responsible for 3D information
        lines = cg_string.split('\n')
        for line in lines:
            line = line.strip()
            parts = line.split()
            if len(parts) == 0:
                continue
            if parts[0] == 'coord':
                name = parts[1]
                self.coords[name] = np.array([list(map(float, parts[2:5])), list(map(float, parts[5:8]))])
            if parts[0] == 'twist':
                name = parts[1]
                self.twists[name] = np.array([list(map(float, parts[2:5])), list(map(float, parts[5:8]))])
            if parts[0] == 'longrange':
                self.longrange[parts[1]].add(parts[2])
                self.longrange[parts[2]].add(parts[1])

            if parts[0] == 'sampled':
                self.sampled[parts[1]] = [parts[2]] + list(map(int, parts[3:]))
            if parts[0] == 'project':
                self.project_from=np.array(parts[1:], dtype=float)

    def to_cg_file(self, filename):
        '''
        Save this structure as a string in a file.
average_stem_vres_atom_positions_new.py
        :param filename: The filename to save it to.
        '''
        with open(filename, 'w') as f:
            s = self.to_cg_string()

            f.write(s)

    def radius_of_gyration(self):
        '''
        Calculate the radius of gyration of this structure.

        :return: A number with the radius of gyration of this structure.
        '''
        coords = []

        for s in self.sorted_stem_iterator():
            coords += list(self.coords[s])

        rmsd = ftur.radius_of_gyration(coords)

        return rmsd

    def get_coordinates_list(self):
        '''
        Get all of the coordinates as a list of vectors.

        The vectors are sorted according to the name of the
        element they correspond to.

        :return: A list of vectors.
        '''
        all_coords = []

        for d in sorted(self.defines):
            all_coords += [self.coords[d][0]]
            all_coords += [self.coords[d][1]]

        return all_coords

    def get_coordinates_array(self):
        '''
        Get all of the coordinates in one large array.

        The coordinates are sorted in the order of the keys
        in coordinates dictionary.

        :return: One large array containing a sequential list of coordinates
        '''
        all_coords = []
        for key in sorted(self.coords.keys()):
            for i in range(len(self.coords[key])):
                all_coords += list(self.coords[key][i])
        return all_coords

    def load_coordinates_array(self, coords):
        '''
        Read in an array of coordinates (as may be produced by get_coordinates_array)
        and replace the coordinates of this structure with it.

        :param coords: A 1D array of coordinates
        :return: self
        '''
        counter = 0

        for key in sorted(self.coords.keys()):
            for i in range(len(self.coords[key])):
                self.coords[key][i] = coords[counter:counter+3]
                counter += 3
        return self

    def get_twists(self, node):
        ''' 
        Get the array of twists for this node. If the node is a stem,
        then the twists will simply those stored in the array. 
        If the node is an interior loop or a junction segment, 
        then the twists will be the ones that are adjacent to it. 
        If the node is a hairpin loop or a free end, then the same twist
        will be duplicated and returned twice.

        :param node: The name of the node
        '''                                                                                                           
        if node[0] == 's':
            return self.twists[node]                                                                                  

        connections = list(self.edges[node])
        (s1b, s1e) = self.get_sides(connections[0], node)

        if len(connections) == 1:
            vec = ftuv.normalize(ftuv.vector_rejection( 
                                  self.twists[connections[0]][s1b],
                                  self.coords[connections[0]][1] -  
                                  self.coords[connections[0]][0]))
            return (vec,vec)                                                  

        if len(connections) == 2: 
            # interior loop or junction segment                                                                  
            (s2b, s2e) = self.get_sides(connections[1], node) 
            bulge_vec = (self.coords[connections[0]][s1b] - 
                         self.coords[connections[1]][s2b])                                                            
            return (ftuv.normalize(ftuv.vector_rejection( 
                    self.twists[connections[0]][s1b], bulge_vec)),
                    ftuv.normalize(ftuv.vector_rejection(self.twists[connections[1]][s2b], bulge_vec)))  

        # uh oh, this shouldn't happen since every node                 
        # should have either one or two edges 
        return None                                                                                                   
    def element_physical_distance(self, element1, element2):
        '''
        Calculate the physical distance between two coarse grain elements.

        :param element1: The name of the first element (e.g. 's1')
        :param element2: The name of the first element (e.g. 's2')
        :return: The closest distance between the two elements.
        '''
        (i1, i2) = ftuv.line_segment_distance(self.coords[element1][0],
                                              self.coords[element1][1],
                                              self.coords[element2][0],
                                              self.coords[element2][1])

        return ftuv.vec_distance(i1, i2)
    

    def longrange_iterator(self, filter_connected=False):
        '''
        Iterate over all long range interactions in this molecule.
        
        :param filter_connected: Filter interactions that are between elements
                                 which are connected (mostly meaning multiloops
                                 which connect to the same end of the same stem)
        :return: A generator yielding long-range interaction tuples (i.e. ('s7', 'i2'))
        '''
        seen = set()

        for partner1 in self.longrange.keys():
            for partner2 in self.longrange[partner1]:
                if filter_connected:
                    if self.connected(partner1, partner2):
                        continue

                interaction = tuple(sorted([partner1, partner2]))

                # check if we've already seen this interaction
                if interaction in seen:
                    continue

                seen.add(interaction)

                yield interaction

    def total_length(self):
        '''
        Calculate the combined length of all the elements.
        '''
        return sum([len(list(self.define_residue_num_iterator(d))) for d in self.defines])

    def virtual_atoms(self, key):
        """
        Get virtual atoms for a key.

        .. note:: In the current implementation, the virtual atoms 
                  do not rely on the virtual residues.

        :param key: An INTEGER: The number of the base in the RNA. 
                    Returns a dict {"C8":np.array([x,y,z]), ...}
        """
        if isinstance(key, int):
             if key not in self._virtual_atom_cache:
                try:
                    self._virtual_atom_cache[key]=ftug.virtual_atoms(self)[key]
                except KeyError:
                    self.add_all_virtual_residues()
                    self._virtual_atom_cache[key]=ftug.virtual_atoms(self)[key]
             return self._virtual_atom_cache[key]
        else:
            raise ValueError("Expected an int, found {}".format(key))
    def reset_vatom_cache(self, key):
        """
        Used as on_call function for the observing of the self.coords dictionary.

        :param key: A coarse grain element name, e.g. "s1" or "m15"
        """
        try:
            if not self._virtual_atom_cache:
                return
        except AttributeError: #Happens during deepcopy
            return 

        define=self.defines[key]

        #Delete virtual residues
        try: del self.vposs[key]
        except KeyError: pass
        try: del self.vbases[key]
        except KeyError: pass
        try: del self.vvecs[key]
        except KeyError: pass
        try: del self.v3dposs[key]
        except KeyError: pass
        try: del self.vinvs[key]
        except KeyError: pass

        #Delete virtual atoms
        if len(define)>1:
            for i in range(define[0], define[1]+1):
                if i in self._virtual_atom_cache:
                    del self._virtual_atom_cache[i]
        if len(define)>3:
            for i in range(define[2], define[3]+1):
                if i in self._virtual_atom_cache:
                    del self._virtual_atom_cache[i]
    #def __deepcopy__(self, memo):
        
def cg_from_sg(cg, sg):
    '''
    Create a coarse-grain structure from a subgraph.
    
    :param cg: The original structure
    :param sg: The list of elements that are in the subgraph
    '''
    new_cg = CoarseGrainRNA()
    
    for d in sg:
        new_cg.defines[d] = cg.defines[d]
        new_cg.coords[d] = cg.coords[d]
        if d in cg.twists:
            new_cg.twists[d] = cg.twists[d]
        new_cg.longrange[d] = cg.longrange[d]
        
        for x in cg.edges[d]:
            if x in new_cg.defines.keys():
                new_cg.edges[d].add(x)
                new_cg.edges[x].add(d)
    
    return new_cg

