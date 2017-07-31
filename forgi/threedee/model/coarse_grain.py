from __future__ import absolute_import, unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)
from past.types import basestring

import collections as c
import contextlib
import os
import os.path as op
import shutil
import subprocess as sp
import sys
import tempfile as tf
import time
import math
import warnings
import itertools as it
try:
    import io
except ImportError: #Python 2
    import StringIO as io
import logging
from pprint import pprint


import numpy as np
import scipy.spatial
import scipy.stats

import Bio.PDB as bpdb

from logging_exceptions import log_to_exception, log_exception


from ...graph import bulge_graph as fgb
from ..utilities import graph_pdb as ftug
from ..model import stats as ftms

from ...aux.k2n_standalone import knotted2nested as cak
from ..utilities import mcannotate as ftum
from ..utilities import pdb as ftup
from . import descriptors as ftud
from ..utilities import vector as ftuv
from ...utilities import debug as fud
from ...utilities import stuff as fus
from ...utilities.observedDict import observedDict
from ...utilities.exceptions import CgConstructionError, CgIntegrityError, GraphConstructionError
from .linecloud import CoordinateStorage, LineSegmentStorage




log = logging.getLogger(__name__)



try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x):
    return x

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
        try:
            seq_id1 = cg.seq_ids.index(fgb.RESID(from_chain, ftum.parse_resid(from_base))) + 1
        except ValueError as e:
            with log_to_exception(log, e):
                log.error("seq_ids are {}".format(cg.seq_ids))
            raise
        seq_id2 = cg.seq_ids.index(fgb.RESID(to_chain, ftum.parse_resid(to_base))) + 1

        node1 = cg.get_node_from_residue_num(seq_id1)
        node2 = cg.get_node_from_residue_num(seq_id2)

        if abs(seq_id2 - seq_id1) > 1 and node1 != node2 and not cg.has_connection(node1, node2):
            cg.longrange[node1].add(node2)
            cg.longrange[node2].add(node1)

def breakpoints_from_residuemap(residue_map):
    """
    Create the list of cofold cutpoints from the list of MC-Annotate identifiers in the
    residue map.
    """
    breakpoints = []
    old_chain = None
    log.debug(residue_map)
    for i, r in enumerate(residue_map):
        (from_chain, _) = ftum.parse_chain_base(r)
        if old_chain is not None and from_chain != old_chain:
            breakpoints.append(i)
        old_chain = from_chain
    return breakpoints

def _are_adjacent_basepairs(cg, edge1, edge2):
    """
    Helper function used by connected_cgs_from_pdb.

    :param cg: A (potentially) inconsistent cg. Only seq_ids are needed.
    :param edge1, edge2: A dictionary with the key "basepair" mapping to a two-tuple of RESIDS
    """
    fromA, toA = edge1["basepair"]
    fromB, toB = edge2["basepair"]
    if fromA.chain != fromB.chain:
        if fromA.chain==toB.chain:
            fromB, toB = toB, fromB
            assert toA.chain==toB.chain
        else:
            assert False
    if (abs(cg.seq_ids.index(fromA)-cg.seq_ids.index(fromB))==1 and
        abs(cg.seq_ids.index(toA)-cg.seq_ids.index(toB))==1):
        log.debug("Basepairs %s and %s are adjacent. No length 1 stem", edge1, edge2)
        return True
    log.debug("Basepairs %s and %s are NOT adjacent.", edge1, edge2)
    return False


def connected_cgs_from_pdb(pdb_filename, remove_pseudoknots=False, dissolve_length_one_stems = True):
    with fus.make_temp_directory() as output_dir:
        chains = ftup.get_all_chains(pdb_filename)
        new_chains = []
        for chain in chains:
            log.debug("Loaded Chain %s", chain.id)
            chain = ftup.clean_chain(chain)
            new_chains.append(chain)

        rna_pdb_fn = op.join(output_dir, 'temp.pdb')
        with open(rna_pdb_fn, 'w') as f:
            #We need to output in pdb format for MC-Annotate
            ftup.output_multiple_chains(new_chains, f.name)
            f.flush()

        # first we annotate the 3D structure
        log.info("Starting MC-Annotate for all chains")
        out = sp.check_output(['MC-Annotate', rna_pdb_fn], universal_newlines=True)
        lines = out.strip().split('\n')

        # convert the mcannotate output into bpseq format
        try:
            (dotplot, residue_map) = ftum.get_dotplot(lines)
        except Exception as e:
            print (e, file=sys.stderr)
            return []

        # f2 will store the dotbracket notation
        #with open(op.join(output_dir, 'temp.bpseq'), 'w') as f2:
        #    f2.write(dotplot)
        #    f2.flush()

        # remove pseudoknots
        if remove_pseudoknots:
            log.info("Removing pseudoknots")
            out = cak.k2n_main(io.StringIO(dotplot), input_format='bpseq',
                               #output_format = 'vienna',
                               output_format = 'bpseq',
                               method = cak.DEFAULT_METHOD,
                               opt_method = cak.DEFAULT_OPT_METHOD,
                               verbose = cak.DEFAULT_VERBOSE,
                               removed= cak.DEFAULT_REMOVED)
        else:
            out = dotplot
        #(out, residue_map) = add_missing_nucleotides(out, residue_map)

        breakpoints = breakpoints_from_residuemap(residue_map)


        import networkx as nx
        chain_connections_multigraph = nx.MultiGraph()
        cg = CoarseGrainRNA()
        cg.seqids_from_residue_map(residue_map)

        for line in out.splitlines():
            try:
                from_, res, to_ = line.split()
            except:
                continue
            from_seqid = cg.seq_ids[int(from_)-1]
            from_chain = from_seqid.chain
            chain_connections_multigraph.add_node(from_chain)
            if int(to_) != 0 and int(to_)>int(from_):
                to_seqid = cg.seq_ids[int(to_)-1]
                from_seqid = cg.seq_ids[int(from_)-1]
                to_seqid = cg.seq_ids[int(to_)-1]
                from_chain = from_seqid.chain
                to_chain   = to_seqid.chain
                if from_chain != to_chain:
                    log.debug("Adding {} - {}".format(from_chain, to_chain))
                    chain_connections_multigraph.add_edge(from_chain, to_chain, basepair=(from_seqid, to_seqid))
        chain_connections = nx.Graph(chain_connections_multigraph)
        if dissolve_length_one_stems:
            for chain1, chain2 in it.combinations(chain_connections_multigraph.nodes(), 2):
                if chain2 in chain_connections_multigraph.edge[chain1]:
                    for edge1, edge2 in it.combinations(chain_connections_multigraph.edge[chain1][chain2].values(), 2):
                        if _are_adjacent_basepairs(cg, edge1, edge2):
                            break
                    else: #break NOT encountered.
                        log.debug("No connection remaining between chains %s and %s", chain1, chain2)
                        chain_connections.remove_edge(chain1, chain2)

        # Now remove multiple edges. We don't care what nx does with the edge attributes.

        log.debug("CONNECTIONS: {}, nodes {}".format(chain_connections, chain_connections.nodes()))
        log.debug("Edges {}".format(chain_connections.edges()))
        cgs = []
        for component in nx.connected_components(chain_connections):
            #print(component, type(component))
            log.info("Loading PDB: Connected component with chains %s", str(list(component)))
            try:
                cgs.append(load_cg_from_pdb(pdb_filename, remove_pseudoknots=remove_pseudoknots, chain_id = list(component)))
            except GraphConstructionError as e:
                log_exception(e, logging.ERROR, with_stacktrace=False)
                log.error("Could not load chains %s, due to the above mentioned error.", list(component))
        cgs.sort(key=lambda x: x.name)
        if dissolve_length_one_stems:
            for cg in cgs:
                cg.dissolve_length_one_stems()
        return cgs

def load_cg_from_pdb_in_dir(pdb_filename, output_dir, secondary_structure='',
                            chain_id=None, remove_pseudoknots=True, parser=None,
                            dissolve_length_one_stems=True):
    '''
    Create the coarse grain model from a pdb file and store all
    of the intermediate files in the given directory.

    :param pdb_filename: The name of the pdb file to be coarseified
    :param output_dir: The name of the output directory
    :param secondary_structure: Specify a particular secondary structure
                                for this coarsification.
    :param chain_id: The id of the chain to create the CG model from.
                    This can be one of the following:

                    * The string "all" in lowercase, to extract all chains.
                    * None, to extract the biggest chain
                    * A (single-letter) string containing the chain id
                    * A list of chain-ids. In that case only chains that match this id will be extracted.
    '''

    #chain = ftup.load_structure(pdb_filename)
    chains = []
    if chain_id == "all":
        chains = ftup.get_all_chains(pdb_filename, parser=parser)
    elif chain_id is None:
        chains = [ftup.get_biggest_chain(pdb_filename, parser=parser)]
    elif isinstance(chain_id, basestring):
        chains = [ftup.get_particular_chain(pdb_filename, chain_id, parser=parser)]
    else:
        chains = ftup.get_all_chains(pdb_filename, parser=parser)
        chains = [ chain for chain in chains if chain.id in chain_id ]
        if len(chain_id) != len(chains):
            raise CgConstructionError("Bad chain-id given. "
                                      "{} not present (or not an RNA)".format( set(chain_id) -
                                                                                       set([chain.id for chain in chains])))

    new_chains = []
    for chain in chains:
        chain = ftup.clean_chain(chain)
        new_chains.append(chain)

    pdb_base = op.splitext(op.basename(pdb_filename))[0]

    if len(new_chains)==1:
        pdb_base += "_" + new_chains[-1].id
    else:
        pdb_base += "_" + "-".join(chain.id for chain in new_chains)

    cg = CoarseGrainRNA()
    if sum(len(chain.get_list()) for chain in new_chains) == 0:
        return cg

    # output a pdb with RNA only (for MCAnnotate):
    output_dir = op.join(output_dir, pdb_base )

    if not op.exists(output_dir):
        os.makedirs(output_dir)
    rna_pdb_fn = op.join(output_dir, 'temp.pdb')
    with open(rna_pdb_fn, 'w') as f:
        #We need to output in pdb format for MC-Annotate
        ftup.output_multiple_chains(new_chains, f.name)
        f.flush()

    # first we annotate the 3D structure
    log.info("Starting MC-Annotate for chain(s) %s", chain_id)
    if hasattr(sp, "DEVNULL"): #python>=3.3
        kwargs = {"stderr":sp.DEVNULL}
    else:
        kwargs = {}
    out = sp.check_output(['MC-Annotate', rna_pdb_fn], universal_newlines=True, **kwargs)

    lines = out.strip().split('\n')

    # convert the mcannotate output into bpseq format
    try:
        (dotplot, residue_map) = ftum.get_dotplot(lines)
    except Exception as e:
        log.exception("Could not convert MC-Annotate output to dotplot")
        raise CgConstructionError("Could not convert MC-Annotate output to dotplot") #from e

    # f2 will store the dotbracket notation
    #with open(op.join(output_dir, 'temp.bpseq'), 'w') as f2:
    #    f2.write(dotplot)
    #    f2.flush()


    if secondary_structure:
        warnings.warn("Not adding any longrange interactions because secondary structure is given.")
        if remove_pseudoknots:
            warnings.warn("Option 'remove_pseudoknots ignored, because secondary structure is given.")
        cg.from_dotbracket(secondary_structure, dissolve_length_one_stems=dissolve_length_one_stems)
        breakpoints = breakpoints_from_residuemap(residue_map)
        cg.seqids_from_residue_map(residue_map)

    else:
        if remove_pseudoknots:
            log.info("Removing pseudoknots")
            out = cak.k2n_main(io.StringIO(dotplot), input_format='bpseq',
                               #output_format = 'vienna',
                               output_format = 'bpseq',
                               method = cak.DEFAULT_METHOD,
                               opt_method = cak.DEFAULT_OPT_METHOD,
                               verbose = cak.DEFAULT_VERBOSE,
                               removed= cak.DEFAULT_REMOVED)

            out = out.replace(' Nested structure', pdb_base)
        else:
            out = dotplot
        breakpoints = breakpoints_from_residuemap(residue_map)
        log.info("Breakpoints in %s are %s",pdb_base, breakpoints)
        cg.from_bpseq_str(out, breakpoints = breakpoints, dissolve_length_one_stems=dissolve_length_one_stems) #Sets the seq without cutpoints
        cg.seqids_from_residue_map(residue_map)
        add_longrange_interactions(cg, lines)


    cg.name = pdb_base
    # Add the 3D information about the starts and ends of the stems
    # and loops
    cg.chains = { chain.id : chain for chain in new_chains }
    cg.chain_ids = [ chain.id for chain in new_chains ]

    log.debug("First 10 seq-IDs of loaded structure are %s", cg.seq_ids[:10])
    log.debug("Elements %s", cg.defines.keys())
    #Stems can span 2 chains.
    ftug.add_stem_information_from_pdb_chains(cg)
    cg.add_bulge_coords_from_stems()
    ftug.add_loop_information_from_pdb_chains(cg)
    cg.incomplete_elements = ftug.get_incomplete_elements(cg)

    assert len(cg.defines)==len(cg.coords), cg.defines.keys()^cg.coords.keys()
    #with open(op.join(output_dir, 'temp.cg'), 'w') as f3:
    #    f3.write(cg.to_cg_string())
    #    f3.flush()

    log.debug("Sequence {}".format(cg.seq))

    return cg

def load_cg_from_pdb(pdb_filename, secondary_structure='',
                     intermediate_file_dir=None, chain_id=None,
                    remove_pseudoknots=True, parser=None,
                    dissolve_length_one_stems=True):
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
                                    remove_pseudoknots=remove_pseudoknots, parser=parser,
                                    dissolve_length_one_stems=dissolve_length_one_stems)
    else:
        with fus.make_temp_directory() as output_dir:
            cg = load_cg_from_pdb_in_dir(pdb_filename, output_dir,
                                         secondary_structure, chain_id = chain_id,
                                        remove_pseudoknots=remove_pseudoknots, parser=parser,
                                        dissolve_length_one_stems=dissolve_length_one_stems)

    return cg

def from_pdb(pdb_filename, secondary_structure='', intermediate_file_dir=None,
             chain_id=None, remove_pseudoknots=True, parser=None, dissolve_length_one_stems=True):
    cg = load_cg_from_pdb(pdb_filename, secondary_structure,
                          intermediate_file_dir, chain_id=chain_id,
                         remove_pseudoknots=remove_pseudoknots, parser=parser,
                         dissolve_length_one_stems=dissolve_length_one_stems)
    return cg

def from_bulge_graph(bulge_graph):
    """
    Create a CoarseGrainRNA from a BulgeGraph
    """
    if not bulge_graph.seq:
        raise CgConstructionError("Sequence missing in BulgeGraph. Cannot create CoarseGrainRNA.")
    bg_str = bulge_graph.to_bg_string()
    cg = CoarseGrainRNA()
    cg.from_cg_string(bg_str)
    return  cg

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
        #: The first value of stem coordinates corresponds to the start of the stem
        #: (the one with the lowest nucleotide number),
        #: The second value to the end of the stem.
        #: If the coordinates for an element change, the virtual atom and virtual residue
        #: coordinates are automatically invalidated.
        self.coords = None #We can only initialize this, when we know defines.keys()
        self.twists = None
        self.sampled = dict()

        if self.defines:
            self._init_coords()

        #:The following 5 defaultdicts are cleared when coords or twists change.
        #: Global (carthesian) position of the virtual residue
        #: (=offset of the residue's coordinate-system)
        #: generated by self.add_all_virtual_residues()
        self.vposs = c.defaultdict( dict )
        #: The coordinate system specific to each virtual residue
        #: (3x3 matrix, carthesian coordiantes)
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
        self.chains = {} #the PDB chain if loaded from a PDB file

        self.incomplete_elements = [] # An estimated list of cg-elements with missing residues.

        if cg_file is not None:
            self.from_file(cg_file)
            if not self.defines:
                raise CgConstructionError("Empty CoarseGrainRNA created. Was '{}' in *.cg/ *.coord file format?".format(cg_file))

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
            if key[0] in ["m", "i"]: continue #Bulge coordinates are redundant. They can be deduced from the stem coordinates.
            [p, n] = self.coords[key]
            out_str += ("coord {k} {x[0]:.16f} {x[1]:.16f} {x[2]:.16f} "
                        "{y[0]:.16f} {y[1]:.16f} {y[2]:.16f}".format(k=key, x=p, y=n))
            out_str += '\n'
        return out_str
    def add_bulge_coords_from_stems(self):
        '''
        Add the information about the starts and ends of the bulges (i and m elements).
        The stems have to be created beforehand.

        This is called during loading of the RNA structure from pdb and from cg files.
        '''
        for d in self.defines.keys():
            if d[0] != 's':
                edges = list(self.edges[d])
                if len(edges) == 2:
                    (s1b, _) = self.get_sides(edges[0], d)
                    (s2b, _) = self.get_sides(edges[1], d)

                    mids1 = self.coords[edges[0]]
                    mids2 = self.coords[edges[1]]

                    #Save coordinates in direction of the strand.
                    if self.get_link_direction(edges[0], edges[1],d)==1:
                        self.coords[d] = (mids1[s1b], mids2[s2b])
                    else:
                        self.coords[d] = (mids2[s2b], mids1[s1b])

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
            except (KeyError, ValueError):
                if np.all(np.isnan(self.coords[stem])):
                    raise RnaMissing3dError("No 3D coordinates available for stem {}".format(stem))
                elif np.all(np.isnan(self.twists[stem])):
                    raise RnaMissing3dError("No twists available for stem {}".format(stem))
                else:
                    log.warning("Reraising in add_all_virtual_residues")
                    raise
    def get_virtual_residue(self, pos, allow_single_stranded = False):
        """
        Get the virtual residue for the nmucleotide at position pos (1-based)

        :param allow_single_stranded: If True and pos is not in a stem, return a
              rough estimate for the residue position instead of raising an error.
              Currenly, for non-stem elements, these positions are on the axis of the cg-element.
        """
        elem = self.get_node_from_residue_num(pos)
        if elem[0]=="s":
            if elem not in self.v3dposs or not self.v3dposs[elem]:
                ftug.add_virtual_residues(self, elem)
            for i in range(self.stem_length(elem)):
                if self.defines[elem][0] + i == pos:
                    vres = self.v3dposs[elem][i]
                    return vres[0] + vres [2]
                elif self.defines[elem][3] - i == pos:
                    vres = self.v3dposs[elem][i]
                    return vres[0] + vres [3]
            else:
                assert False
        else:
            if not allow_single_stranded:
                raise ValueError("Position {} is not in a stem! It is in {}.".format(pos, elem))
            if elem[0]=="h":
                #We estimate the vres position along the axis of the hairpin.
                h_length = self.element_length(elem) / 2
                if pos - self.defines[elem][0] > h_length:
                    l = 2 * h_length - (pos - self.defines[elem][0])
                else:
                    l = pos - self.defines[elem][0]
                perc = l / h_length
            elif elem[0]=="i":
                if pos<=self.defines[elem][1]:
                    l = pos - self.defines[elem][0]
                    perc = l / (self.defines[elem][1]-self.defines[elem][0]+1)
                else:
                    l = pos - self.defines[elem][2]
                    tl = (self.defines[elem][3]-self.defines[elem][2]+1)
                    perc = (tl - l)/tl
            else:
                l = pos - self.defines[elem][0]
                perc = l / self.element_length(elem)
            return self.coords[elem][0] + (self.coords[elem][1]-self.coords[elem][0]+1) * perc

    def get_ordered_stem_poss(self):
        points = []
        for s in self.sorted_stem_iterator():
            points += list(self.coords[s])
        return np.array(points)

    def get_ordered_virtual_residue_poss(self, return_elements = False):
        """
        Get the coordinates of all stem's virtual residues in a consistent order.

        This is used for RMSD calculation and is ment to replace ftug.bg_virtual_residues
        If no virtual_residue_positions are known, self.add_all_virtual_residues() is called
        automatically.

        :param return_elements: In addition to the positions, return a list with
                                the cg-elements these coordinates belong to
        :returns: A numpy array.
        """
        if not self.v3dposs:
            self.add_all_virtual_residues()
        vress = []
        elems = []
        for s in self.sorted_stem_iterator():
            for i in range(self.stem_length(s)):
                vres=self.v3dposs[s][i]
                vress.append(vres[0] + vres[2])
                vress.append(vres[0] + vres[3])
                elems.extend([s]*2)
        if return_elements:
            return np.array(vress), elems
        return np.array(vress)

    def get_poss_for_domain(self, elements, mode="vres"):
        """
        Get an array of coordinates only for the elements specified.

        ..note::

            This code is still experimental in the current version of forgi.

        :param elements: A list of coarse grain element names.
        """
        points=[]
        if mode=="fast":
            for e in sorted(elements):
                points += list(self.coords[e])
            return np.array(points)
        elif mode=="vres":
            if not self.v3dposs:
                self.add_all_virtual_residues()
            vress = []
            for s in sorted(elements):
                if s[0]!="s": continue
                for i in range(self.stem_length(s)):
                    vres=self.v3dposs[s][i]
                    vress += [vres[0] + vres[2], vres[0] + vres[3]]
            return np.array(vress)
        assert False


    def steric_value(self, elements, method = "r**-2"):
        """
        Estimate, how difficult a set of elements was to build,
        by counting the atom density around the center of these elements
        """
        try:
            if isinstance(elements, list) or isinstance(elements, tuple):
                center = ftuv.get_vector_centroid(self.coords[elements])
            elif elements.shape==(3,):
                center = elements
            else: assert False, repr(elements)
        except:
            print(elements, repr(elements))
            raise
        if method == "kde":
        #print(center)
            all_vas = []
            for pos in range(1,self.seq_length+1):
                for va in self.virtual_atoms(pos).values():
                    all_vas.append(va)

            all_vas = np.array(all_vas).T
            log.debug("Shape of all atoms {}".format(all_vas.shape))
            kde = scipy.stats.gaussian_kde(all_vas, 50) #randomly take 50 Angstrom bandwidth
            return kde(center)
        elif method.startswith("r**"):
            power=-int(method[3:5])
            exclude=method[5:]
            if exclude and exclude!="e":
                raise ValueError("Not supported method")
            value = 0
            for pos in range(1,self.seq_length+1):
                if exclude and self.get_node_from_residue_num(pos) in elements:
                    continue
                for va in self.virtual_atoms(pos).values():
                    value+=1/(1+ftuv.vec_distance(va, center))**power
            return value
        elif method.startswith("cutoff"):
            cutoff = float(method.split()[1])
            value = 0
            for pos in range(1,self.seq_length+1):
                for va in self.virtual_atoms(pos).values():
                    if ftuv.vec_distance(va, center)< cutoff:
                        value+=1
            return value
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
            out_str += ("twist {k} {x[0]:.16f} {x[1]:.16f} {x[2]:.16f} "
                        "{y[0]:.16f} {y[1]:.16f} {y[2]:.16f}".format(k=key, x=p, y=n))
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
        log.debug("Define %s connections %s",define, connections)
        (stem1, twist1, stem2, twist2, bulge) = ftug.get_stem_twist_and_bulge_vecs(self, define, connections)
        log.debug("stem1 %s, twist1 %s, stem2 %s, twist2 %s, bulge %s", stem1, twist1, stem2, twist2, bulge)

        if round(np.dot(stem1, twist1),10)!=0 or round(np.dot(stem2, twist2),10)!=0:
            err = CgIntegrityError("The twists are inconsistent. "
                               "They should be orthogonal to the corresponding stem vectors."
                               "Inconsistency found for {},{}".format(define, connections))
            with log_to_exception(log,err):
                log.error("Angle stem1-twist1 %s dot_product=%s, Angle stem2-twist2 %s degrees dot_product=%s",
                          math.degrees(ftuv.vec_angle(stem1, twist1)), np.dot(stem1, twist1),
                          math.degrees(ftuv.vec_angle(stem2, twist2)), np.dot(stem2, twist2),)
            raise err

        try:
            # Get the orientations for orienting these two stems
            (r, u, v, t) = ftug.get_stem_orientation_parameters(stem1, twist1,
                                                                stem2, twist2)
            (r1, u1, v1) = ftug.get_stem_separation_parameters(stem1, twist1, bulge)
        except ZeroDivisionError:
            print ("Cannot get stats for {}. The 3D coodinates are probably wrong.".format(define), file=sys.stderr)
            raise
        dims =self.get_bulge_dimensions(define)
        ang_type = self.connection_type(define, connections)
        seqs = self.get_define_seq_str(define, adjacent=True)
        log.debug("u %s, v %s", u, v)

        if define[0]=="m":
            mls = self.find_mlonly_multiloops()
            for ml in mls:
                if define in ml:
                    descr = self.describe_multiloop([x for x in ml if x[0] !="s"])
                    if "pseudoknot" in descr:
                        stat_type = "pseudo"
                    elif "open" in descr:
                        stat_type = "open"
                    else:
                        stat_type = "angle" #ML
                    break
            else:
                assert False
        else:
            stat_type="angle" #IL


        angle_stat = ftms.AngleStat(stat_type,
                                    self.name, dims[0], dims[1], u, v, t, r1,
                                    u1, v1, ang_type, self.defines[define],
                                    seqs)

        return angle_stat

    def get_stats(self, d):
        '''
        Calls get_loop_stat/ get_bulge_angle_stats or get_stem_stats, depending on the element d.

        :returns: A 1- or 2 tuple of stats (2 in case of bulges. One for each direction)
        '''
        if d[0]=="s":
            return (self.get_stem_stats(d),)
        elif d[0] in "mi":
            return self.get_bulge_angle_stats(d)
        else:
            try:
                stat = self.get_loop_stat(d)
            except ValueError:
                if len(self.defines)==1:
                    return tuple() #A structure without any stem has no stats.
                else:
                    raise
            if d[0] == "f":
                stat.stat_type="5prime"
            elif d[0] == "t":
                stat.stat_type="3prime"
            return (stat,)

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
        stem1, = self.edges[d] # Make sure there is only one edge

        (s1b, s1e) = self.get_sides(stem1, d)

        stem1_vec = self.coords[stem1][s1b] - self.coords[stem1][s1e]
        twist1_vec = self.twists[stem1][s1b]
        bulge_vec = self.coords[d][1] - self.coords[d][0]

        if ftuv.magnitude(bulge_vec)<10**-3: #To avoid loops with 0 physical length. (If disconnects in the structure are modelled as loop)
            bulge_vec += (10**-3) * (stem1_vec / ftuv.magnitude(stem1_vec))

        (r,u,v) = ftug.get_stem_separation_parameters(stem1_vec, twist1_vec,
                                                      bulge_vec)
        (loop_stat.r, loop_stat.u, loop_stat.v) = (r, u, v)
        loop_stat.r = loop_stat.phys_length # Will this cause problems in other parts of the code base???
        loop_stat.define=self.defines[d]
        loop_stat.seqs = self.get_define_seq_str(d, adjacent=True)
        if d[0]=="f":
            loop_stat.stat_type = "5prime"
        elif d[0]=="t":
            loop_stat.stat_type = "3prime"
        return loop_stat

    def get_bulge_angle_stats(self, bulge):
        '''
        Return the angle stats for a particular bulge. These stats describe
        the relative orientation of the two stems that it connects.

        :param bulge: The name of the bulge.
        :param connections: The two stems that are connected by it.
        :return: The angle statistics in one direction and angle statistics in
                 the other direction
        '''
        if bulge == 'start':
            return (ftms.AngleStat(), ftms.AngleStat())

        connections = self.connections(bulge)

        angle_stat1 = self.get_bulge_angle_stats_core(bulge, connections)
        angle_stat2 = self.get_bulge_angle_stats_core(bulge, list(reversed(connections)))

        assert round(angle_stat1.get_angle(),5) == round(angle_stat2.get_angle(),5), ("{}!={}".format(angle_stat1.get_angle(), angle_stat2.get_angle()))
        return (angle_stat1, angle_stat2)

    def get_stacking_helices(self, method="Tyagi"):
        """
        Return all helices (longer stacking regions) as sets.

        Two stems and one bulge are in a stacking relation, if self.is_stacking(bulge) is true and the stems are connected to the bulge.
        Further more, a stem is in a stacking relation with itself.
        A helix is the transitive closure this stacking relation.

        :returns: A list of sets of element names.
        """
        helices=[]
        for d in self.defines:
            if d[0] in "mi" and self.is_stacking(d, method):
                s1, s2 = self.connections(d)
                helices.append(set([d, s1, s2]))
            if d[0]=="s":
                helices.append(set([d]))
        while True:
            for i,j in it.combinations(range(len(helices)),2):
                stack_bag1 = helices[i]
                stack_bag2 = helices[j]
                if stack_bag1 & stack_bag2:
                    stack_bag1|=stack_bag2
                    del helices[j]
                    break
            else:
                break
        return helices
    def is_stacking(self, bulge, method="Tyagi", verbose=False):
        """
        Reports, whether the stems connected by the given bulge are coaxially stacking.



        :param bulge: STRING. Name of a interior loop or multiloop (e.g. "m3")
        :param method": STRING. "Tyagi": Use cutoffs from doi:10.1261/rna.305307, PMCID: PMC1894924.
        :returns: A BOOLEAN.
        """

        assert method in ["Tyagi", "CG" ]
        if method=="Tyagi":
            return self._is_stacking_tyagi(bulge, verbose)
        return self._is_stacking_CG(bulge, verbose)

    def _is_stacking_CG(self, bulge, verbose=False):
        """"""
        stem1, stem2 = self.connections(bulge)
        angle = ftuv.vec_angle(self.coords[stem1][1]-self.coords[stem1][0],
                               self.coords[stem2][1]-self.coords[stem2][0])
        if angle>math.pi/2:
            angle=math.pi-angle
        if angle>math.radians(45):
            if verbose: print("Angle {}>45".format(math.degrees(angle)))
            return False
        shear_angle1 = ftuv.vec_angle(self.coords[stem1][1]-self.coords[stem1][0],
                               self.coords[bulge][1]-self.coords[bulge][0])
        if shear_angle1>math.pi/2:
            shear_angle1=math.pi-shear_angle1
        shear_angle2 = ftuv.vec_angle(self.coords[stem2][1]-self.coords[stem2][0],
                               self.coords[bulge][1]-self.coords[bulge][0])
        if shear_angle2>math.pi/2:
            shear_angle2=math.pi-shear_angle2
        if shear_angle1>math.radians(60) or shear_angle2>math.radians(60):
            if verbose: print ("Shear angle 1 {}>60 or shear angle 2 {}>60".format(math.degrees(shear_angle1), math.degrees(shear_angle1)))
            return False
        return True

    def _is_stacking_tyagi(self, bulge, verbose=False):
        """
        Implementation of the method described in doi:10.1261/rna.305307 (Tyagi and Matthews) for
        the detection of coaxial stacking.

        Called by self.is_stacking(bulge, "Tyagi")

        ..note::
            This does NOT implement the method for coaxial stacking prediction which is the main
            focus of the paper, only the method for the detection of stacking in pdb files.
        """
        assert bulge[0] in "mi"
        DISTANCE_CUTOFF = [ 14, 6 ]
        ANGLE_CUTOFF    = [  math.acos(0.75), math.acos(0.8) ]
        SHEAR_ANGLE_CUTOFF = math.radians(60) #Relaxed compared to 60 in the paper, because we use
                                              #virtual atom positions
        SHEAR_OFFSET_CUTOFF = 10
        if bulge[0]=="m" and self.get_length(bulge) == 0:
            is_flush = True #flush-stack vs. mismatch-mediated stack
        else:
            is_flush = False

        stem1, stem2 = self.connections(bulge)
        side_nts = self.get_connected_residues(stem1, stem2, bulge)[0]
        #Distance
        bp_center1 = ftug.get_basepair_center(self, side_nts[0])
        bp_center2 = ftug.get_basepair_center(self, side_nts[1])
        if ftuv.vec_distance(bp_center1, bp_center2)>DISTANCE_CUTOFF[is_flush]:
            if verbose: print ("Distance {} > {}".format(ftuv.vec_distance(bp_center1, bp_center2), DISTANCE_CUTOFF[is_flush]))
            return False
        normalvec1 = ftug.get_basepair_plane(self, side_nts[0])
        normalvec2 = ftug.get_basepair_plane(self, side_nts[1])
        #Coaxial
        angle = ftuv.vec_angle(normalvec1, normalvec2)
        if angle>math.pi/2:
            #Triggered frequently
            #warnings.warn("Angle > 90 degrees: {} ({})".format(angle, math.degrees(angle)))
            angle=math.pi-angle
        if angle>ANGLE_CUTOFF[is_flush]:
            if verbose: print ("Angle {} > {}".format(angle, ANGLE_CUTOFF[is_flush]))
            return False
        #Shear Angle
        shear_angle1 = ftuv.vec_angle(normalvec1, bp_center2-bp_center1)
        if shear_angle1>math.pi/2:
            shear_angle1=math.pi-shear_angle1
        if shear_angle1>SHEAR_ANGLE_CUTOFF:
            if verbose: print ("Shear angle 1 {} > {}".format(shear_angle1, SHEAR_ANGLE_CUTOFF))
            return False
        shear_angle2 = ftuv.vec_angle(normalvec2, bp_center1-bp_center2)
        if shear_angle2>math.pi/2:
            shear_angle2=math.pi-shear_angle2
        if shear_angle2>SHEAR_ANGLE_CUTOFF:
            if verbose: print ("Shear angle 2 {} > {}".format(shear_angle2, SHEAR_ANGLE_CUTOFF))
            return False
        #Shear Offset
        #Formula for distance between a point and a line
        #from http://onlinemschool.com/math/library/analytic_geometry/p_line/
        if (ftuv.magnitude(np.cross((bp_center1-bp_center2), normalvec2))/
                            ftuv.magnitude(normalvec2))>SHEAR_OFFSET_CUTOFF:
            if verbose: print ("Shear offset 1 wrong:", (ftuv.magnitude(np.cross((bp_center1-bp_center2), normalvec2))/
                                                         ftuv.magnitude(normalvec2)), ">" , SHEAR_OFFSET_CUTOFF)
            return False
        if (ftuv.magnitude(np.cross((bp_center1-bp_center2), normalvec1))/
                            ftuv.magnitude(normalvec1))>SHEAR_OFFSET_CUTOFF:
            if verbose: print ("Shear offset 2 wrong")
            return False
        return True

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
        # Stems may or may not have adjacent nucleotides, so we do not include the adj nucleotides
        ss.seqs = self.get_define_seq_str(stem, adjacent=False)
        return ss

    #def get_loop_from_residue(self, residue) ->  use BulgeGraph.get_node_from_residue_num()!
    def _init_coords(self):
        self.coords = LineSegmentStorage(self.defines.keys(), on_change = self.reset_vatom_cache)
        self.twists = CoordinateStorage([x for x in self.defines if x[0] =="s"], on_change = self.reset_vatom_cache)
    def from_bpseq_str(self, bpseq_str, dissolve_length_one_stems=False, breakpoints = []):
        super(CoarseGrainRNA, self).from_bpseq_str(bpseq_str, dissolve_length_one_stems, breakpoints)
        self._init_coords()
    def from_dotbracket(self,  dotbracket_str, dissolve_length_one_stems=False):
        super(CoarseGrainRNA, self).from_dotbracket(dotbracket_str, dissolve_length_one_stems)
        self._init_coords()
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
        self._init_coords()
        log.debug("BG read. Now reading 3D information")
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
        self.add_bulge_coords_from_stems() #Old versions of the file may contain bulge coordinates in the wrong order.
        self.incomplete_elements = ftug.get_incomplete_elements(self)


    def to_cg_file(self, filename):
        '''
        Save this structure as a string in a file.

        :param filename: The filename to save it to.
        '''
        warnings.warn("to_cg_file is deprecated. Use to_file!", stacklevel=2)
        with open(filename, 'w') as f:
            s = self.to_cg_string()

            f.write(s)

    def radius_of_gyration(self, method = "vres"):
        '''
        Calculate the radius of gyration of this structure.

        :param method: A STRING. one of
                       "fast" (use only coordinates of coarse grained stems) or
                       "vres" (use virtual residue coordinates of stems)

        :return: A number with the radius of gyration of this structure.
        '''
        if len(list(self.stem_iterator()))==0:
            log.warning("Cannnot calculate ROG (%s) for structure without stems", method)
            return float("nan")
        if method=="fast":
            coords=self.get_ordered_stem_poss()
        elif method=="vres":
            coords = self.get_ordered_virtual_residue_poss()
        else:
            raise ValueError("Wrong method {}. Choose one of 'fast' and 'vres'".format(method))

        rog = ftud.radius_of_gyration(coords)
        return rog

    def get_coordinates_list(self):
        warnings.warn("CoarseGrainRNA.get_coordinates_list is deprecated and being "
                      "replaced by get_coordinates_array!")
        return self.get_coordinates_array()

    def get_coordinates_array(self):
        '''
        Get all of the coordinates in one large array.

        The coordinates are sorted in the order of the keys
        in coordinates dictionary.

        :return: A 2D numpy array containing all coordinates
        '''
        all_coords = []
        assert len(self.coords) == len(self.defines), self.coords.keys()^self.defines.keys()
        for key in sorted(self.coords.keys()):
            for i in range(len(self.coords[key])):
                all_coords.append(self.coords[key][i])
        return np.array(all_coords)

    def load_coordinates_array(self, coords):
        '''
        Read in an array of coordinates (as may be produced by get_coordinates_array)
        and replace the coordinates of this structure with it.

        :param coords: A 2D array of coordinates
        :return: self
        '''

        for j, key in enumerate(sorted(self.coords.keys())):
            self.coords[key] = coords[ 2*j ], coords[ 2*j+1 ]
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
                                  self.coords[node][1] -
                                  self.coords[node][0]))

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
        assert False

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
        total_length = sum([len(list(self.define_residue_num_iterator(d))) for d in self.defines])
        assert total_length == self.seq_length
        return self.seq_length

    def sorted_edges_for_mst(self):
        """
        Keep track of all linked nodes. Used for the generation of the minimal spanning tree.

        This overrides the function in bulge graph and adds an additional sorting criterion
        with lowest priority.
        Elements that have no entry in self.sampled should be preferedly broken.
        This should ensure that the minimal spanning tree is the same after saving
        and loading an RNA to/from a file, if changes of the minimal spanning tree
        were performed by ernwin.
        """
        priority = {'s': 1, 'i': 2, 'm': 3, 'f': 4, 't': 5}
        edges = sorted(it.chain(self.mloop_iterator(),
                                self.iloop_iterator()),
                       key=lambda x: (priority[x[0]], min(self.get_node_dimensions(x)),not x in self.sampled))
        return edges

    def coords_to_directions(self):
        """
        The directions of each coarse grain element. One line per cg-element.

        The array is sorted by the corresponding element names alphabetically (`sorted(defines.keys()`)
        The directions usually point away from the elemnt's lowest nucleotide.
        However h,t and f elements always point away from the connected stem.
        """
        coords = self.get_coordinates_array()
        directions = coords[1::2]-coords[0::2]
        return directions
    def coords_from_directions(self, directions):
        """
        Generate coordinates from direction vectors (using also their lengths)

        Currently ignores the twists!

        :param directions: An array of vectors from the side of a cg-element with lower nucleotide number to the side with higher number
                           The array is sorted by the corresponding element names alphabetically (`sorted(defines.keys()`)

        """
        sorted_defines = sorted(self.defines.keys())
        assert len(sorted_defines)==len(directions), "{} != {}".format(len(sorted_defines), len(directions))
        if self.build_order is None:
            self.traverse_graph()
        self.coords["s0"]=np.array([0,0,0]), directions[sorted_defines.index("s0")]

        for stem1, link, stem2 in self.build_order: #Bulges and stems
            conn = self.connection_ends(self.connection_type(link, [stem1,stem2]))
            link_dir = self.get_link_direction(stem1, stem2, link)
            if link_dir==1:
                self.coords[link] = self.coords[stem1][conn[0]], self.coords[stem1][conn[0]]+directions[sorted_defines.index(link)]
                if conn[1]==0:
                    self.coords[stem2] = self.coords[link][1], self.coords[link][1]+directions[sorted_defines.index(stem2)]
                else:
                    self.coords[stem2] = self.coords[link][1] - directions[sorted_defines.index(stem2)], self.coords[link][1]
            else:
                self.coords[link] = self.coords[stem1][conn[0]] - directions[sorted_defines.index(link)], self.coords[stem1][conn[0]]
                if conn[1]==0:
                    self.coords[stem2] = self.coords[link][0], self.coords[link][0]+directions[sorted_defines.index(stem2)]
                else:
                    self.coords[stem2] = self.coords[link][0] - directions[sorted_defines.index(stem2)], self.coords[link][0]
        for d in self.defines:
            if d[0] == "m" and d not in self.mst:
                edges = list(self.edges[d])
                (s1b, _) = self.get_sides(edges[0], d)
                (s2b, _) = self.get_sides(edges[1], d)
                mids1 = self.coords[edges[0]]
                mids2 = self.coords[edges[1]]
                #Save coordinates in direction of the strand.
                if self.get_link_direction(edges[0], edges[1],d)==1:
                    self.coords[d] = (mids1[s1b], mids2[s2b])
                else:
                    self.coords[d] = (mids2[s2b], mids1[s1b])

            if d[0] in "hft": #Loops
                stem, = self.edges[d]
                (s1b, _) = self.get_sides(stem, d)
                self.coords[d] = self.coords[stem][s1b], self.coords[stem][s1b] + directions[sorted_defines.index(d)]

    def virtual_atoms(self, key):
        """
        Get virtual atoms for a key.

        :param key: An INTEGER: The number of the base in the RNA.

        :returns: A dict {atom:coords}, e.g. {"C8":np.array([x,y,z]), ...}
        """
        if isinstance(key, int):
             if key not in self._virtual_atom_cache:
                try:
                    self._virtual_atom_cache[key]=ftug.virtual_atoms(self)[key]
                except KeyError as e:
                    log.info("Key {} not present. Need to recreate all virtual residues".format(e))
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

    def rotate(self, angle, axis="x", unit="radians"):
        if unit=="degrees":
            angle=math.radians(angle)
        elif unit!="radians":
            raise ValueError("Unit {} not understood. Use 'degrees' or 'radians'".format(unit))
        s = math.sin(angle)
        cosi = math.cos(angle)
        rotation_matrix = np.zeros((3,3))
        if axis=="x":
            rotation_matrix[0,0]=1
            rotation_matrix[1,1]=rotation_matrix[2,2]=cosi
            rotation_matrix[1,2]=-s
            rotation_matrix[2,1]=s
        elif axis=="y":
            rotation_matrix[0,0]=1
            rotation_matrix[0,0]=rotation_matrix[2,2]=cosi
            rotation_matrix[0,2]=-s
            rotation_matrix[2,0]=s
        elif axis=="z":
            rotation_matrix[2,2]=1
            rotation_matrix[0,0]=rotation_matrix[1,1]=cosi
            rotation_matrix[0,1]=-s
            rotation_matrix[1,0]=s
        self.coords.rotate(rotation_matrix)
        self.twists.rotate(rotation_matrix)

        #Caching for virtual residues
        self.vposs = c.defaultdict( dict )
        self.vbases = c.defaultdict( dict )
        self.vvecs = c.defaultdict( dict )
        self.v3dposs = c.defaultdict( dict )
        self.vinvs = c.defaultdict( dict )
