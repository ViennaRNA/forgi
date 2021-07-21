from __future__ import absolute_import, unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

import collections as c
import contextlib
import os
import os.path as op
import shutil
import subprocess as sp
import sys
import time
import math
import warnings
import itertools as it
import json
try:
    from json import JSONDecodeError #py3k
except ImportError: #py2k
    JSONDecodeError = ValueError
try:
    import io
except ImportError:  # Python 2
    import StringIO as io
try:
    import shutil
    which = shutil.which
except (ImportError, AttributeError):
    import distutils.spawn
    which = distutils.spawn.find_executable

from io import open
import logging
from pprint import pprint


import numpy as np

import scipy.spatial
import scipy.stats

import Bio.PDB as bpdb

from logging_exceptions import log_to_exception, log_exception

from ...graph import bulge_graph as fgb
import forgi.graph.residue as fgr
from ...graph.sequence import Sequence, _insert_breakpoints_simple
from ...graph._graph_construction import _BulgeGraphConstruction
from ..utilities import graph_pdb as ftug
from . import stats as ftms
from . import transform_cg as ftmt
from ..._k2n_standalone import knotted2nested as cak
from ..utilities import mcannotate as ftum
from ..utilities import pdb as ftup
from . import descriptors as ftmd
from ..utilities import _dssr as ftud
from ..utilities import vector as ftuv
from ...utilities import debug as fud
from ...utilities import stuff as fus
import forgi.threedee.utilities.virtual_residues as ftuvres
from ...utilities.observedDict import observedDict
from ...utilities.exceptions import CgConstructionError, CgIntegrityError, GraphConstructionError
from .linecloud import CoordinateStorage, LineSegmentStorage
from ...utilities.stuff import is_string_type
from ... import config

log = logging.getLogger(__name__)

class AnnotationToolNotInstalled(ValueError):
    pass


try:
    profile  # The @profile decorator from line_profiler (kernprof)
except:
    def profile(x):
        return x


def add_longrange_interactions(cg, lines):
    '''
    Iterate over the lines in an MC-Annotate file and add information
    about interactions between non-adjacent elements.

    :param cg: A CoarseGrainRNA structure
    :param lines: All the lines in an MC-Annotate file
    '''
    for line in ftum.iterate_over_interactions(lines):
        (from_chain, from_base, to_chain,
         to_base) = ftum.get_interacting_base_pairs(line)
        try:
            seq_id1 = cg.seq_ids.index(
                fgr.RESID(from_chain, ftum.parse_resid(from_base))) + 1
        except ValueError as e:
            with log_to_exception(log, e):
                log.error("seq_ids are {}".format(cg.seq_ids))
            raise
        seq_id2 = cg.seq_ids.index(
            fgr.RESID(to_chain, ftum.parse_resid(to_base))) + 1

        node1 = cg.get_node_from_residue_num(seq_id1)
        node2 = cg.get_node_from_residue_num(seq_id2)

        if abs(seq_id2 - seq_id1) > 1 and node1 != node2 and not cg.has_connection(node1, node2):
            cg.longrange[node1].add(node2)
            cg.longrange[node2].add(node1)


def breakpoints_from_seqids(seqids):
    """
    Create the list of cofold cutpoints from the seq_ids.
    Return 1-based breakpoints
    """
    breakpoints = []
    old_chain = None
    for i, seqid in enumerate(seqids):
        if old_chain is not None and seqid.chain != old_chain:
            assert i >= 1
            breakpoints.append(i)  # Breakpoint BEFORE current seqid!
        old_chain = seqid.chain
    return breakpoints


def _are_adjacent_basepairs(seq_ids, edge1, edge2):
    """
    Helper function used by connected_cgs_from_pdb.

    :param cg: A (potentially) inconsistent cg. Only seq_ids are needed.
    :param edge1, edge2: A dictionary with the key "basepair" mapping to a two-tuple of RESIDS
    """
    fromA, toA = edge1["basepair"]
    fromB, toB = edge2["basepair"]
    if fromA.chain != fromB.chain:
        if fromA.chain == toB.chain:
            fromB, toB = toB, fromB
            assert toA.chain == toB.chain
        else:
            assert False
    if (abs(seq_ids.index(fromA) - seq_ids.index(fromB)) == 1 and
            abs(seq_ids.index(toA) - seq_ids.index(toB)) == 1):
        log.debug(
            "Basepairs %s and %s are adjacent. No length 1 stem", edge1, edge2)
        return True
    log.debug("Basepairs %s and %s are NOT adjacent.", edge1, edge2)
    return False


def _annotate_pdb(filename, annotation_tool, filetype, subprocess_kwargs={}):
    """
    Get the secondary structure of the pdb by using an external tool.

    :param filename: The name of the temporary pdb file
    :param subprocess_kwargs: Will be passed as keyword arguments to subprocess.call/ subprcess.check_output
    :param annotation_tool: See docstring of CoarseGrainRNA.from_pdb
    """
    not_installed_msg = ("{} was requested as annotation tool, but {} "
                         "is not installed or not in the PATH. (Hint: run "
                         "forgi_config.py to set the preferred annotation tool)")
    program = annotation_tool
    if program is None:
        c = config.read_config()
        if "PDB_ANNOTATION_TOOL" in c:
            program = c["PDB_ANNOTATION_TOOL"]
    log.info("Requested annotation program is %s", program)
    if program == "DSSR" or program is None:
        if which("x3dna-dssr"):
            return _run_dssr(filename, subprocess_kwargs)
        else:
            log.info("x3dna-dssr is not installed or not in the PATH.")
            if program is not None:
                raise AnnotationToolNotInstalled(
                    not_installed_msg.format(program, "x3dna-dssr"))
    if program == "MC-Annotate" or (program is None and filetype == "pdb"):
        if filetype != "pdb":
            raise ValueError("MC-Annotate does not support MMCIF")
        if which("MC-Annotate"):
            lines = _run_mc_annotate(filename, subprocess_kwargs)
            try:
                bpseq, seq_ids = ftum.get_dotplot(lines)
                return bpseq, seq_ids, {}
            except Exception as e:
                log.exception(
                    "Could not convert MC-Annotate output to dotplot")
                raise CgConstructionError(
                    "Could not convert MC-Annotate output to dotplot")  # from e
        else:
            log.info("MC-Annotate is not installed or not in the PATH.")
            if program is not None:
                raise AnnotationToolNotInstalled(not_installed_msg.format(
                    program, "MC-Annotate"))
    elif program == "forgi" or program is None:
        return None
    else:
        raise ValueError("Supported programs for annotating the pdb are: "
                         "'MC-Annotate', 'DSSR' and 'forgi', not '{}'".format(program))


def _run_dssr(filename, subprocess_kwargs={}):
    # See http://forum.x3dna.org/rna-structures/redirect-auxiliary-file-output/
    with fus.make_temp_directory() as dssr_output_dir:
        dssr_out = os.path.join(dssr_output_dir, "out")

        arguments = ['x3dna-dssr', "-i=" + filename,
                     "--prefix=" + os.path.join(dssr_output_dir, "d"),
                     "-o=" + dssr_out, "--json"]
        # https://stackoverflow.com/a/14837250/5069869
        log.info("Running DSSR: %s", sp.list2cmdline(arguments))
        with open(os.path.join(dssr_output_dir, "stderror"), "w+") as errfile:
            try:
                ret_code = sp.call(
                    arguments, universal_newlines=True, stderr=errfile, **subprocess_kwargs)
                try:
                    with open(os.path.join(dssr_output_dir, "d-2ndstrs.bpseq"), encoding='ascii') as f:
                        bpseq = f.read()
                    with open(dssr_out) as f:
                        try:
                            dssr_dict = json.load(f)
                        except JSONDecodeError as e: # DSSR does not escape backslashes (see "http://forum.x3dna.org/bug-reports/json-output-should-escape-backslashes/"), which can be in windows paths. Escape them manually.
                            f.seek(0)
                            json_stri = f.read().replace('\\','\\\\')
                            log.error(json_stri)
                            dssr_dict = json.loads(json_stri)
                        nts = dssr_dict["nts"]
                        seq_ids = list(map(ftud.dssr_to_pdb_resid, [
                                       nt["nt_id"] for nt in nts]))
                except (OSError, IOError) as e:
                    with log_to_exception(log, e):
                        log.error("Content of the directory %s is %s",
                                  dssr_output_dir, os.listdir(dssr_output_dir))
                    raise
            except (OSError, IOError) as e:
                assert op.isfile(
                    filename), "File {} (created by forgi) no longer exists".format(filename)
                errfile.seek(0)
                err_msg = errfile.readlines()
                if err_msg:
                    if len(err_msg) >= 3:
                        e = CgConstructionError(
                            "DSSR could not process the file: " + err_msg[-3])
                    with log_to_exception(log, e):
                        log.error("Captured Stderr is:\n%s", "".join(
                            ["... " + line for line in err_msg]))
                    raise e
                else:
                    # On Py2.7, the filename of the executable (x3dna-dssr) is not part of the
                    # error message raised by sp.call,
                    # For this reason, we add the following hint (unneeded on Py3k)
                    e.strerror += ". (Hint: Did you install x3dna-dssr?)"
                raise e
    return bpseq, seq_ids, dssr_dict


def _run_mc_annotate(filename, subprocess_kwargs={}):
    """
    Returns: A tuple (bpseq, seq_ids`)
    """
    log.info("Running MC-Annotate")
    try:
        out = sp.check_output(['MC-Annotate', filename],
                              universal_newlines=True, **subprocess_kwargs)
    except OSError as e:
        assert op.isfile(
            filename), "File {} (created by forgi) no longer exists".format(filename)
        e.strerror += ". Hint: Did you install MC-Annotate?"
        raise e

    lines = out.strip().split('\n')
    return lines


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

    def __init__(self, graph_construction, sequence, name=None, infos=None, _dont_split=False):
        '''
        Initialize the new structure.
        '''
        super(CoarseGrainRNA, self).__init__(
            graph_construction, sequence, name, infos, _dont_split)

        self.dssr = None

        self._virtual_atom_cache = {}
        #: Keys are element identifiers (e.g.: "s1" or "i3"), values are 2-tuples of vectors
        #: The first value of stem coordinates corresponds to the start of the stem
        #: (the one with the lowest nucleotide number),
        #: The second value to the end of the stem.
        #: If the coordinates for an element change, the virtual atom and virtual residue
        #: coordinates are automatically invalidated.
        self.coords = None  # We can only initialize this, when we know defines.keys()
        self.twists = None
        self.sampled = dict()

        if self.defines:
            self._init_coords()

        #: A dictionary storing the stem-bases (not the vres basis)
        self.bases = {}
        self.stem_invs = {}


        #: generated by self.add_all_virtual_residues()
        self.vposs = c.defaultdict(dict)
        self.vbase = c.defaultdict(dict)
        self.vsugar = c.defaultdict(dict)
        self.vbackbone = c.defaultdict(dict)
        #:The following defaultdicts are cleared when coords or twists change.
        #: The coordinate system specific to each virtual residue
        #: (3x3 matrix, carthesian coordiantes)
        #: generated by self.add_all_virtual_residues()
        self.vbases = c.defaultdict(dict)
        #: generated by self.add_all_virtual_residues()
        self.vvecs = c.defaultdict(dict)
        #: generated by self.add_all_virtual_residues()
        self.v3dposs = c.defaultdict(dict)
        #: generated by self.add_all_virtual_residues()
        #: The inverse of the transposed basis.
        self.vinvs = c.defaultdict(dict)

        #: A 3D vector. Used as hint, from what direction the Projection2D object
        #: should be generated in the default case.
        self.project_from = None

        self.longrange = c.defaultdict(set)
        self.chains = {}  # the PDB chains if loaded from a PDB file

        self.interacting_residues = []

        # Lazily loaded:
        self._incomplete_elements = None
        # Lazily calculated from interacting_residues
        self._interacting_elements = None
    ############################################################################
    # Factory functions
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    @classmethod
    def from_bg_string(cls, cg_string):
        '''
        Populate this structure from the string
        representation of a graph.
        '''
        # pylint: disable=no-member  # See https://github.com/PyCQA/pylint/issues/981
        # Reading the bulge_graph-part of the file
        cg = super(CoarseGrainRNA, cls).from_bg_string(cg_string)
        cg._init_coords()
        log.debug("BG read. Now reading 3D information")
        # Reading the part of the file responsible for 3D information
        lines = cg_string.split('\n')
        for line in lines:
            line = line.strip()
            parts = line.split()
            if len(parts) == 0:
                continue
            if parts[0] == 'coord':
                name = parts[1]
                cg.coords[name] = np.array(
                    [list(map(float, parts[2:5])), list(map(float, parts[5:8]))])
            if parts[0] == 'twist':
                name = parts[1]
                cg.twists[name] = np.array(
                    [list(map(float, parts[2:5])), list(map(float, parts[5:8]))])
            if parts[0] == 'longrange':
                cg.longrange[parts[1]].add(parts[2])
                cg.longrange[parts[2]].add(parts[1])

            if parts[0] == 'sampled':
                cg.sampled[parts[1]] = [parts[2]] + list(map(int, parts[3:]))
            if parts[0] == 'project':
                cg.project_from = np.array(parts[1:], dtype=float)
            if parts[0] == "interacting":
                cg.interacting_residues.append(fgr.resid_from_str(parts[1]))
            if parts[0] == "vres":
                elem = parts[1]
                cg.vposs[elem] = ftuvres.parse_vres(parts[2:])[0]
            if parts[0] == "vbase":
                elem = parts[1]
                cg.vbase[elem] = ftuvres.parse_vres(parts[2:])[0]
            if parts[0] == "vsugar":
                elem = parts[1]
                cg.vsugar[elem] = ftuvres.parse_vres(parts[2:])[0]
            if parts[0] == "vbackbone":
                elem = parts[1]
                cg.vbackbone[elem] = ftuvres.parse_vres(parts[2:])[0]
        # Old versions of the file may contain bulge coordinates in the wrong order.
        cg.add_bulge_coords_from_stems()
        return cg

    @classmethod
    def from_pdb(cls, pdb_filename, load_chains=None, remove_pseudoknots=False,
                 dissolve_length_one_stems=True, secondary_structure=None,
                 filetype="pdb", annotation_tool=None, query_PDBeChem=False):
        """
        :param load_chains: A list of chain_ids or None (all chains)
        :param secondary_structure: Only useful if we load only 1 component
        :param filetype: One of 'pdb' or 'cif'
        :param query_PDBeChem: If true, query the PDBeChem database whenever a
                        modified residue with unknown 3-letter code
                        is encountered.
        :param annotation_tool: One of "DSSR", "MC-Annotate", "forgi" or None.
                        If this is None, we take the value of the configuration
                        file (run forgi_config.py to create a config file).
                        If no config file is given either, we see what tools are
                        installed, preferring the newer DSSR over MC-Annotate and
                        falling back to the fogi implementation if neither is in
                        the PATH variable.
                        If a string is given or the configuration file set,
                        we never fall back to a different option but raise an
                        error, if the requested tool is unavailable.
        """
        #warnings.warn("We currently do not load any long-range interactions")
        # We need to create files, so we can interface with
        # the annotation program.
        if remove_pseudoknots and secondary_structure:
            warnings.warn(
                "Option 'remove_pseudoknots ignored for user-supplied secondary structure.")

        with fus.make_temp_directory() as output_dir:
            if load_chains == "biggest":
                chain, missing_res, ir = ftup.get_biggest_chain(pdb_filename)
                chains = [chain]
            else:
                chains, missing_res, ir = ftup.get_all_chains(pdb_filename)
            new_chains = []
            for chain in chains:
                if load_chains in [None, "biggest"] or chain.id in load_chains:
                    log.debug("Loaded Chain %s", chain.id)
                    chain, modifications = ftup.clean_chain(chain, query_PDBeChem)
                    new_chains.append(chain)
            log.debug("%s, %s", pdb_filename, os.path.split(pdb_filename))
            fn_basename = os.path.split(pdb_filename)[1]
            if load_chains is None:
                chain_string = "None"
            else:
                chain_string = "-".join(map(str, load_chains))

            rna_pdb_fn = op.join(output_dir, fn_basename +
                                 "_" + chain_string + '.temp.' + filetype)
            # Output cleaned version for annotating
            ftup.output_multiple_chains(new_chains, rna_pdb_fn, filetype)

            # first we annotate the 3D structure
            log.info("Starting annotation program for all chains")
            annotation = _annotate_pdb(rna_pdb_fn, annotation_tool, filetype)
            if annotation is None:
                # Fallback-annotation using forgi
                bpseq, seq_ids = ftup.annotate_fallback(new_chains)
                dssr_dict = {}
            else:
                bpseq, seq_ids, dssr_dict = annotation
            breakpoints = breakpoints_from_seqids(seq_ids)

        assert bpseq is not None, "{} {} {}".format(bpseq, seq_ids, dssr_dict)
        # Here, we remove Pseudoknots on the bp-seq level. not the BulgeGraph
        # level, because there may be multiple components, so we cannot create
        # a BulgeGraph.
        if remove_pseudoknots:
            #log.info("Removing pseudoknots")
            #log.info(bpseq)
            bpseq = cak.k2n_main(io.StringIO(bpseq), input_format='bpseq',
                                 #output_format = 'vienna',
                                 output_format='bpseq',
                                 method=cak.DEFAULT_METHOD,
                                 opt_method=cak.DEFAULT_OPT_METHOD,
                                 verbose=cak.DEFAULT_VERBOSE,
                                 removed=cak.DEFAULT_REMOVED)

        import networkx as nx
        chain_connections_multigraph = nx.MultiGraph()
        bpseq_lines = bpseq.splitlines()
        for line in bpseq_lines:
            try:
                from_, res, to_ = line.split()
            except:
                continue
            from_seqid = seq_ids[int(from_) - 1]
            from_chain = from_seqid.chain
            chain_connections_multigraph.add_node(from_chain)
            if int(to_) != 0 and int(to_) > int(from_):
                to_seqid = seq_ids[int(to_) - 1]
                to_chain = to_seqid.chain
                if from_chain != to_chain:
                    log.debug("Adding {} - {}".format(from_chain, to_chain))
                    chain_connections_multigraph.add_edge(
                        from_chain, to_chain, basepair=(from_seqid, to_seqid))
        chain_connections = nx.Graph(chain_connections_multigraph)
        if dissolve_length_one_stems:
            for chain1, chain2 in it.combinations(chain_connections_multigraph.nodes(), 2):
                if chain2 in chain_connections_multigraph.adj[chain1]:
                    for edge1, edge2 in it.combinations(chain_connections_multigraph.adj[chain1][chain2].values(), 2):
                        if _are_adjacent_basepairs(seq_ids, edge1, edge2):
                            break
                    else:  # break NOT encountered.
                        log.debug(
                            "No connection remaining between chains %s and %s", chain1, chain2)
                        chain_connections.remove_edge(chain1, chain2)

        log.debug("CONNECTIONS: {}, nodes {}".format(
            chain_connections, chain_connections.nodes()))
        log.debug("Edges {}".format(chain_connections.edges()))

        pdb_base = op.splitext(op.basename(pdb_filename))[0]

        components = list(nx.connected_components(chain_connections))
        cgs = []
        if len(components) != 1 and secondary_structure:
            raise GraphConstructionError("User-provided secondary structure is "
                                         "ambigouse if more than 1 connected "
                                         "component loaded.")
        for component in nx.connected_components(chain_connections):
            try:
                cgs.append(cls._load_pdb_component(bpseq_lines, pdb_base, new_chains,
                                                   component, missing_res, modifications,
                                                   seq_ids, secondary_structure,
                                                   dissolve_length_one_stems,
                                                   ir))
            except GraphConstructionError as e:
                log_exception(e, logging.ERROR, with_stacktrace=False)
                log.error(
                    "Could not load chains %s, due to the above mentioned error.", list(component))
                raise

        cgs.sort(key=lambda x: x.name)
        if dssr_dict:
            for cg in cgs:
                cg.dssr = ftud.DSSRAnnotation(dssr_dict, cg)
                cg.infos["dssr_stacks"] = cg.dssr.stacking_loops()
        log.debug("Returning %s", cgs)
        return cgs

    @classmethod
    def _load_pdb_component(cls, original_bpseq_lines, name, chains, chain_ids,
                            missing_res, modifications, seq_ids,
                            secondary_structure="", dissolve_length_one_stems=False,
                            interacting_residues=[]):
        """
        :param original_bpseq_lines: List of strings. Will be filtered for chains.
        """
        #print(component, type(component))
        log.info("Loading PDB: Connected component with chains %s",
                 str(list(chain_ids)))
        log.debug("missing residues %s", missing_res)
        # Since the external annotation program can take some time,
        # we do not re-annotate, but filter the bpseq_string instead.
        new_bpseq = []
        new_seqids = []

        # Filter the bpseq_lines
        for line in original_bpseq_lines:
            from_, res, to_ = line.split()
            from_seqid = seq_ids[int(from_) - 1]
            if from_seqid.chain in chain_ids:
                new_seqids.append(from_seqid)
                if to_ != '0':
                    to_seqid = seq_ids[int(to_) - 1]
                    if to_seqid.chain not in chain_ids:  # Because of length-1-basepair
                        to_ = 0
                new_bpseq.append((from_, res, to_))

        new_bpseq_str = fus.renumber_bpseq(new_bpseq)
        tuples, seq_str = fus.bpseq_to_tuples_and_seq(new_bpseq_str)
        breakpoints = breakpoints_from_seqids(new_seqids)
        seq_str = _insert_breakpoints_simple(seq_str, breakpoints, 1)

        if secondary_structure:
            pt = fus.dotbracket_to_pairtable(secondary_structure)
            # Overwrite tuples!
            tuples = fus.pairtable_to_tuples(pt)
            stru_fragments = secondary_structure.split("&")
            seq_fragments = seq_str.split("&")
            seq_lengths = list(map(len, seq_fragments))
            if list(map(len, stru_fragments)) != seq_lengths:
                raise ValueError("The given secondary structure is inconsistent "
                                 "with the pdb-chain-lengths of %s", seq_lengths)
        sequence = Sequence(seq_str, new_seqids,
                            [r for r in missing_res if r["chain"] in chain_ids],
                            {k: v for k, v in modifications.items() if k.chain in chain_ids})

        name = name + "_" + "-".join(c for c in sorted(chain_ids))
        graph_constr = _BulgeGraphConstruction(tuples)
        cg = cls(graph_constr, sequence, name)
        cg = fgb._cleaned_bg(cg, dissolve_length_one_stems)
        cg.chains = {
            chain.id: chain for chain in chains if chain.id in chain_ids}

        ftug.add_stem_information_from_pdb_chains(cg)
        cg.add_bulge_coords_from_stems()
        ftug.add_loop_information_from_pdb_chains(cg)
        ftug._add_loop_vres(cg)
        assert len(cg.defines) == len(
            cg.coords), cg.defines.keys() ^ cg.coords.keys()
        cg.interacting_residues = list(r for r in map(fgr.resid_from_biopython, interacting_residues)
                                       if r.chain in chain_ids)
        return cg

    ############################################################################

    @property
    def transformed(self):
        return ftmt.CGTransformer(self)

    def _get_coord_str(self):
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
            if key[0] in ["m", "i"]:
                # Bulge coordinates are redundant. They can be deduced from the stem coordinates.
                continue
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

                    # Save coordinates in direction of the strand.
                    if self.get_link_direction(edges[0], edges[1], d) == 1:
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
           residue positions for loops usually do not make sense. If the file
           was loaded from the PDB, residue positions from the PDB file are
           stored already.
        """
        for stem in self.stem_iterator():
            try:
                log.debug(
                    "Adding virtual residues for stem %s with coords %s", stem, self.coords[stem])
                ftug.add_virtual_residues(self, stem)
            except (KeyError, ValueError, AssertionError):
                if np.all(np.isnan(self.coords[stem])):
                    raise RnaMissing3dError(
                        "No 3D coordinates available for stem {}".format(stem))
                elif stem[0] == "s" and np.all(np.isnan(self.twists[stem])):
                    raise RnaMissing3dError(
                        "No twists available for stem {}".format(stem))
                else:
                    log.info("Reraising an ERROR in add_all_virtual_residues")
                    raise

    def iter_three_points(self, pos):
      elem = self.get_node_from_residue_num(pos)
      d = self.defines[elem]
      if d[0] <= pos <= d[1]:
          i = pos - d[0]
      else:
          i = d[1] - d[0] + 1 + pos - d[2]
      for coord_type in self.vbase, self.vsugar, self.vbackbone:
          elem_coords = coord_type[elem][i]
          origin, basis = ftug.element_coord_system(self, elem)
          new_coords = ftuv.change_basis(
              elem_coords, ftuv.standard_basis, basis) + origin
          log.debug("%s %s: coords %s mapped to %s", elem, i, elem_coords, new_coords)
          yield new_coords

    def get_virtual_residue(self, pos, allow_single_stranded=False):
        """
        Get the virtual residue position in the global coordinate system
        for the nucleotide at position pos (1-based)

        :param pos: A 1-based nucleotide number
        :param allow_single_stranded: If True and pos is not in a stem, return a
              rough estimate for the residue position instead of raising an error.
              Currenly, for non-stem elements, these positions are on the axis of the cg-element.
        """
        mult = 5
        if isinstance(pos, fgr.RESID):
            pos = self.seq.to_integer(pos)
        elem = self.get_node_from_residue_num(pos)
        if elem[0] == "s" and elem not in self.v3dposs or not self.v3dposs[elem]:
            ftug.add_virtual_residues(self, elem)
        if elem[0] == "s":
            for i in range(self.stem_length(elem)):
                if self.defines[elem][0] + i == pos:
                    vres = self.v3dposs[elem][i]
                    return vres[0] + vres[3] * mult
                elif self.defines[elem][3] - i == pos:
                    vres = self.v3dposs[elem][i]
                    return vres[0] + vres[2] * mult
            else:
                assert False
        else:
            d = self.defines[elem]
            if d[0] <= pos <= d[1]:
                i = pos - d[0]
            else:
                i = d[1] - d[0] + 1 + pos - d[2]
            try:
                elem_coords = self.vposs[elem][i]
                origin, basis = ftug.element_coord_system(self, elem)
                new_coords = ftuv.change_basis(
                    elem_coords, ftuv.standard_basis, basis) + origin
                log.debug("%s %s: coords %s mapped to %s", elem, i, elem_coords, new_coords)
                return new_coords
            except KeyError as e:
                try:
                    self._has_warned_old_vres
                except AttributeError:
                    self._has_warned_old_vres = set()
                if elem not in self._has_warned_old_vres:
                    log.warning(
                        "No virtual residues have been loaded for loops: %s (for pos %s) in elem %s."
                        "Using inaccurate position along the cylinder instead.", e, pos, elem)
                    self._has_warned_old_vres.add(elem)
                if not allow_single_stranded:
                    raise ValueError(
                        "Position {} is not in a stem! It is in {}.".format(pos, elem))
                if elem[0] == "h":
                    # We estimate the vres position along the axis of the hairpin.
                    h_length = self.element_length(elem) / 2
                    pos_in_h = pos - self.defines[elem][0] +1
                    if pos_in_h > math.ceil(h_length):
                        l = self.defines[elem][1] - pos+1
                    elif pos_in_h <=int(h_length):
                        l = pos - self.defines[elem][0]+1
                    else:
                        l = math.ceil(h_length)
                    perc = (l) / math.ceil(h_length)
                elif elem[0] == "i":
                    if pos <= self.defines[elem][1]:
                        l = pos - self.defines[elem][0]
                        tl = (self.defines[elem][1] - self.defines[elem][0])
                    else:
                        l = self.defines[elem][3] - pos
                        tl = (self.defines[elem][3] -
                              self.defines[elem][2] )
                    perc = (l+1) / (tl+2)
                else:
                    l = pos - self.defines[elem][0]
                    perc = (l+1) / (self.element_length(elem)+1)
                log.debug("Calculating vres: %s,%s, %s: length=%s, perc=%s",elem, pos, l,(self.coords[elem][1] - self.coords[elem][0]), perc)
                return self.coords[elem][0] + (self.coords[elem][1] - self.coords[elem][0]) * perc

    def get_ordered_stem_poss(self):
        points = []
        for s in self.sorted_stem_iterator():
            points += list(self.coords[s])
        points = np.array(points)
        if np.any(np.isnan(points)):
            raise RnaMissing3dError("No 3D coordinates present.")
        return points

    def get_ordered_virtual_residue_poss(self, return_elements=False):
        """
        Get the coordinates of all stem's virtual residues in a consistent order.

        This is used for RMSD calculation.
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
        for i in range(1, len(self.seq)+1):
            pos = self.get_virtual_residue(i, allow_single_stranded = True)
            vress.append(pos)
            elems.append(self.get_elem(i))
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
        points = []
        if mode == "fast":
            for e in sorted(elements):
                points += list(self.coords[e])
            return np.array(points)
        elif mode == "vres":
            if not self.v3dposs:
                self.add_all_virtual_residues()
            vress = []
            for s in sorted(elements):
                if s[0] != "s":
                    continue
                for i in range(self.stem_length(s)):
                    vres = self.v3dposs[s][i]
                    vress += [vres[0] + vres[2], vres[0] + vres[3]]
            return np.array(vress)
        assert False

    def steric_value(self, elements, method="r**-2"):
        """
        Estimate, how difficult a set of elements was to build,
        by counting the atom density around the center of these elements
        """
        try:
            if isinstance(elements, list) or isinstance(elements, tuple):
                center = ftuv.get_vector_centroid(self.coords[elements])
            elif elements.shape == (3,):
                center = elements
            else:
                assert False, repr(elements)
        except:
            print(elements, repr(elements))
            raise
        if method == "kde":
            # print(center)
            all_vas = []
            for pos in range(1, self.seq_length + 1):
                for va in self.virtual_atoms(pos).values():
                    all_vas.append(va)

            all_vas = np.array(all_vas).T
            log.debug("Shape of all atoms {}".format(all_vas.shape))
            # randomly take 50 Angstrom bandwidth
            kde = scipy.stats.gaussian_kde(all_vas, 50)
            return kde(center)
        elif method.startswith("r**"):
            power = -int(method[3:5])
            exclude = method[5:]
            if exclude and exclude != "e":
                raise ValueError("Not supported method")
            value = 0
            for pos in range(1, self.seq_length + 1):
                if exclude and self.get_node_from_residue_num(pos) in elements:
                    continue
                for va in self.virtual_atoms(pos).values():
                    value += 1 / (1 + ftuv.vec_distance(va, center))**power
            return value
        elif method.startswith("cutoff"):
            cutoff = float(method.split()[1])
            value = 0
            for pos in range(1, self.seq_length + 1):
                for va in self.virtual_atoms(pos).values():
                    if ftuv.vec_distance(va, center) < cutoff:
                        value += 1
            return value

    def _get_twist_str(self):
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

    def _get_long_range_str(self):

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

    def _get_interacting_str(self):
        out = []
        for seqid in self.interacting_residues:
            out.append("interacting\t{}".format(fgr.resid_to_str(seqid)))
        return "\n".join(out) + "\n"

    def _get_sampled_stems_str(self):
        out_str = ''
        for key in self.sampled.keys():
            out_str += 'sampled %s %s\n' % (key,
                                            " ".join(map(str, self.sampled[key])))
        return out_str

    def _get_vres_str(self):
        out_str = ""
        for elem in self.vposs:
            if elem[0] != "s":
                out_str += "vres {} ".format(elem)
                out_str += ftuvres.serialize_vres(self.vposs[elem])
                out_str += "\n"
        for elem in self.vbase: 
            out_str += "vbase {} ".format(elem)
            out_str += ftuvres.serialize_vres(self.vbase[elem])
            out_str += "\n"
            out_str += "vsugar {} ".format(elem)
            out_str += ftuvres.serialize_vres(self.vsugar[elem])
            out_str += "\n"
            out_str += "vbackbone {} ".format(elem)
            out_str += ftuvres.serialize_vres(self.vbackbone[elem])
            out_str += "\n"
        return out_str

    def to_cg_string(self):
        '''
        Output this structure in string form.
        '''
        curr_str = self.to_bg_string()
        curr_str += self._get_coord_str()
        curr_str += self._get_twist_str()
        curr_str += self._get_sampled_stems_str()
        curr_str += self._get_long_range_str()
        curr_str += self._get_interacting_str()
        curr_str += self._get_vres_str()
        if self.project_from is not None:
            curr_str += "project {} {} {}\n".format(*self.project_from)
        return curr_str

    def to_file(self, filename):
        with open(filename, 'w') as f:
            cg_str = self.to_cg_string()
            f.write(cg_str)

    def get_bulge_angle_stats_core(self, elem, forward=True):
        '''
        Return the angle stats for a particular bulge. These stats describe the
        relative orientation of the two stems that it connects.

        :param elem: The name of the bulge.
        :param connections: The two stems that are connected by it.
        :return: ftms.AngleStat object
        '''
        connections = self.connections(elem)
        if not forward:
            connections = connections[::-1]

        log.debug("elem %s connections %s", elem, connections)
        (stem1, twist1, stem2, twist2, bulge) = ftug.get_stem_twist_and_bulge_vecs(
            self, elem, connections)
        log.debug("stem1 %s, twist1 %s, stem2 %s, twist2 %s, bulge %s",
                  stem1, twist1, stem2, twist2, bulge)

        if round(np.dot(stem1, twist1), 10) != 0 or round(np.dot(stem2, twist2), 10) != 0:
            err = CgIntegrityError("The twists are inconsistent. "
                                   "They should be orthogonal to the corresponding stem vectors."
                                   "Inconsistency found for {},{}".format(elem, connections))
            with log_to_exception(log, err):
                log.error("Angle stem1-twist1 %s dot_product=%s, Angle stem2-twist2 %s degrees dot_product=%s",
                          math.degrees(ftuv.vec_angle(stem1, twist1)
                                       ), np.dot(stem1, twist1),
                          math.degrees(ftuv.vec_angle(stem2, twist2)), np.dot(stem2, twist2),)
            raise err

        try:
            u, v, t, r1, u1, v1 = ftug.get_angle_stat_geometry(
                stem1, twist1, stem2, twist2, bulge)
        except ZeroDivisionError as e:
            with log_to_exception(log, e):
                log.error("Error getting stats for elem %s", elem)
            raise
        dims = self.get_bulge_dimensions(elem)
        ang_type = self.connection_type(elem, connections)
        seq = "&".join(self.get_define_seq_str(elem, adjacent=True))
        log.debug("u %s, v %s", u, v)

        if elem[0] == "m":
            mls = self.find_mlonly_multiloops()
            for ml in mls:
                if elem in ml:
                    descr = self.describe_multiloop(
                        [x for x in ml if x[0] != "s"])
                    if "pseudoknot" in descr:
                        stat_type = "pseudo"
                    elif "open" in descr:
                        stat_type = "open"
                    else:
                        stat_type = "angle"  # ML
                    break
            else:
                assert False
        else:
            stat_type = "angle"  # IL

        angle_stat = ftms.AngleStat(stat_type,
                                    self.name, dims[0], dims[1], u, v, t, r1,
                                    u1, v1, ang_type, self.defines[elem],
                                    seq, self.vposs[elem], self.vbase[elem],
                                    self.vsugar[elem], self.vbackbone[elem])

        return angle_stat

    def get_stats(self, d):
        '''
        Calls get_loop_stat/ get_bulge_angle_stats or get_stem_stats, depending on the element d.

        :returns: A 1- or 2 tuple of stats (2 in case of bulges. One for each direction)
        '''
        if d[0] == "s":
            return (self.get_stem_stats(d),)
        elif d[0] in "mi":
            return self.get_bulge_angle_stats(d)
        else:
            try:
                stat = self.get_loop_stat(d)
            except ValueError:
                if len(self.defines) == 1:
                    # A structure without any stem has no stats.
                    return tuple()
                else:
                    raise
            if d[0] == "f":
                stat.stat_type = "5prime"
            elif d[0] == "t":
                stat.stat_type = "3prime"
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
        loop_stat.phys_length = ftuv.magnitude(
            self.coords[d][1] - self.coords[d][0])
        stem1, = self.edges[d]  # Make sure there is only one edge

        (s1b, s1e) = self.get_sides(stem1, d)

        stem1_vec = self.coords[stem1][s1b] - self.coords[stem1][s1e]
        twist1_vec = self.twists[stem1][s1b]
        bulge_vec = self.coords[d][1] - self.coords[d][0]

        # To avoid loops with 0 physical length. (If disconnects in the structure are modelled as loop)
        if ftuv.magnitude(bulge_vec) < 10**-3:
            bulge_vec += (10**-3) * (stem1_vec / ftuv.magnitude(stem1_vec))

        (r, u, v) = ftug.get_stem_separation_parameters(stem1_vec, twist1_vec,
                                                        bulge_vec)
        (loop_stat.r, loop_stat.u, loop_stat.v) = (r, u, v)
        # Will this cause problems in other parts of the code base???
        loop_stat.r = loop_stat.phys_length
        loop_stat.define = self.defines[d]
        loop_stat.seq, = self.get_define_seq_str(d, adjacent=True)
        loop_stat.vres = self.vposs[d]
        loop_stat.vbase = self.vbase[d]
        loop_stat.vsugar = self.vsugar[d]
        loop_stat.vbackbone = self.vbackbone[d]

        if d[0] == "f":
            loop_stat.stat_type = "5prime"
        elif d[0] == "t":
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

        angle_stat1 = self.get_bulge_angle_stats_core(bulge, True)
        angle_stat2 = self.get_bulge_angle_stats_core(bulge, False)
        # If we go into the reverse direction, the first vector of the element basis is inverted.
        # and thus also the last one (cross product of first and second)
        # for k,v in angle_stat2.vres.items():
        #    angle_stat2.vres[k]= v*[-1,1,-1]
        assert round(angle_stat1.get_angle(), 5) == round(angle_stat2.get_angle(
        ), 5), ("{}!={}".format(angle_stat1.get_angle(), angle_stat2.get_angle()))
        return (angle_stat1, angle_stat2)

    def get_stacking_helices(self, method="Tyagi"):
        """
        EXPERIMENTAL

        Return all helices (longer stacking regions) as sets.

        Two stems and one bulge are in a stacking relation, if self.is_stacking(bulge) is true and the stems are connected to the bulge.
        Further more, a stem is in a stacking relation with itself.
        A helix is the transitive closure this stacking relation.

        :returns: A list of sets of element names.
        """
        helices = []
        for d in self.defines:
            if d[0] in "mi" and self.is_stacking(d, method):
                s1, s2 = self.connections(d)
                helices.append(set([d, s1, s2]))
            if d[0] == "s":
                helices.append(set([d]))
        while True:
            for i, j in it.combinations(range(len(helices)), 2):
                stack_bag1 = helices[i]
                stack_bag2 = helices[j]
                if stack_bag1 & stack_bag2:
                    stack_bag1 |= stack_bag2
                    del helices[j]
                    break
            else:
                break
        return helices

    def is_stacking(self, bulge, method="Tyagi", verbose=False):
        """
        EXPERIMENTAL

        Reports, whether the stems connected by the given bulge are coaxially stacking.

        :param bulge: STRING. Name of a interior loop or multiloop (e.g. "m3")
        :param method": STRING. "Tyagi": Use cutoffs from doi:10.1261/rna.305307, PMCID: PMC1894924.
        :returns: A BOOLEAN.
        """

        assert method in ["Tyagi", "CG"]
        if method == "Tyagi":
            return self._is_stacking_tyagi(bulge, verbose)
        return self._is_stacking_CG(bulge, verbose)

    def _is_stacking_CG(self, bulge, verbose=False):
        """"""
        stem1, stem2 = self.connections(bulge)
        angle = ftuv.vec_angle(self.coords[stem1][1] - self.coords[stem1][0],
                               self.coords[stem2][1] - self.coords[stem2][0])
        if angle > math.pi / 2:
            angle = math.pi - angle
        if angle > math.radians(45):
            if verbose:
                print("Angle {}>45".format(math.degrees(angle)))
            return False
        shear_angle1 = ftuv.vec_angle(self.coords[stem1][1] - self.coords[stem1][0],
                                      self.coords[bulge][1] - self.coords[bulge][0])
        if shear_angle1 > math.pi / 2:
            shear_angle1 = math.pi - shear_angle1
        shear_angle2 = ftuv.vec_angle(self.coords[stem2][1] - self.coords[stem2][0],
                                      self.coords[bulge][1] - self.coords[bulge][0])
        if shear_angle2 > math.pi / 2:
            shear_angle2 = math.pi - shear_angle2
        if shear_angle1 > math.radians(60) or shear_angle2 > math.radians(60):
            if verbose:
                print ("Shear angle 1 {}>60 or shear angle 2 {}>60".format(
                    math.degrees(shear_angle1), math.degrees(shear_angle1)))
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
        DISTANCE_CUTOFF = [14, 6]
        ANGLE_CUTOFF = [math.acos(0.75), math.acos(0.8)]
        # Relaxed compared to 60 in the paper, because we use
        SHEAR_ANGLE_CUTOFF = math.radians(60)
        # virtual atom positions
        SHEAR_OFFSET_CUTOFF = 10
        if bulge[0] == "m" and self.get_length(bulge) == 0:
            is_flush = True  # flush-stack vs. mismatch-mediated stack
        else:
            is_flush = False

        stem1, stem2 = self.connections(bulge)
        side_nts = self.get_connected_residues(stem1, stem2, bulge)[0]
        # Distance
        try:
            bp_center1 = ftug.get_basepair_center(self, side_nts[0])
            bp_center2 = ftug.get_basepair_center(self, side_nts[1])
        except KeyError:
            return False
        if ftuv.vec_distance(bp_center1, bp_center2) > DISTANCE_CUTOFF[is_flush]:
            if verbose:
                print ("Distance {} > {}".format(ftuv.vec_distance(
                    bp_center1, bp_center2), DISTANCE_CUTOFF[is_flush]))
            return False
        normalvec1 = ftug.get_basepair_plane(self, side_nts[0])
        normalvec2 = ftug.get_basepair_plane(self, side_nts[1])
        # Coaxial
        angle = ftuv.vec_angle(normalvec1, normalvec2)
        if angle > math.pi / 2:
            # Triggered frequently
            #warnings.warn("Angle > 90 degrees: {} ({})".format(angle, math.degrees(angle)))
            angle = math.pi - angle
        if angle > ANGLE_CUTOFF[is_flush]:
            if verbose:
                print ("Angle {} > {}".format(angle, ANGLE_CUTOFF[is_flush]))
            return False
        # Shear Angle
        shear_angle1 = ftuv.vec_angle(normalvec1, bp_center2 - bp_center1)
        if shear_angle1 > math.pi / 2:
            shear_angle1 = math.pi - shear_angle1
        if shear_angle1 > SHEAR_ANGLE_CUTOFF:
            if verbose:
                print ("Shear angle 1 {} > {}".format(
                    shear_angle1, SHEAR_ANGLE_CUTOFF))
            return False
        shear_angle2 = ftuv.vec_angle(normalvec2, bp_center1 - bp_center2)
        if shear_angle2 > math.pi / 2:
            shear_angle2 = math.pi - shear_angle2
        if shear_angle2 > SHEAR_ANGLE_CUTOFF:
            if verbose:
                print ("Shear angle 2 {} > {}".format(
                    shear_angle2, SHEAR_ANGLE_CUTOFF))
            return False
        # Shear Offset
        # Formula for distance between a point and a line
        # from http://onlinemschool.com/math/library/analytic_geometry/p_line/
        if (ftuv.magnitude(np.cross((bp_center1 - bp_center2), normalvec2)) /
                ftuv.magnitude(normalvec2)) > SHEAR_OFFSET_CUTOFF:
            if verbose:
                print ("Shear offset 1 wrong:", (ftuv.magnitude(np.cross((bp_center1 - bp_center2), normalvec2)) /
                                                 ftuv.magnitude(normalvec2)), ">", SHEAR_OFFSET_CUTOFF)
            return False
        if (ftuv.magnitude(np.cross((bp_center1 - bp_center2), normalvec1)) /
                ftuv.magnitude(normalvec1)) > SHEAR_OFFSET_CUTOFF:
            if verbose:
                print ("Shear offset 2 wrong")
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
        ss.phys_length = ftuv.magnitude(
            self.coords[stem][0] - self.coords[stem][1])
        ss.twist_angle = ftug.get_twist_angle(
            self.coords[stem], self.twists[stem])
        ss.define = self.defines[stem]
        # Stems may or may not have adjacent nucleotides, so we do not include the adj nucleotides
        ss.seqs = self.get_define_seq_str(stem, adjacent=False)
        ss.vbase = self.vbase[stem]
        ss.vsugar = self.vsugar[stem]
        ss.vbackbone = self.vbackbone[stem]
        return ss

    # def get_loop_from_residue(self, residue) ->  use BulgeGraph.get_node_from_residue_num()!
    def _init_coords(self):
        self.coords = LineSegmentStorage(
            self.defines.keys(), on_change=self.reset_vatom_cache)
        self.twists = CoordinateStorage(
            [x for x in self.defines if x[0] == "s"], on_change=self.reset_vatom_cache)

    @property
    def incomplete_elements(self):
        if self._incomplete_elements is None:
            self._incomplete_elements = ftug.get_incomplete_elements(self)
        return set(self._incomplete_elements)

    @property
    def interacting_elements(self):
        if self._interacting_elements is None:
            interacting_nts = []
            for r in self.interacting_residues:
                try:
                    interacting_nts.append(self.seq_id_to_pos(r))
                except ValueError as e:  # interacting missing residue
                    assert "not in list" in str(e)
                    pass
            self._interacting_elements = set(
                self.nucleotides_to_elements(interacting_nts))
        return self._interacting_elements

    def radius_of_gyration(self, method="vres"):
        '''
        Calculate the radius of gyration of this structure.

        :param method: A STRING. one of
                       "fast" (use only coordinates of coarse grained stems) or
                       "vres" (use virtual residue coordinates of stems)

        :return: A number with the radius of gyration of this structure.
        '''
        if len(list(self.stem_iterator())) == 0:
            log.warning(
                "Cannnot calculate ROG (%s) for structure %s without stems. Returning 'nan'", method, self.name)
            return float("nan")
        if method == "fast":
            coords = self.get_ordered_stem_poss()
        elif method == "vres":
            coords = self.get_ordered_virtual_residue_poss()
        else:
            raise ValueError(
                "Wrong method {}. Choose one of 'fast' and 'vres'".format(method))

        rog = ftmd.radius_of_gyration(coords)
        return rog

    def _assign_loop_roles(self, loop):
        """
        EXPERIMENTAL

        For 3-way junctions: Assign the roles P1, P2 and P3
        and J12, J31, J23 according to Lescoute and Westhof,
        RNA 2006 (doi: 10.1261/rna.2208106)

        We try to detect which helices stack and assign the roles in a way
        that J12 is always the stacking segment.

        :returns: None if no clear classification was possible,
                  a dictionary otherwise.
        """
        if len(loop) != 3:
            raise ValueError(
                "junction_family is currentlty only implemented for 3-way junctions.")
        # Find out which helices stack.
        # Version 1: colinearity of helix axes:
        stems = self.edges[loop[0]] | self.edges[loop[1]] | self.edges[loop[2]]
        if len(stems) != 3:
            e = ValueError(
                "The multiloop segments do not form a 3-way junction.")
            with log_to_exception(log, e):
                log.error("stems are %s", stems)
            raise e
        collinearities = {}
        for stem1, stem2 in it.combinations(stems, 2):
            c = ftuv.line_segments_collinearity(
                self.coords[stem1], self.coords[stem2])
            collinearities[(stem1, stem2)] = c
        # The two helices that are most likely to stack.
        stacking_stems = max(collinearities, key=lambda k: collinearities[k])
        if collinearities[stacking_stems] < 0.8:
            return None
        # The second best pair
        second_best_pair, second_best_score = list(
            sorted(collinearities.items(), key=lambda x: x[1]))[1]
        if collinearities[stacking_stems] - second_best_score < 0.5:
            # Two are possible. We look now which one has fewer offset to the 3rd stem.
            # The reference stem could stack with both.
            reference = set(stacking_stems) & set(second_best_pair)
            other = set(stems) - reference
            reference, = reference
            stem1, stem2 = other
            o1 = self.stem_offset(reference, stem1)
            o2 = self.stem_offset(reference, stem2)
            if abs(o1 - o2) < 0.2:  # Angstrom
                log.warning("Cannot assign loop roles. Stacking is ambiguous.")
                return None
            elif o1 < o2:
                stacking_stems = (reference, stem1)
            else:
                stacking_stems = (reference, stem2)
        # We have clearly identified the stacking helices.
        # j12 is the ML-segment that connects the stacking helices.
        j12, = self.edges[stacking_stems[0]] & self.edges[stacking_stems[1]]
        # stem P1 is at the 5' side of the connecting multiloop...
        p1 = self.get_node_from_residue_num(self.define_a(j12)[0])
        # ...and P2 at the 3' side
        p2 = self.get_node_from_residue_num(self.define_a(j12)[1])
        assert set([p1, p2]) == set(stacking_stems)
        assert p1 != p2
        # p3 is the third stem,
        p3, = stems - set([p1, p2])
        # j31 and j23 are the remaining ml segments
        j31, = self.edges[p3] & self.edges[p1]
        j23, = self.edges[p2] & self.edges[p3]
        return {"P1": p1, "P2": p2, "P3": p3, "J12": j12, "J31": j31, "J23": j23}

    def _junction_family_westhof1(self, loop):
        """
        EXPERIMENTAL

        For 3-way junctions: Return the junction family according to
        Lescoute and Westhof, RNA 2006 (doi: 10.1261/rna.2208106).

        In this method, the junction family is defined via the relative
        lengths of fragments, not the orientation of the elements.

        :param loop: Either a dictionary, as returned by self._assign_loop_roles
                     or a list of 3 loop segments.
        """
        if isinstance(loop, dict):
            j_roles = loop
        else:
            j_roles = self._assign_loop_roles(loop)
        # Now compare the lengths of J31 and J23 to determine the junction family
        if self.get_length(j_roles["J31"]) > self.get_length(j_roles["J23"]):
            return "A"
        elif self.get_length(j_roles["J31"]) == self.get_length(j_roles["J23"]):
            return "B"
        else:
            return "C"

    def _junction_family_3d(self, loop):
        """
        EXPERIMENTAL

        For 3-way junctions: Return the junction family depending
        on the 3D orientation.

        :param loop: Either a dictionary, as returned by self._assign_loop_roles
                     or a list of 3 loop segments.
        """
        if isinstance(loop, dict):
            j_roles = loop
        else:
            j_roles = self._assign_loop_roles(loop)
        a31 = self.stem_angle(j_roles["P3"], j_roles["P1"])
        a23 = self.stem_angle(j_roles["P2"], j_roles["P3"])
        # If the stack is perfectly straight, a31+a23=180 degrees.
        if a31 > a23:
            return 2
        else:
            return 1

    def _junction_family_is_perpenticular(self, loop):
        """
        EXPERIMENTAL

        For 3-way junctions Return whether or not the not-stacking loop
        is very roughly perpenticular to the stack.

        :param loop: Either a dictionary, as returned by self._assign_loop_roles
                     or a list of 3 loop segments.
        """
        if isinstance(loop, dict):
            j_roles = loop
        else:
            j_roles = self._assign_loop_roles(loop)
        a31 = self.stem_angle(j_roles["P3"], j_roles["P1"])
        a23 = self.stem_angle(j_roles["P2"], j_roles["P3"])
        # If the stack is perfectly straight, a31+a23=180 degrees.
        if math.radians(60) < a31 < math.radians(120) or math.radians(60) < a23 < math.radians(120):
            return 1
        else:
            return 0

    def stem_offset(self, ref_stem, stem2):
        """
        How much is the offset between the start of stem 2
        and the axis of stem1.

        Assumes that stem1 and stem 2 are connected by a single bulge.
        Then the start of stem2 is defined to be the stem side
        closer to the bulge.
        """
        common_edges = self.edges[ref_stem] & self.edges[stem2]
        if len(common_edges) != 1:
            raise ValueError("Stem1 and stem2 must be connected by a single bulge "
                             "to calculate stem_offset.")
        bulge, = common_edges
        side = self.get_sides(stem2, bulge)[0]
        return ftuv.point_line_distance(self.coords[stem2][side],
                                        self.coords[ref_stem][0],
                                        self.coords.get_direction(ref_stem))

    def stem_angle(self, stem1, stem2):
        """
        Returns the angle between two stems.

        If they are connected via a single element,
        use the direction pointing away from this element for both stems.
        Otherwise, use the direction from start to end.
        """
        vec1 = self.coords.get_direction(stem1)
        vec2 = self.coords.get_direction(stem2)

        common_edges = self.edges[stem1] & self.edges[stem2]
        if len(common_edges) == 1:
            bulge, = common_edges
            if self.get_sides(stem1, bulge) == (1, 0):
                vec1 = -vec1
            else:
                assert self.get_sides(stem1, bulge) == (0, 1)
            if self.get_sides(stem2, bulge) == (1, 0):
                vec2 = -vec2
            else:
                assert self.get_sides(stem2, bulge) == (0, 1)
        return ftuv.vec_angle(vec1, vec2)

    def get_coordinates_array(self):
        '''
        Get all of the coordinates in one large array.

        The coordinates are sorted in the order of the keys
        in coordinates dictionary.

        :return: A 2D numpy array containing all coordinates
        '''
        all_coords = []
        assert len(self.coords) == len(
            self.defines), self.coords.keys() ^ self.defines.keys()
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
            self.coords[key] = coords[2 * j], coords[2 * j + 1]
        return self

    def get_twists(self, node):
        '''
        Get the array of twists for this node. If the node is a stem,
        then the twists will simply those stored in the array.
        If the node is an interior loop or a junction segment,
        then the twists will be the ones that are adjacent to it,
        projected to the plane normal to the element vector.
        If the node is a hairpin loop or a free end, then the same twist
        will be duplicated and returned twice.

        :param node: The name of the node
        '''
        if node[0] == 's':
            return self.twists[node]

        connections = list(self.edges[node])
        if not connections:
            assert node[0] == "f"
            log.warning("Twists for structures without any stems are 'nan'")
            return np.zeros(3) * float("nan"), np.zeros(3) * float("nan")
        (s1b, s1e) = self.get_sides(connections[0], node)

        if len(connections) == 1:
            vec = ftuv.normalize(ftuv.vector_rejection(
                self.twists[connections[0]][s1b],
                self.coords[node][1] -
                self.coords[node][0]))

            return (vec, vec)

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
        total_length = sum(
            [len(list(self.define_residue_num_iterator(d))) for d in self.defines])
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
                       key=lambda x: (priority[x[0]], min(self.get_node_dimensions(x)), not x in self.sampled, x))
        return edges

    def coords_to_directions(self):
        """
        The directions of each coarse grain element. One line per cg-element.

        The array is sorted by the corresponding element names alphabetically (`sorted(defines.keys()`)
        The directions usually point away from the elemnt's lowest nucleotide.
        However h,t and f elements always point away from the connected stem.
        """
        coords = self.get_coordinates_array()
        directions = coords[1::2] - coords[0::2]
        return directions

    def coords_from_directions(self, directions):
        """
        Generate coordinates from direction vectors (using also their lengths)

        Currently ignores the twists!

        :param directions: An array of vectors from the side of a cg-element with lower nucleotide number to the side with higher number
                           The array is sorted by the corresponding element names alphabetically (`sorted(defines.keys()`)

        """
        sorted_defines = sorted(self.defines.keys())
        assert len(sorted_defines) == len(directions), "{} != {}".format(
            len(sorted_defines), len(directions))
        if self.build_order is None:
            self.traverse_graph()
        self.coords["s0"] = np.array(
            [0, 0, 0]), directions[sorted_defines.index("s0")]

        for stem1, link, stem2 in self.build_order:  # Bulges and stems
            conn = self.connection_ends(
                self.connection_type(link, [stem1, stem2]))
            link_dir = self.get_link_direction(stem1, stem2, link)
            if link_dir == 1:
                self.coords[link] = self.coords[stem1][conn[0]
                                                       ], self.coords[stem1][conn[0]] + directions[sorted_defines.index(link)]
                if conn[1] == 0:
                    self.coords[stem2] = self.coords[link][1], self.coords[link][1] + \
                        directions[sorted_defines.index(stem2)]
                else:
                    self.coords[stem2] = self.coords[link][1] - \
                        directions[sorted_defines.index(
                            stem2)], self.coords[link][1]
            else:
                self.coords[link] = self.coords[stem1][conn[0]] - \
                    directions[sorted_defines.index(
                        link)], self.coords[stem1][conn[0]]
                if conn[1] == 0:
                    self.coords[stem2] = self.coords[link][0], self.coords[link][0] + \
                        directions[sorted_defines.index(stem2)]
                else:
                    self.coords[stem2] = self.coords[link][0] - \
                        directions[sorted_defines.index(
                            stem2)], self.coords[link][0]
        for d in self.defines:
            if d[0] == "m" and d not in self.mst:
                edges = list(self.edges[d])
                (s1b, _) = self.get_sides(edges[0], d)
                (s2b, _) = self.get_sides(edges[1], d)
                mids1 = self.coords[edges[0]]
                mids2 = self.coords[edges[1]]
                # Save coordinates in direction of the strand.
                if self.get_link_direction(edges[0], edges[1], d) == 1:
                    self.coords[d] = (mids1[s1b], mids2[s2b])
                else:
                    self.coords[d] = (mids2[s2b], mids1[s1b])

            if d[0] in "hft":  # Loops
                stem, = self.edges[d]
                (s1b, _) = self.get_sides(stem, d)
                self.coords[d] = self.coords[stem][s1b], self.coords[stem][s1b] + \
                    directions[sorted_defines.index(d)]

    def virtual_atoms(self, key):
        """
        Get virtual atoms for a key.

        :param key: An INTEGER: The number of the base in the RNA.

        :returns: A dict {atom:coords}, e.g. {"C8":np.array([x,y,z]), ...}
        """
        if isinstance(key, int):
            if key not in self._virtual_atom_cache:
                try:
                    self._virtual_atom_cache[key] = ftug.virtual_atoms(self)[
                        key]
                except KeyError as e:
                    log.info(
                        "Key {} not present. Need to recreate all virtual residues".format(e))
                    self.add_all_virtual_residues()
                    self._virtual_atom_cache[key] = ftug.virtual_atoms(self)[
                        key]
            return self._virtual_atom_cache[key]
        else:
            raise ValueError("Expected an int, found {}".format(key))

    def reset_vatom_cache(self, key):
        """
        Delete all cached information about virtual residues and virtual atoms.
        Used as on_call function for the observing of the self.coords dictionary.

        :param key: A coarse grain element name, e.g. "s1" or "m15"
        """
        try:
            if not self._virtual_atom_cache:
                return
        except AttributeError:  # Happens during deepcopy
            return

        define = self.defines[key]

        # Do not delete self.vposs, it's in the element's coordinate system

        # Delete virtual residues
        try:
            del self.vbases[key]
        except KeyError:
            pass
        try:
            del self.vvecs[key]
        except KeyError:
            pass
        try:
            del self.v3dposs[key]
        except KeyError:
            pass
        try:
            del self.vinvs[key]
        except KeyError:
            pass

        # Delete virtual atoms
        if len(define) > 1:
            for i in range(define[0], define[1] + 1):
                if i in self._virtual_atom_cache:
                    del self._virtual_atom_cache[i]
        if len(define) > 3:
            for i in range(define[2], define[3] + 1):
                if i in self._virtual_atom_cache:
                    del self._virtual_atom_cache[i]
    # def __deepcopy__(self, memo):

    def rotate(self, angle, axis="x", unit="radians"):
        if unit == "degrees":
            angle = math.radians(angle)
        elif unit != "radians":
            raise ValueError(
                "Unit {} not understood. Use 'degrees' or 'radians'".format(unit))
        s = math.sin(angle)
        cosi = math.cos(angle)
        rotation_matrix = ftuv.rotation_matrix(axis, angle)
        self.coords.rotate(rotation_matrix)
        self.twists.rotate(rotation_matrix)
        self.after_coordinates_changed()
        for chain in self.chains.values():
            chain.transform(rotation_matrix.T, [0, 0, 0])

    def rotate_translate(self, offset, rotation_matrix):
        """
        First translate the RNA by offset, then rotate by rotation matrix
        """
        self.coords._coordinates -= offset
        self.coords.rotate(rotation_matrix)
        self.twists.rotate(rotation_matrix)
        self.after_coordinates_changed()
        for chain in self.chains.values():
            chain.transform([[1, 0, 0], [0, 1, 0], [0, 0, 1]], -offset)
            chain.transform(rotation_matrix.T, [0, 0, 0])

    def after_coordinates_changed(self):
        # vposs is in the element coordinate system ==> does not have to be changed.
        #self.vposs = c.defaultdict( dict )
        # Caching for virtual residues
        self.vbases = c.defaultdict(dict)
        self.vvecs = c.defaultdict(dict)
        self.v3dposs= c.defaultdict(dict)
        self.vinvs = c.defaultdict(dict)
