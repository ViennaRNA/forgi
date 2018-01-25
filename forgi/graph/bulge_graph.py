#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip, object)
"""bulge_graph.py: A graph representation of RNA secondary structure based
   on its decomposition into primitive structure types: stems, hairpins,
   interior loops, multiloops, etc..."""


import sys
import collections as col
import random
import re
import itertools as it
import os
import warnings
import operator as oper
import functools
import traceback
import math
from string import ascii_lowercase, ascii_uppercase
import logging
from pprint import pprint, pformat

import numpy as np

from logging_exceptions import log_to_exception, log_at_caller

from ..utilities import debug as fud
from ..utilities import stuff as fus
from ..utilities.exceptions import GraphConstructionError, GraphIntegrityError
from ..threedee.utilities import mcannotate as ftum
from .sequence import Sequence
from .residue import RESID, resid_to_str, resid_from_str
from ._basegraph import BaseGraph
from ._graph_construction import _BulgeGraphConstruction

log = logging.getLogger(__name__)

VALID_CHAINIDS = ascii_uppercase+ascii_lowercase

try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x):
    return x

def add_bulge(bulges, bulge, context, message):
    """
    A wrapper for a simple dictionary addition
    Added so that debugging can be made easier

    :param bulges:
    :param bulge:
    :param context:
    :param message:
    :return:
    """
    # bulge = (context, bulge)
    bulges[context] = bulges.get(context, []) + [bulge]
    return bulges


def from_id_seq_struct(id_str, seq, struct):
    """
    Return a new BulgeGraph with the given id,
    sequence and structure.

    :param id_str: The id (i.e. >1y26)
    :param seq: the sequence (i.e. 'ACCGGG')
    :param struct: The dotplot secondary structure (i.e. '((..))')
    """
    if seq is not None and len(seq)!=len(struct):
        raise GraphConstructionError("Sequence and structure length are not equal for id {}".format(id_str))
    bg = BulgeGraph()
    bg.from_dotbracket(struct)
    if id_str is not None:
        bg.name = id_str
    if seq is not None:
        bg.seq = seq
    bg.seq_ids_from_seq()

    return bg


def from_fasta_text(fasta_text, dissolve_length_one_stems=False):
    """
    Create a bulge graph or multiple bulge
    graphs from some fasta text.
    """
    # compile searches for the fasta id, sequence and
    # secondary structure respectively
    id_search = re.compile('>(.+)')
    seq_search = re.compile('^([acgutACGUT&]+)$')
    stru_search = re.compile('^([(){}<>.A-Za-z&\[\]]+)$')

    prev_id = None
    prev_seq = None
    prev_struct = None
    curr_id=None

    bgs = []

    for i, line in enumerate(fasta_text.split('\n')):
        # newlines suck
        line = line.strip()
        # We allow comments
        if line.startswith("#"):
            continue
        # find out what this line contains
        id_match = id_search.match(line)
        seq_match = seq_search.match(line)
        stru_match = stru_search.match(line)

        if id_match is not None:
            prev_id=curr_id
            # we found an id, check if there's a previous
            # sequence and structure, and create a BG
            curr_id = id_match.group(1)

            if prev_seq is None and prev_struct is None:
                # must be the first sequence/structure
                continue

            if prev_struct is None:
                raise GraphConstructionError("No structure for id: {}", prev_id)

            bg = from_id_seq_struct(prev_id, prev_seq, prev_struct)
            if dissolve_length_one_stems:
                bg.dissolve_length_one_stems()
            bgs.append(bg)

            prev_seq = None
            prev_struct = None
            prev_id=None

        if seq_match is not None:
            curr_seq = seq_match.group(1)
            if "t" in curr_seq or "T" in curr_seq:
                warnings.warn("Original sequence contained T. All occurrences of T/t were replaced by U/u respectively!")
                curr_seq=curr_seq.replace("T", "U")
                curr_seq=curr_seq.replace("t", "u")
            if prev_seq:
                prev_seq+=curr_seq
            else:
                prev_seq=curr_seq

        if id_match is None and seq_match is None:
            if stru_match:
                if prev_struct:
                    prev_struct += line
                else:
                    prev_struct = line
            elif line:
                raise GraphConstructionError("Cannot parse line {}: '{}' is neither sequence, nor structure, nor name (starting with '>'), nor comment (starting with '#').".format(i, line))

    if prev_struct is None:
        raise GraphConstructionError("Error during parsing of fasta file. No structure found for id {} and sequence {}".format(prev_id, prev_seq))

    bg = from_id_seq_struct(curr_id, prev_seq, prev_struct)
    if dissolve_length_one_stems:
        bg.dissolve_length_one_stems()
    bgs.append(bg)

    if len(bgs) == 1:
        return bgs[0]
    else:
        return bgs

def from_fasta(filename, dissolve_length_one_stems=False):
    with open(filename) as f:
        fasta_text = f.read()
    return from_fasta_text(fasta_text, dissolve_length_one_stems)

def any_difference_of_one(stem, bulge):
    """
    See if there's any difference of one between the two
    ends of the stem [(a,b),(c,d)] and a bulge (e,f)

    :param stem: A couple of couples (2 x 2-tuple) indicating the start and end
                 nucleotides of the stem in the form ((s1, e1), (s2, e2))
    :param bulge: A couple (2-tuple) indicating the first and last position
                  of the bulge.
    :return: True if there is an overlap between the stem nucleotides and the
                  bulge nucleotides. False otherwise
    """
    for stem_part in stem:
        for part in stem_part:
            for bulge_part in bulge:
                if abs(bulge_part - part) == 1:
                    return True
    return False


def print_bulges(bulges):
    """
    Print the names and definitions of the bulges.

    :param bulges: A list of tuples of the form [(s, e)] where s and e are the
                   numbers of the nucleotides at the start and end of the bulge.
    """
    for i in range(len(bulges)):
        # print "bulge:", bulge
        bulge_str = "define b{} 1".format(i)
        bulge = bulges[i]
        bulge_str += " {} {}".format(bulge[0] + 1, bulge[1] + 1)
        print (bulge_str)


def condense_stem_pairs(stem_pairs):
    """
    Given a list of stem pairs, condense them into stem definitions

    I.e. the pairs (0,10),(1,9),(2,8),(3,7) can be condensed into
    just the ends of the stem: [(0,10),(3,7)]

    :param stem_pairs: A list of tuples containing paired base numbers.

    :returns: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    """
    stem_pairs.sort()

    prev_pair = (-10, -10)

    stems = []
    start_pair = None

    for pair in stem_pairs:
        # There's a potential bug here since we don't check the direction
        # but hopefully it won't bite us in the ass later
        if abs(pair[0] - prev_pair[0]) != 1 or abs(pair[1] - prev_pair[1]) != 1:
            if start_pair is not None:
                stems += [(start_pair, prev_pair)]
            start_pair = pair

        prev_pair = pair

    if start_pair is not None:
        stems += [(start_pair, prev_pair)]

    return stems


def print_brackets(brackets):
    """
    Print the brackets and a numbering, for debugging purposes

    :param brackets: A string with the dotplot passed as input to this script.
    """
    numbers = [chr(ord('0') + i % 10) for i in range(len(brackets))]
    tens = [chr(ord('0') + i // 10) for i in range(len(brackets))]
    print ("brackets:\n", brackets, "\n", "".join(tens), "\n", "".join(numbers))

# @Coverage: Seems to be unused. If this is removed, condense_stem_pairs can be removed as well.
def find_bulges_and_stems(brackets):
    """
    Iterate through the structure and enumerate the bulges and the stems that are
    present.

    The returned stems are of the form [[(s1, s2), (e1,e2)], [(s1,s2),(e1,e2)],...]
    where (s1,s2) are the residue numbers of one end of the stem and (e1,e2) are the
    residue numbers at the other end of the stem
    (see condense_stem_pairs)

    The returned bulges are of the form [(s,e), (s,e),...] where s is the start of a bulge
    and e is the end of a bulge

    :param brackets: A string with the dotbracket passed as input to this script.
    """
    prev = 'x'
    context = 0

    bulges = dict()
    finished_bulges = []
    context_depths = dict()

    opens = []
    stem_pairs = []

    dots_start = 0
    context_depths[0] = 0

    i = 0
    for i in range(len(brackets)):
        if brackets[i] == '(':
            opens.append(i)

            if prev == '(':
                context_depths[context] = context_depths.get(context, 0) + 1
                continue
            else:
                context += 1
                context_depths[context] = 1

            if prev == '.':
                dots_end = i - 1
                bulges = add_bulge(bulges, (dots_start, dots_end), context, "4")

        if brackets[i] == ')':
            if len(opens) == 0:
                raise Exception("Unmatched close bracket")

            stem_pairs.append((opens.pop(), i))

            context_depths[context] -= 1

            if context_depths[context] == 0:
                if context in bulges:
                    finished_bulges += bulges[context]
                bulges[context] = []
                context -= 1

            if prev == '.':
                dots_end = i - 1
                bulges = add_bulge(bulges, (dots_start, dots_end), context, "2")

        if brackets[i] == '.':
            if prev == '.':
                continue

            dots_start = i

        prev = brackets[i]

    if prev == '.':
        dots_end = i
        bulges = add_bulge(bulges, (dots_start, dots_end), context, "7")
    elif prev == '(':
        print ("Unmatched bracket at the end", file=sys.stderr)
        sys.exit(1)
    """
    elif prev == ')':
        bulges = add_bulge(bulges, (i+1, i+1), context, "8")
    """

    if context in bulges.keys():
        finished_bulges += bulges[context]

    if len(opens) > 0:
        raise Exception("Unmatched open bracket")

    stem_pairs.sort()
    stems = condense_stem_pairs(stem_pairs)

    return finished_bulges, stems

class BulgeGraph(BaseGraph):
    def __init__(self, bg_file=None, dotbracket_str='', seq=''):
        """
        A bulge graph object.

        :var self.defines: The coarse grain element definitions: Keys are for example 's1'/ 'm2'/ 'h3'/ 'f1'/ 't1'
                       Values are the positions in the sequence (1D-coordinate) of start , end, ...
        """
        self.ang_types = None
        self.mst = None
        self.build_order = None
        self.name = "untitled"
        #: The coarse grain element definitions: Keys are for example 's1'/ 'm2'/ 'h3'/ 'f1'/ 't1'
        #: Values are the positions in the sequence (1D-coordinate) of start , end, ...
        self.defines = dict()
        self.edges = col.defaultdict(set)
        self.longrange = col.defaultdict(set)

        # Some cached values:
        self.nx_graph = None
        self.nuc_bp_dists = None
        self._elem_bp_dists = {}

        # Additional infos as key-value pairs are stored here.
        self.infos = col.defaultdict(list)


        #Consistency check, only if both dotbracket and sequence are present.
        if dotbracket_str and seq:
            db_strs = dotbracket_str.split('&')
            seq_strs = seq.split('&')
            if not len(seq_strs)==len(db_strs) or any(len(db_strs[i])!=len(seq_strs[i])
                                                      for i in range(len(db_strs))):
                raise GraphConstructionError("Sequence and dotbracket string are not consistent!")
        elif dotbracket_str:
            seq = []
            for substr in dotbracket_str.split('&'):
                seq.append("N"*len(substr))
            seq="&".join(seq)

        seq_strs = seq.split('&')
        seq_ids = []
        for i, seq_str in enumerate(seq_strs):
            for j, s in enumerate(seq):
                seq_ids += [resid_from_str("{}:{}".format(VALID_CHAINIDS[i], j+1))]

        self._seq = Sequence(seq, seq_ids)

        if dotbracket_str:
            self._from_dotbracket(dotbracket_str)
        if bg_file is not None:
            self.from_bg_file(bg_file)

    @property
    def seq_length(self):
        return len(self.seq)

    @property
    def seq(self):
        return self._seq



    def element_length(self, key):
        """
        Get the number of residues that are contained within this element.

        :param key: The name of the element.
        """
        d = self.defines[key]
        length = 0

        for i in range(0, len(d), 2):
            length += d[i + 1] - d[i] + 1

        return length

    def stem_length(self, key):
        """
        Get the length of a particular element. If it's a stem, it's equal to
        the number of paired bases. If it's an interior loop, it's equal to the
        number of unpaired bases on the strand with less unpaired bases. If
        it's a multiloop, then it's the number of unpaired bases.
        """
        d = self.defines[key]
        if key[0] == 's' or key[0] == 'y':
            return (d[1] - d[0]) + 1
        elif key[0] == 'f':
            return self.get_bulge_dimensions(key)[0]
        elif key[0] == 't':
            return self.get_bulge_dimensions(key)[1]
        elif key[0] == 'h':
            return self.get_bulge_dimensions(key)[0]
        else:
            return min(self.get_bulge_dimensions(key))

    def add_info(self, key, value):
        self.infos[key].append(value)

    def get_single_define_str(self, key):
        """
        Get a define string for a single key.
        """
        return "define {} {}".format(key, " ".join([str(d) for d in self.defines[key]]))

    def get_define_str(self):
        """
        Convert the defines into a string.

        Format:

        define [name] [start_res1] [end_res1] [start_res2] [end_res2]
        """
        defines_str = ''

        # a method for sorting the defines
        def define_sorter(k):
            drni = self.define_residue_num_iterator(k, adjacent=True)
            return next(drni) #.next()

        for key in sorted(self.defines.keys(), key=define_sorter):
            defines_str += self.get_single_define_str(key)
            # defines_str += "define %s %s" % ( key, " ".join([str(d) for d in self.defines[key]]))
            defines_str += '\n'
        return defines_str

    def get_length_str(self):
        return "length " + str(self.seq_length) + '\n'

    def get_info_str(self):
        out=""
        for info in self.infos:
            for value in self.infos[info]:
                out+="info "+info+" "+value+"\n"
        return out

    def get_connect_str(self):
        """
        Get the connections of the bulges in the graph.

        Format:

        connect [from] [to1] [to2] [to3]
        """

        whole_str = ''
        for key in sorted(self.edges):
            if len(self.edges[key]) == 0:
                continue

            # Our graph will be defined by the stems and the bulges they connect to
            name = key
            if name[0] == 's':
                out_str = "connect {}".format(name)

                for dest in self.connections(key):
                    out_str += " {}".format(dest)

                whole_str += out_str
                whole_str += '\n'

        return whole_str

    def get_sequence_str(self):
        """
        Return the sequence along with its keyword. I.e.

            seq ACGGGCC
        """
        if len(self.seq) > 0:
            return "seq {}\n".format(self.seq)
        else:
            return ""

    def get_seq_ids_str(self):
        """
        Return the sequence id string

        seq_ids 1 2 2.A 17
        """
        out_str = "seq_ids "
        out_str += " ".join(map(resid_to_str, self.seq_ids))
        out_str += "\n"

        return out_str

    def get_name_str(self):
        """
        Return the name of this structure along with its keyword:

            name 1y26
        """
        return "name {}\n".format(self.name)

    def to_bg_string(self):
        """
        Output a string representation that can be stored and reloaded.

        """
        out_str = ''
        out_str += self.get_name_str()
        out_str += self.get_length_str()
        out_str += self.get_sequence_str()
        out_str += self.get_seq_ids_str()
        out_str += self.get_define_str()
        out_str += self.get_connect_str()
        out_str += self.get_info_str()
        return out_str

    def to_file(self, filename):
        with open(filename, 'w') as f:
            out_str = self.to_bg_string()

            f.write(out_str)

    def to_element_string(self, with_numbers=False):
        """
        Create a string similar to dotbracket notation that identifies what
        type of element is present at each location.

        For example the following dotbracket:

        ..((..))..

        Should yield the following element string:

        ffsshhsstt

        Indicating that it begins with a fiveprime region, continues with a
        stem, has a hairpin after the stem, the stem continues and it is terminated
        by a threeprime region.

        :param with_numbers: show the last digit of the element id in a second line.::

                                 (((.(((...))))))

                             Could result in::

                                 sssissshhhssssss
                                 0000111000111000

                             Indicating that the first stem is named 's0', followed by 'i0','
                             s1', 'h0', the second strand of 's1' and the second strand of 's0'
        """
        output_str = [' '] * (self.seq_length + 1)
        output_nr = [' '] * (self.seq_length + 1)
        for d in self.defines.keys():
            for resi in self.define_residue_num_iterator(d, adjacent=False):
                output_str[resi] = d[0]
                output_nr[resi] = d[-1]
        if with_numbers:
            return "".join(output_str).strip()+"\n"+"".join(output_nr).strip()
        else:
            return "".join(output_str).strip()
    def to_neato_string(bg):

        # The different nodes for different types of bulges
        node_lines = dict()

        fontsize=20
        out = []
        out.append("graph G {")
        out.append("\tgraph [overlap=false,splines=true];")
        out.append("\tnode [shape=box];")

        for key2 in bg.defines.keys():
            # Create the nodes with a different color for each type of element
            if key2[0] == 's':
                out.append( '\t{node [style=filled,fillcolor="#B3E2CD",fontsize=%d,label=\"%s\\n(%d)\"] %s};' % (fontsize,key2, bg.stem_length(key2), key2) )
                continue
            elif key2[0] == 'i':
                out.append(  '\t{node [style=filled,shape=circle,fillcolor="#FFF2AE",fontsize=%d' % (fontsize) )
            elif key2[0] == 'm':
                out.append(  '\t{node [style=filled,shape=circle,fillcolor="#F4CAE4",fontsize=%d' % (fontsize) )
            elif key2[0] == 'f':
                out.append(  '\t{node [style=filled,shape=circle,fillcolor="#FDCDAC",fontsize=%d' % (fontsize) )
            elif key2[0] == 't':
                out.append(  '\t{node [style=filled,shape=circle,fillcolor="#E6F5C9",fontsize=%d' % (fontsize) )
            else:
                out.append(  '\t{node [style=filled,shape=circle,fillcolor="#CBD5E8",fontsize=%d' % (fontsize) )

            out[-1] += ',label=\"%s \\n' % (key2)

            # figure out the size of the node and use that as a label
            node_dims = bg.get_node_dimensions(key2)
            total_bulge = sum(node_dims)

            if node_dims[0] == -1 or node_dims[0] == 1000:
                node_dims = "({})".format(node_dims[1])
            elif node_dims[1] == -1 or node_dims[1] == 1000:
                node_dims = "({})".format(node_dims[0])
            log.info("Dims of node %s are %r", key2, node_dims)
            out[-1] += str(node_dims)

            # make bigger interior loops visually bigger
            width = math.sqrt(1.5 * total_bulge / 10.0)
            height = width

            if key2[0] == 'i':
                out[-1] += "\",width=%f,heigh=%f] %s};" % (width, height, key2)
            else:
                out[-1] += "\"] %s};" % (key2)


        for key1 in bg.edges:
            if key1[0] == 's':
                for key2 in bg.edges[key1]:
                    out.append(  "\t%s -- %s;" % (key1, key2) )

        #for key1 in bg.longrange.keys():
    #        for key2 in bg.longrange[key1]:
#                out.append( "\t{%s -- %s [style=dashed]};" % (key1, key2))

        out.append("}")
        return "\n".join(out)

    def define_a(self, elem):
        # Special case, because interior loops can have
        # defines of length 2 or 4
        if elem[0]=="i":
            conns = self.connections(elem)
            s1 = self.defines[conns[0]]
            s2 = self.defines[conns[1]]
            return [s1[1] , s2[0], s2[3] , s1[2]]
        elif:
            # The following may call this classes
            # _define_a_nonzero implementation.
            return super(BulgeGraph, self).define_a(elem)

    def _define_a_nonzero(self, elem):
        """
        Get a define including the adjacent nucleotides.
        """
        define = self.defines[elem]
        new_def = []
        for i in range(0,len(define),2):
            if define[i]-1 in self.backbone_breaks_after:
                new_def.append(define[i])
            else:
                new_def.append(max(define[i]-1, 1))
            if define[i+1] in self.backbone_breaks_after:
                new_def.append(define[i+1])
            else:
                new_def.append(min(define[i+1]+1, self.seq_length))
        return new_def

    def define_residue_num_iterator(self, node, adjacent=False, seq_ids=False):
        """
        Iterate over the residue numbers that belong to this node.

        :param node: The name of the node
        """
        visited=set()

        for r in self.define_range_iterator(node, adjacent):
            for i in range(r[0], r[1] + 1):
                if seq_ids:
                    if self.seq_ids[i-1] not in visited:
                        visited.add(self.seq_ids[i-1])
                        yield self.seq_ids[i - 1]
                else:
                    if i not in visited:
                        visited.add(i)
                        yield i

    def iterate_over_seqid_range(self, start_id, end_id):
        """
        Iterate over the seq_ids between the start_id and end_id.
        """
        i1 = self.seq_ids.index(start_id)
        i2 = self.seq_ids.index(end_id)

        for i in range(i1, i2 + 1):
            yield self.seq_ids[i]

    def seq_id_to_pos(self, seq_id):
        """
        Convert a pdb seq_id to a 1-based nucleotide position

        :param seq_id: An instance of RESID
        """
        assert isinstance(seq_id, RESID)
        #if isinstance(seq_id, int):
        #    seq_id=(" ", seq_id, " ")
        try:
            return self.seq_ids.index(seq_id)+1
        except ValueError as e:
            log.debug("seq_id is {}, self.seq_ids is {}".format(seq_id, self.seq_ids))
            raise

    def shortest_bg_loop(self, vertex):
        """
        Find a shortest loop containing this node. The vertex should
        be a multiloop.

        :param vertex: The name of the vertex to find the loop.
        :return: A list containing the elements in the shortest cycle.
        """
        log.debug("Starting shortest BG loop for {}".format(vertex))
        G = self.to_networkx()
        log.debug("nx graph  %r with edges %r ", G, G.adj)

        # use the nucleotide in the middle of this element as the starting point
        residues = sorted(list(self.define_residue_num_iterator(vertex, adjacent=True)))
        mid_res = residues[len(residues) // 2]
        log.debug("mid_residue %s", mid_res)
        if len(residues) == 2:
            # no-residue multiloop
            # find the neighbor which isn't part of the multiloop
            neighbors = [n for n in G.neighbors(mid_res) if n != residues[0]]

            if len(neighbors) == 2:
                # break the chain so that we don't get cycles within a stem
                for n in neighbors:
                    if abs(n - mid_res) == 1:
                        G.remove_edge(n, mid_res)
                        break

        import forgi.utilities.graph as fug

        path = fug.shortest_cycle(G, mid_res)
        log.debug("Shortest cycle is %s", path)
        return path

    def _chain_start_end(self, pos):
      if pos not in self.backbone_breaks_after:
        if pos == self.seq_length:
          if self.backbone_breaks_after:
            return self.backbone_breaks_after[-1]
          else:
            return 1
        else:
          raise ValueError("Pos {} is not at the end of a chain.".format(pos))
      else:
        i = self.backbone_breaks_after.index(pos)
        if i==0:
          return 1
        else:
          return self.backbone_breaks_after[i-1]+1
    @profile
    def get_next_ml_segment(self, ml_segment):
        """
        Get the adjacent multiloop-segment (or 3' loop) next to the 3' side of ml_segment.

        If there is no other single stranded RNA after the stem, the backbone must end there.
        In that case return None.
        """
        log.debug("get_next_ml_segment called for {}".format(ml_segment))
        if ml_segment.startswith("t"):
            return None
        else:
            if ml_segment[0] in "mf":
                f = self.define_a(ml_segment)[-1]
            else:
                raise ValueError("{} is not a multiloop".format(ml_segment))

            # The stem following the ml-segment
            s = self.get_node_from_residue_num(f)
            if s == ml_segment:
                # The Cg consists of only a single f-element.
                assert len(self.defines)==1
                return None
            side_stem, _ = self._get_sides_plus(s, ml_segment)
            # Get the stem-side where we expect to find the next ML-segment
            if side_stem == 0:
                side_stem = 3
            elif side_stem == 3:
                side_stem = 0
            elif side_stem ==1:
                side_stem = 2
            elif side_stem == 2:
                side_stem = 1
            else:
                assert False

        log.debug("flanking_nuc_at_stem_side called for %s, side %s with defines %s.", s, side_stem, self.defines[s])
        ml_nuc = self.flanking_nuc_at_stem_side(s, side_stem)
        log.debug("ml_nucleotide is %s (sequence length is %s).", ml_nuc, self.seq_length)
        # End of the backbone
        if ml_nuc>self.seq_length or ml_nuc-1 in self.backbone_breaks_after:
            return None
        elem =  self.get_node_from_residue_num(ml_nuc)
        log.debug("side now %s, ml_nuc %s, ml %s", side_stem, ml_nuc, elem)
        if elem[0]=="s":
            #0-length multiloop
            elems=self.edges[elem]&self.edges[s]
            for elem in elems:
                if self.element_length(elem)==0:
                    return elem
            assert False
        if elem[0] not in "mft":
            self.log()
            log.error("%s is not a multiloop node", elem)
            return None
        return elem

    def shortest_mlonly_multiloop(self, ml_segment):
        loops = self.find_mlonly_multiloops()
        for loop in loops:
            if ml_segment in loop:
                return loop

    def flanking_nuc_at_stem_side(self, s, side):
        """
        Return the nucleotide number that is next to the stem at the given stem side.

        :param side: 0, 1, 2 or 3, as returned by self._get_sides_plus
        :returns: The nucleotide position. If the stem has no neighbor at that side,
                  0 or self.seq_length+1 is returned instead.
        """
        assert s[0]=="s", "{} is not a stem".format(s)
        stem_nuc = self.defines[s][side]
        if side == 0 or side ==2:
            return stem_nuc - 1
        else:
            return stem_nuc + 1

    def all_connections(self, elem):
        """
        Return the connected elements in order along the backbone.

        The difference to self.connections is that the returned list
        always contains as many elements, as the define of elem has numbers.
        If there is no connected element at this side,the returned list contains None.
        If elem is a stem connected to a hairpin or interior loop,
        this loop will be contained twice in the resulting output list.
        """
        connections = []
        # To correctly account for 0-length elements, we have to treat stems seperately.
        if elem[0]!="s":
            for nt in self.define_a(elem):
                neighbor = self.get_node_from_residue_num(nt)
                if neighbor == elem:
                    connections.append(None)
                else:
                    connections.append(neighbor)
        else:
            connections = [ None, None, None, None ]
            for neighbor in self.edges[elem]:
                for side in [0,1,2,3]:
                    if any(x==self.defines[elem][side] for x in self.define_a(neighbor)):
                        connections[side]=neighbor
        return connections

    def nucleotides_to_elements(self, nucleotides):
        """
        Convert a list of nucleotides (nucleotide numbers) to element names.

        Remove redundant entries and return a set.

        ..note::
            Use `self.get_node_from_residue_num` if you have only a single nucleotide number.
        """
        return set([self.get_node_from_residue_num(n) for n in nucleotides])



    def elements_to_nucleotides(self, elements):
        """
        Convert a list of element names to a list of nucleotide numbers.

        Remove redundant entries.

        """
        nucs = set()
        for elem in elements:
            for def_range in self.define_range_iterator(elem, adjacent = False):
                for nuc in range(def_range[0], def_range[1]+1):
                    nucs.add(nuc)
        return sorted(nucs)
    def find_bulge_loop(self, vertex, max_length=4):
        """
        Find a set of nodes that form a loop containing the
        given vertex and being no greater than max_length nodes long.

        :param vertex: The vertex to start the search from.
        :param max_length: Only fond loops that contain no more then this many elements
        :returns: A list of the nodes in the loop.
        """
        visited = set()
        to_visit = [(key, 1) for key in self.edges[vertex]]
        visited.add(vertex)
        in_path = [vertex]

        while len(to_visit) > 0:
            (current, depth) = to_visit.pop()
            visited.add(current)

            in_path = in_path[:depth]
            in_path.append(current)

            for key in self.edges[current]:
                if key == vertex and depth > 1:
                    if len(in_path[:depth + 1]) > max_length:
                        continue
                    else:
                        return in_path[:depth + 1]

                if key not in visited:
                    to_visit.append((key, depth + 1))
        return []

    def add_node(self, name, edges, define, weight=1):
        self.defines[name] = define
        self.edges[name] = edges
        self.weights[name] = weight

        for edge in self.edges[name]:
            self.edges[edge].add(name)
    def _remove_node(self, name):
        """
        Remove a node and all its (in- and out-) edges.

        .. warning::

            This does NOT update the defines/ cponnections of other nodes, and
            may result in the graph falling apart into several connected components.

        """
        for edge in list(self.edges[node]):
            self.remove_edge(*edge)
        del self.edges[node]
        del self.defines[node]

    def length_one_stem_basepairs(self):
        """
        Return a list of basepairs that correspond to length-1 stems.
        """
        # dissolve all stems which have a length of one
        stems_to_dissolve = [ s for s in self.stem_iterator() if self.stem_length(s)==1 ]
        if log.isEnabledFor(logging.WARNING):
            if len(stems_to_dissolve) == len(list(self.stem_iterator))!=0:
                log.warning("All stems of the structure have length 1!")
        log.info("Stems with length 1: %s", stems_to_dissolve)
        bps_to_dissolve = []
        for s in stems_to_dissolve:
            bps_to_dissolve.extend(self.stem_bp_iterator(s))
        return bps_to_dissolve

    def compare_stems(self, b):
        """
        A function that can be passed in as the key to a sort.
        """
        return (self.defines[b][0], 0)


    def compare_bulges(self, b, flank_nucs = False):
        """
        A function that can be passed in as the key to a sort.

        Compares based on the nucleotide number
        (using define_a to allow for sorting 0-length MLs)
        """
        return self.define_a(b)

    def compare_hairpins(self, b):
        connections = self.connections(b)

        return (self.defines[connections[0]][1], sys.maxsize)


    def has_connection(self, v1, v2):
        """ Is there an edge between these two nodes """

        if v2 in self.edges[v1]:
            return True
        else:
            # two multiloops can be connected at the end of a stem
            for e in self.edges[v1]:
                if e[0] != 's':
                    continue

                if v2 in self.edges[e]:
                    (s1b, s1e) = self.get_sides(e, v1)
                    (s2b, s2e) = self.get_sides(e, v2)

                    if s1b == s2b:
                        return True

            return False

    def connection_type(self, define, connections):
        """
        Classify the way that two stems are connected according to the type
        of bulge that separates them.

        Potential angle types for single stranded segments, and the ends of
        the stems they connect:

        =   = ======  ===========
        1   2 (1, 1)  #pseudoknot
        1   0 (1, 0)
        3   2 (0, 1)
        3   0 (0, 0)
        =   = ======  ===========

        :param define: The name of the bulge separating the two stems
        :param connections: The two stems and their separation

        :returns: INT connection type

                  =   ======================================================================
                  +   positive values mean forward (from the connected stem starting at the
                      lower nucleotide number to the one starting at the higher nuc. number)
                  -   negative values mean backwards.
                  1   interior loop
                  2   first multi-loop segment of normal multiloops and most pseudoknots
                  3   middle segment of a normal multiloop
                  4   last segment of normal multiloops and most pseudoknots
                  5   middle segments of pseudoknots
                  =   ======================================================================

        """
        if define[0] == 'i':
            # interior loop, we just have to check if
            # connections[0] < connections[1]
            if self.defines[connections[0]][0] < self.defines[connections[1]][0]:
                return 1
            else:
                return -1
        elif define[0] == 'm':
            (s1c, b1c) = self._get_sides_plus(connections[0], define)
            (s2c, b2c) = self._get_sides_plus(connections[1], define)

            if (s1c, s2c) == (1, 0):
                return 2
            elif (s1c, s2c) == (0, 1):
                return -2
            elif (s1c, s2c) == (3, 0):
                return 3
            elif (s1c, s2c) == (0, 3):
                return -3
            elif (s1c, s2c) == (2, 3):
                return 4
            elif (s1c, s2c) == (3, 2):
                return -4

            # the next two refer to pseudoknots
            elif (s1c, s2c) == (2, 1):
                return 5
            elif (s1c, s2c) == (1, 2):
                return -5
            else:
                raise GraphIntegrityError("Weird angle type: (s1c, s2c) = (%d, %d)" %
                                (s1c, s2c))
        else:
            raise ValueError("connection_type called on non-interior loop/multiloop")

    def connection_ends(self, connection_type):
        """
        Find out which ends of the stems are connected by a particular angle
        type.

        :param connection_type: The angle type, as determined by which corners
                                of a stem are connected
        :return: (s1e, s2b) 0 means the side of the stem with the lowest nucleotide, 1 the other side
        """
        ends = ()

        if abs(connection_type) == 1:
            ends = (1, 0)
        elif abs(connection_type) == 2:
            ends = (1, 0)
        elif abs(connection_type) == 3:
            ends = (0, 0)
        elif abs(connection_type) == 4:
            ends = (1, 0)
        elif abs(connection_type) == 5:
            ends = (1, 1)
        else:
            raise GraphIntegrityError('Unknown connection type: %d' % (connection_type))

        if connection_type < 0:
            return ends[::-1]
        else:
            return ends

    # This function seems to be unused. Consider deprecation...
    def get_multiloop_nucleotides(self, multiloop_loop):
        """
        Return a list of nucleotides which make up a particular
        multiloop.

        :param multiloop_loop: The elements which make up this multiloop
        :return: A list of nucleotides
        """
        stems = [d for d in multiloop_loop if d[0] == 's']
        multis = [d for d in multiloop_loop if d[0] == 'm']
        residues = []

        for s in stems:
            relevant_edges = [c for c in self.edges[s] if c in multiloop_loop]
            sides = [self._get_sides_plus(s, c)[0] for c in relevant_edges]
            sides.sort()

            # the whole stem is part of this multiloop
            if sides == [2, 3] or sides == [0, 1]:
                residues += list(range(self.defines[s][sides[0]], self.defines[s][sides[1]] + 1))
            else:
                residues += [self.defines[s][sides[0]], self.defines[s][sides[1]]]

        for m in multis:
            residues += self.define_residue_num_iterator(m, adjacent=False)
        return residues

    # This function seems to be unused. Consider deprecation...
    def find_external_loops(self):
        '''
        Return all of the elements which are part of
        an external loop.

        :return: A list containing the external loops in this molecule
                 (i.e. ['f0, m3, m5, t0'])
        '''
        ext_loop = []

        for d in it.chain(self.floop_iterator(),
                          self.tloop_iterator(),
                          self.mloop_iterator()):
            loop_nts = self.shortest_bg_loop(d)
            if len(loop_nts) == 0:
                ext_loop += [d]

        return ext_loop

    @profile
    def find_mlonly_multiloops(self):
        import networkx as nx
        ml_graph = nx.Graph()
        for d in it.chain(self.mloop_iterator(), self.floop_iterator(),self.tloop_iterator()):
            next_ml = self.get_next_ml_segment(d)
            if next_ml is not None:
                ml_graph.add_edge(d, next_ml)
            else:
                ml_graph.add_node(d)
        loops = []
        for comp in nx.connected_components(ml_graph):
            #Order along the cycle, in arbitrary direction.

            #We need to start at a node with only 1 connection, if present
            for x in comp:
                if len(ml_graph.edges(x))==1:
                    st_node=x
                    break
            else:
                st_node=x #Just take any node
            #Sort nodes along the cycle
            loop = list(nx.dfs_preorder_nodes(ml_graph.subgraph(comp), st_node))
            #See if we need to reverse the order
            for i,l in enumerate(loop):
                next_l = self.get_next_ml_segment(l)
                if i+1<len(loop):
                    if loop[i+1]==next_l:
                        break
                else:
                    if loop[0]==next_l:
                        break
            else:
                loop.reverse()
            #Find first node
            first = min(loop, key=lambda x: sorted(self.flanking_nucleotides(x)))
            first_i = loop.index(first)
            loop=loop[first_i:]+loop[:first_i]
            loops.append(tuple(loop))
        return loops


    def describe_multiloop(self, multiloop):
        """
        :param multiloop: An iterable of nodes (only "m", "t" and "f" elements)
        """
        descriptors = set()
        all_stems = col.Counter()
        angle_types = col.Counter()
        for elem in multiloop:
            if elem[0] in "ft":
                descriptors.add("open")
                continue
            elif elem[0] !="m":
                raise ValueError("Non multiloop element '{}' encountered in describe_multiloop.".foemat(elem))
            conn = self.connections(elem)
            ctype = abs(self.connection_type(elem, conn))
            angle_types[ctype] += 1
            if ctype == 5:
                descriptors.add("pseudoknot")
            all_stems.update(self.edges[elem])
        if sum(v % 2 for v in all_stems.values())==2: #Odd number of occurrences for 2 stems.
            descriptors.add("open")
        elif "open" not in descriptors:
            if any(v!=2 for v in all_stems.values()):
                print(all_stems)
                print(multiloop)
                print(self.to_dotbracket_string())
                print (self.to_element_string(True))
            assert sum(v % 2 for v in all_stems.values())==0
        if angle_types[2]==1 and angle_types[4]==1 and "pseudoknot" not in descriptors:
            descriptors.add("regular_multiloop")
        elif angle_types[2]>1 or angle_types[4]>1:
            descriptors.add("pseudoknot")
        elif "open" in descriptors and (angle_types[2]>0 or angle_types[4]>0):
            descriptors.add("pseudoknot")
        return descriptors

    def _from_dotbracket(self, dotbracket_str, dissolve_length_one_stems=False, remove_pseudoknots=False):
        """
        Called from the BulgeGraphs __init__ function.
        """
        pt = fus.dotbracket_to_pairtable(dotbracket_str)
        tuples = fus.pairtable_to_tuples(pt)
        self._from_tuples(tuples, dissolve_length_one_stems, remove_pseudoknots)

    def _from_tuples(self, tuples, dissolve_length_one_stems=False,
                    remove_pseudoknots=False):
        c = _BulgeGraphConstruction()
        c.from_tuples(tuples, self.seq.backbone_breaks_after)
        self._from_graph_construction(c)
        bps_to_remove = []
        if dissolve_length_one_stems:
            bps_to_remove.extend(self.length_one_stem_basepairs())
        if remove_pseudoknot:
            bps_to_remove.extend(self.pseudoknotted_basepairs(ignore_basepairs=bps_to_remove))
        if bps_to_remove:
            log.info("Recreating without the following "
                     "basepairs: %s", bps_to_remove)
            new_tuples = self.to_pair_tuples(bps_to_remove)
            self._from_tuples(new_tuples, False, False)

    def _from_graph_construction(self, c):
        """
        Copies defines and edges from A BulgeGraphConstruction instance
        to this graph.

        :param c: A BulgeGraphConstruction object after completed calculation.
        """
        self.defines = c.defines
        self.edges = d.edges

    def to_pair_table(self):
        """
        Create a pair table from the list of elements.

        The first element in the returned list indicates the number of
        nucleotides in the structure.

        i.e. [5,5,4,0,2,1]
        """
        pair_tuples = self.to_pair_tuples()

        return fus.tuples_to_pairtable(pair_tuples)

    def to_pair_tuples(self, remove_basepairs=[]):
        """
        Create a list of tuples corresponding to all of the base pairs in the
        structure. Unpaired bases will be shown as being paired with a
        nucleotide numbered 0.

        i.e. [(1,5),(2,4),(3,0),(4,2),(5,1)]
        """
        # iterate over each element
        table = []

        for d in self.defines:
            # iterate over each nucleotide in each element
            for b in self.define_residue_num_iterator(d):
                p = self.pairing_partner(b)
                if p is None:
                    p = 0
                table.append((b, p))


        if remove_basepairs:
            nt = []
            for p in table:
                to_add = p
                for s in remove_basepairs:
                    if sorted(p) == sorted(s):
                        to_add = (p[0], 0)
                        break
                nt += [to_add]
            table = nt
        return table

    def to_bpseq_string(self):
        """
        Create a bpseq string from this structure.
        """
        out_str = ''
        for i in range(1, self.seq_length + 1):
            pp = self.pairing_partner(i)
            if pp is None:
                pp = 0
            out_str += "{} {} {}\n".format(i, self.seq[i - 1], pp)

        return out_str

    def bpseq_to_tuples_and_seq(self, bpseq_str):
        """
        Convert a bpseq string to a list of pair tuples and a sequence
        dictionary. The return value is a tuple of the list of pair tuples
        and a sequence string.

        :param bpseq_str: The bpseq string
        :return: ([(1,5),(2,4),(3,0),(4,2),(5,1)], 'ACCAA')
        """
        lines = bpseq_str.split('\n')
        seq = []
        tuples = []
        pairing_partner={}
        for line in lines:
            parts = line.split()

            if len(parts) == 0:
                continue

            (t1, s, t2) = (int(parts[0]), parts[1], int(parts[2]))
            if t2 in pairing_partner and t1!=pairing_partner[t2]:
                raise GraphConstructionError("Faulty bpseq string. {} pairs with {}, "
                                 "but {} pairs with {}".format(t2, pairing_partner[t2], t1, t2))
            if t1 in pairing_partner and t2!=pairing_partner[t1]:
                raise GraphConstructionError("Faulty bpseq string. {} pairs with {}, "
                                 "but {} pairs with {}".format(pairing_partner[t1], t1,  t1, t2))

            pairing_partner[t1]=t2
            if t2!=0:
                pairing_partner[t2]=t1
            tuples += [(t1, t2)]
            seq += [s]

        seq = "".join(seq).upper().replace('T', 'U')

        return (tuples, seq)

    def from_bpseq_str(self, bpseq_str, dissolve_length_one_stems=False, breakpoints=[]):
        """
        Create the graph from a string listing the base pairs.

        The string should be formatted like so:

            1 G 115
            2 A 0
            3 A 0
            4 U 0
            5 U 112
            6 G 111

        :param bpseq_str: The string, containing newline characters.
        :return: Nothing, but fill out this structure.
        """
        self.__init__()
        log.debug(bpseq_str)
        #: This stores backbone breaks before they have been implemented!
        self._backbone_will_break_after = breakpoints
        tuples, seq = self.bpseq_to_tuples_and_seq(bpseq_str)

        self.seq = seq
        self.from_tuples(tuples)

        log.info("From bpseq_str: Secondary structure: %s", self.to_dotbracket_string())
        if dissolve_length_one_stems:
            self.dissolve_length_one_stems()

        self.seq_ids_from_seq()

    def  _zerolen_defines_a_between(self, stem1, stem2):
        zl_coordinates = set()
        for k, l in it.product(range(4), repeat=2):
            if abs(self.defines[stem1][k]-self.defines[stem2][l])==1:
                d = [self.defines[stem1][k], self.defines[stem2][l]]
                d.sort()
                log.debug("Zero-length element found: %s", d)
                if d[0] not in self.backbone_breaks_after:
                    zl_coordinates.add(tuple(d))
                else:
                    log.debug("But backbone-break encountered!")
        return zl_coordinates

    def to_dotbracket_string(self):
        """
        Convert the BulgeGraph representation to a dot-bracket string
        and return it.

        :return: A dot-bracket representation of this BulgeGraph
        """
        pt = self.to_pair_table()
        db_string = fus.pairtable_to_dotbracket(pt)
        for breakpoint in reversed(sorted(self.backbone_breaks_after)):
            db_string = db_string[:breakpoint]+"&"+db_string[breakpoint:]
        return db_string

    def log(self, level=logging.DEBUG):
        with log_at_caller(log):
            log.log(level, self.seq)
            log.log(level, self.to_dotbracket_string())
            es = self.to_element_string(with_numbers=True).split("\n")
            log.log(level, es[0])
            log.log(level, es[1])
            log.log(level, "DEFINES: %s", pformat(self.defines))
            log.log(level, "EDGES: %s", pformat(self.edges))

    def to_fasta_string(self):
        """
        Output the BulgeGraph representation as a fast string of the
        format::

            >id
            AACCCAA
            ((...))

        """
        output_string = ''

        output_string += ">%s\n" % (self.name)
        output_string += "%s\n" % (self.seq)
        output_string += "%s" % (self.to_dotbracket_string())

        return output_string

    def from_bg_file(self, bg_file):
        """
        Load from a file containing a text-based representation
        of this BulgeGraph.

        :param bg_file: The filename.
        :return: No return value since the current structure is the one
                 being loaded.
        """
        with open(bg_file, 'r') as f:
            bg_string = "".join(f.readlines())
            self.from_bg_string(bg_string)

    def from_bg_string(self, bg_str):
        """
        Populate this BulgeGraph from the string created by the method
        to_bg_string.

        :param bg_str: The string representation of this BugleGraph.
        """
        lines = bg_str.split('\n')
        for line in lines:
            line = line.strip()
            parts = line.split()
            if len(parts) == 0:
                # blank line
                continue
            if parts[0] == 'length':
                self.seq_length = int(parts[1])
            elif parts[0] == 'define':
                self.defines[parts[1]] = list(map(int, parts[2:]))
            elif parts[0] == 'connect':
                for p in parts[2:]:
                    self.edges[parts[1]].add(p)
                    self.edges[p].add(parts[1])
            elif parts[0] == 'seq':
                self.seq = parts[1]
                log.debug("from_bg_string: seq {}, breakpoints {}".format(self.seq, self.seq.backbone_breaks_after))
                self.backbone_breaks_after = self.seq.backbone_breaks_after
            elif parts[0] == 'seq_ids':
                self.seq_ids = list(map(resid_from_str, parts[1:]))
                self.chain_ids = []
                for res in self.seq_ids:
                    if res.chain not in self.chain_ids:
                        self.chain_ids.append(res.chain)
            elif parts[0] == 'name':
                self.name = parts[1].strip()
            elif parts[0] == 'info':
                self.infos[parts[1]].append(" ".join(parts[2:]))
        if not self.seq_ids:
            self.seq_ids_from_seq()

    def sorted_stem_iterator(self):
        """
        Iterate over a list of the stems sorted by the lowest numbered
        nucleotide in each stem.
        """
        stems = [d for d in self.defines if d[0] == 's']
        stems.sort(key=lambda s: self.defines[s][0])

        for s in stems:
            yield s

    def sorted_element_iterator(self):
        """
        Iterate over a list of the coarse grained elements sorted by the lowest numbered
        nucleotide in each stem. Multiloops with no nucleotide coordinates come last.
        """
        elements = [d for d in self.defines ]
        elements.sort(key=lambda s: self.defines[s][0] if self.defines[s] else 10000+int(s[1:]))

        for e in elements:
            yield e

    def is_single_stranded(self, node):
        """
        Does this node represent a single-stranded region?

        Single stranded regions are five-prime and three-prime unpaired
        regions, multiloops, and hairpins

        .. warning::
            Interior loops are never considered single stranded by this function.

        :param node: The name of the node
        :return: True if yes, False if no
        """
        if node[0] == 'f' or node[0] == 't' or node[0] == 'm' or node[0] == 'h':
            return True
        else:
            return False

    def get_node_dimensions(self, node):
        """
        Return the dimensions of a node.

        If the node is a stem, then the dimensions will be l where l is
        the length of the stem.

        Otherwise, see get_bulge_dimensions(node)

        :param node: The name of the node
        :return: A pair containing its dimensions
        """
        if node not in self.defines:
            self.log()
        if node[0] == 's':
            return (self.stem_length(node), self.stem_length(node))
            """
            return (self.defines[node][1] - self.defines[node][0] + 1,
                    self.defines[node][1] - self.defines[node][0] + 1)
            """
        else:
            return self.get_bulge_dimensions(node)

    def adjacent_stem_pairs_iterator(self):
        """
        Iterate over all pairs of stems which are separated by some element.

        This will always yield triples of the form (s1, e1, s2) where s1 and
        s2 are the stem identifiers and e1 denotes the element that separates
        them.
        """
        for d in self.defines.keys():
            if len(self.edges[d]) == 2:
                edges = list(self.edges[d])

                if edges[0][0] == 's' and edges[1][0] == 's':
                    yield (edges[0], d, edges[1])

    def stem_bp_iterator(self, stem):
        """
        Iterate over all the base pairs in the stem.
        """
        d = self.defines[stem]
        stem_length = self.stem_length(stem)

        for i in range(stem_length):
            yield (d[0] + i, d[3] - i)

    def get_connected_residues(self, s1, s2, bulge=None):
        """
        Get the nucleotides which are connected by the element separating
        s1 and s2. They should be adjacent stems.

        :param s1, s2: 2 adjacent stems
        :param bulge: Optional: The bulge seperating the two stems.
                      If s1 and s2 are connected by more than one element,
                      this has to be given, or a ValueError will be raised.
                      (useful for pseudoknots)

        The connected nucleotides are those which are spanned by a single
        interior loop or multiloop. In the case of an interior loop, this
        function will return a list of two tuples and in the case of multiloops
        if it will be a list of one tuple.

        If the two stems are not separated by a single element, then return
        an empty list.
        """

        c1 = self.edges[s1]
        c2 = self.edges[s2]

        # find out which edges they share
        common_edges = c1.intersection(c2)

        if len(common_edges) == 0:
            # not connected
            return []

        if len(common_edges) > 1 and bulge is None:
            raise ValueError("Too many connections between the stems. "
                            "Please provide the connectiong bulge you are interested in.")

        if bulge and bulge not in common_edges:
            raise ValueError("{}  does not connecty the stems {} and {}.".format(bulge, s1, s2))

        # the element linking the two stems
        if bulge:
            conn = bulge
        else:
            conn = list(common_edges)[0]

        # find out the sides of the stems that face the bulge
        (s1b, s1e) = self.get_sides(s1, conn)
        (s2b, s2e) = self.get_sides(s2, conn)

        # get the nucleotides on the side facing the stem
        s1_nucleotides = self.get_side_nucleotides(s1, s1b)
        s2_nucleotides = self.get_side_nucleotides(s2, s2b)

        # find out the distances between all the nucleotides flanking
        # the bulge
        dists = []
        for n1 in s1_nucleotides:
            for n2 in s2_nucleotides:
                dists += [(abs(n2 - n1), n1, n2)]
        dists.sort()

        # return the ones which are closest to each other first
        if conn[0] == 'i':
            return sorted([list(dists[0][1:]), list(dists[1][1:])])
        else:
            return [list(dists[0][1:])]

    def get_link_direction(self, stem1, stem2, bulge = None):
        """
        Get the direction in which stem1 and stem2 are linked (by the bulge)

        :returns: 1 if the bulge connects stem1 with stem2 in forward direction (5' to 3')
                  -1 otherwise
        """
        linked = self.get_connected_residues(stem1, stem2, bulge)
        if linked[0][0]<linked[0][1]:
            return 1
        else:
            return -1

    def get_side_nucleotides(self, stem, side):
        """
        Get the nucleotide numbers on the given side of
        them stem. Side 0 corresponds to the 5' end of the
        stem whereas as side 1 corresponds to the 3' side
        of the stem.

        :param stem: The name of the stem
        :param side: Either 0 or 1, indicating the 5' or 3' end of the stem
        :return: A tuple of the nucleotide numbers on the given side of
                 the stem.
        """
        if side == 0:
            return (self.defines[stem][0], self.defines[stem][3])
        elif side == 1:
            return (self.defines[stem][1], self.defines[stem][2])

        raise ValueError("Invalid side (%d) for the stem (%s)." % (stem, side))

    def get_stem_edge(self, stem, pos):
        """
        Returns the side of the stem that position is on.
        Side 0 corresponds to the 5' pairing residues in the
        stem whereas as side 1 corresponds to the 3' pairing
        residues in the stem.
        :param stem: The name of the stem
        :param pos: A position in the stem
        :return: 0 if pos on 5' edge of stem
        """
        fp_side = self.get_side_nucleotides(stem, 0)
        tp_side = self.get_side_nucleotides(stem, 1)

        fp_edge = range(fp_side[0],tp_side[0]+1)
        tp_edge = range(tp_side[1],fp_side[1]+1)

        if pos in fp_edge:
            return 0
        elif pos in tp_edge:
            return 1

        raise ValueError("Position (%d) not in stem (%s)." % (pos, stem))

    # Seems to be unused. Consider deprecation
    def get_any_sides(self, e1, e2):
        """
        Get the side of e1 that e2 is on. The only difference from the get_sides
        method is the fact that e1 does not have to be a stem.

        0 indicates that e2 is on the side with lower numbered
        nucleotides and 1 indicates that e2 is on the side with
        greater nucleotide numbers.

        :param e1: The name of the first element.
        :param e2: The name of the second element.
        :return: A tuple indicating the side of e1 adjacent to e2 and the side of e2
                 adjacent to e1
        """
        if e1[0] == 's':
            return self.get_sides(e1, e2)
        elif e2[0] == 's':
            return self.get_sides(e2, e1)[::-1]

        return None

    def get_sides(self, s1, b):
        """
        Get the side of s1 that is next to b.

        s1e -> s1b -> b

        :param s1: The stem.
        :param b: The bulge.
        :return: A tuple indicating which side is the one next to the bulge
                 and which is away from the bulge.
        """
        s1d = self.defines[s1]
        bd = self.defines[b]

        # Special case if the bulge is a length 0 multiloop
        if len(bd) == 0:
            bd = self._zero_length_element_adj_position(b)
            bd[0]+=1
            bd[1]-=1

        if bd[1]+1==s1d[0]:
            return (0,1)
        elif bd[1]+1==s1d[2]:
            return (1,0)
        elif s1d[1]+1==bd[0]:
            return (1,0)
        elif s1d[3]+1==bd[0]:
            return (0,1)
        else:
            raise GraphIntegrityError("Faulty bulge {}:{} connected to {}:{}".format(b, bd, s1, s1d))


    def stem_resn_to_stem_vres_side(self, stem, res):
        d = self.defines[stem]
        if res<=d[1]:
            assert res>= d[0]
            pos=res-d[0]
            side = 0
        elif res<=d[3]:
            assert res>=d[2]
            pos=d[3]-res
            side = 1
        else:
            raise ValueError("Residue {} not in stem {} with define {}".format(res, stem, d))
        return pos, side

    def stem_side_vres_to_resn(self, stem, side, vres):
        """
        Return the residue number given the stem name, the strand (side) it's on
        and the virtual residue number.
        """
        d = self.defines[stem]

        if side == 0:
            return d[0] + vres
        else:
            return d[3] - vres

    def stem_iterator(self):
        """
        Iterator over all of the stems in the structure.
        """
        for d in self.defines.keys():
            if d[0] == 's':
                yield d

    def hloop_iterator(self):
        """
        Iterator over all of the hairpin in the structure.
        """
        for d in self.defines.keys():
            if d[0] == 'h':
                yield d

    def mloop_iterator(self):
        """
        Iterator over all of the multiloops in the structure.
        """
        for d in self.defines.keys():
            if d[0] == 'm':
                yield d

    def iloop_iterator(self):
        """
        Iterator over all of the interior loops in the structure.
        """
        for d in self.defines.keys():
            if d[0] == 'i':
                yield d

    def floop_iterator(self):
        """
        Yield the name of the 5' prime unpaired region if it is
        present in the structure.
        """
        for d in self.defines.keys():
            if d[0] == 'f':
                yield d

    def tloop_iterator(self):
        """
        Yield the name of the 3' prime unpaired region if it is
        present in the structure.
        """
        for d in self.defines.keys():
            if d[0] == 't':
                yield d
    def pairing_partner(self, nucleotide_number):
        """
        Return the base pairing partner of the nucleotide at position
        nucleotide_number. If this nucleotide is unpaired, return None.

        :param nucleotide_number: The position of the query nucleotide in the
                                  sequence.
        :return: The number of the nucleotide base paired with the one at
                 position nucleotide_number.
        """
        for d in self.stem_iterator():
            for (r1, r2) in self.stem_bp_iterator(d):
                if r1 == nucleotide_number:
                    return r2
                elif r2 == nucleotide_number:
                    return r1
        return None

    def get_define_seq_str(self, elem, adjacent=False):
        """
        Get a list containing the sequences for the given define.

        :param d: The element name for which to get the sequences
        :param adjacent: Boolean. Include adjacent nucleotides (for single stranded RNA only)
        :return: A list containing the sequence(s) corresponding to the defines
        """
        if adjacent:
            define = self.define_a(elem)
        else:
            define = self.defines[elem]
        seqs=[]
        for i in range(0,len(define), 2):
            seqs.append(self.seq[define[i]:define[i+1]+1]) #seq is 1-based!
            if elem[0]=="i" and not adjacent:
                def_a = self.define_a(elem)
                if define[0]<def_a[1]:
                    seqs.append("")
                else:
                    seqs.insert(0,"")
        return seqs

    def get_stem_direction(self, s1, s2):
        """
        Return 0 if the lowest numbered residue in s1
        is lower than the lowest numbered residue in s2.
        """
        warnings.warn("get_stem_direction is deprecated and will be removed in the future!", stacklevel=2)
        if self.defines[s1][0] < self.defines[s2][0]:
            return 0
        return 1

    def get_multiloop_side(self, m):
        """
        Find out which strand a multiloop is on. An example of a situation in
        which the loop can be on both sides can be seen in the three-stemmed
        structure below:

            (.().().)

        In this case, the first multiloop section comes off of the 5' strand of
        the first stem (the prior stem is always the one with a lower numbered
        first residue). The second multiloop section comess of the 3' strand of
        the second stem and the third loop comes off the 3' strand of the third
        stem.
        """
        c = self.connections(m)

        p1 = self._get_sides_plus(c[0], m)
        p2 = self._get_sides_plus(c[1], m)

        return (p1[0], p2[0])

    # This function seems to be unused here. This code is possible duplicated somewhere.
    # Requires cleanup.
    def get_strand(self, multiloop):
        """
        Get the strand on which this multiloop is located.

        :param multiloop: The name of the multiloop
        :return: 0 for being on the lower numbered strand and 1 for
                 being on the higher numbered strand.
        """
        conn = self.connections(multiloop)
        t = self.connection_type(multiloop, conn)

        if abs(t) == 2:
            return 1
        elif abs(t) == 3:
            return 0
        else:
            return 2
        pass

    def get_bulge_dimensions(self, bulge):
        """
        Return the dimensions of the bulge.

        If it is single stranded it will be (x, -1) for h,t,f or (x, 1000) for m.
        Otherwise it will be (x, y).

        :param bulge: The name of the bulge.
        :return: A pair containing its dimensions
        """

        bd = self.defines[bulge]
        c = self.connections(bulge)

        if bulge[0] == 'i':
            # if this interior loop only has one unpaired region
            # then we have to find out if it's on the 5' strand or
            # the 3' strand
            # Example:
            # s1 1 3
            # 23 25
            # s2 5 10
            # 15 20
            s1 = self.defines[c[0]]
            s2 = self.defines[c[1]]

            dims = (s2[0] - s1[1] - 1, s1[2] - s2[3] - 1)

        if bulge[0] == 'm':
            # Multiloops are also pretty easy
            if len(bd) == 2:
                dims = (bd[1] - bd[0] + 1, 1000)
            else:
                dims = (0, 1000)
        if bulge[0] == 'f' or bulge[0] == 't':
            dims = (bd[1] - bd[0] + 1, -1)

        if bulge[0] == 'h':
            dims = (bd[1] - bd[0] + 1, -1)

        return dims


    def get_length(self, vertex):
        """
        Get the minimum length of a vertex.

        If it's a stem, then the result is its length (in base pairs).

        If it's a bulge, then the length is the smaller of it's dimensions.

        :param vertex: The name of the vertex.
        """
        if vertex[0] == 's':
            return abs(self.defines[vertex][1] - self.defines[vertex][0]) + 1
        else:
            if len(self.edges[vertex]) == 1:
                return self.defines[vertex][1] - self.defines[vertex][0] + 1
            else:
                dims = list(self.get_bulge_dimensions(vertex))
                dims.sort()

                if vertex[0] == 'i':
                    return sum(dims) / float(len(dims))

                else:
                    return min(dims)

    def get_flanking_region(self, bulge_name, side=0):
        """
        If a bulge is flanked by stems, return the lowest residue number
        of the previous stem and the highest residue number of the next
        stem.

        :param bulge_name: The name of the bulge
        :param side: The side of the bulge (indicating the strand)
        """
        c = self.connections(bulge_name)

        if bulge_name[0] == 'h':
            s1 = self.defines[c[0]]
            return (s1[0], s1[3])

        s1 = self.defines[c[0]]
        s2 = self.defines[c[1]]

        if bulge_name[0] == 'i':
            # interior loop
            if side == 0:
                return (s1[0], s2[1])
            else:
                return (s2[2], s1[3])

        elif bulge_name[0] == 'm':
            ss = self.get_multiloop_side(bulge_name)
            st = [s1, s2]

            ends = []

            # go through the two sides and stems and pick
            # the other end of the same strand
            for i, s in enumerate(ss):
                if s == 0:
                    ends += [st[i][1]]
                elif s == 1:
                    ends += [st[i][0]]
                elif s == 2:
                    ends += [st[i][3]]
                elif s == 3:
                    ends += [st[i][2]]
                else:
                    raise GraphIntegrityError("Weird multiloop sides: %s" %
                                    bulge_name)

            ends.sort()

            return tuple(ends)
            # multiloop

        return (None, None)


    def get_flanking_sequence(self, bulge_name, side=0):
        if len(self.seq) == 0:
            raise ValueError("No sequence present in the bulge_graph: %s" % (self.name))

        (m1, m2) = self.get_flanking_region(bulge_name, side)
        return self.seq[m1:m2+1] #1 based indexing

    def get_flanking_handles(self, bulge_name, side=0):
        """
        Get the indices of the residues for fitting bulge regions.

        So if there is a loop like so (between residues 7 and 16)::

          (((...))))
          7890123456
            ^   ^

        Then residues 9 and 13 will be used as the handles against which
        to align the fitted region.

        In the fitted region, the residues (2,6) will be the ones that will
        be aligned to the handles.

        :return: (orig_chain_res1, orig_chain_res1, flanking_res1, flanking_res2)
        """
        f1 = self.get_flanking_region(bulge_name, side)
        c = self.connections(bulge_name)

        if bulge_name[0] == 'h':
            s1 = self.defines[c[0]]
            ab = [s1[1], s1[2]]
            return (ab[0], ab[1], ab[0] - f1[0], ab[1] - f1[0])

        s1 = self.defines[c[0]]
        s2 = self.defines[c[1]]

        if bulge_name[0] == 'm':
            sides = self.get_multiloop_side(bulge_name)
            ab = [s1[sides[0]], s2[sides[1]]]
            ab.sort()

            return (ab[0], ab[1], ab[0] - f1[0], ab[1] - f1[0])

        if bulge_name[0] == 'i':
            if side == 0:
                ab = [s1[1], s2[0]]
            else:
                ab = [s2[3], s1[2]]

            return (ab[0], ab[1], ab[0] - f1[0], ab[1] - f1[0])

        # probably still have to include the 5' and 3' regions, but that
        # will come a little later
        return None

    def are_adjacent_stems(self, s1, s2, multiloops_count=True):
        """
        Are two stems separated by only one element. If multiloops should not
        count as edges, then the appropriate parameter should be set.

        :param s1: The name of the first stem
        :param s2: The name of the second stem
        :param multiloops_count: Whether to count multiloops as an edge linking
                                 two stems
        """
        for e in self.edges[s1]:
            if not multiloops_count and e[0] == 'm':
                continue
            if s2 in self.edges[e]:
                return True

        return False

    def random_subgraph(self, subgraph_length=None):
        """
        Return a random subgraph of this graph.

        :return: A list containing a the nodes comprising a random subgraph
        """
        if subgraph_length == None:
            subgraph_length = random.randint(1, len(self.defines.keys()))

        start_node = random.choice(list(self.defines.keys()))
        curr_length = 0
        visited = set()
        next_nodes = [start_node]
        new_graph = []

        while curr_length < subgraph_length:
            curr_node = random.choice(next_nodes)
            if curr_node[0] == 'i' or curr_node[0] == 'm':
                # if it's an interior loop or a multiloop, then we have to
                # add the adjacent stems
                for e in self.edges[curr_node]:
                    if e in new_graph:
                        continue
                    visited.add(e)
                    new_graph += [e]
                    next_nodes += list(self.edges[e])
                    curr_length += 1

            visited.add(curr_node)
            next_nodes += list(self.edges[curr_node])
            next_nodes = [n for n in next_nodes if n not in visited]
            new_graph += [curr_node]
            curr_length += 1  # self.element_length(curr_node)

        return new_graph

    def get_resseqs(self, define, seq_ids=True):
        """
        Return the pdb ids of the nucleotides in this define.

        :param define: The name of this element.
        :param: Return a tuple of two arrays containing the residue ids
                on each strand
        """
        resnames = []
        ranges = zip(*[iter(self.defines[define])] * 2)

        for r in ranges:
            strand_resnames = []
            for x in range(r[0], r[1] + 1):
                if seq_ids:
                    try:
                        res_id = self.seq_ids[x - 1]
                    except IndexError as e:
                        with log_to_exception(log, e):
                            log.error("Index %s not in seq_ids.", (x-1))
                        raise
                    if hasattr(self, "chain") and self.chain is not None:
                        assert res_id in self.chain
                    strand_resnames.append(res_id)
                else:
                    strand_resnames += [x]

            resnames += [strand_resnames]

        return resnames

    def _insert_cutpoints_into_seq(self):
        if self.seq:
            for breakpoint in self.backbone_breaks_after:
                log.debug("Inserting breakpoint into seq '{}'".format(self.seq))
                self.seq = self.seq.subseq_with_cutpoints(1,breakpoint+1)+"&"+self.seq.subseq_with_cutpoints(breakpoint+1, None)
                log.info("seq now has {} cutpoints".format(self.seq.count('&')))



    # This function seems to be dead code, but might be useful in the future.
    # Consider adding this to whitelist.py
    def connected_stem_iterator(self):
        """
        Iterate over all pairs of connected stems.
        """
        for l in it.chain(self.mloop_iterator(), self.iloop_iterator()):
            edge_list = list(self.edges[l])
            yield (edge_list[0], l, edge_list[1])

    def sorted_edges_for_mst(self):
        """
        Keep track of all linked nodes. Used for the generation of the minimal spanning tree.
        """
        priority = {'s': 1, 'i': 2, 'm': 3, 'f': 4, 't': 5}
        edges = sorted(it.chain(self.mloop_iterator(),
                                self.iloop_iterator()),
                       key=lambda x: (priority[x[0]], min(self.get_node_dimensions(x)),x) )
        return edges
    def get_mst(self):
        """
        Create a minimum spanning tree from this BulgeGraph. This is useful
        for constructing a structure where each section of a multiloop is
        sampled independently and we want to introduce a break at the largest
        multiloop section.
        """
        # keep track of all linked nodes
        edges = self.sorted_edges_for_mst()

        mst = set(it.chain(self.stem_iterator(),
                           self.floop_iterator(),
                           self.tloop_iterator(),
                           self.hloop_iterator()))

        # store all of the disconnected trees
        forest = [set([m]) for m in mst]

        # get the tree containing a particular element
        def get_tree(elem):
            for t in forest:
                if elem in t:
                    return t

        while len(edges) > 0:
            conn = edges.pop(0)
            neighbors = list(self.edges[conn])

            # get the trees containing the neighbors of this node
            # the node should be an interior loop or multiloop so
            # the neighbors should necessarily be stems, 5' or 3'
            t1 = get_tree(neighbors[0])
            t2 = get_tree(neighbors[1])

            if len(set.intersection(t1, t2)) == 0:
                # if this node connects two disparate trees, then add it to the mst
                new_tree = t1.union(t2)
                forest.remove(t1)
                forest.remove(t2)
                forest.append(new_tree)

                mst.add(conn)

        return mst

    def traverse_graph(self):
        """
        Traverse the graph to get the angle types. The angle type depends on
        which corners of the stem are connected by the multiloop or internal
        loop.

        :returns: A list of triples (stem, loop, stem)
        """
        if self.mst is None:
            self.mst = self.get_mst()

        build_order = []
        to_visit = [('s0', 'start')]
        visited = set(['s0'])
        build_paths = col.defaultdict(list)

        while len(to_visit) > 0:
            to_visit.sort(key=lambda x: min(self.get_node_dimensions(x[0])))
            (current, prev) = to_visit.pop(0)

            for e in self.edges[current]:
                if e not in visited and e in self.mst:
                    # make sure the node hasn't been visited
                    # and is in the minimum spanning tree
                    to_visit.append((e, current))

                    build_paths[e] += [e]
                    build_paths[e] += build_paths[current]

                    visited.add(e)

            if current[0] != 's' and len(self.edges[current]) == 2:
                # multiloop or interior loop

                # overkill method of getting the stem that isn't
                # equal to prev
                next_stem = set.difference(self.edges[current],
                                           set([prev]))
                build_order += [(prev, current, list(next_stem)[0])]
                # If pseudoknots exist, the direction is not always 0!
                # assert self.get_stem_direction(prev, build_order[-1][2])==0
		# does not hold for pseudoknots!
        self.build_order = build_order
        self.ang_type = None

        return build_order

    def set_angle_types(self):
        """
        Fill in the angle types based on the build order
        """
        if self.build_order is None:
            self.traverse_graph()

        self.ang_types = dict()
        for (s1, b, s2) in self.build_order:
            self.ang_types[b] = self.connection_type(b, [s1, s2])

    def get_angle_type(self, bulge, allow_broken = False):
        """
        Return what type of angle this bulge is, based on the way this
        would be built using a breadth-first traversal along the minimum
        spanning tree.

        :param allow_broken: How to treat broken multiloop segments.
                             * False (default): Return None
                             * The string "bo" or "build_order": Return the
                               angle type according to the build-order
                               (i.e. from the first built stem to the last-built stem)
                             * True: Return the angle_type from the stem with
                               lower nt number to the stem with higher nt number.
                               In forgi 2.0 this will be removed and the behavior of "bo"
                               will be used instead.
        """
        if self.ang_types is None:
            self.set_angle_types()

        if bulge in self.ang_types:
            return self.ang_types[bulge]
        else:
            if allow_broken == "bo" or allow_broken ==" build_order":
                stems = self.connections[bulge]
                s1, s2 = sorted(stems, key=lambda x: self.buildorder_of(x))
                return self.connection_type(bulge, [s1, s2])
            elif allow_broken:
                warnings.warn("The behavior of 'allow_broken=True' will change "
                              "to reflect the behavior of 'allow_broken=\"bo\"' "
                              "with forgi version 2.0", DeprecationWarning)
                s1, s2 = self.connections(bulge) #Ordered by nucleotide number
                return self.connection_type(bulge, [s1, s2])
            else:
                return None

    def buildorder_of(self, element):
        """
        Returns the index into build_order where the element FIRST appears.

        :param element: Element name, a string. e.g. "m0" or "s0"
        :returns: An index into self.build_order or None, if the element is not
                  part of the build_order (e.g. hairpin loops)
        """
        if self.build_order is None:
            self.traverse_graph()
        for i, elements in enumerate(self.build_order):
            if element in elements:
                return i
        return None


    def is_loop_pseudoknot(self, loop):
        """
        Is a particular loop a pseudoknot?

        :param loop: A list of elements that are part of the loop (only m,f and t elements).

        :return: Either True or false
        """
        return "pseudoknot" in self.describe_multiloop(loop)

    def iter_elements_along_backbone(self, startpos = 1):
        """
        Iterate all coarse grained elements along the backbone.

        Note that stems are yielded twice (for forward and backward strand).
        Interior loops may be yielded twice or once (if one side has no nucleotide)

        :param startpos: The nucleotide position at which tio start
        :yields: Coarse grained element names, like "s0", "i0"
        """
        nuc = startpos
        try:
            node = self.get_node_from_residue_num(nuc)
        except:
            assert len(self.defines)==0
            raise StopIteration("Empty Graph")
        while True:
            yield node
            if node[0] in "si": #The strand matters
                if node[0]=="s":
                    strand = self.get_stem_edge(node, nuc)
                    if strand == 0: #forward
                        nuc = self.flanking_nuc_at_stem_side(node, 1)
                    else:
                        nuc = self.flanking_nuc_at_stem_side(node, 3)
                else:
                    f1,f2,f3,f4 = sorted(self.flanking_nucleotides(node))
                    if f1<nuc<f2:
                        nuc = f2
                    elif f3<nuc<f4:
                        nuc = f4
                    else:
                        assert False
                try:
                    next_node = self.get_node_from_residue_num(nuc)
                except LookupError:
                    raise StopIteration("End of chain reached")
                else:
                    #We need to make sure there is no 0-length multiloop between the two stems.
                    intersect = self.edges[node] & self.edges[next_node]
                    for el in intersect:
                        if el[0]=="m" and self.defines[el]==[]:
                            #In case of this structuire ([)], there are 2 0-length multiloops between the two stems.
                            prev_nuc = min(self.flanking_nucleotides(el))
                            if self.get_node_from_residue_num(prev_nuc) == node:
                                node = el
                                break
                    else:
                        node = next_node

            else:
                try:
                    f1, f2 = self.flanking_nucleotides(node)
                except ValueError as e: #Too few values to unpack
                    if node[0]=="f":
                        if len(self.defines)==1:
                            raise StopIteration("Single stranded-only RNA")
                        nuc, = self.flanking_nucleotides(node)
                    else:
                        raise StopIteration("End of chain reached")
                else:
                    if f1>f2:
                        nuc=f1
                    else:
                        nuc=f2
                log.debug("Next nuc is {} ({})".format(nuc, repr(nuc)))
                node =  self.get_node_from_residue_num(nuc)



    def to_networkx(self):
        """
        Convert this graph to a networkx representation. This representation
        will contain all of the nucleotides as nodes and all of the base pairs
        as edges as well as the adjacent nucleotides.
        """
        import networkx as nx

        G = nx.Graph()

        residues = []
        for d in self.defines:
            prev = None

            for r in self.define_residue_num_iterator(d):
                G.add_node(r)
                residues += [r]

        #Add links along the backbone
        residues.sort()
        prev = None
        for r in residues:
            if prev is not None:
                G.add_edge(prev, r)
            prev = r

        #Add links along basepairs
        for s in self.stem_iterator():
            for (f, t) in self.stem_bp_iterator(s):
                G.add_edge(f, t)

        return G

    def pseudoknotted_basepairs(self, ignore_basepairs=[]):
        """
        Return a list of base-pairs that will be removed to
        remove pseudoknots using the knotted2nested.py script.

        :param ignore_basepairs: An optional list of basepairs that
                                 knested2knotted will not consider present
                                 in the structure.
        :return: A list of base-pairs that can be removed.
        """
        # remove unpaired bases and redundant pairs (i.e. (2,3) and (3,2))
        pairs = sorted([tuple(sorted(p)) for p in self.to_pair_tuples() if p[1] != 0])
        pairs = set(pairs)
        pairs = [ p  for p in pairs if p not in ignore_basepairs
                                       and (p[1], p[0] not in ignore_basepairs)]

        import forgi._k2n_standalone.knots as fakk

        pk_function = fakk.eg
        nested_pairs, removed_pairs = pk_function(pairs, return_removed=True)

        return removed_pairs

    #Seems to be unused...
    def ss_distance(self, e1, e2):
        '''
        Calculate the distance between two elements (e1, e2)
        along the secondary structure. The distance only starts
        at the edge of each element, and is the closest distance
        between the two elements.

        :param e1: The name of the first element
        :param e2: The name of the second element
        :return: The integer distance between the two along the secondary
                 structure.
        '''
        # get the edge nucleotides
        # thanks to:
        # http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
        # we get the edges, except that they might be one too close because we use adjacent
        # nucleotides, nevertheless we'll take care of that later
        d1_corners = []
        d2_corners = []

        for key, group in it.groupby(enumerate(self.define_residue_num_iterator(e1, adjacent=True)),
                                  lambda index_item: index_item[0] - index_item[1]):
            group = list(map(oper.itemgetter(1), group))
            d1_corners += group

        for key, group in it.groupby(enumerate(self.define_residue_num_iterator(e2, adjacent=True)),
                                  lambda index_item: index_item[0] - index_item[1]):
            group = list(map(oper.itemgetter(1), group))
            d2_corners += group

        import networkx as nx

        G = self.to_networkx()
        path_lengths = []
        for c1, c2 in it.product(d1_corners, d2_corners):
            path_lengths += [nx.shortest_path_length(G, c1, c2)]

        if e1 == e2:
            return 0

        if e1 in self.edges[e2]:
            return min(path_lengths) + 1

        # make some exceptions for edges which have length 0
        common_edges = set.intersection(self.edges[e1], self.edges[e2])
        for e in common_edges:
            if e[0] == 'i' and len(self.defines[e]) < 4:
                return min(path_lengths) + 1
            elif e[0] == 'm' and len(self.defines[e]) < 2:
                return min(path_lengths) + 1

        return min(path_lengths) + 2

    def shortest_path(self, e1, e2):
        '''
        Determine the shortest path between two elements (e1, e2)
        along the secondary structure.

        :param e1: The name of the first element
        :param e2: The name of the second element
        :return: A list of the element names along the shortest path

        '''

        import networkx as nx

        # Get residue numbers of source and targets, for shortest_path in nx
        source = min( [res for res in self.define_residue_num_iterator(e1)] )
        target = min( [res for res in self.define_residue_num_iterator(e2)] )

        # Get nx graph, and the shortest path
        G = self.to_networkx()
        nx_sp = nx.shortest_path(G, source=source, target=target)

        # Convert shortest path of residue numbers to a shortest path of node names
        sp, sp_set = [], set() # Use set to keep track of additions for faster lookup
        for res in nx_sp:
            node = self.get_node_from_residue_num(res)
            if node not in sp_set:
                sp_set.add(node)
                sp.append(node)

        # assymetric bulges with a length of 0 on 1 side are missed,
        # two adjacent stems indicate a bulge with length 0 along the path
        shortest_path, sp_set = [], set()
        traversal = self.traverse_graph() # Connections are ordered compared to connected_stem_iterator()

        for n1, n2 in zip(sp, sp[1:]): # Iterate through adjacent pairs of elements in the list
            if n1.startswith('s') and n2.startswith('s'): # If two elements are both stems
                connection = list([conn for conn in traversal if n1 in conn and n2 in conn][0]) # Find their connection in graph traversal
                if connection.index(n1) > connection.index(n2): #If we're moving 'backwards' on the traversal
                    connection.reverse()
                for node in connection:
                    if node not in sp_set:
                        sp_set.add(node)
                        shortest_path.append(node)
            else:
                if n1 not in sp_set:
                    sp_set.add(n1)
                    shortest_path.append(n1)
        if n2 not in sp_set:
            shortest_path.append(n2) # Append last item in path

        return shortest_path

    def get_position_in_element(self, resnum):
        """
        Return the position of the residue in the cg-element and the length of the element.

        :param resnum: An integer. The 1-based position in the total sequence.
        :returns: A tuple (p,l) where p is the position of the residue in the cg-element
                  (0-based for stems, 1-based for loops) and p/l gives a measure for the position
                  of the residue along the cg-element's axis (0 means at cg.coords[elem][0],
                  1 at cg.coords[elem][1] and 0.5 exactely in the middle of these two. )
        """
        node = self.get_node_from_residue_num(resnum)

        if node[0] == 's':
            if self.defines[node][0] <= resnum <= self.defines[node][1]:
                return resnum - self.defines[node][0], self.defines[node][1] - self.defines[node][0]
            else:
                return abs(resnum - self.defines[node][3]), self.defines[node][1] - self.defines[node][0]
        elif node[0] == 'i':
            s0,s1 = self.connections(node)
            if self.defines[s0][1] <= resnum <= self.defines[s1][0]:
                return resnum - self.defines[s0][1], self.defines[s1][0] - self.defines[s0][1]
            else:
                return abs(resnum - self.defines[s0][2]) - 1, self.defines[s0][2] - self.defines[s1][3]
        elif node[0] == 'h':
            pos1 = resnum - self.defines[node][0]
            pos2 = abs(resnum - self.defines[node][1])

            return min(pos1, pos2)+1, (self.defines[node][1] - self.defines[node][0] + 2) // 2


        i = 0
        while i < len(self.defines[node]):
            s = self.defines[node][i]
            e = self.defines[node][i+1]

            if s <= resnum <= e:
                return resnum - s+1, e - s + 2

            i += 2

        return None

    def connected(self, n1, n2):
        '''
        Are the nucleotides n1 and n2 connected?

        :param n1: A node in the BulgeGraph
        :param n2: Another node in the BulgeGraph
        :return: True or False indicating whether they are connected.
        '''
        if n1 in self.edges[n2] or n2 in self.edges[n1]:
            return True

        # two multiloops can be considered connected if they both
        # link to the same side of the same stem
        if n1[0] == 'm' and n2[0] == 'm':
            common_stems = list(set.intersection(self.edges[n1], self.edges[n2]))
            if len(common_stems) == 0:
                return False

            common_stem = common_stems[0]

            (s1c, b1c) = self._get_sides_plus(common_stem, n1)
            (s2c, b1c) = self._get_sides_plus(common_stem, n2)

            if sorted([s1c, s2c]) == [0,3] or sorted([s1c, s2c]) == [1,2]:
                return True

        return False

    def flanking_nucleotides(self, d):
        '''
        Return the nucleotides directly flanking an element.

        :param d: the name of the element
        :return: a list of nucleotides
        '''
        set_adjacent = set(self.define_residue_num_iterator(d, adjacent=True))
        set_not_adjacent = set(self.define_residue_num_iterator(d, adjacent=False))

        return list(sorted(set_adjacent - set_not_adjacent))

    def min_max_bp_distance(self, e1, e2):
        '''
        Get the minimum and maximum base pair distance between
        these two elements.

        If they are connected, the minimum distance will be 1.
        The maximum will be 1 + length(e1) + length(e1)

        :param e1: The name of the first element
        :param e2: The name of the second element
        :return:   A tuple containing the minimum and maximum distance between
                   the two elements.
        '''

        if (e1,e2) in self._elem_bp_dists: #Shortcut if cached.
            return self._elem_bp_dists[(e1,e2)]

        min_bp = sys.maxsize
        max_bp = 0

        if self.nx_graph is None:
            self.nx_graph = self.to_networkx()

        if self.nuc_bp_dists is None:
            import networkx as nx
            self.nuc_bp_dists = np.zeros((self.seq_length+1, self.seq_length+1))#col.defaultdict(dict)
            self.nuc_bp_dists[0,:]=np.nan
            self.nuc_bp_dists[:,0]=np.nan
            dist_matrix = np.array(nx.floyd_warshall_numpy(self.nx_graph))
            for (i1,n1), (i2,n2) in it.product(enumerate(self.nx_graph.nodes()),
                                               enumerate(self.nx_graph.nodes())):
                self.nuc_bp_dists[n1,n2] = dist_matrix[i1][i2]

        for f1, f2 in it.product(set(self.defines[e1]), set(self.defines[e2])):
            #d =  nx.dijkstra_path_length(self.nx_graph, f1, f2)
            d = self.nuc_bp_dists[f1,f2]

            if d < min_bp:
                min_bp = d
            if d > max_bp:
                max_bp = d
        self._elem_bp_dists[(e1,e2)] = (min_bp, max_bp)
        self._elem_bp_dists[(e2,e1)] = (min_bp, max_bp)
        return (min_bp, max_bp)

    def nd_define_iterator(self):
        '''
        Iterate over defines which contain some nucleotides.

        :return: An iterator over all defines which contain some
                 nucleotides.
        '''
        for d in self.defines:
            if len(self.defines[d]) > 0:
                yield d

    def get_domains(self):
        """
        Get secondary structure domains.

        Currently domains found are:
          * multiloops (without any connected stems)
          * rods: stretches of stems + interior loops (without branching), with trailing hairpins
          * pseudoknots
        """
        domains = col.defaultdict(list)
        multiloops = self.find_mlonly_multiloops()
        for ml in multiloops:
            ml = sorted(ml)
            if self.is_loop_pseudoknot(ml):
                domains["pseudoknots"].append(ml)
            else:
                domains["multiloops"].append(ml)

        doublestr = []
        for s in self.stem_iterator():
            neighbors = self.edges[s]
            for region in doublestr:
                if s in region or any(n in region for n in neighbors):
                    curr_region = region
                    curr_region.add(s)
                    break
            else:
                doublestr.append(set([s]))
                curr_region = doublestr[-1]

            for n in neighbors:
                if n[0] in "sih":
                    curr_region.add(n)
        #print(doublestr)
        while True:
            for reg1, reg2 in it.combinations(doublestr,2):
                if reg1 & reg2:
                    doublestr.remove(reg1)
                    doublestr.remove(reg2)
                    doublestr.append(reg1|reg2)
                    break
            else:
                break

        for region in doublestr:
            domains["rods"].append(sorted(region))
        domains["pseudoknots"].sort()
        domains["multiloops"].sort()
        domains["rods"].sort()
        #print(domains)
        return domains
