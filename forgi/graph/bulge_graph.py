#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)
"""bulge_graph.py: A graph representation of RNA secondary structure based
   on its decomposition into primitive structure types: stems, hairpins,
   interior loops, multiloops, etc..."""

__author__ = "Peter Kerpedjiev, Bernhard Thiel"
__copyright__ = "Copyright 2012 - 2016"
__version__ = "0.3"
__maintainer__ = "Peter Kerpedjiev, Bernhard Thiel"
__email__ = "pkerp@tbi.univie.ac.at, thiel@tbi.univie.ac.at"

import sys
import collections as col
import random
import re
import itertools as it
from ..aux.k2n_standalone import knotted2nested as fak
from ..utilities import debug as fud
from ..utilities import stuff as fus
from ..threedee.utilities import mcannotate as ftum
import os, warnings
import operator as oper
import numpy as np

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
    bg = BulgeGraph()
    bg.from_dotbracket(struct)
    bg.name = id_str
    bg.seq = seq

    return bg


def from_fasta_text(fasta_text):
    """
    Create a bulge graph or multiple bulge
    graphs from some fasta text.
    """
    # compile searches for the fasta id, sequence and 
    # secondary structure respectively
    id_search = re.compile('>(.+)')
    seq_search = re.compile('^([acguACGU]+)$') #BT: This does not allow for "t"/"T"(DNA)/"N". Is it on purpose?

    prev_id = None
    prev_seq = None
    prev_struct = None
    curr_id=None

    bgs = []

    for line in fasta_text.split('\n'):
        # newlines suck
        line = line.strip()

        # find out what this line contains
        id_match = id_search.match(line)
        seq_match = seq_search.match(line)

        if id_match is not None:
            prev_id=curr_id
            # we found an id, check if there's a previous
            # sequence and structure, and create a BG
            curr_id = id_match.group(0).strip('>')

            if prev_seq is None and prev_struct is None:
                # must be the first sequence/structure
                continue

            # make sure we have
            if prev_seq is None:
                raise Exception("No sequence for id: {}", prev_id)
            if prev_struct is None:
                raise Exception("No sequence for id: {}", prev_id) 
                #BT: This message not very helpful, if wrong character ("N"/..) in sequence
            if prev_id is None:
                raise Exception("No previous id")

            bgs += [from_id_seq_struct(prev_id, prev_seq, prev_struct)]

        if seq_match is not None:
            prev_seq = seq_match.group(0)
        if id_match is None and seq_match is None:
            if len(line) > 0:
                prev_struct = line

    bgs += [from_id_seq_struct(curr_id, prev_seq, prev_struct)]

    if len(bgs) == 1:
        return bgs[0]
    else:
        return bgs


def from_fasta(filename):
    """
    Load a bulge graph from a fasta file. The format of the fasta
    file is roughly:

        >1
        AACCCAA
        ((...))
    """
    with open(filename, 'r') as f:
        text = f.read()
        bg = BulgeGraph()
        bg.from_fasta(text)
        return bg


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

    :param brackets: A string with the dotplot passed as input to this script.
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

class BulgeGraph(object):
    def __init__(self, bg_file=None, dotbracket_str='', seq=''):
        """
        A bulge graph object.

        :var self.defines: The coarse grain element definitions: Keys are for example 's1'/ 'm2'/ 'h3'/ 'f1'/ 't1'
                       Values are the positions in the sequence (1D-coordinate) of start , end, ...
        """
        self.seq_length = 0
        self.ang_types = None
        self.mst = None
        self.build_order = None
        self.name = "untitled"
        #: The coarse grain element definitions: Keys are for example 's1'/ 'm2'/ 'h3'/ 'f1'/ 't1'
        #: Values are the positions in the sequence (1D-coordinate) of start , end, ...
        self.defines = dict()
        self.edges = col.defaultdict(set)
        self.longrange = col.defaultdict(set)
        self.weights = dict()
        self.nx_graph = None
        self.nuc_bp_dists = None
        self._elem_bp_dists = {}

        # store the coordinate basis for each stem
        self.bases = dict()
        self.stem_invs = dict()
        self.seq_ids = []

        # Additional infos as key-value pairs are stored here.
        self.infos = col.defaultdict(list)

        self.name_counter = 0

        if dotbracket_str != '':
            self.from_dotbracket(dotbracket_str)

        self.seq = seq
        for i, s in enumerate(seq):
            self.seq_ids += [(' ', str(i + 1), ' ')]

        if bg_file is not None:
            self.from_bg_file(bg_file)

    # get an internal index for a named vertex
    # this applies to both stems and edges
    def get_vertex(self, name=None):
        """
        Return a new unique vertex name.
        """

        if name is None:
            name = "x{}".format(self.name_counter)
            self.name_counter += 1

        return name

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

                for dest in self.edges[key]:
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
        out_str += " ".join(map(ftum.format_resid, self.seq_ids))
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

    def to_element_string(self):
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
        """
        output_str = [' '] * (self.seq_length + 1)

        for d in self.defines.keys():
            for resi in self.define_residue_num_iterator(d, adjacent=False):
                output_str[resi] = d[0]

        return "".join(output_str).strip()

    def define_range_iterator(self, node, adjacent=False, seq_ids=False):
        """
        Return the ranges of the nucleotides in the define.

        In other words, if a define contains the following: [1,2,7,8]
        The ranges will be [1,2] and [7,8].

        :param adjacent: Use the nucleotides in the neighboring element which
                         connect to this element as the range starts and ends.
        :return: A list of two-element lists
        """
        a = iter(self.defines[node])
        ranges = zip(a, a)

        if node[0] == 'i':
            # interior loops have to be treated specially because
            # they might have a bulge that has no unpaired nucleotides on one strand

            if adjacent:
                conns = self.connections(node)
                s1 = self.defines[conns[0]]
                s2 = self.defines[conns[1]]

                # offset by one, which will be reversed in the yield step
                # below
                ranges = [[s1[1] + 1, s2[0] - 1], [s2[3] + 1, s1[2] - 1]]

        if node[0] == 'm':
            if adjacent:
                conns = self.connections(node)
                s1 = self.get_sides_plus(conns[0], node)[0]
                s2 = self.get_sides_plus(conns[1], node)[0]

                rnge = sorted([self.defines[conns[0]][s1],
                               self.defines[conns[1]][s2]])
                ranges = [[rnge[0] + 1, rnge[1] - 1]]

        for (ds1, ds2) in ranges:
            if adjacent:
                if ds1 > 1:
                    ds1 -= 1
                if ds2 < self.seq_length:
                    ds2 += 1

            if seq_ids:
                # this will cause problems if the nucleotide has insertion codes
                yield [self.seq_ids[ds1 - 1], self.seq_ids[ds2 - 1]]
            else:
                yield [ds1, ds2]

    @profile
    def define_residue_num_iterator(self, node, adjacent=False, seq_ids=False):
        """
        Iterate over the residue numbers that belong to this node.

        :param node: The name of the node
        """
        visited=set()

        for r in self.define_range_iterator(node, adjacent, seq_ids=False):
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

        :param seq_id: A triple or an integer.
        """
        if isinstance(seq_id, int):
            seq_id=(" ", seq_id, " ")
        return self.seq_ids.index(seq_id)+1

    def create_bulge_graph(self, stems, bulges):
        """
        Find out which stems connect to which bulges

        Stems and bulges which share a nucleotide are considered connected.

        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.

        :param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                       numbers of the nucleotides at the start and end of the bulge.
        """
        for i in range(len(stems)):
            stem = stems[i]
            for j in range(len(bulges)):
                bulge = bulges[j]
                if any_difference_of_one(stem, bulge):
                    self.edges['y{}'.format(i)].add('b{}'.format(j))
                    self.edges['b{}'.format(j)].add('y{}'.format(i))

    def create_stem_graph(self, stems, bulge_counter):
        """
        Determine which stems are connected to each other. A stem can be connected to
        another stem when there is an interior loop with an unpaired nucleotide on
        one side. In this case, a bulge will be created on the other side, but it
        will only consist of the two paired bases around where the unpaired base 
        would be if it existed.

        The defines for these bulges will be printed as well as the connection strings
        for the stems they are connected to.

        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.
        :param bulge_counter: The number of bulges that have been encountered so far.

        :returns: A dictionary indexed by the number of a stem, containing a set of the 
                 other stems that the index is connected to.
        """
        # print "stems:", stems
        stem_stems = dict()
        for i in range(len(stems)):
            for j in range(i + 1, len(stems)):
                for k1 in range(2):
                    # don't fear the for loop
                    for k2 in range(2):
                        for l1 in range(2):
                            for l2 in range(2):
                                s1 = stems[i][k1][l1]
                                s2 = stems[j][k2][l2]
                                if abs(s1 - s2) == 1:
                                    stem_stems_set = stem_stems.get(i, set())
                                    if j not in stem_stems_set:
                                        bn = 'b{}'.format(bulge_counter)
                                        # self.defines[bn] = [min(s1, s2)+1, max(s1, s2)+1]
                                        self.defines[bn] = []
                                        self.weights[bn] = 1

                                        self.edges['y{}'.format(i)].add(bn)
                                        self.edges[bn].add('y{}'.format(i))

                                        self.edges['y{}'.format(j)].add(bn)
                                        self.edges[bn].add('y{}'.format(j))

                                        bulge_counter += 1
                                        stem_stems_set.add(j)
                                    stem_stems[i] = stem_stems_set

        for d in list(self.defines.keys()): #list needed for 4W92.pdb??? TODO: Need to check this!
            if d[0] != 'y':
                continue

            (s1, e1, s2, e2) = self.defines[d]
            if abs(s2 - e1) == 1:
                bn = 'b{}'.format(bulge_counter)

                self.defines[bn] = []
                self.weights[bn] = 1

                self.edges[bn].add(d)
                self.edges[d].add(bn)

                bulge_counter += 1

        return stem_stems

    def remove_vertex(self, v):
        """
        Delete a node after merging it with another

        :param v: The name of the node
        """
        # delete all edges to this node
        for key in self.edges[v]:
            self.edges[key].remove(v)

        for edge in self.edges:
            if v in self.edges[edge]:
                self.edges[edge].remove(v)

        # delete all edges from this node
        del self.edges[v]
        del self.defines[v]

    def reduce_defines(self):
        """
        Make defines like this:

        define x0 2 124 124 3 4 125 127 5 5
        
        Into this:

        define x0 2 3 5 124 127

        That is, consolidate contiguous bulge region defines.
        """
        for key in self.defines.keys():
            if key[0] != 's':
                assert (len(self.defines[key]) % 2 == 0)
                new_j = 0

                while new_j < len(self.defines[key]):

                    j = new_j
                    new_j += j + 2

                    (f1, t1) = (int(self.defines[key][j]), int(self.defines[key][j + 1]))

                    # remove bulges of length 0
                    if f1 == -1 and t1 == -2:
                        del self.defines[key][j]
                        del self.defines[key][j]

                        new_j = 0
                        continue

                    # merge contiguous bulge regions
                    for k in range(j + 2, len(self.defines[key]), 2):
                        if key[0] == 'y':
                            # we can have stems with defines like: [1,2,3,4]
                            # which would imply a non-existant loop at its end
                            continue

                        (f2, t2) = (int(self.defines[key][k]), int(self.defines[key][k + 1]))

                        if t2 + 1 != f1 and t1 + 1 != f2:
                            continue

                        if t2 + 1 == f1:
                            self.defines[key][j] = str(f2)
                            self.defines[key][j + 1] = str(t1)
                        elif t1 + 1 == f2:
                            self.defines[key][j] = str(f1)
                            self.defines[key][j + 1] = str(t2)

                        del self.defines[key][k]
                        del self.defines[key][k]

                        new_j = 0

                        break

    def merge_vertices(self, vertices):
        """
        This is done when two of the outgoing strands of a stem
        go to different bulges
        It is assumed that the two ends are on the same sides because
        at least one vertex has a weight of 2, implying that it accounts
        for all of the edges going out of one side of the stem

        :param vertices: A list of vertex names to combine into one.
        """
        new_vertex = self.get_vertex()
        self.weights[new_vertex] = 0

        # assert(len(vertices) == 2)

        connections = set()

        for v in vertices:

            # what are we gonna merge?
            for item in self.edges[v]:
                connections.add(item)

            # Add the definition of this vertex to the new vertex
            # self.merge_defs[new_vertex] = self.merge_defs.get(new_vertex, []) + [v]

            if v[0] == 's':
                self.defines[new_vertex] = self.defines.get(new_vertex, []) + [self.defines[v][0],
                                                            self.defines[v][2]] + [
                                                            self.defines[v][1], self.defines[v][3]]
            else:
                self.defines[new_vertex] = self.defines.get(new_vertex, []) + self.defines[v]

            self.weights[new_vertex] += 1

            # remove the old vertex, since it's been replaced by new_vertex
            self.remove_vertex(v)
            self.reduce_defines()

        # self.weights[new_vertex] = 2
        for connection in connections:
            self.edges[new_vertex].add(connection)
            self.edges[connection].add(new_vertex)

        return new_vertex

    def shortest_bg_loop(self, vertex):
        """
        Find the shortest loop containing this node. The vertex should
        be a multiloop.

        :param vertex: The name of the vertex to find the loop.
        :return: A list containing the elements in the shortest cycle.
        """
        G = self.to_networkx()

        # use the nucleotide in the middle of this element as the starting point
        residues = sorted(list(self.define_residue_num_iterator(vertex, adjacent=True)))
        mid_res = residues[len(residues) // 2]

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
        return path

    def nucleotides_to_elements(self, nucleotides):
        """
        Convert a list of nucleotides (nucleotide numbers) to element names.

        Remove redundant entries and return a set.
        """
        return set([self.get_node_from_residue_num(n) for n in nucleotides])

    def find_bulge_loop(self, vertex, max_length=4):
        """
        Find a set of nodes that form a loop containing the
        given vertex and being no greater than 4 nodes long.

        :param vertex: The vertex to start the search from.
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

    def dissolve_stem(self, key):
        """
        Remove a stem. This means that we need
        to reconfigure all of the adjacent elements in such a manner
        that they now include the nucleotides that were formerly 
        in this stem.
        """
        st = list(self.stem_bp_iterator(key))

        self.remove_base_pairs(st)

    def remove_base_pairs(self, to_remove):
        """
        Remove all of the base pairs which are in pair_list.

        :param to_remove: A list of tuples containing the names of the base pairs.
        :return: nothing
        """
        pt = self.to_pair_tuples()

        nt = []
        for p in pt:
            to_add = p
            for s in to_remove:
                if sorted(p) == sorted(s):
                    to_add = (p[0], 0)
                    break
            nt += [to_add]

        self.defines = dict()
        # self.edges = dict()

        self.from_tuples(nt)

    def collapse(self):
        """
        If any vertices form a loop, then they are either a bulge region of 
        a fork region. The bulge (interior loop) regions will be condensed 
        into one node.
        """

        new_vertex = True
        while new_vertex:
            new_vertex = False
            bulges = [k for k in self.defines if k[0] != 'y']

            for (b1, b2) in it.combinations(bulges, r=2):
                if self.edges[b1] == self.edges[b2] and len(self.edges[b1]) > 1:
                    connections = self.connections(b1)

                    all_connections = [sorted((self.get_sides_plus(connections[0], b1)[0],
                                               self.get_sides_plus(connections[0], b2)[0])),
                                       sorted((self.get_sides_plus(connections[1], b1)[0],
                                               self.get_sides_plus(connections[1], b2)[0]))]

                    if all_connections == [[1, 2], [0, 3]]:
                        # interior loop
                        self.merge_vertices([b1, b2])
                        new_vertex = True
                        break

    def interior_loop_iterator(self):
        """
        Iterate over all of the interior loops.

        An interior loop can only have two connections: to the two stems which it links. 
        """

        for key in self.defines.keys():
            if key[0] == 'i':
                yield key

    def relabel_node(self, old_name, new_name):
        """
        Change the name of a node.

        param old_name: The previous name of the node
        param new_name: The new name of the node
        """
        # replace the define name
        define = self.defines[old_name]

        del self.defines[old_name]
        self.defines[new_name] = define

        # replace the index into the edges array
        edge = self.edges[old_name]
        del self.edges[old_name]
        self.edges[new_name] = edge

        #replace the name of any edge that pointed to old_name
        for k in self.edges.keys():
            new_edges = set()
            for e in self.edges[k]:
                if e == old_name:
                    new_edges.add(new_name)
                else:
                    new_edges.add(e)
            self.edges[k] = new_edges

    def compare_stems(self, b):
        """
        A function that can be passed in as the key to a sort.
        """
        return (self.defines[b][0], 0)


    def compare_bulges(self, b):
        connections = self.connections(b)

        return (self.defines[connections[0]][0],
                self.defines[connections[1]][0])

    def compare_hairpins(self, b):
        connections = self.connections(b)

        return (self.defines[connections[0]][1], sys.maxint)

    def relabel_nodes(self):
        """
        Change the labels of the nodes to be more indicative of their nature.

        s: stem
        h: hairpin
        i: interior loop
        m: multiloop
        f: five-prime unpaired
        t: three-prime unpaired
        """
        stems = []
        hairpins = []
        interior_loops = []
        multiloops = []
        fiveprimes = []
        threeprimes = []

        for d in self.defines.keys():
            if d[0] == 'y' or d[0] == 's':
                stems += [d]

                stems.sort(key=self.compare_stems)
                continue

            if len(self.defines[d]) == 0 and len(self.edges[d]) == 1:
                hairpins += [d]
                continue

            if len(self.defines[d]) == 0 and len(self.edges[d]) == 2:
                multiloops += [d]
                continue

            if len(self.edges[d]) <= 1 and self.defines[d][0] == 1:
                fiveprimes += [d]
                continue

            if len(self.edges[d]) == 1 and self.defines[d][1] == self.seq_length:
                threeprimes += [d]
                continue

            if (len(self.edges[d]) == 1 and
                        self.defines[d][0] != 1 and
                        self.defines[d][1] != self.seq_length):
                hairpins += [d]

                hairpins.sort(key=self.compare_hairpins)
                continue

            if d[0] == 'm' or (d[0] != 'i' and len(self.edges[d]) == 2 and
                                       self.weights[d] == 1 and
                                       self.defines[d][0] != 1 and
                                       self.defines[d][1] != self.seq_length):
                multiloops += [d]

                multiloops.sort(key=self.compare_bulges)
                continue

            if d[0] == 'i' or self.weights[d] == 2:
                interior_loops += [d]
                interior_loops.sort(key=self.compare_stems)

        for d in fiveprimes:
            self.relabel_node(d, 'f1')
        for d in threeprimes:
            self.relabel_node(d, 't1')

        for i, d in enumerate(stems):
            self.relabel_node(d, 's%d' % (i))
        for i, d in enumerate(interior_loops):
            self.relabel_node(d, 'i%d' % (i))
        for i, d in enumerate(multiloops):
            self.relabel_node(d, 'm%d' % (i))
        for i, d in enumerate(hairpins):
            self.relabel_node(d, 'h%d' % (i))

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

            1   2 (1, 1) #pseudoknot
            1   0 (1, 0)
            3   2 (0, 1)
            3   0 (0, 0)

        :param define: The name of the bulge separating the two stems
        :param connections: The two stems and their separation

        :returns: INT connection type
                  *   positive values mean forward (from the connected stem starting at the 
                      lower nucleotide number to the one starting at the higher nuc. number)
                  *   negative values mean backwards.
                  1   interior loop
                  2   first multi-loop segment of normal multiloops and most pseudoknots
                  3   middle segment of a normal multiloop
                  4   last segment of normal multiloops and most pseudoknots
                  5   middle segments of pseudoknots
                  
        """

        if define[0] == 'i':
            # interior loop, we just have to check if 
            # connections[0] < connections[1]
            if self.defines[connections[0]][0] < self.defines[connections[1]][0]:
                return 1
            else:
                return -1
        elif define[0] == 'm':
            (s1c, b1c) = self.get_sides_plus(connections[0], define)
            (s2c, b2c) = self.get_sides_plus(connections[1], define)

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
                raise Exception("Weird angle type: (s1c, s2c) = (%d, %d)" %
                                (s1c, s2c))
        else:
            raise Exception("connection_type called on non-interior loop/multiloop")

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
            raise Exception('Unknown connection type: %d' % (connection_type))

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
            sides = [self.get_sides_plus(s, c)[0] for c in relevant_edges]
            sides.sort()

            # the whole stem is part of this multiloop
            if sides == [2, 3] or sides == [0, 1]:
                residues += list(range(self.defines[s][sides[0]], self.defines[s][sides[1]] + 1))
            else:
                residues += [self.defines[s][sides[0]], self.defines[s][sides[1]]]

        for m in multis:
            residues += self.define_residue_num_iterator(m, adjacent=False)
        return residues

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

    def find_multiloop_loops(self):
        """
        Find out which defines are connected in a multiloop.

        :return: Two lists, one containing the sets of nucleotides comprising the shortest loops
                 and the other containing sets of nucleotides comprising the shortest loops.
        """
        loops = set()

        for d in self.mloop_iterator():
            loop_nts = self.shortest_bg_loop(d)

            if len(loop_nts) > 0:
                loops.add(tuple(sorted(loop_nts)))

        loops = list(loops)
        loop_elems = []

        for loop in loops:
            all_loops = set([self.get_node_from_residue_num(n) for n in loop])

            # some multiloops might not contain any nucleotides, so we
            # have to explicitly add these
            for a, b in it.combinations(all_loops, r=2):
                common_edges = set.intersection(self.edges[a], self.edges[b])
                for e in common_edges:
                    all_loops.add(e)

            loop_elems += [all_loops]

        return loop_elems, loops

    def from_fasta(self, fasta_str, dissolve_length_one_stems=False):
        """
        Create a bulge graph from a fasta-type file containing the following
        format:

            > id
            ACCGGGG
            ((...))
        """
        lines = fasta_str.split('\n')
        self.from_dotbracket(lines[2].strip(), dissolve_length_one_stems)
        self.name = lines[0].strip('>')
        self.seq = lines[1].strip()

        self.seq_ids_from_seq()

    def seq_ids_from_seq(self):
        """
        Get the sequence ids of the string.
        """
        self.seq_ids = []

        # when provided with just a sequence, we presume that the
        # residue ids are numbered from 1-up
        for i, s in enumerate(self.seq):
            self.seq_ids += [(' ', i + 1, ' ')]

    def remove_degenerate_nodes(self):
        """
        For now just remove all hairpins that have no length.
        """
        to_remove = []
        for d in self.defines:
            if d[0] == 'h' and len(self.defines[d]) == 0:
                to_remove += [d]

        for r in to_remove:
            self.remove_vertex(r)

    def from_stems_and_bulges(self, stems, bulges):
        """
        Create the graph from the list of stems and bulges.

        :param stems: A list of tuples of two two-tuples, each containing the start
                      and end nucleotides of each strand of the stem.
        :param bulges: A list of tuples containing the starts and ends of the 
                       of the bulge regions.
        :return: Nothing, just make the bulgegraph
        """
        for i in range(len(stems)):
            # one is added to each coordinate to make up for the fact that residues are 1-based
            ss1 = stems[i][0][0] + 1
            ss2 = stems[i][0][1] + 1
            se1 = stems[i][1][0] + 1
            se2 = stems[i][1][1] + 1

            self.defines['y%d' % (i)] = [min(ss1, se1), max(ss1, se1),
                                         min(ss2, se2), max(ss2, se2)]
            self.weights['y%d' % (i)] = 1

        for i in range(len(bulges)):
            bulge = bulges[i]
            self.defines['b%d' % (i)] = sorted([bulge[0] + 1, bulge[1] + 1])
            self.weights['b%d' % (i)] = 1

        self.create_bulge_graph(stems, bulges)
        self.create_stem_graph(stems, len(bulges))
        self.collapse()
        self.relabel_nodes()
        self.remove_degenerate_nodes()
        self.sort_defines()

    def dissolve_length_one_stems(self):
        # dissolve all stems which have a length of one
        repeat = True
        while repeat:
            repeat = False
            for k in self.defines:
                if k[0] == 's' and self.stem_length(k) == 1:
                    self.dissolve_stem(k)
                    repeat = True
                    break

    def from_dotbracket(self, dotbracket_str, dissolve_length_one_stems=False):
        """
        Populate the BulgeGraph structure from a dotbracket representation.

        ie: ..((..))..

        :param dotbracket_str: A string containing the dotbracket representation
                               of the structure
        """
        self.__init__()
        self.dotbracket_str = dotbracket_str
        self.seq_length = len(dotbracket_str)

        if len(dotbracket_str) == 0:
            return

        pt = fus.dotbracket_to_pairtable(dotbracket_str)
        tuples = fus.pairtable_to_tuples(pt)
        self.from_tuples(tuples)

        if dissolve_length_one_stems:
            self.dissolve_length_one_stems()

    def to_pair_table(self):
        """
        Create a pair table from the list of elements.

        The first element in the returned list indicates the number of
        nucleotides in the structure.

        i.e. [5,5,4,0,2,1]
        """
        pair_tuples = self.to_pair_tuples()

        return fus.tuples_to_pairtable(pair_tuples)

    def to_pair_tuples(self):
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
                table += [(b, p)]

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
        for line in lines:
            parts = line.split()

            if len(parts) == 0:
                continue

            (t1, s, t2) = (int(parts[0]), parts[1], int(parts[2]))
            tuples += [(t1, t2)]
            seq += [s]

        seq = "".join(seq).upper().replace('T', 'U')

        return (tuples, seq)

    def from_bpseq_str(self, bpseq_str, dissolve_length_one_stems=False):
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

        tuples, seq = self.bpseq_to_tuples_and_seq(bpseq_str)

        self.seq = seq
        self.seq_length = len(seq)
        self.from_tuples(tuples)

        if dissolve_length_one_stems:
            self.dissolve_length_one_stems()

    def from_tuples(self, tuples):
        """
        Create a bulge_graph from a list of pair tuples. Unpaired
        nucleotides have a pairing partner of 0.
        """
        stems = []
        bulges = []

        tuples.sort()
        tuples = iter(tuples)
        (t1, t2) = next(tuples) #.next()

        prev_from = t1
        prev_to = t2

        start_from = prev_from
        start_to = prev_to
        last_paired = prev_from

        for t1, t2 in tuples:
            (from_bp, to_bp) = (t1, t2)

            if abs(to_bp - prev_to) == 1 and prev_to != 0:
                # stem
                if (((prev_to - prev_from > 0 and to_bp - from_bp > 0) or
                         (prev_to - prev_from < 0 and to_bp - from_bp < 0)) and
                            (to_bp - prev_to) == -(from_bp - prev_from)):
                    (prev_from, prev_to) = (from_bp, to_bp)
                    last_paired = from_bp
                    continue

            if to_bp == 0 and prev_to == 0:
                # bulge
                (prev_from, prev_to) = (from_bp, to_bp)
                continue
            else:
                if prev_to != 0:
                    new_stem = tuple(sorted([tuple(sorted([start_from - 1, start_to - 1])),
                                             tuple(sorted([prev_from - 1, prev_to - 1]))]))
                    if new_stem not in stems:
                        stems += [new_stem]

                    last_paired = from_bp
                    start_from = from_bp
                    start_to = to_bp
                else:
                    new_bulge = ((last_paired - 1, prev_from - 1))
                    bulges += [new_bulge]

                    start_from = from_bp
                    start_to = to_bp

            prev_from = from_bp
            prev_to = to_bp

        # Take care of the last element
        if prev_to != 0:
            new_stem = tuple(sorted([tuple(sorted([start_from - 1, start_to - 1])),
                                     tuple(sorted([prev_from - 1, prev_to - 1]))]))
            if new_stem not in stems:
                stems += [new_stem]
        if prev_to == 0:
            new_bulge = ((last_paired - 1, prev_from - 1))
            bulges += [new_bulge]
        self.from_stems_and_bulges(stems, bulges)

    def sort_defines(self):
        """
        Sort the defines of interior loops and stems so that the 5' region
        is always first.
        """
        for k in self.defines.keys():
            d = self.defines[k]

            if len(d) == 4:
                if d[0] > d[2]:
                    new_d = [d[2], d[3], d[0], d[1]]
                    self.defines[k] = new_d

    def to_dotbracket_string(self):
        """
        Convert the BulgeGraph representation to a dot-bracket string
        and return it.

        :return: A dot-bracket representation of this BulgeGraph
        """
        pt = self.to_pair_table()
        return fus.pairtable_to_dotbracket(pt)

    def to_fasta_string(self):
        """
        Output the BulgeGraph representation as a fast string of the
        format:

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
            elif parts[0] == 'seq_ids':
                self.seq_ids = list(map(ftum.parse_resid, parts[1:]))
            elif parts[0] == 'name':
                self.name = parts[1].strip()
            elif parts[0] == 'info':
                self.infos[parts[1]].append(" ".join(parts[2:]))

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
                      (useful for pseudeknots)
        The connected nucleotides are those which are spanned by a single
        interior loop or multiloop. In the case of an interior loop, this
        function will return a list of two tuples and in the case of multiloops
        if it will be a list of one tuple.

        If the two stems are not separated by a single element, then return
        an empty list.
        """
        ## sort the stems according to the number of their first nucleotide
        #stems = [s1, s2]
        #stems.sort(key=lambda x: self.defines[x][0])

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

        # return the ones which are closest to each other
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
        #print ("LINKED FOR {}: {}".format(bulge, linked))
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

        raise Exception("Invalid side (%d) for the stem (%s)." % (stem, side))
    
    def get_stem_edge(self, stem, pos):
        """
        Returns the side of the stem that position is on
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
        
        raise Exception("Position (%d) not in stem (%s)." % (pos, stem))
    
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

        # if the bulge is a length 0 multiloop then use the adjacent
        # stem to determine its side
        if len(bd) == 0:
            edges = self.edges[b]

            for e in edges:
                if e != s1:
                    bd = self.defines[e]
                    break

        for i in range(4):
            for k in range(len(bd)):
                if s1d[i] - bd[k] == 1:
                    if i == 0:
                        s1b = 0
                        break
                    if i == 2:
                        s1b = 1
                        break
                elif s1d[i] - bd[k] == -1:
                    if i == 1:
                        s1b = 1
                        break
                    if i == 3:
                        s1b = 0
                        break
        if s1b == 0:
            s1e = 1
        else:
            s1e = 0

        return (s1b, s1e)

    def get_sides_plus(self, s1, b):
        """
        Get the side of s1 that is next to b.

        s1e -> s1b -> b

        :param s1: The stem.
        :param b: The bulge.
        :return: A tuple indicating the corner of the stem that connects
                 to the bulge as well as the corner of the bulge that connects
                 to the stem.
                 These sides are equivalent to the indices of the define.
        """
        s1d = self.defines[s1]
        bd = self.defines[b]

        if len(bd) == 0:
            edges = self.edges[b]

            for e in edges:
                if e != s1:
                    bd = self.defines[e] #For bulges of length 0, use the next stem
                    break

        for k in range(len(bd)):
            # before the stem on the 5' strand
            if s1d[0] - bd[k] == 1:
                return (0, k)
            # after the stem on the 5' strand
            elif bd[k] - s1d[1] == 1:
                return (1, k)
            # before the stem on the 3' strand
            elif s1d[2] - bd[k] == 1:
                return (2, k)
            # after the stem on the 3' strand
            elif bd[k] - s1d[3] == 1:
                return (3, k)

        raise Exception("Faulty multiloop %s connecting %s"
                        % (" ".join(map(str, bd)),
                           " ".join(map(str, s1d))))

    def stem_resn_to_stem_vres_side(self, stem, res):
        d = self.defines[stem]
        if res<=d[1]:
            pos=res-d[0]
            side = 0
        elif res<=d[3]:
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
        if 'f1' in self.defines.keys():
            yield 'f1'

    def tloop_iterator(self):
        """
        Yield the name of the 3' prime unpaired region if it is
        present in the structure.
        """
        if 't1' in self.defines.keys():
            yield 't1'

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

    @profile
    def connections(self, bulge):
        """
        Return the edges that connect to a bulge in a list form,
        sorted by lowest res number of the connection.
        """
        def sort_key(x):
            if len(self.defines[x]) > 0:
                if self.defines[x][0] == 1:
                    # special case for stems at the beginning since there is no
                    # adjacent nucleotide 0
                    return 0
            return list(self.define_residue_num_iterator(x, adjacent=True))[0]


        connections = list(self.edges[bulge])
        connections.sort(key=sort_key)

        return connections

    def get_define_seq_str(self, d, adjacent=False):
        """
        Get an array containing the sequences for the given define.
        Non-stem sequences will contain the sequence without the overlapping
        stem residues that are part of the define.

        :param d: The define for which to get the sequences
        :return: An array containing the sequences corresponding to the defines
        """
        define = self.defines[d]
        ranges = zip(*[iter(define)] * 2)
        c = self.connections(d)

        if d[0] == 'i':
            s1 = self.defines[c[0]]
            s2 = self.defines[c[1]]
            if adjacent:
                return [self.seq[s1[1] - 1:s2[0]],
                        self.seq[s2[3] - 1:s1[2]]]
            else:
                return [self.seq[s1[1]:s2[0] - 1],
                        self.seq[s2[3]:s1[2] - 1]]
        if d[0] == 'm':
            s1 = self.defines[c[0]]
            s2 = self.defines[c[1]]

            i1 = s1[self.get_sides_plus(c[0], d)[0]]
            i2 = s2[self.get_sides_plus(c[1], d)[0]]

            (i1, i2) = (min(i1, i2), max(i1, i2))

            if adjacent:
                return [self.seq[i1 - 1:i2]]
            else:
                return [self.seq[i1:i2 - 1]]
        else:
            seqs = []
            for r in ranges:
                if d[0] == 's':
                    seqs += [self.seq[r[0] - 1:r[1]]]
                else:
                    if adjacent:
                        if r[0] > 1:
                            seqs += [self.seq[r[0] - 2:r[1] + 1]]
                        else:
                            seqs += [self.seq[r[0] - 1:r[1] + 1]]
                    else:
                        seqs += [self.seq[r[0] - 1:r[1]]]

            return seqs


    
    #Seems to be unused!
    def get_stem_direction(self, s1, s2):
        """
        Return 0 if the lowest numbered residue in s1
        is lower than the lowest numbered residue in s2.
        """
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

        p1 = self.get_sides_plus(c[0], m)
        p2 = self.get_sides_plus(c[1], m)

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

        If it is single stranded it will be (0, x). Otherwise it will be (x, y).

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

    def get_node_from_residue_num(self, base_num, seq_id=False):
        """
        Iterate over the defines and see which one encompasses this base.
        """
        for key in self.defines.keys():
            define = self.defines[key]

            for i in range(0, len(define), 2):
                a = [int(define[i]), int(define[i + 1])]
                a.sort()

                if seq_id:
                    for i in range(a[0], a[1] + 1):
                        if self.seq_ids[i - 1][1] == base_num:
                            return key
                else:
                    if base_num >= a[0] and base_num <= a[1]:
                        return key

        raise Exception("Base number %d not found in the defines." % (base_num))

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
                    raise Exception("Weird multiloop sides: %s" %
                                    bulge_name)

            ends.sort()

            return tuple(ends)
            # multiloop

        return (None, None)


    def get_flanking_sequence(self, bulge_name, side=0):
        if len(self.seq) == 0:
            raise Exception("No sequence present in the bulge_graph: %s" % (self.name))

        (m1, m2) = self.get_flanking_region(bulge_name, side)

        return self.seq[m1 - 1:m2]

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

    def same_stem_end(self, sd):
        """
        Return the index of the define that is on the same end of the
        stem as the index sd.

        :param sd: An index into a define.
        :return: The index pointing to the nucleotide on the other strand 
                 on the same side as the stem.
        """
        warnings.warn("BulgeGraph.same_stem_end() is deprecated and will be removed in the future."
                       "Use ??? instead." , DeprecationWarning, stacklevel=2)
        if sd == 0:
            return 3
        elif sd == 1:
            return 2
        elif sd == 2:
            return 1
        else:
            return 0

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
                    strand_resnames += [self.seq_ids[x - 1]]
                else:
                    strand_resnames += [x]

            resnames += [strand_resnames]

        return resnames

    def seqids_from_residue_map(self, residue_map):
        """
        Create the list of seq_ids from the list of MC-Annotate identifiers in the
        residue map.
        """
        self.seq_ids = []

        for r in residue_map:
            (from_chain, from_base) = ftum.parse_chain_base(r)

            self.seq_ids += [ftum.parse_resid(from_base)]

    # This function seems to be dead code, but might be useful in the future.
    # Consider adding this to whitelist.oy
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
                       key=lambda x: (priority[x[0]], min(self.get_node_dimensions(x))))
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
                           self.tloop_iterator()))

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
                # assert  self.get_stem_direction(prev, build_order[-1][2])==0 does not hold!
        self.build_order = build_order

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

    def get_angle_type(self, bulge):
        """
        Return what type of angle this bulge is, based on the way this
        would be built using a breadth-first traversal along the minimum
        spanning tree.
        """
        if self.ang_types is None:
            self.set_angle_types()

        if bulge in self.ang_types:
            return self.ang_types[bulge]
        else:
            return None

    def is_node_pseudoknot(self, d):
        """
        Is a particular multiloop part of a pseudoknot?
        """
        conn = self.connections(d)
        ct = self.connection_type(d, conn)
        if abs(ct) == 5:
            return True
        return False

    def is_loop_pseudoknot(self, loop):
        """
        Is a particular loop a pseudoknot?

        :param loop: A list of elements that are part of the loop.
        :return: Either True or false
        """
        allowed_ang_types = [2, 3, 4]
        found_ang_types = set()

        for l in loop:
            if l[0] != 'm':
                continue

            conn = self.connections(l)
            ctype = self.connection_type(l, conn)

            if ctype not in allowed_ang_types:
                return True

            found_ang_types.add(ctype)

        if len(found_ang_types) == 3:
            return False

        return True


    def is_pseudoknot(self):
        """
        Is this bulge part of a pseudoknot?
        """
        for d in self.mloop_iterator():
            if self.is_node_pseudoknot(d):
                return True

        return False

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

        residues.sort()
        prev = None
        for r in residues:
            if prev is not None:
                G.add_edge(prev, r)
            prev = r

        for s in self.stem_iterator():
            for (f, t) in self.stem_bp_iterator(s):
                G.add_edge(f, t)

        return G

    def remove_pseudoknots(self):
        """
        Remove all of the pseudoknots using the knotted2nested.py script.

        :return: A list of base-pairs that were removed.
        """
        # remove unpaired bases and redundant pairs (i.e. (2,3) and (3,2))
        pairs = sorted([tuple(sorted(p)) for p in self.to_pair_tuples() if p[1] != 0])
        pairs = list(set(pairs))

        # knotted_struct = fak.KnottedStructure(pairs, Seq=self.seq, Header=[])

        import forgi.aux.k2n_standalone.knots as fakk

        pk_function = fakk.eg
        nested_pairs, removed_pairs = pk_function(pairs, return_removed=True)

        self.remove_base_pairs(removed_pairs)
        return removed_pairs

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

            (s1c, b1c) = self.get_sides_plus(common_stem, n1)
            (s2c, b1c) = self.get_sides_plus(common_stem, n2)

            if sorted([s1c, s2c]) == [0,3] or sorted([s1c, s2c]) == [1,2]:
                return True

        return False

    # This function seems to be unused. Consider deprecation...
    def flanking_nucleotides(self, d):
        '''
        Return the nucleotides directly flanking an element.
    
        :param d: the name of the element
        :return: a set of nucleotides
        '''
        set_adjacent = set(self.define_residue_num_iterator(d, adjacent=True))
        set_not_adjacent = set(self.define_residue_num_iterator(d, adjacent=False))

        return set_adjacent - set_not_adjacent




    @profile
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

        #flanking1 = self.flanking_nucleotides(e1)
        #flanking2 = self.flanking_nucleotides(e2)

        if (e1,e2) in self._elem_bp_dists: #Shortcut if cached.
            return self._elem_bp_dists[(e1,e2)]

        min_bp = sys.maxint
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
          * multiloops with connected stems
          * rods: stretches of stems + interior loops (without branching), with trailing hairpins
          * pseudoknots
        """
        domains = col.defaultdict(list)
        multiloops, nucleotides = self.find_multiloop_loops()
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

def bg_from_subgraph(bg, sg):
    """
    Create a BulgeGraph from a list containing the nodes
    to take from the original.

    WARNING: The sequence information is not copied
    """
    nbg = BulgeGraph()
    nbg.seq_length = 0

    for d in sg:
        # copy the define
        nbg.defines[d] = bg.defines[d][::]

    # copy edges only if they connect elements which 
    # are also in the new structure
    for e in bg.edges.keys():
        for conn in bg.edges[e]:
            if conn in sg:
                nbg.edges[e].add(conn)

    return nbg

