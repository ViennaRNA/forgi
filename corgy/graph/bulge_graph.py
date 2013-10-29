#!/usr/bin/env python

"""bulge_graph.py: A graph representation of RNA secondary structure based
   on its decomposition into primitive structure types: stems, hairpins,
   interior loops, multiloops, etc..."""

__author__      = "Peter Kerpedjiev"
__copyright__   = "Copyright 2012, 2013"
__license__     = "GPL"
__version__     = "0.1"
__maintainer__  = "Peter Kerpedjiev"
__email__       = "pkerp@tbi.univie.ac.at"

import sys
import collections as c
import math
import random
import itertools as it
import corgy.utilities.debug as cud

def error_exit(message):
    print >> sys.stderr, message
    sys.exit(1)

# A wrapper for a simple dictionary addition
# Added so that debugging can be made easier
def add_bulge(bulges, bulge, context, message):
    #print >>sys.stderr,"Adding bulge", context, bulge, message
    #bulge = (context, bulge)
    bulges[context] = bulges.get(context, []) + [bulge]
    return bulges

def any_difference_of_one(stem, bulge):
    '''
    See if there's any difference of one between the two
    ends of the stem [(a,b),(c,d)] and a bulge (e,f)

    :param stem: A couple of couples (2 x 2-tuple) indicating the start and end
                 nucleotides of the stem in the form ((s1, e1), (s2, e2))
    :param bulge: A couple (2-tuple) indicating the first and last position
                  of the bulge.
    :return: True if there is an overlap between the stem nucleotides and the 
                  bulge nucleotides. False otherwise
    '''
    for stem_part in stem:
        for part in stem_part:
            for bulge_part in bulge:
                if abs(bulge_part - part) == 0:
                    return True
    return False


def print_bulge_graph(graph):
    '''
    Print out the connections in the graph.

    :param graph: A dictionary indexed by stem number containing a set
                  of the bulge numbers that it is connected to.
    '''
    for key in graph.keys():
        stem_str = "connect s%d" % (key)
        for item in graph[key]:
            stem_str += " b%d" % (item)
        print stem_str

def print_stems(stems):
    '''
    Print the names and definitions of the stems.

    :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    '''
    for i in range(len(stems)):
        # one is added to each coordinate to make up for the fact that residues are 1-based
        ss1 = stems[i][0][0]+1
        ss2 = stems[i][0][1]+1
        se1 = stems[i][1][0]+1
        se2 = stems[i][1][1]+1

        print "define s%d 0 %d %d %d %d" % (i, min(ss1,se1), max(ss1,se1), min(ss2,se2), max(ss2,se2))

def print_bulges(bulges):
    '''
    Print the names and definitions of the bulges.

    :param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                   numbers of the nucleotides at the start and end of the bulge.
    '''
    for i in range(len(bulges)):
            #print "bulge:", bulge
        bulge_str = "define b%d 1" % (i)
        bulge = bulges[i]
        #print >>sys.stderr, "bulge:", bulge
        bulge_str += " %d %d" % (bulge[0]+1, bulge[1]+1)
        print bulge_str

def condense_stem_pairs(stem_pairs):
    '''
    Given a list of stem pairs, condense them into stem definitions

    I.e. the pairs (0,10),(1,9),(2,8),(3,7) can be condensed into
    just the ends of the stem: [(0,10),(3,7)]

    :param stem_pairs: A list of tuples containing paired base numbers.

    :returns: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    '''
    stem_pairs.sort()

    prev_pair = (-10, -10)

    stems = []
    start_pair = None

    for pair in stem_pairs:
        # There's a potential bug here since we don't check the direction
        # but hopefully it won't bite us in the ass later
        if abs(pair[0] - prev_pair[0]) != 1 or abs(pair[1] - prev_pair[1]) != 1:
            if start_pair != None:
                stems += [(start_pair, prev_pair)]
            start_pair = pair
    
        prev_pair = pair

    if start_pair != None:
        stems += [(start_pair, prev_pair)]

    return stems

def print_brackets(brackets):
    '''
    Print the brackets and a numbering, for debugging purposes

    :param brackets: A string with the dotplot passed as input to this script.
    '''
    numbers = [chr(ord('0') + i % 10) for i in range(len(brackets))]
    tens = [chr(ord('0') + i / 10) for i in range(len(brackets))]
    print "brackets:\n", brackets, "\n", "".join(tens), "\n" ,"".join(numbers)

def find_bulges_and_stems(brackets):
    '''
    Iterate through the structure and enumerate the bulges and the stems that are
    present.

    The returned stems are of the form [[(s1, s2), (e1,e2)], [(s1,s2),(e1,e2)],...]
    where (s1,s2) are the residue numbers of one end of the stem and (e1,e2) are the
    residue numbers at the other end of the stem
    (see condense_stem_pairs)

    The returned bulges are of the form [(s,e), (s,e),...] where s is the start of a bulge
    and e is the end of a bulge

    :param brackets: A string with the dotplot passed as input to this script.
    '''
    prev = '.'
    context = 0

    bulges = dict()
    finished_bulges = []
    context_depths = dict()

    opens = []
    stem_pairs = []

    stems = dict()

    dots_start = 0
    dots_end = 0

    context_depths[0] = 0

    i = 0
    for i in range(len(brackets)):
        #print >> sys.stderr, "bracket:", brackets[i]
        if brackets[i] == '(':
            opens.append(i)

            if prev == '(':
                context_depths[context] = context_depths.get(context, 0) + 1
                continue
            else:
                context += 1
                context_depths[context] = 1

            if prev == '.':
                dots_end = i-1
                bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "4")
            """
            if prev == ')':
                bulges = add_bulge(bulges, (-1,-1), context, "3")
            """

        if brackets[i] == ')':
            if len(opens) == 0:
                error_exit("ERROR: Unmatched close bracket")

            stem_pairs.append((opens.pop(), i))

            context_depths[context] -= 1

            if context_depths[context] == 0:
                finished_bulges += bulges[context]
                bulges[context] = []
                context -= 1
 

            if prev == '.':
                dots_end = i-1
                bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "2")
  
        if brackets[i] == '.':
            if prev == '.':
                continue

            dots_start = i

        prev = brackets[i]
    if prev == '.':
        dots_end = i
        bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "7")
    elif prev == '(':
        print >>sys.stderr, "Unmatched bracket at the end"
        sys.exit(1)
    elif prev == ')':
        bulges = add_bulge(bulges, (i, i+1), context, "8")
    
    if context in bulges.keys():
        finished_bulges += bulges[context]

    if len(opens) > 0:
        error_exit("ERROR: Unmatched open bracket")

    stem_pairs.sort()
    stems = condense_stem_pairs(stem_pairs)
    
    return (finished_bulges, stems)

def print_name(filename):
    print "name", os.path.splitext(filename)[0]


class BulgeGraph:
    def __init__(self, dotbracket_str=''):
        self.name = "untitled"
        self.seq = ""
        self.defines = dict()
        self.edges = c.defaultdict(set)
        self.longrange = c.defaultdict(set)
        self.weights = dict()
        self.seq_length = 0

        self.name_counter = 0

        if dotbracket_str != '':
            self.from_dotbracket(dotbracket_str)

    # get an internal index for a named vertex
    # this applies to both stems and edges
    def get_vertex(self, name = None):
        '''
        Return a new unique vertex name.
        '''

        if name == None:
            name = "x%d" % (self.name_counter)
            self.name_counter += 1

        return name

    def stem_length(self, key):
        d = self.defines[key]
        return (d[1] - d[0]) + 1

    def get_define_str(self):
        '''
        Convert the defines into a string. 

        Format:

        define [name] [weight] [start_res1] [end_res1] [start_res2] [end_res2]
        '''
        defines_str = ''
        for key in self.defines.keys():
            defines_str += "define %s %s" % ( key, " ".join([str(d) for d in self.defines[key]]))
            defines_str += '\n'
        return defines_str

    def get_length_str(self):
        return "length " + str(self.seq_length) + '\n'

    def get_connect_str(self):
        '''
        Get the connections of the bulges in the graph.

        Format:

        connect [from] [to1] [to2] [to3]
        '''

        whole_str = ''
        for key in self.edges:
            # Our graph will be defined by the stems and the bulges they connect to
            name = key
            if name[0] == 's':
                out_str = "connect %s" % (name)

                for dest in self.edges[key]:
                    out_str += " %s" % (dest)

                whole_str += out_str
                whole_str += '\n'
        return whole_str

    def get_sequence_str(self):
        '''
        Return the sequence along with its keyword. I.e.

            seq ACGGGCC
        '''
        if len(self.seq) > 0:
            return "seq %s\n" % (self.seq)
        else:
            return ""

    def get_name_str(self):
        '''
        Return the name of this structure along with its keyword:

            name 1y26
        '''
        return "name %s\n" % (self.name)



    def to_bg_string(self):
        '''
        Output a string representation that can be stored and reloaded.

        '''
        out_str = ''
        out_str += self.get_name_str()
        out_str += self.get_length_str()
        out_str += self.get_sequence_str()
        out_str += self.get_define_str()
        out_str += self.get_connect_str()

        return out_str

    def to_element_string(self):
        '''
        Create a string similar to dotbracket notation that identifies what
        type of element is present at each location.

        For example the following dotbracket:

        ..((..))..

        Should yield the following element string:

        ffsshhsstt

        Indicating that it begins with a fiveprime region, continues with a
        stem, has a hairpin after the stem, the stem continues and it is terminated
        by a threeprime region.
        '''
        output_str = [' ' for i in range(self.seq_length+1)]

        for d in self.defines.keys():
            for resi in self.define_residue_num_iterator(d, adjacent=False):
                output_str[resi] = d[0]

        return "".join(output_str).strip()

    def define_residue_num_iterator(self, node, adjacent=True):
        '''
        Iterate over the residue numbers that belong to this node.

        Note that everything except stems starts and ends one node
        before and after its actual start. So a multiloop section
        from residue 10 to 13 will only actually contain residues
        11 and 12, since residues 10 and 13 will belong to the adjacent
        stems.

        :param node: The name of the node
        '''
        if node[0] == 's':  
            # iterate in sets of two
            a = iter(self.defines[node])
            for (ds1, ds2) in it.izip(a,a):
                for i in range(ds1, ds2+1):
                    yield i
        else:
            # iterate in sets of two elements
            a = iter(self.defines[node])
            for (ds1, ds2) in it.izip(a,a):
                if ds1 == 0 and ds2 == 1:
                    continue

                if ds1 == self.seq_length:
                    continue

                if adjacent:
                    for i in range(ds1, ds2+1):
                        if i <= 0 or i > self.seq_length:
                            continue
                        yield i
                else:
                    for i in range(ds1+1, ds2):
                        yield i

    def create_bulge_graph(self, stems, bulges):
        '''
        Find out which stems connect to which bulges

        Stems and bulges which share a nucleotide are considered connected.

        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.

        :param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                       numbers of the nucleotides at the start and end of the bulge.
        '''
        for i in range(len(stems)):
            stem = stems[i]
            for j in range(len(bulges)):
                bulge = bulges[j]
                if any_difference_of_one(stem, bulge):
                    self.edges['s%d' % (i)].add('b%d' % (j))
                    self.edges['b%d' % (j)].add('s%d' % (i))

    def create_stem_graph(self, stems, bulge_counter):
        '''
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
        '''
        #print "stems:", stems
        stem_stems = dict()
        define_text = ""
        connect_text = ""
        for i in range(len(stems)):
            stem = stems[i]
            for j in range(i+1, len(stems)):
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
                                        bn = 'b%d' % (bulge_counter)
                                        self.defines[bn] = [min(s1, s2)+1, max(s1, s2)+1]
                                        self.weights[bn] = 1

                                        self.edges['s%d' % i].add(bn)
                                        self.edges[bn].add('s%d' % i)

                                        self.edges['s%d' % j].add(bn)
                                        self.edges[bn].add('s%d' % j)

                                        bulge_counter += 1
                                        stem_stems_set.add(j)
                                    stem_stems[i] = stem_stems_set
        return stem_stems

    def remove_vertex(self, v):
        '''
        Delete a node after merging it with another

        :param v: The name of the node
        '''
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

        new_defines = []

        for key in self.defines.keys():
            if key[0] != 's':
                assert(len(self.defines[key]) % 2 == 0)

                j = 0
                new_j = 0
                #print >>sys.stderr,"self.defines[key]:", self.defines[key]

                while new_j < len(self.defines[key]):

                    j = new_j
                    new_j += j + 2

                    #print >>sys.stderr, "reduce_defines", i, j
                    (f1, t1) = (int(self.defines[key][j]), int(self.defines[key][j+1]))

                    # remove bulges of length 0
                    if f1 == -1 and t1 == -2:
                        del self.defines[key][j]
                        del self.defines[key][j]
                        
                        new_j = 0
                        continue

                    # merge contiguous bulge regions
                    for k in range(j+2, len(self.defines[key]), 2):
                        (f2, t2) = (int(self.defines[key][k]), int(self.defines[key][k+1]))

                        #print >>sys.stderr, "j: %d f1: %d, t1: %d, f2: %d, t2: %d" % (j, f1, t1, f2, t2)

                        if t2 + 1 != f1 and t1 + 1 != f2:
                            continue

                        #print >>sys.stderr, "pre self.defines[key]:", self.defines[key]

                        if t2 + 1 == f1:
                            self.defines[key][j] = str(f2)
                            self.defines[key][j+1] = str(t1)
                        elif t1 + 1 == f2:
                            self.defines[key][j] = str(f1)
                            self.defines[key][j+1] = str(t2)

                        del self.defines[key][k]
                        del self.defines[key][k]

                        #print >>sys.stderr, "post self.defines[key]:", self.defines[key]
                        
                        new_j = 0

                        break

    def merge_vertices(self, vertices):
        '''
        This is done when two of the outgoing strands of a stem
        go to different bulges
        It is assumed that the two ends are on the same sides because
        at least one vertex has a weight of 2, implying that it accounts
        for all of the edges going out of one side of the stem

        :param vertices: A list of vertex names to combine into one.
        '''
        merge_str = ""
        new_vertex = self.get_vertex()
        self.weights[new_vertex] = 0

        #assert(len(vertices) == 2)

        connections = set()
        needs_merging = set()

        for v in vertices:
            merge_str += " %s" % (v)

            # what are we gonna merge?
            for item in self.edges[v]:
                connections.add(item)

            # Add the definition of this vertex to the new vertex
            #self.merge_defs[new_vertex] = self.merge_defs.get(new_vertex, []) + [v]

            if v[0] == 's':
                self.defines[new_vertex] = self.defines.get(new_vertex, []) + [self.defines[v][0], self.defines[v][2]] + [self.defines[v][1], self.defines[v][3]]
            else:
                self.defines[new_vertex] = self.defines.get(new_vertex,[]) + self.defines[v]


            self.weights[new_vertex] += 1

            # remove the old vertex, since it's been replaced by new_vertex
            self.remove_vertex(v)
            self.reduce_defines()

        #self.weights[new_vertex] = 2
        for connection in connections:
            self.edges[new_vertex].add(connection)
            self.edges[connection].add(new_vertex)

        return new_vertex

    def find_bulge_loop(self, vertex, max_length=4):
        '''
        Find a set of nodes that form a loop containing the
        given vertex and being no greater than 4 nodes long.

        :param vertex: The vertex to start the search from.
        :returns: A list of the nodes in the loop.
        '''
        visited = set()
        to_visit = [(key, 0) for key in self.edges[vertex]]
        visited.add(vertex)
        in_path = [vertex]

        while len(to_visit) > 0:
            (current, depth) = to_visit.pop()
            visited.add(current)

            in_path = in_path[:depth]
            in_path.append(current)

            for key in self.edges[current]:
                if key == vertex and depth > 1:
                    if len(in_path[:depth+1]) > max_length:
                        continue
                    else:
                        return in_path[:depth+1]
                
                if key not in visited:
                    to_visit.append((key, depth+1))
        return []

    def add_node(self, name, edges, define, weight=1):
        self.defines[name] = define
        self.edges[name] = edges
        self.weights[name] = weight

        for edge in self.edges[name]:
            self.edges[edge].add(name)

    def dissolve_stem(self, key):
        #print >>sys.stderr, "hi"
        d = self.defines[key]

        bulge_sides = dict()
        bulge_sides[0] = set()
        bulge_sides[1] = set()

        #print >>sys.stderr,"dissolving:", key
        #print >>sys.stderr, "edges:", self.edges[key]
        for edge in self.edges[key]:
            be = self.defines[edge]

            for i in range(0, 4):
                for j in range(0, len(be)):
                    if d[i] == be[j]:
                        #print >>sys.stderr, key, edge
                        bulge_sides[i / 2].add(edge)
        
        #print >>sys.stderr, "bulge_sides:", bulge_sides
        
        new_nodes = [0,0]
        for i in range(2):
            new_node = self.get_vertex()
            edges = set()

            mins = 10000
            maxs = -10000

            for bulge in bulge_sides[i]:
                if bulge in self.defines:
                    bd = self.defines[bulge]
                else:
                    continue
                
                for edge in self.edges[bulge]:
                    edges.add(edge)

                mins = min(mins, min(bd))
                maxs = max(maxs, max(bd))
                self.remove_vertex(bulge)

            edges.remove(key)
            #print >> sys.stderr, "new_node", new_node, "edges:", edges, "mins:", mins, "maxs:", maxs
           
            self.add_node(new_node, edges, [mins, maxs], self.weights[bulge])
            new_nodes[i] = new_node

        if len(self.edges[new_nodes[0]]) == 1 and len(self.edges[new_nodes[1]]) == 1:
            dnn0 = self.defines[new_nodes[0]]
            dnn1 = self.defines[new_nodes[1]]
            newer_node = self.get_vertex()

            define = [min(dnn0[0], dnn1[0]), max(dnn0[1], dnn1[1])]
            edges = self.edges[new_nodes[0]].union(self.edges[new_nodes[1]])

            self.add_node(newer_node, edges, define, self.weights[new_nodes[0]])

            self.remove_vertex(new_nodes[0])
            self.remove_vertex(new_nodes[1])

        self.remove_vertex(key)

    def collapse(self):
        '''
        If any vertices form a loop, then they are either a bulge region of 
        a fork region. The bulge (interior loop) regions will be condensed 
        into one node.
        '''
        # for all stems of length 1, merge with adjacent bulges
        for key in self.edges.keys():
            if self.dissolve_length_one_stems:
                if key[0] == 's':
                    if self.stem_length(key) == 1:
                        self.dissolve_stem(key)

        for key in self.edges.keys():
            if key[0] == 's':
                while True:
                    loop = self.find_bulge_loop(key)
                    bulge_loop = [v for v in loop if v[0] == 'b' or v[0] == 'x']

                    if len(bulge_loop) > 0:
                        assert(len(bulge_loop) != 1)
                        self.merge_vertices(bulge_loop)

                    if len(bulge_loop) == 0:
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
        '''
        Change the name of a node.

        param old_name: The previous name of the node
        param new_name: The new name of the node
        '''
        #replace the define name
        define = self.defines[old_name]
        del self.defines[old_name]
        self.defines[new_name] = define

        #replace the index into the edges array
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

    def relabel_nodes(self):
        '''
        Change the labels of the nodes to be more indicative of their nature.

        s: stem
        h: hairpin
        i: interior loop
        m: multiloop
        f: five-prime unpaired
        t: three-prime unpaired
        '''
        hairpins = []
        interior_loops = []
        multiloops = []
        fiveprimes = []
        threeprimes = []

        for d in self.defines.keys():
            if d[0] == 's':
                continue

            if len(self.edges[d]) == 1 and self.defines[d][0] == 0:
                fiveprimes += [d]

            if len(self.edges[d]) == 1 and self.defines[d][1] == self.seq_length+1:
                threeprimes += [d]

            if (len(self.edges[d]) == 1 and 
                self.defines[d][0] != 0 and 
                self.defines[d][1] != self.seq_length+1):
                hairpins += [d]

            if (len(self.edges[d]) == 2 and 
                self.weights[d] == 1 and 
                self.defines[d][0] != 0 and
                self.defines[d][1] != self.seq_length+1):
                multiloops += [d]

            if self.weights[d] == 2:
                interior_loops += [d]

        for d in fiveprimes:
            self.relabel_node(d, 'f1')
        for d in threeprimes:
            self.relabel_node(d, 't1')
        for i,d in enumerate(interior_loops):
            self.relabel_node(d, 'i%d' % (i))
        for i,d in enumerate(multiloops):
            self.relabel_node(d, 'm%d' % (i))
        for i,d in enumerate(hairpins):
            self.relabel_node(d, 'h%d' % (i))

    def find_multiloop_loops(self):
        '''
        Find out which defines are connected in a multiloop.

        @return: A list of sets, where each set contains the names of
                 the elements in a particular multiloop.
        '''
        multis = []
        visited = set()
        for d in self.mloop_iterator():
            if d in visited:
                continue
            # use a really high loop length
            v = self.find_bulge_loop(d , max_length=400)
            v += [d]
            multis += [set(v)]

            for d in v:
                visited.add(d)

        return multis

    def from_fasta(self, fasta_str, dissolve_length_one_stems):
        '''
        Create a bulge graph from a fasta-type file containing the following
        format:

            > id
            ACCGGGG
            ((...))
        '''
        lines = fasta_str.split('\n')
        self.from_dotbracket(lines[2].strip(), dissolve_length_one_stems)
        self.name = lines[0].strip('>')
        self.seq = lines[1].strip()

    def from_dotbracket(self, dotbracket_str, dissolve_length_one_stems = False):
        '''
        Populate the BulgeGraph structure from a dotbracket representation.

        ie: ..((..))..

        :param dotbracket_str: A string containing the dotbracket representation
                               of the structure
        '''
        self.__init__()
        self.dissolve_length_one_stems = dissolve_length_one_stems
        (bulges, stems) = find_bulges_and_stems(dotbracket_str)

        for i in range(len(stems)):
            # one is added to each coordinate to make up for the fact that residues are 1-based
            ss1 = stems[i][0][0]+1
            ss2 = stems[i][0][1]+1
            se1 = stems[i][1][0]+1
            se2 = stems[i][1][1]+1

            self.defines['s%d' % (i)] = [min(ss1,se1), max(ss1,se1), 
                                         min(ss2,se2), max(ss2,se2)]
            self.weights['s%d' % (i)] = 1
        for i in range(len(bulges)):
            bulge = bulges[i]
            self.defines['b%d' % (i)] = [bulge[0]+1, bulge[1]+1]
            self.weights['b%d' % (i)] = 1

        self.dotbracket_str = dotbracket_str
        self.seq_length = len(dotbracket_str)
        self.create_bulge_graph(stems, bulges)
        self.create_stem_graph(stems, len(bulges))
        self.collapse()
        self.relabel_nodes()

    def to_dotbracket(self):
        '''
        Convert the BulgeGraph representation to a dot-bracket string
        and return it.

        @return: A dot-bracket representation of this BulgeGraph
        '''
        out = ['.' for i in xrange(self.seq_length)]
        for s in self.stem_iterator():
            for i in xrange(self.defines[s][0], self.defines[s][1]+1):
                out[i-1] = '('
            for i in xrange(self.defines[s][2], self.defines[s][3]+1):
                out[i-1] = ')'

        return "".join(out)

    def from_bg_string(self, bg_str):
        '''
        Populate this BulgeGraph from the string created by the method
        to_bg_string.

        @param bg_str: The string representation of this BugleGraph.
        '''
        lines = bg_str.split('\n')
        for line in lines:
            line = line.strip()
            parts = line.split()
            if len(parts) == 0:
                #blank line
                continue
            if parts[0] == 'length':
                self.seq_length = int(parts[1])
            elif parts[0] == 'define':
                self.defines[parts[1]] = map(int, parts[2:])
            elif parts[0] == 'connect':
                for p in parts[2:]:
                    self.edges[parts[1]].add(p)
                    self.edges[p].add(parts[1])
            elif parts[0] == 'seq':
                self.seq = parts[1]
            elif parts[0] == 'name':
                self.name = parts[1].strip()

    def stem_iterator(self):
        '''
        Iterate over all of the stem elements.
        '''
        for d in self.defines.keys():
            if d[0] == 's':
                yield d

    def is_single_stranded(self, node):
        '''
        Does this node represent a single-stranded region?

        Single stranded regions are five-prime and three-prime unpaired
        regions, multiloops, and hairpins

        :param node: The name of the node
        :return: True if yes, False if no
        '''
        if node[0] == 'f' or node[0] == 't' or node[0] == 'm' or node[0] == 'h':
            return True
        else:
            return False

    def get_node_dimensions(self, node):
        '''
        Return the dimensions of a node.

        If the node is a stem, then the dimensions will be (l, l) where l is
        the length of the stem.

        If it is single stranded it will be (0, x). Otherwise it will be (x, y).

        :param node: The name of the node
        :return: A pair containing its dimensions
        '''
        if self.is_single_stranded(node):
            return (0, self.defines[node][1] - self.defines[node][0])
        else:
            return (self.defines[node][1] - self.defines[node][0],
                    self.defines[node][3] - self.defines[node][2])

    def stem_length(self, s):
        '''
        Return the number of base pairs that comprise
        one strand of the stem.

        :param s: The name of the stem.
        '''
        return self.defines[s][1] - self.defines[s][0] + 1

    def adjacent_stem_pairs_iterator(self):
        '''
        Iterate over all pairs of stems which are separated by some element.

        This will always yield triples of the form (s1, e1, s2) where s1 and
        s2 are the stem identifiers and e1 denotes the element that separates
        them.
        '''
        for d in self.defines.keys():
            if len(self.edges[d]) == 2:
                edges = list(self.edges[d])

                if edges[0][0] == 's' and edges[1][0] == 's':
                    yield (edges[0], d, edges[1])

    def stem_bp_iterator(self, stem):
        '''
        Iterate over all the base pairs in the stem.
        '''
        d = self.defines[stem]
        stem_length = self.stem_length(stem)

        for i in range(stem_length):
            yield (d[0] + i, d[3] - i)

    def get_side_nucleotides(self, stem, side):
        '''
        Get the nucleotide numbers on the given side of
        them stem. Side 0 corresponds to the 5' end of the
        stem whereas as side 1 corresponds to the 3' side
        of the stem.

        @param stem: The name of the stem
        @param side: Either 0 or 1, indicating the 5' or 3' end of the stem
        @return: A tuple of the nucleotide numbers on the given side of
                 the stem.
        '''
        if side == 0:
            return (self.defines[stem][0], self.defines[stem][3])
        elif side == 1:
            return (self.defines[stem][1], self.defines[stem][2])

        raise Exception("Invalid side (%d) for the stem (%s)." % (stem, side))

    def get_sides(self, s1, b):
        '''
        Get the side of s1 that is next to b.

        s1e -> s1b -> b

        @param s1: The stem.
        @param b: The bulge.
        @return: A tuple indicating which side is the one next to the bulge
                 and which is away from the bulge.
        '''
        s1d = self.defines[s1]
        bd = self.defines[b]

        #print >>sys.stderr, "s1: %s b: %s" % (s1, b)

        for i in xrange(4):
            for k in xrange(len(bd)):
                if s1d[i] == bd[k]:
                    if i == 0 or i == 3:
                        s1b = 0
                    else:
                        s1b = 1

        if s1b == 0:
            s1e = 1
        else:
            s1e = 0

        return (s1b, s1e)

    def stem_iterator(self):
        '''
        Iterator over all of the stems in the structure.
        '''
        for d in self.defines.keys():
            if d[0] == 's':
                yield d
    
    def hloop_iterator(self):
        '''
        Iterator over all of the hairpin in the structure.
        '''
        for d in self.defines.keys():
            if d[0] == 'h':
                yield d

    def mloop_iterator(self):
        '''
        Iterator over all of the multiloops in the structure.
        '''
        for d in self.defines.keys():
            if d[0] == 'm':
                yield d

    def iloop_iterator(self):
        '''
        Iterator over all of the interior loops in the structure.
        '''
        for d in self.defines.keys():
            if d[0] == 'i':
                yield d

    def pairing_partner(self, nucleotide_number):
        '''
        Return the base pairing partner of the nucleotide at position
        nucleotide_number. If this nucleotide is unpaired, return None.

        @param nucleotide_number: The position of the query nucleotide in the
                                  sequence.
        @return: The number of the nucleotide base paired with the one at
                 position nucleotide_number.
        '''
        for d in self.stem_iterator():
            for (r1, r2) in self.stem_bp_iterator(d):
                if r1 == nucleotide_number:
                    return r2
                elif r2 == nucleotide_number:
                    return r1
        return None

    def get_length(self, vertex):
        '''
        Get the minimum length of a vertex.

        If it's a stem, then the result is its length (in base pairs).

        If it's a bulge, then the length is the smaller of it's dimensions.

        @param vertex: The name of the vertex.
        '''
        if vertex[0] == 's':
            return abs(self.defines[vertex][1] - self.defines[vertex][0]) + 1
        else:
            if len(self.edges[vertex]) == 1:
                return self.defines[vertex][1] - self.defines[vertex][0]
            else:
                dims = list(self.get_bulge_dimensions(vertex))
                dims.sort()

                if dims[0] == 0:
                    return dims[1]
                else:
                    return dims[0]
