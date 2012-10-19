#!/usr/bin/python

import sys
import collections as c
import math
import random

def error_exit(message):
    print >> sys.stderr, message
    sys.exit(1)

# A wrapper for a simple dictionary addition
# Added so that debugging can be amde easier
def add_bulge(bulges, bulge, context, message):
    #print >>sys.stderr,"Adding bulge", context, bulge, message
    #bulge = (context, bulge)
    bulges[context] = bulges.get(context, []) + [bulge]
    return bulges

def any_difference_of_one(stem, bulge):
    '''
    See if there's any difference of one between the two
    ends of the stem [(a,b),(c,d)] and a bulge (e,f)

    @param stem: A couple of couples (2 x 2-tuple) indicating the start and end
                 nucleotides of the stem in the form ((s1, e1), (s2, e2))
    @param bulge: A couple (2-tuple) indicating the first and last position
                  of the bulge.

    @return: True if there is an overlap between the stem nucleotides and the 
                  bulge nucleotides
             False otherwise
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

    @param graph: A dictionary indexed by stem number containing a set
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

    @param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
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

    @param bulges: A list of tuples of the form [(s, e)] where s and e are the 
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

    @param stem_pairs: A list of tuples containing paired base numbers.

    @return: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
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

    @param brackets: A string with the dotplot passed as input to this script.
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

    @param brackets: A string with the dotplot passed as input to this script.
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
        bulges = add_bulge(bulges, (dots_start-1, dots_end), context, "7")
    elif prev == '(':
        print >>sys.stderr, "Unmatched bracket at the end"
        sys.exit(1)
    """
    elif prev == ')':
        bulges = add_bulge(bulges, (i+1, i), context, "8")
    """
    
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
    def __init__(self):
        self.defines = dict()
        self.edges = c.defaultdict(set)
        self.weights = dict()
        self.seq_length = 0

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
            defines_str += "define %s %d %s" % ( key, self.weights[key], " ".join([str(d) for d in self.defines[key]]))
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

    def to_bg_string(self):
        '''
        Output a string representation that can be stored and reloaded.

        '''
        out_str = ''
        out_str += self.get_length_str()
        out_str += self.get_define_str()
        out_str += self.get_connect_str()

        return out_str

    def create_bulge_graph(self, stems, bulges):
        '''
        Find out which stems connect to which bulges

        Stems and bulges which share a nucleotide are considered connected.

        @param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.

        @param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                       numbers of the nucleotides at the start and end of the bulge.
        '''
        for i in range(len(stems)):
            stem = stems[i]
            for j in range(len(bulges)):
                bulge = bulges[j]
                if any_difference_of_one(stem, bulge):
                    self.edges['s%d' % (i)].add('b%d' % (j))

    def create_stem_graph(self, stems, bulge_counter):
        '''
        Determine which stems are connected to each other. A stem can be connected to
        another stem when there is an interior loop with an unpaired nucleotide on
        one side. In this case, a bulge will be created on the other side, but it
        will only consist of the two paired bases around where the unpaired base 
        would be if it existed.

        The defines for these bulges will be printed as well as the connection strings
        for the stems they are connected to.

        @param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.
        @param bulge_counter: The number of bulges that have been encountered so far.

        @return: A dictionary indexed by the number of a stem, containing a set of the 
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
                                        self.defines[bn] = (min(s1, s2)+1, max(s1, s2)+1)
                                        self.weights[bn] = 1

                                        self.edges['s%d' % i].add(bn)
                                        self.edges['s%d' % j].add(bn)

                                        bulge_counter += 1
                                        stem_stems_set.add(j)
                                    stem_stems[i] = stem_stems_set
        return stem_stems

    def find_bulge_loop(self, vertex):
        '''
        Find a set of nodes that form a loop containing the
        given vertex and being no greater than 4 nodes long.

        @param vertex: The vertex to start the search from.
        @return: A list of the nodes in the loop.
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
                    if len(in_path[:depth+1]) > 4:
                        continue
                    else:
                        return in_path[:depth+1]
                
                if key not in visited:
                    to_visit.append((key, depth+1))
        return []

    def collapse(self):
        '''
        If any vertices form a loop, then they are either a bulge region of a fork region. The bulge (interior loop) regions will be condensed into one node.
        '''
        # for all stems of length 1, merge with adjacent bulges
        for key in self.edges.keys():
            if key[0] == 's':
                if self.stem_length(key) == 1:
                    self.dissolve_stem(key)

        for key in self.edges.keys():
            if key[0] == 's':
                while True:
                    loop = self.find_bulge_loop(key)
                    bulge_loop = [v for v in loop if v[0] == 'b' or v[0] == 'x']

                    if len(bulge_loop) > 0:
                        #assert(len(bulge_loop) != 1)
                        self.merge_vertices(bulge_loop)

                    if len(bulge_loop) == 0:
                        break

    def from_dotbracket(self, dotbracket_str):
        (bulges, stems) = find_bulges_and_stems(dotbracket_str)

        for i in range(len(stems)):
            # one is added to each coordinate to make up for the fact that residues are 1-based
            ss1 = stems[i][0][0]+1
            ss2 = stems[i][0][1]+1
            se1 = stems[i][1][0]+1
            se2 = stems[i][1][1]+1

            self.defines['s%d' % (i)] = (min(ss1,se1), max(ss1,se1), 
                                         min(ss2,se2), max(ss2,se2))
            self.weights['s%d' % (i)] = 1
        for i in range(len(bulges)):
            bulge = bulges[i]
            self.defines['b%d' % (i)] = (bulge[0]+1, bulge[1]+1)
            self.weights['b%d' % (i)] = 1

        self.dotbracket_str = dotbracket_str
        self.seq_length = len(dotbracket_str)
        self.create_bulge_graph(stems, bulges)
        self.create_stem_graph(stems, len(bulges))
        self.collapse()

if __name__ == "__main__":
    main()
