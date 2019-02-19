import sys
import logging
import itertools as it
from collections import defaultdict
from pprint import pformat

from ._basegraph import BaseGraph

log = logging.getLogger(__name__)


class _BulgeGraphConstruction(BaseGraph):
    """
    An intermediate object that is used during BulgeGraphConstruction.
    It is responsible ONLY for the structure without cofold cutpoints,
    holds no sequence information and uses arbitrary graph labels during
    it's intialization!
    """

    def __init__(self, tuples):  # pylint: disable=W0231
        self.defines = {}
        self.edges = defaultdict(set)
        self.weights = {}
        self._name_counter = 0
        self.from_tuples(tuples)

    def from_tuples(self, tuples):
        """
        Create a bulge_graph from a list of pair tuples. Unpaired
        nucleotides have a pairing partner of 0.
        """
        stems = []
        bulges = []

        tuples.sort()  # We move along the backbone
        tuples = iter(tuples)
        (t1, t2) = next(tuples)

        prev_from = t1
        prev_to = t2

        start_from = prev_from
        start_to = prev_to
        last_paired = prev_from

        for t1, t2 in tuples:
            (from_bp, to_bp) = (t1, t2)

            if abs(to_bp - prev_to) == 1 and prev_to != 0:  # adjacent basepairs on 3' strand
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

        log.debug("from_tuples: stems %s, bulges %s", stems, bulges)
        self.from_stems_and_bulges(stems, bulges)

    def from_stems_and_bulges(self, stems, bulges):
        """
        Create the graph from the list of stems and bulges.

        :param stems: A list of tuples of two two-tuples, each containing the start
                      and end nucleotides of each strand of the stem.
        :param bulges: A list of tuples containing the starts and ends of the
                       of the bulge regions.
        :return: Nothing, just make the bulgegraph
        """
        assert self.defines == {}
        assert self.edges == defaultdict(set)
        for i in range(len(stems)):
            # one is added to each coordinate to make up for the fact that residues are 1-based
            ss1 = stems[i][0][0] + 1
            ss2 = stems[i][0][1] + 1
            se1 = stems[i][1][0] + 1
            se2 = stems[i][1][1] + 1
            log.debug("stem define not sorted: %s %s %s %s",
                      ss1, ss2, se1, se2)
            log.debug("self.defines %s", self.defines)

            self.defines['y%d' % (i)] = [min(ss1, se1), max(ss1, se1),
                                         min(ss2, se2), max(ss2, se2)]
            self.weights['y%d' % (i)] = 1

        for i in range(len(bulges)):
            bulge = bulges[i]
            self.defines['b%d' % (i)] = sorted([bulge[0] + 1, bulge[1] + 1])
            self.weights['b%d' % (i)] = 1

        log.debug("from_stems_and_bulges: %s; %s", self.defines, self.edges)
        self.create_bulge_graph(stems, bulges)
        log.debug("after create_bulge_graph: DEFINES:\n %s;\n EDGES:\n %s", pformat(
            self.defines), pformat(self.edges))
        self.create_stem_graph(stems, len(bulges))
        log.debug("after create_stem_graph: DEFINES \n%s;\nEDGES \n%s",
                  pformat(self.defines), pformat(self.edges))
        self.collapse()
        log.debug("after collapse: DEFINES:\n %s;\n EDGES:\n %s",
                  pformat(self.defines), pformat(self.edges))
        self.sort_defines()
        log.debug("after _sort_defines: DEFINES:\n%s;\n EDGES:\n%s",
                  pformat(self.defines), pformat(self.edges))
        self.relabel_nodes()
        log.debug("after relabel_nodes: DEFINES:\n %s;\n EDGES:\n %s",
                  pformat(self.defines), pformat(self.edges))
        self.remove_degenerate_nodes()

    def remove_degenerate_nodes(self):
        """
        For now just remove all hairpins that have no length.
        """
        to_remove = []
        for d in self.defines:
            if d[0] == 'h' and len(self.defines[d]) == 0:
                to_remove += [d]

        for r in to_remove:
            remove_vertex(self, r)

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

                    (f1, t1) = (int(self.defines[key][j]), int(
                        self.defines[key][j + 1]))

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

                        (f2, t2) = (int(self.defines[key][k]), int(
                            self.defines[key][k + 1]))

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
                self.defines[new_vertex] = self.defines.get(
                    new_vertex, []) + self.defines[v]

            self.weights[new_vertex] += 1

            # remove the old vertex, since it's been replaced by new_vertex
            remove_vertex(self, v)
            self.reduce_defines()

        # self.weights[new_vertex] = 2
        for connection in connections:
            self.edges[new_vertex].add(connection)
            self.edges[connection].add(new_vertex)

        return new_vertex

    def get_vertex(self, name=None):
        """
        Return a new unique vertex name starting with `x`.

        At this stage stems and bulges are not distinguished
        """

        if name is None:
            name = "x{}".format(self._name_counter)
            self._name_counter += 1

        return name

    def collapse(self):
        """
        If any vertices form a loop, then they are either a bulge region or
        a fork region. The bulge (interior loop) regions will be condensed
        into one node.
        """

        new_vertex = True
        while new_vertex:
            new_vertex = False
            bulges = [k for k in self.defines if k[0] != 'y']

            for (b1, b2) in it.combinations(bulges, r=2):
                if self.edges[b1] == self.edges[b2] and len(self.edges[b1]) == 2:
                    connections = self.connections(b1)

                    all_connections = [sorted((self._get_sides_plus(connections[0], b1)[0],
                                               self._get_sides_plus(connections[0], b2)[0])),
                                       sorted((self._get_sides_plus(connections[1], b1)[0],
                                               self._get_sides_plus(connections[1], b2)[0]))]

                    if all_connections == [[1, 2], [0, 3]]:
                        # interior loop
                        log.debug("Collapsing %s and %s", b1, b2)
                        self.merge_vertices([b1, b2])
                        new_vertex = True
                        break

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

        :returns: None
        """
        # print "stems:", stems
        for i, j in it.combinations(range(len(stems)), 2):
            for k1, k2, l1, l2 in it.product(range(2), repeat=4):
                s1 = stems[i][k1][l1]
                s2 = stems[j][k2][l2]
                if k1 == 1 and stems[i][0][l1] == stems[i][1][l1]:
                    continue
                if k2 == 1 and stems[j][0][l2] == stems[j][1][l2]:
                    continue
                if abs(s1 - s2) == 1:
                    bn = 'b{}'.format(bulge_counter)
                    log.debug("Adding bulge %s between %s and %s. (%s is next to %s ) k1 %s, k2 %s, l1 %s, l2 %s",
                              bn, stems[i], stems[j], s1, s2, k1, k2, l1, l2)
                    # self.defines[bn] = [min(s1, s2)+1, max(s1, s2)+1]
                    self.defines[bn] = []
                    self.weights[bn] = 1

                    self.edges['y{}'.format(i)].add(bn)
                    self.edges[bn].add('y{}'.format(i))

                    self.edges['y{}'.format(j)].add(bn)
                    self.edges[bn].add('y{}'.format(j))

                    bulge_counter += 1

        for d in list(self.defines.keys()):  # 0-nt Hairpins
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

        return

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

        seq_length = 0
        for d in self.defines:
            for r in self.define_range_iterator(d):
                seq_length += r[1] - r[0] + 1

        for d in self.defines.keys():
            if d[0] == 'y' or d[0] == 's':
                stems += [d]
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

            if len(self.edges[d]) == 1 and self.defines[d][1] == seq_length:
                threeprimes += [d]
                continue

            if (len(self.edges[d]) == 1 and
                self.defines[d][0] != 1 and
                    self.defines[d][1] != seq_length):
                hairpins += [d]
                continue

            if d[0] == 'm' or (d[0] != 'i' and len(self.edges[d]) == 2 and
                                       self.weights[d] == 1 and
                                       self.defines[d][0] != 1 and
                                       self.defines[d][1] != seq_length):
                multiloops += [d]
                continue

            if d[0] == 'i' or self.weights[d] == 2:
                interior_loops += [d]

        stems.sort(key=self.compare_stems)
        hairpins.sort(key=self.compare_hairpins)
        multiloops.sort(key=self.compare_bulges)
        interior_loops.sort(key=self.compare_stems)

        if fiveprimes:
            d, = fiveprimes
            relabel_node(self, d, 'f0')
        if threeprimes:
            d, = threeprimes
            relabel_node(self, d, 't0')
        for i, d in enumerate(stems):
            relabel_node(self, d, 's%d' % (i))
        for i, d in enumerate(interior_loops):
            relabel_node(self, d, 'i%d' % (i))
        for i, d in enumerate(multiloops):
            relabel_node(self, d, 'm%d' % (i))
        for i, d in enumerate(hairpins):
            relabel_node(self, d, 'h%d' % (i))

    def compare_stems(self, b):
        """
        A function that can be passed in as the key to a sort.
        """
        return (self.defines[b][0], 0)

    def compare_bulges(self, b):
        """
        A function that can be passed in as the key to a sort.

        Compares based on the nucleotide number
        (using define_a to allow for sorting 0-length MLs)
        """
        return self.define_a(b)

    def compare_hairpins(self, b):
        connections = self.connections(b)

        return (self.defines[connections[0]][1], sys.maxsize)

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


def remove_vertex(bg, v):
    """
    Delete a node after merging it with another

    :param v: The name of the node
    """
    # delete all edges to this node
    for key in bg.edges[v]:
        bg.edges[key].remove(v)

    for edge in bg.edges:
        if v in bg.edges[edge]:
            bg.edges[edge].remove(v)

    # delete all edges from this node
    del bg.edges[v]
    del bg.defines[v]


def relabel_node(bg, old_name, new_name):
    """
    Change the name of a node.

    param old_name: The previous name of the node
    param new_name: The new name of the node
    """
    #log.debug("Relabelling node {} to {}".format(old_name, new_name))
    # replace the define name
    define = bg.defines[old_name]

    del bg.defines[old_name]
    bg.defines[new_name] = define

    # replace the index into the edges array
    edge = bg.edges[old_name]
    del bg.edges[old_name]
    bg.edges[new_name] = edge

    # replace the name of any edge that pointed to old_name
    for k in bg.edges.keys():
        new_edges = set()
        for e in bg.edges[k]:
            if e == old_name:
                new_edges.add(new_name)
            else:
                new_edges.add(e)
        bg.edges[k] = new_edges
