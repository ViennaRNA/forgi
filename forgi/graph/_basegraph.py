from collections import defaultdict
from pprint import pformat
import logging
import itertools as it

from ..utilities.exceptions import GraphConstructionError, GraphIntegrityError

log = logging.getLogger(__name__)

class BaseGraph(object):
    """
    A Base-class for the BulgeGraph and BulgeGraphConstruction.

    In contrast to BulgeGraph, this class does not enforce any nomenclature
    on the nodes.

    It does not support backbone breaks!

    It has no sequence.
    """
    def __init__(self):
        self.defines = []
        self.edges = defaultdict(set)

    def connections(self, bulge):
        """
        :param g: Graph-like: A BulgeGraph or BulgeGraphConstruction.
        """
        def sort_key(x):
            if len(self.defines[x]) > 0 and self.defines[x][0] == 1:
                # special case for stems at the beginning since there is no
                # adjacent nucleotide 0
                return 0
            return self.define_a(x)[0]
        connections = list(self.edges[bulge])
        connections.sort(key=sort_key)
        return connections

    def define_a(self, elem):
        if self.defines[elem] == []:
            return self._define_a_zerolength(elem)
        else:
            return self._define_a_nonzero(elem)

    def _define_a_nonzero(self, elem):
        define = self.defines[elem]
        log.debug("Define_a nonzero of BaseGraph called for %s with define %s", elem, define)
        new_def = []
        for i in range(0,len(define),2):
            new_def.append(max(define[i]-1, 1))
            new_def.append(define[i+1]+1)
        log.debug("Define_a nonzero returning %s", new_def)
        return new_def


    def flanking_nucleotides(self, d):
        '''
        Return the nucleotides directly flanking an element.

        :param d: the name of the element
        :return: a list of nucleotides
        '''
        set_adjacent = set(self.define_a(d))
        set_not_adjacent = set(self.defines[d])
        flanking = list(sorted(set_adjacent - set_not_adjacent))
        log.debug("Flanking nts are %s - %s = %s", set_adjacent, set_not_adjacent, flanking)
        return flanking

    def _define_a_zerolength(self, elem): #TODO speed-up by caching
        """
        Return the define with adjacent nucleotides for a zero-length element.

        Hereby we define that in cases of ambiuigity, the alphabetically first
        zero-length element comes at the lowest nucleotide position etc.

        :param elem: An element, e.g. "m0
        """
        if self.defines[elem]!=[]:
            raise ValueError("{} does not have zero length".format(elem))
        edges = self.edges[elem]
        if len(edges)==1: #Hairpin
            stem, = edges
            define = self.defines[stem]
            if define[2]==define[1]+1:
                return [define[1], define[2]]
            else:
                raise GraphIntegrityError("Very strange zero-length hairpin {} "
                                          "(not?) connected to {}".format(elem, stem))
        elif len(edges)==2:
            stem1, stem2 = edges
            #See if this is the only element connecting the two stems.
            connections = self.edges[stem1] & self.edges[stem2]
            log.debug("Stems %s and %s, connected by %s have the following common edges: %s with defines %s",
                      stem1, stem2, elem, connections, list(map(lambda x: self.defines[x], connections)))
            zl_connections = []
            for conn in connections:
                if self.defines[conn]==[]:
                    zl_connections.append(conn)
            assert elem in zl_connections
            #We DEFINE the 0-length connections to be sorted alphabetically by position
            zl_connections.sort()
            zl_coordinates = self._zerolen_defines_a_between(stem1, stem2)

            if len(zl_connections)!=len(zl_coordinates):
                raise GraphIntegrityError("Expecting stems {} and {} to have {} zero-length "
                                          "connections at nucleotide positions {}, however, "
                                          "found {} elements: {}".format(stem1, stem2,
                                                            len(zl_coordinates),
                                                            zl_coordinates,
                                                            len(zl_connections),
                                                            zl_connections))
            zl_coordinates = list(zl_coordinates)
            zl_coordinates.sort()
            i = zl_connections.index(elem)
            return list(zl_coordinates[i])
        else:
            raise GraphIntegrityError("Very strange zero length bulge {} with more than 2 adjacent "
                                      "elements: {}.".format(elem, edges))

    def _zerolen_defines_a_between(self, stem1, stem2):
        zl_coordinates = set()
        for k, l in it.product(range(4), repeat=2):
            if abs(self.defines[stem1][k]-self.defines[stem2][l])==1:
                d = [self.defines[stem1][k], self.defines[stem2][l]]
                d.sort()
                log.debug("Zero-length element found: %s", d)
                zl_coordinates.add(tuple(d))
        return zl_coordinates


    def _get_sides_plus(self, s1, b):
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
        if b not in self.edges[s1]:
            raise ValueError("_get_sides_plus expects stem to be connected to bulge!")

        s1d = self.defines[s1]
        bd = self.defines[b]

        if len(bd) == 0:
            bd = self._define_a_zerolength(b)
            bd[0]+=1
            bd[1]-=1

        # before the stem on the 5' strand
        if s1d[0] - bd[1] == 1:
            return (0, 1)
        # after the stem on the 5' strand
        elif bd[0] - s1d[1] == 1:
            return (1, 0)
        # before the stem on the 3' strand
        elif s1d[2] - bd[1] == 1:
            return (2, 1)
        # after the stem on the 3' strand
        elif bd[0] - s1d[3] == 1:
            return (3, 0)

        raise GraphIntegrityError("Faulty bulge {}:{} connected to {}:{}".format(b, bd, s1, s1d))

    def get_node_from_residue_num(self, base_num):
        """
        Iterate over the defines and see which one encompasses this base.
        """
        seq_id=False
        for key in self.defines.keys():
            for r in self.define_range_iterator(key):
                if base_num >= r[0] and base_num <= r[1]:
                    return key
        raise LookupError("Base number {} not found in the defines {}.".format(base_num, self.defines))

    def define_range_iterator(self, node, adjacent=False):
        """
        Return the ranges of the nucleotides in the define.

        In other words, if a define contains the following: [1,2,7,8]
        The ranges will be [1,2] and [7,8].

        :param adjacent: Use the nucleotides in the neighboring element which
                         connect to this element as the range starts and ends.
        :return: A list of two-element lists
        """
        if adjacent:
            define = self.define_a(node)
        else:
            define = self.defines[node]

        if define:
            yield [ define[0], define[1] ]
            if len(define)>2:
                yield [ define[2], define[3] ]
