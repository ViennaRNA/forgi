"""
Take a BulgeGraph and return a copy of it with cofold splitpoints inserted.
"""
import logging

from logging_exceptions import log_to_exception


from ._graph_construction import remove_vertex, relabel_node
from ..utilities.exceptions import GraphConstructionError


log = logging.getLogger(__name__)


# This module belongs to the BulgeGraph creation machinenery and is allowed to access its private members
# pylint: disable=protected-access


def split_at_cofold_cutpoints(bg, cutpoints):
    """
    Multiple sequences should not be connected along the backbone.

    We have constructed the bulge graph, as if they were connected along the backbone, so
    now we have to split it.
    """

    for splitpoint in cutpoints:
        element_left = bg.get_node_from_residue_num(splitpoint)
        element_right = bg.get_node_from_residue_num(splitpoint + 1)
        if element_left[0] in "ft" or element_right[0] in "ft":
            if element_left[0] == "t" and element_left[0] != "t":
                continue  # Splitpoint already implemented
            elif element_right[0] == "f" and element_left[0] != "f":
                continue  # Splitpoint already implemented
            else:
                # No cofold structure. First sequence is disconnected from rest
                e = GraphConstructionError("Cannot create BulgeGraph. Found two sequences not "
                                           "connected by any base-pair.")
                with log_to_exception(log, e):
                    log.error("Trying to split between %s and %s",
                              element_left, element_right)
                raise e
            return
        elif element_left[0] == "i" or element_right[0] == "i":
            _split_interior_loop(bg, splitpoint, element_left, element_right)
        elif element_left != element_right:
            _split_between_elements(
                        bg, splitpoint, element_left, element_right)
        elif element_left[0] == "s":
            _split_inside_stem(bg, splitpoint, element_left)
        else:
            _split_inside_loop(bg, splitpoint, element_left)
        bg._node_to_resnum = {}
    if not _is_connected(bg):
        raise GraphConstructionError("Cannot create BulgeGraph. Found two sequences not connected by any "
                                     " base-pair.")


def _is_connected(bg):
    if not bg.defines:
        return True  # We define an empty Graph as connected.
    start_node = list(bg.defines.keys())[0]
    known_nodes = set([start_node])
    pending = list(bg.edges[start_node])
    while pending:
        next_node = pending.pop()
        if next_node in known_nodes:
            continue
        pending.extend(bg.edges[next_node])
        known_nodes.add(next_node)
    log.info("Testing connectivity: connected component =?= all nodes:\n%s =?= %s",
             list(sorted(known_nodes)), list(sorted(set(bg.defines.keys()))))
    return known_nodes == set(bg.defines.keys())


def _split_between_elements(bg, splitpoint, element_left, element_right):
    log.debug("Before splitting between %s and %s at %s: %s", element_left,
              element_right, splitpoint, bg.defines)
    if element_left[0] in "mh":
        next3 = _next_available_element_name(bg, "t")
        relabel_node(bg, element_left, next3)
        if element_left[0] != "h":
            _remove_edge(bg, next3, element_right)
    elif element_right[0] in "mh":
        next5 = _next_available_element_name(bg, "f")
        relabel_node(bg, element_right, next5)
        if element_right[0] != "h":
            _remove_edge(bg, next5, element_left)
    else:
        assert element_left[0] == "s" and element_right[0] == "s"
        # Zero-length i or m element!
        connections = bg.edges[element_left] & bg.edges[element_right]
        # If multiple ml connections exist, find the one which is 0-len
        # Remove the highest numbered element first
        # We know that adjacent stems are connected
        # by either an i-loop or a 0-len ml.
        for connection in reversed(sorted(connections)):
            if connection[0] == "i":
                break
            if not bg.defines[connection]:
                assert splitpoint in bg.defines[element_left]
                break
        else:
            raise GraphConstructionError(
                        "Cannot split at cofold cutpoint. "
                        "No (suitable) connection between "
                        "{} and {}.".format(element_left, element_right))
        # pylint: disable=undefined-loop-variable
        if connection[0] == "m":
            # Just remove it without replacement
            remove_vertex(bg, connection)
        else:
            assert connection[0] == "i"
            # Replace i by ml (this is then located on the other strand than the splitpoint)
            next_ml = _next_available_element_name(bg, "m")
            assert next_ml not in bg.defines
            relabel_node(bg, connection, next_ml)


def _split_inside_loop(bg, splitpoint, element):
    if element[0] in "hm":
        from_, to_ = bg.defines[element]
        stem_left = bg.get_node_from_residue_num(from_ - 1)
        stem_right = bg.get_node_from_residue_num(to_ + 1)

        next3 = _next_available_element_name(bg, "t")
        next5 = _next_available_element_name(bg, "f")
        bg.defines[next3] = [from_, splitpoint]
        bg.defines[next5] = [splitpoint + 1, to_]
        _add_edge(bg, stem_left, next3)
        _add_edge(bg, next5, stem_right)
        remove_vertex(bg, element)
    else:
        assert False


def _split_inside_stem(bg, splitpoint, element):
    assert element[0] == "s"
    log.debug("Split inside stem %s at %s", element, splitpoint)
    if splitpoint == bg.defines[element][1]:
        # Nothing needs to be done. 2 strands split at end
        log.debug("Nothing to do")
        return
    if splitpoint < bg.defines[element][1]:
        # Splitpoint in forward strand:
        define1 = [bg.defines[element][0], splitpoint,
                   bg.pairing_partner(splitpoint), bg.defines[element][3]]
        define2 = [splitpoint + 1, bg.defines[element][1],
                   bg.defines[element][2], bg.pairing_partner(splitpoint + 1)]
        log.debug("Split in forward strand")
    else:
        # Splitpoint in backwards strand:
        define1 = [bg.defines[element][0], bg.pairing_partner(splitpoint + 1),
                   splitpoint + 1, bg.defines[element][3]]
        define2 = [bg.pairing_partner(splitpoint), bg.defines[element][1],
                   bg.defines[element][2], splitpoint]
        log.debug("Split in backwards strand")
    edges1 = []
    edges2 = []

    for edge in bg.edges[element]:
        log.debug("Checking edge %s with define %s connected to %s",
                  edge, bg.defines[edge], bg.edges[edge])
        flank = bg.flanking_nucleotides(edge)
        found = 0
        if define1[0] in flank or define1[3] in flank:
            edges1.append(edge)
            found += 1
        if define2[1] in flank or define2[2] in flank:
            edges2.append(edge)
            found += 1
        if found != 1:
            log.error("For stem %s with define %s and cutpoint %s:",
                      element, bg.defines[element], splitpoint)
            log.error("Edge %s, with flanking nts %s, define1 %s, "
                      "define2 %s", edge, bg.flanking_nucleotides(edge),
                      define1, define2)
            assert False
    remove_vertex(bg, element)
    next_s1 = _next_available_element_name(bg, "s")
    bg.defines[next_s1] = define1
    next_m = _next_available_element_name(bg, "m")
    bg.defines[next_m] = []
    next_s2 = _next_available_element_name(bg, "s")
    bg.defines[next_s2] = define2

    for e1 in edges1:
        bg.edges[e1].add(next_s1)
    for e2 in edges2:
        bg.edges[e2].add(next_s2)
    edges1.append(next_m)
    edges2.append(next_m)
    bg.edges[next_s1] = set(edges1)
    bg.edges[next_s2] = set(edges2)
    bg.edges[next_m] = set([next_s1, next_s2])


def _next_available_element_name(bg, element_type):
    """
    :param element_type: A single letter ("t", "f", "s"...)
    """
    i = 0
    while True:
        name = "{}{}".format(element_type, i)
        if name not in bg.defines:
            return name
        i += 1


def _remove_edge(bg, from_element, to_element):
    bg.edges[from_element].remove(to_element)
    bg.edges[to_element].remove(from_element)


def _add_edge(bg, from_element, to_element):
    bg.edges[from_element].add(to_element)
    bg.edges[to_element].add(from_element)


def _split_interior_loop_at_side(bg, splitpoint, strand, other_strand, stems):
    """
    Called by _split_at_cofold_cutpoints
    """
    next_ml = _next_available_element_name(bg, "m")
    next_a = _next_available_element_name(bg, "t")
    next_b = _next_available_element_name(bg, "f")

    if other_strand[0] > other_strand[1]:
        bg.defines[next_ml] = []
    else:
        bg.defines[next_ml] = other_strand
    _add_edge(bg, next_ml, stems[0])
    _add_edge(bg, next_ml, stems[1])

    if splitpoint >= strand[0]:
        bg.defines[next_a] = [strand[0], splitpoint]
        _add_edge(bg, next_a, stems[0])
    if splitpoint < strand[1]:
        bg.defines[next_b] = [splitpoint + 1, strand[1]]
        _add_edge(bg, next_b, stems[1])


def _split_interior_loop(bg, splitpoint, element_left, element_right):
    if element_left[0] == "i":
        iloop = element_left
    elif element_right[0] == "i":
        iloop = element_right
    else:
        assert False
    c = bg.connections(iloop)
    s1 = bg.defines[c[0]]
    s2 = bg.defines[c[1]]
    forward_strand = [s1[1] + 1, s2[0] - 1]
    back_strand = [s2[3] + 1, s1[2] - 1]
    if forward_strand[0] - 1 <= splitpoint <= forward_strand[1]:
        # Split forward strand, relabel backwards strand to multiloop.
        _split_interior_loop_at_side(bg, splitpoint,
                                     forward_strand, back_strand, c)
    elif back_strand[0] - 1 <= splitpoint <= back_strand[1]:
        _split_interior_loop_at_side(bg, splitpoint,
                                     back_strand, forward_strand, [c[1], c[0]])
    else:
        assert False
    remove_vertex(bg, iloop)
