from __future__ import print_function

import networkx as nx


def shortest_cycle(G, v):
    '''
    Calculate one of the shortest cycles which includes the node v.

    @param G: A networkx graph
    @param v: The node to include in the cycle
    '''
    lengths = []

    for n in list(G.neighbors(v)):
        G.remove_edge(v, n)
        try:
            p = nx.shortest_path(G, v, n)
            l = nx.shortest_path_length(G, v, n)
        except nx.exception.NetworkXNoPath:
            G.add_edge(v, n)
            continue

        lengths += [(l, p)]
        G.add_edge(v, n)

    if len(lengths) > 0:
        return min(lengths)[1]
    else:
        return []


if __name__ == '__main__':
    edges = [(1, 2), (2, 3), (2, 5), (1, 4), (4, 5), (1, 3)]
    G = nx.Graph()
    G.add_edges_from(edges)

    print(shortest_cycle(G, 1))
    print(shortest_cycle(G, 3))
    print(shortest_cycle(G, 4))
