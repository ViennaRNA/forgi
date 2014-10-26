import sys
import networkx as nx
import unittest
import forgi.utilities.graph as fug

class GraphTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_shortest_cycle(self):
        edges = [(1,2),(2,3),(2,5),(1,4),(4,5)]
        G = nx.Graph()
        G.add_edges_from(edges)

        sc = fug.shortest_cycle(G, 3)
        self.assertEqual(sorted(sc), [1,2,3])
        sc = fug.shortest_cycle(G, 1)
        self.assertEqual(sorted(sc), [1,2,3])
        sc = fug.shortest_cycle(G, 2)
        self.assertEqual(sorted(sc), [1,2,3])

        sc = fug.shortest_cycle(G, 4)
        self.assertEqual(sorted(sc), [1,2,4,5])
        sc = fug.shortest_cycle(G, 5)
        self.assertEqual(sorted(sc), [1,2,4,5])

