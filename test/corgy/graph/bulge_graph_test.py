import unittest, os

import corgy.graph.bulge_graph as cgb

import copy, time

class TestBulgeGraph(unittest.TestCase):
    def setUp(self):
        self.dotbracket = '....((((((....((.......((((.((((.(((...(((((..........)))))...((.......))....)))......))))))))......))...)).))))......(((....((((((((...))))))))...)))........'

    def test_from_dotplot(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket(self.dotbracket)
        print bg.to_bg_string()


