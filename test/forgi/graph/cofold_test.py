import unittest
import logging

import forgi.graph.bulge_graph as fgb
import forgi.graph._cofold as fgc
from forgi.graph._graph_construction import remove_vertex
log = logging.getLogger(__name__)


class CofoldPrivateMemberTest(unittest.TestCase):
    def test_split_interior_loop_at_side(self):
        # Normal, forward strand
        db = "(...(...)..)"
        bg = fgb.BulgeGraph.from_dotbracket(db)
        fgc._split_interior_loop_at_side(bg, 2, [2, 4], [10, 11], ["s0", "s1"])
        self.assertEqual(bg.defines["t0"], [2, 2])
        self.assertEqual(bg.defines["f0"], [3, 4])
        self.assertEqual(bg.defines["m0"], [10, 11])

        # No nt at back, forward strand
        db = "(...(...))"
        bg = fgb.BulgeGraph.from_dotbracket(db)
        fgc._split_interior_loop_at_side(bg, 2, [2, 4], [10, 9], ["s0", "s1"])
        self.assertEqual(bg.defines["t0"], [2, 2])
        self.assertEqual(bg.defines["f0"], [3, 4])
        self.assertEqual(bg.defines["m0"], [])
        # Normal, backwards strand
        db = "(...(...)..)"
        bg = fgb.BulgeGraph.from_dotbracket(db)
        fgc._split_interior_loop_at_side(
            bg, 10, [10, 11], [2, 4], ["s1", "s0"])
        self.assertEqual(bg.defines["t0"], [10, 10])
        self.assertEqual(bg.defines["f0"], [11, 11])
        self.assertEqual(bg.defines["m0"], [2, 4])

    def test_is_connected(self):
        db = "(...)(...)"
        bg = fgb.BulgeGraph.from_dotbracket(db)
        self.assertTrue(fgc._is_connected(bg))
        remove_vertex(bg, "m0")
        log.error(bg.edges)
        self.assertFalse(fgc._is_connected(bg))

    def test_split_inside_stem(self):
        db = "(.((...)).)"
        bg = fgb.BulgeGraph.from_dotbracket(db)
        fgc._split_inside_stem(bg, 8, "s1")
        self.assertEqual(bg.defines["s0"], [1, 1, 11, 11])
        self.assertEqual(bg.defines["s1"], [3, 3, 9, 9])
        self.assertEqual(bg.defines["s2"], [4, 4, 8, 8])
        self.assertEqual(bg.defines["m0"], [])
