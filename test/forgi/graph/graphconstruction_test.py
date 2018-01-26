import unittest
import logging

import forgi.utilities.stuff as fus
from forgi.graph._graph_construction import _BulgeGraphConstruction

log=logging.getLogger(__name__)

class ConstructionSimple(unittest.TestCase):
    def test_construction_rod(self):
        db = "(((.(((...)))..)))"
        tuples = fus.pairtable_to_tuples(fus.dotbracket_to_pairtable(db))
        c = _BulgeGraphConstruction(tuples)
        self.assertEqual(c.defines, {"s0":[1,3,16,18], "i0":[4,4,14,15],
                                     "s1":[5,7,11,13], "h0":[8,10]})
        self.assertEqual(c.edges, {"s0":{"i0"}, "i0":{"s0", "s1"},
                                   "s1":{"i0","h0"},"h0":{"s1"}})
    def test_construction_ml(self):
        db = "(.(.).(.).).."
        tuples = fus.pairtable_to_tuples(fus.dotbracket_to_pairtable(db))
        c = _BulgeGraphConstruction(tuples)
        self.assertEqual(c.defines, {"s0":[1,1,11,11], "s1":[3,3,5,5],
                                     "s2":[7,7,9,9], "m0":[2,2], "m1":[6,6],
                                     "m2":[10,10], "h0":[4,4], "h1":[8,8],
                                     "t0":[12,13]})
        self.assertEqual(c.edges["s1"], {"m0", "m1", "h0"})

class ConstructionZeroLength(unittest.TestCase):
    def test_constr_0_len_ml(self):
        db = ".((.).(.)).."
        tuples = fus.pairtable_to_tuples(fus.dotbracket_to_pairtable(db))
        c = _BulgeGraphConstruction(tuples)
        self.assertEqual(c.defines, {"f0":[1,1], "s0":[2,2,10,10], "s1":[3,3,5,5],
                                     "s2":[7,7,9,9], "m0":[], "m1":[6,6],
                                     "m2":[], "h0":[4,4], "h1":[8,8],
                                     "t0":[11,12]})
        self.assertEqual(c.edges["m0"], {"s0", "s1"})

class PrivateMemberTests(unittest.TestCase):
    def test_get_sides_plus(self):
        db = ".((.).(.)).."
        tuples = fus.pairtable_to_tuples(fus.dotbracket_to_pairtable(db))
        c = _BulgeGraphConstruction(tuples)

        p1 = c._get_sides_plus('s0', 'm0')
        self.assertEquals(p1[0], 1)

        p1 = c._get_sides_plus('s0', 'm2')
        self.assertEquals(p1[0], 2)

        p1 = c._get_sides_plus('s1', 'm0')
        self.assertEquals(p1[0], 0)
