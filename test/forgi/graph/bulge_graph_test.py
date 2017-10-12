from __future__ import absolute_import
from __future__ import print_function

from builtins import zip
from builtins import next
from builtins import range

import unittest
import itertools as it
from pprint import pprint
import collections as col
import sys
import os
import logging

import forgi.graph.bulge_graph as fgb
import forgi.utilities.debug as fud
from forgi.utilities.exceptions import GraphConstructionError
import forgi.utilities.stuff as fus

log=logging.getLogger(__name__)


class GraphVerification(unittest.TestCase):
    def check_for_overlapping_defines(self, bg):
        """
        Check to make sure none of the defines overlap.
        """
        for d1, d2 in it.combinations(bg.defines.keys(), 2):
            for dx in bg.defines[d1]:
                for dy in bg.defines[d2]:
                    self.assertNotEqual(dx, dy)

    def check_for_all_nucleotides(self, bg):
        """
        Check to make sure that the bulge_graph covers each nucleotide
        in the structure.
        """
        nucs = [False] * bg.seq_length
        for d in bg.defines.keys():
            for r in bg.define_residue_num_iterator(d):
                nucs[r-1] = True

        for i, n in enumerate(nucs):
            self.assertTrue(n)

    def check_node_labels(self, bg):
        """
        There should be only six types of nodes in the graph. The internal
        representation sometimes uses nodes that start with 'x' or 'y' as
        intermediates, but these should always be removed.
        """
        for k in bg.defines:
            self.assertTrue(k[0] in ['s', 'h', 'i', 'm', 't', 'f'])

    def check_graph_integrity(self, bg):
        self.check_node_labels(bg)
        self.check_for_all_nucleotides(bg)
        self.check_for_overlapping_defines(bg)


class BulgeGraphCofoldPrivateMemberTest(GraphVerification):
    def test_split_interior_loop_at_side(self):
        #Normal, forward strand
        db = "(...(...)..)"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        cg._split_interior_loop_at_side(2, [2,4], [10,11], ["s0", "s1"])
        self.assertEqual(cg.defines["t0"], [2,2])
        self.assertEqual(cg.defines["f0"], [3,4])
        self.assertEqual(cg.defines["m0"], [10,11])
        #No nt at back, forward strand
        db = "(...(...))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        cg._split_interior_loop_at_side( 2, [2,4], [10,9], ["s0", "s1"])
        self.assertEqual(cg.defines["t0"], [2,2])
        self.assertEqual(cg.defines["f0"], [3,4])
        self.assertEqual(cg.defines["m0"], [])
        #Normal, backwards strand
        db = "(...(...)..)"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        cg._split_interior_loop_at_side( 10, [10,11], [2,4], ["s1", "s0"])
        self.assertEqual(cg.defines["t0"], [10,10])
        self.assertEqual(cg.defines["f0"], [11,11])
        self.assertEqual(cg.defines["m0"], [2,4])
    def test_is_connected(self):
        db = "(...)(...)"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertTrue(cg._is_connected())
        cg.remove_vertex("m0")
        log.error(cg.edges)
        self.assertFalse(cg._is_connected())

class BulgeGraphCofoldOverallTest(GraphVerification):

    def test_cutpoint_in_stem_f(self):
        db = "(((&(((...))))))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "s1", "h0", "m0"]))
        self.assertEqual(cg.defines["s0"], [1,3,13,15])
        self.assertEqual(cg.defines["s1"], [4,6,10,12])
        self.assertEqual(cg.edges["s0"], set(["m0"]))
        self.assertEqual(cg.edges["s1"], set(["m0", "h0"]))
    def test_cutpoint_in_stem_b(self):
        db = "((((((...)))&)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "s1", "h0", "m0"]))
        self.assertEqual(cg.defines["s0"], [1,3,13,15])
        self.assertEqual(cg.defines["s1"], [4,6,10,12])
        self.assertEqual(cg.edges["s0"], set(["m0"]))
        self.assertEqual(cg.edges["s1"], set(["m0", "h0"]))

    def test_cutpoint_in_ml(self):
        db = "(((.(((...)))..&..(((...))).)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "h0", "f0", "t0", "s2", "h1", "m2"]))

    def test_cutpoint_in_il(self):
        db = "(((..&..(((...))))))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "f0", "t0", "s1", "h0", "m0"]))

    def test_cutpoint_in_h(self):
        db = "(((..&..)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "f0", "t0"]))

    def test_cutpoint_in_f(self):
        db = "...&..."
        with self.assertRaises(ValueError):
            cg = fgb.BulgeGraph(dotbracket_str=db)

    def test_cutpoint_in_t(self):
        db = "(((...)))...&..."
        with self.assertRaises(ValueError):
            cg = fgb.BulgeGraph(dotbracket_str=db)

    def test_cutpoint_between_i_s(self):
        db = "(((...&(((...))))))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "t0", "s1", "h0", "m0"]))

    def test_cutpoint_between_m_s(self):
        db = "(((.(((...)))..&(((...))).)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "s1", "s2", "m0", "m2", "t0", "h0", "h1"]))

    def test_cutpoint_between_s_i0_s(self):
        db = "(((&(((...))).)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "h0"]))

    def test_cutpoint_between_s_m0_s(self):
        db = "(((.(((...)))&(((...))).)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "h0", "s2", "h1", "m2"]))

    def test_cutpoint_between_s_pk(self):
        db = "((([[[&)))]]]"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "m2"]))

    def test_cutpoint_between_s_pk2(self):
        db = "((([[[)))&]]]"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "m1"]))
    def test_cutpoint_between_s_pk3(self):
        db = "(((..[[[&)))]]]"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        log.error(cg.defines)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "m2"]))

    def test_cutpoint_between_s_pk4(self):
        db = "((([[[..)))&]]]"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        print(cg.defines)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "m0", "s1", "m1"]))


    def test_cutpoint_between_s_i(self):
        db = "(((&...(((...))))))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "f0", "s1", "h0", "m0"]))

    def test_cutpoint_between_s_m(self):
        db = "(((.(((...)))&..(((...))).)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "s1", "s2", "m0", "m2", "f0", "h0", "h1"]))

    def test_cutpoint_between_s_h(self):
        db = "(((&..)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "f0"]))

    def test_cutpoint_between_h_s(self):
        db = "(((..&)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0", "t0"]))

    def test_cutpoint_instead_of_h(self):
        db="(((&)))"
        cg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(set(cg.defines.keys()), set(["s0"]))

    def test_seq_id_to_pos(self):
        fasta = """>1L2X
                   GCGCG&CGUGC
                   (((((&)))))"""
        bg = fgb.from_fasta_text(fasta)
        bg.seq_ids = [fgb.RESID("A", (" ", 1," ")),fgb.RESID("A", (" ", 2," ")),
                      fgb.RESID("A", (" ", 3," ")),fgb.RESID("A", (" ", 4," ")),
                      fgb.RESID("A", (" ", 5," ")),
                      fgb.RESID("B", (" ", 1," ")),fgb.RESID("B", (" ", 2," ")),
                      fgb.RESID("B", (" ", 3," ")),fgb.RESID("B", (" ", 4," ")),
                      fgb.RESID("B", (" ", 5," "))]
        print(bg.seq_ids)
        self.assertEqual(bg.seq_id_to_pos(fgb.RESID("A", (" ", 1," "))), 1)
        self.assertEqual(bg.seq_id_to_pos(fgb.RESID("B", (" ", 1," "))), 6)
        self.assertEqual(bg.seq[bg.seq_id_to_pos(fgb.RESID("B", (" ", 3," ")))], "U")

    def test_dissolve_length_one_stem_cofold(self):
        db = "(((.(&...).)))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db, dissolve_length_one_stems=True)
        self.assertEqual(len(bg.defines), 4)
        self.assertIn("s0", bg.defines)
        self.assertIn("t0", bg.defines)
        self.assertIn("f0", bg.defines)
        self.assertIn("t1", bg.defines)
    def test_dissolve_length_one_stem_cofold_2(self):
        db = "(((.(...).(.&..).)))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db, dissolve_length_one_stems=True)
        self.assertEqual(len(bg.defines), 4)
        self.assertIn("s0", bg.defines)
        self.assertIn("t0", bg.defines)
        self.assertIn("f0", bg.defines)
        self.assertIn("t1", bg.defines)
        self.assertEqual(bg.edges["s0"], {"f0", "t0", "t1"})
        self.assertEqual(bg.edges["f0"], {"s0"})

    def test_to_dotbracket_string_with_cutpoints(self):
        db = "(((.(.&..).(.&..).))&)..."
        bg = fgb.from_fasta_text(db)
        self.assertEqual(bg.to_dotbracket_string(), db)
        db = "(((.&..)))"
        bg = fgb.from_fasta_text(db)
        self.assertEqual(bg.to_dotbracket_string(), db)
        db = "((([[[..)))&]]]"
        bg = fgb.from_fasta_text(db)
        self.assertEqual(bg.to_dotbracket_string(), db)

class BulgeGraphZeroLengthTest(GraphVerification):
    def test__zero_length_element_adj_position_single_ml(self):
        db="(((...)))(((...)))"
           #123456789012345678
        bg=fgb.from_id_seq_struct("test1", "A"*len(db), db)
        self.assertEqual(bg._zero_length_element_adj_position("m0"), [9,10])
    def test__zero_length_element_adj_position_two_ml(self):
        db="((([[[)))..]]]"
           #12345678901234
        bg=fgb.from_id_seq_struct("test1", "A"*len(db), db)
        zl_elems=[]
        c_l=0
        for elem, d in bg.defines.items():
            if elem[0]=="m":
                if d:
                    with self.assertRaises(ValueError):
                        bg._zero_length_element_adj_position(elem)
                    c_l+=1
                else:
                    zl_elems.append(elem)
        zl_elems.sort()
        self.assertEqual(c_l,1)
        self.assertEqual(len(zl_elems),2)
        self.assertEqual(bg._zero_length_element_adj_position(zl_elems[0]), [3,4])
        self.assertEqual(bg._zero_length_element_adj_position(zl_elems[1]), [6,7])
    def test__zero_length_element_adj_position_single_ml2(self):
        db="(((..[[[)))..]]]"
           #12345678901234
        bg=fgb.from_id_seq_struct("test1", "A"*len(db), db)
        zl_elems=[]
        c_l=0
        for elem, d in bg.defines.items():
            if elem[0]=="m":
                if d:
                    with self.assertRaises(ValueError):
                        bg._zero_length_element_adj_position(elem)
                    c_l+=1
                else:
                    zl_elems.append(elem)
        zl_elems.sort()
        self.assertEqual(c_l,2)
        self.assertEqual(len(zl_elems),1)
        self.assertEqual(bg._zero_length_element_adj_position(zl_elems[0]), [8,9])
    def test__zero_length_element_adj_position_three_ml(self):
        db="((([[[)))]]]"
           #12345678901234
        bg=fgb.from_id_seq_struct("test1", "A"*len(db), db)
        zl_elems=[]
        for elem, d in bg.defines.items():
            if elem[0]=="m":
                self.assertEqual(d, [])
        self.assertEqual(bg._zero_length_element_adj_position("m0"), [3,4])
        self.assertEqual(bg._zero_length_element_adj_position("m1"), [6,7])
        self.assertEqual(bg._zero_length_element_adj_position("m2"), [9,10])

    def test_breakpoint_at_zero_length_element_graph_construction(self):
        # Test needed because of a subtile bug with cofold structures and
        # dissolve_length_one_stems and 0-length elements.
        # It appeared during loading of PDB 1U6B.pdb when keeping pseudoknots
        # and dissolving length 1 stems.
        # and cuased a GraphConstructionError to be raised.
        db = "((([[[.(..)..)))&]]]"
        bg = fgb.from_fasta_text(db, dissolve_length_one_stems=True)
        self.assertEqual(len(bg.defines), 4)
        self.assertEqual(bg.to_dotbracket_string(), "((([[[.......)))&]]]")
class BulgeGraphTest(GraphVerification):
    """
    Simple tests for the BulgeGraph data structure.

    For now the main objective is to make sure that a graph is created
    and nothing crashes in the process. In the future, test cases for
    bugs should be added here.
    """

    def setUp(self):
        self.dotbracket = '....((((((....((.......((((.((((.(((...(((((..........)))))...\
((.......))....)))......))))))))......))...)).))))......(((....((((((((...))))))))...)))........'

        self.bg_string = """
name temp
length 71
seq CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
define f0 0 1
define h1 47 55
define s3 42 47 55 60
define s2 13 19 27 33
define h0 19 27
define s0 1 9 63 71
define t0 71 72
define m1 9 13
define m2 33 42
define m0 60 63
connect s3 h1 m0 m2
connect s2 h0 m1 m2
connect s0 f0 m1 m0 t0
"""
        self.bpseq = dict()
        self.bpseq['1y26'] = """1 G 26
2 A 25
3 G 24
4 C 23
5 U 22
6 G 21
7 C 0
8 A 0
9 G 0
10 C 19
11 A 18
12 C 17
13 G 0
14 A 0
15 A 0
16 A 0
17 G 12
18 U 11
19 G 10
20 A 0
21 C 6
22 G 5
23 G 4
24 C 3
25 U 2
26 C 1
"""

        self.bpseq['pseudoknot'] = """1 A 9
2 A 8
3 A 0
4 A 14
5 A 13
6 A 0
7 A 0
8 A 2
9 A 1
10 A 0
11 A 0
12 A 0
13 A 5
14 A 4
"""

    def test_to_bg_string(self):
        self.fasta = """>1y26
CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
(((((((((...((((((.........))))))........((((((.......))))))..)))))))))
"""
        bg = fgb.from_fasta_text(self.fasta, dissolve_length_one_stems=True)
        stri = bg.to_bg_string()
        bg2 = fgb.BulgeGraph()
        bg2.from_bg_string(stri)
        stri2 = bg2.to_bg_string()
        self.assertEqual(stri, stri2)
        self.assertTrue(bg.defines, bg2.defines)
        self.assertTrue(bg.edges, bg2.edges)

    def test_bg_string_infos(self):
        self.fasta = """>1y26
CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
(((((((((...((((((.........))))))........((((((.......))))))..)))))))))
"""
        bg = fgb.from_fasta_text(self.fasta, dissolve_length_one_stems=True)
        bg.add_info("test", "This is a test info")
        stri1 = bg.to_bg_string()
        bg2 = fgb.BulgeGraph()
        bg2.from_bg_string(stri1)
        self.assertEqual(bg2.infos["test"], ["This is a test info"])

    def test_from_fasta(self):

        with open('test/forgi/threedee/data/3V2F.fa', 'r') as f:
            text = f.read()
            bg = fgb.from_fasta_text(text, dissolve_length_one_stems=False)

        for s in bg.stem_iterator():
            bg.stem_length(s)

    def test_from_fasta2(self):
        fasta_str = (">a\n"
                     "AAAA\n"
                     "([)]\n")
        bg = fgb.from_fasta_text(fasta_str)

        self.assertEqual(bg.defines['s0'], [1, 1, 3, 3])
        self.assertEqual(bg.defines['s1'], [2, 2, 4, 4])

        fasta_str = (">a\n"
                     "AAAAAAA\n"
                     "(.[.).]\n")
        bg = fgb.from_fasta_text(fasta_str)

        self.assertEqual(bg.defines['s0'], [1, 1, 5, 5])
        self.assertEqual(bg.defines['s1'], [3, 3, 7, 7])

    def test_from_fasta_text_with_whitespace(self):
        a = (">a \n"
             "ACGCCA \n"
             "((..)) \n")

        x = fgb.from_fasta_text(a)
        self.assertEqual(x.seq, 'ACGCCA')
        self.assertEqual(x.name, 'a')
        self.assertEqual(x.to_dotbracket_string(), '((..))')


    def test_from_fasta1(self):
        a = (">a\n"
             "ACGCCA\n"
             "((..))\n")

        x = fgb.from_fasta_text(a)
        self.assertEqual(x.seq, 'ACGCCA')

        a = (
"""
>a
ACGCCA
((..))
>b
CCCCCC
((()))
>c
AAAAAA
(....)
""")
        bgs = fgb.from_fasta_text(a)
        self.assertEqual(len(bgs), 3)
        self.assertEqual(bgs[0].seq, 'ACGCCA')
        self.assertEqual(bgs[1].seq, 'CCCCCC')
        self.assertEqual(bgs[2].seq, 'AAAAAA')
        self.assertEqual(bgs[0].name, 'a')
        self.assertEqual(bgs[1].name, 'b')
        self.assertEqual(bgs[2].name, 'c')
        self.assertEqual(bgs[0].dotbracket_str, '((..))')
        self.assertEqual(bgs[1].dotbracket_str, '((()))')
        self.assertEqual(bgs[2].dotbracket_str, '(....)')
        a = """
GGGGGG
......
"""
        bg = fgb.from_fasta_text(a)
        self.assertEqual(bg.seq, 'GGGGGG')

    def test_from_fasta_pseudoknot(self):
        a = """
>3NKB_B
GGUCCGCAGCCUCCUCGCGGCGCAAGCUGGGCAACAUUCCGAAAGGUAAUGGCGAAUGCGGACC
(((((((((((.[[....))).......]]....((((((....)).)))).....))))))))
"""
        bg = fgb.from_fasta_text(a)
        self.assertEqual(len(list(bg.stem_iterator())), 5)


    def test_from_fasta4(self):
        struct = ".(((((((((((((((((.......))))))))(((((...)))))..((((...(((((....)))))...)))).)))))))))..((((.((((((((((...((((..((((((...)))))).(((((.((((((((.....)))))))).((((.....)))))))))....((((.(((((((.(((((((.....)))))........))))).)))).))))(((((((..(((((...((((....(((((((....((.((((((((((((.....((((.......)))).)))).(((((((((((((((((......(((((....)))))))))))))).)))..)))))..(((((((...))))))).)))))))).))....)))))))...))))...)))))..))))))).(((....)))..)))))))))).(((.((((....(((.((((((((((((......)))))).)))))).))).))))...))).)))).)))).......((((((((((..((((((((((.((((((((((((((((((((((((.(((((((((.....))))))))).))))))))))......(((....))).)))))).((((((..(((((.(((......)))(((...((((((.(((..((((((((.....(((....((((((..(((((.(((......)))))))).....))))))........)))...))).))))))))))))))((((((((((((((((...((((..((.(((..((((((.((.(((...))).)).)))))).((.(((((((....))).)))).)).((((((((((..................(((((((((((((((((((...)).)))))).)))))))))))))))(((((...))))).))))))...))).))..)))).....))))).))))))).)))).))).((((((((....))))))))..))))).))))))..))))))..)).)))))))..))).....))))))))))(((((((.(((..((((((((.....))))))))((((((((((....(((....)))......((((.((((((..(((((((...((((((((((((.((.....)))))))))....((((...))))(((((..(((.....))).)))))))))))))))))))).))).)))).)))))))))).........((((((((...(((...((((((((..((((((......((((..(((.(((.(((((.(((.(((((((((((.((((((.((((.(((....))))))).....(((.((...)).)))))))))...((((((((((((.(((((((..((((((((..((..((((((((.......))).)))))..)).(((..(((((((((((....)))))).)).)))..))).))))).))).))))))).....(((..(((....)))..)))........)))).)))))))).)))((((....))))...)))))))).)))))))).))).)))....))))))))))))).)))))....))).))))))))))).)))))))....((((.((((........)))).))))((((((..(((((.((.(((((((....(((((((((((((((.(((((((((...(((((..(((.....(((((.............)))))....))))))))))).(((((((..((.....)))))))))..((((...))))..)))))).).))))))((((.....)))).))))))))....)))))))(((((((((......))))))))))))))))..))))))....((((((..(((.(((((.(((((.((((.((((((((((....(((((.(((((........)))((((((((.(((...((............)).))).))))))))..((((((...)))))).(((((((((((...(((.(((((...(((..(((((((((((((....)))))......((......))..(((((((((((.((((((((((((.((((((.(((.((..(((((...))))))).)))...............((((.....)))))))))).))))))))))))...)))..))))))))..))))))))..(((((((.....(((.(....))))((((....)))).((((..(((((((((((.((....)).)))))))))))..))))))).))))))).)))))..)))..)))))))))))........((((((...(((((...))))))))))))))))))....))))))))))(((((((..(((...)))....)).))).)).))))..)))))))))).))).))))))...............................................((((((.(((((...((.(((((((.(((((...(((((((..(((.(((((((((..(((...)))))))))))).))).....)))))))...))))).(((((...))))).))))).)))).)))))))))))((((((((...((.(((((((..((.((((((...(((((((.((((.....))))))).......))))..))))))))..))..))))).))...))))))))...................................................................................(((((((..((((...)))))))))))(((((((.((((((((((.(((((.((((((..(((........((((((((....(((.(((((((((.............................................((((((.(((......(((..(((.((.....)).)))........(((..(((((((....))))))).)))..(((((....))))).(((((....))))).))))))..))))))........)))))).))).))))))))))))))...))).))))))))((((((....))))))((((....(((((..(((((....)))))...))))))))))))))..))))).)))))))(((.((((((((.(((..(((.....))).))))))))))).)))....((((((((((((...(((.((.(((((.(((......))).))))).)))))......)))))))))))).((....))................................................................................................((((((((((....(((....)))..))))))))))........(((((((((.(((((((((..(((..(((((.(((.......)))))).))..)))))))))))))))))).)))..(((((......)))))((((...(((.(((((((((.(((((........))))))).))))))).)))...((.....((((.((((.((((((((..(((.((((((........))))))....)))))))))))......)))))))))).))))...........((((((((.......)))))))).(((((((....((((((((((.....))))))(((((..((((.....((((.(((.((((((..((((((((......)))))..((((((((((((((((.(((...))))))......((((..((((((..........)))))).((((((...))))))...((((..((((((((((((((((((..(((((.(((..((((......))))...))).)))))............))))..))))))))))).))))))).((....))))))...(((((((((((........)).)))...)))))).))).))))))))))..)))))))))..))).))))....(((........))).(((((...)))))......))))...))))).....(((((((...((.((..(((((((((.(((..(((....)))(((((.((((............)))).)))))....))))))))).))).)).)).))))).))......)))))))))))..((((((((.((((....)))).)))))))).(((.(((((.....))))))))(((((((((((((((((((..((((((((....(((....)))...))))))))...)))).)).))...)))))((((...))))((((((...))))))))))))......................................((((.(.((((((.(((.....................................................((((((((((....((((((..((((.....(((.(((((((((((..((.((.(((..(((.....)))..)))))...)).)))))))))))))))))).)))))).((....))(((((....((((((....(((((..(((........))).)))))(((...)))..))))))))))).((((((((((((((.((((((((((.(((((((.(((((..........)))....)).))))).........)).))...))))))..)).))).)))).)))))))........)))))))))))))))))))).)))).......(((....))).(((...)))....................................((((((.((.........................(((((((....))))))).......((((((.......)))))).......................(((((((((..(((((........(((((((((((((((.(((((.......))))).((((.((.....))))))..))).))))))))))))........)))))..))....)))))))......................................................(((((((((................(((((((...((((.((.(((((.........)))))......((((...(((......)))..)))).(((((((.....))))....(((((....)))))....))).)))))).))))))).......((((....((((...))))))))..)))))))))(((((.((.....))))))).)))))))).((((((.((((((((.((((((((.....((....)).(((.((((((...((((((.....(((((....))))))))))).)))))).)))......)))))))).((((((((((((.((.(((((((((.((((.((((.....)))).)))))))).)))))))..)))))))....(((((((((((.((((.((((((...(((((((..(((..((((...(((...))))))))))....(((((.......))))).......)))))))))))))))))..(((((((....))))))).((((((((.(((....(((((..((((....(((.(((..((((.((((.((((((((........)))))))).....)))).))))))))))........)))))))))......))).))))))))..((((((((..((((((((.(((((............))).)).))))))))..))...))))))......)))))......))))))..))))).((....)).)))))))).))))))((((((........))))))....(((((..(((.((((.((((..(((((..(((((...(((.(((.....((((((..((((((.((((.(((((.((.....(((..((.....))))))).))))).)))).)))))).((((.((((.....))))..)))).(((((((((...)))))))))...))))))))).)))...)))))..)))))...)))).)))))))..))))).((((((((.((......)))))))))).((((((((.....))))))))...((((((((((..((((((((((((((..((((......))))))))))))((((.((..(((.((((((.(((((.....)))))))))))..)))..))))))))))))))))))))))......((((((.....)))))).....(((((..(((((.((((((((...((((((((((((((.((((((((..(((((((((....((((((((.(((((((..((((.....((((.((((...((((...((((((..(((......)))..)))))).)))))))).))))))))..)))))))...(((((.((...((((......))))..)))))))(((((((((((((.(((.(((((((((((...)))))))....))))))).)).....((((....((((...)))))))))))))))..((((((((((((...))))))))...)))).(((((((((......)))))....)))).....)))).))))))))...)))))))))...))))))))))......))).)).(((((.((((...(((((((((((.((((((((.....))))))))..))))..)))))))(((((...)))))((((..(((..((((((((.((....))..))))))))..))).)))).))))))))))))))))....))))))))....))))).)))))((((((.(((((((((...))))))))).))))))"
        seq = """ttaaaactggatccaggttgttcccacctggatttcccacagggagtggtactctgttat
tacggtaactttgtacgccagttttatctcccttcccccatgtaacttagaagtttttca
caaagaccaatagccggtaatcagccagattactgaaggtcaagcacttctgtttccccg
gtcaatgttgatatgctccaacagggcaaaaacaactgcgatcgttaaccgcaaagcgcc
tacgcaaagcttagtagcatctttgaaatcgtttggctggtcgatccgccatttcccctg
gtagacctggcagatgaggctagaaataccccactggcgacagtgttctagcctgcgtgg
ctgcctgcacaccctatgggtgtgaagccaaacaatggacaaggtgtgaagagccccgtg
tgctcgctttgagtcctccggcccctgaatgtggctaaccttaaccctgcagctagagca
cgtaacccaatgtgtatctagtcgtaatgagcaattgcgggatgggaccaactactttgg
gtgtccgtgtttcactttttcctttatatttgcttatggtgacaatatatacaatatata
tattggcaccatgggtgcacaggtttcaagacaaaatgttggaactcactccacgcaaaa
ctctgtatcaaatgggtctagtttaaattattttaacatcaattatttcaaagatgctgc
ttcaaatggtgcatcaaaactggaattcacacaagatcctagtaaatttactgacccagt
taaggatgttttggaaaagggaataccaacactacagtcccccacagtggaggcttgtgg
atactctgataggattatacagattaccagaggagattcaaccataacctcacaagatgt
ggctaatgctatcgttgcgtatggtgtttggccacattatctatcctccaaggatgcctc
tgcaattgataaaccctctcaaccagatacatcttctaatagattttatactctaaggag
tgtgacctggagcagttcctcaaagggttggtggtggaaactacctgatgcactcaagga
catgggtatttttggtgaaaacatgttttatcattacctgggtaggagtggatacacaat
acatgtgcagtgtaatgctagtaaatttcaccagggtacactaattgttgctctgatacc
tgagcatcagattgcaagtgccttacatggcaatgtgaatgttggttacaactacacaca
cccaggtgaaacaggcagggaagttaaagctgagacgagattgaatcctgatctacaacc
tactgaagagtattggctaaactttgatgggacactccttggaaatattaccatattccc
tcatcaatttatcaacttgaggagtaataattctgccacaataattgccccttatgtcaa
tgcagttcctatggattcaatgcggagccacaataattggagtttggtaataataccaat
atgtccccttgagacatcaagtgcaattaacacaatacctattacaatatctataagccc
catgtgtgcagagttttccggcgcgcgtgccaagcgtcaaggattaccagttttcatcac
accaggttcaggacagtttttgacaacagatgatttccaatccccatgtgcacttccctg
gtatcacccaactaaggaaatttctattccaggtgaggttaaaaatttggttgaaatttg
tcaagtagacagcctagtaccaataaataacactgacacctacatcaatagtgaaaatat
gtattctgttgtattgcaatcatcaattaatgcaccagataagatcttctctattcgaac
agatgttgcttcccaacctttagctactactttgattggtgagatatctagctatttcac
ccactggacagggagtctccgtttcagcttcatgttttgtggtactgccaacactactgt
taagcttttgttggcatacacaccacctggtatcgcagaacccaccacaagaaaggatgc
aatgctaggcactcatgttatatgggatgtggggttgcagtctacaatatcaatggtagt
gccatggattagcgctagtcattatagaaacacatcaccaggtagatctacatctgggta
cataacatgctggtatcagactagattagtcattccacctcagaccccaccaacagctag
attgttatgttttgtatctgggtgcaaagacttttgcttgcgcatggcacgagatactaa
cctacacctgcaaagtggtgcaatagcacagaaccctgttgagaattatatagatgaagt
tcttaatgaagttttagttgtcccaaatattaatagtagtaaccccacaacatcaaattc
tgccccagcattagatgctgcagaaacagggcacactagtagtgttcaaccagaggatgt
cattgaaactaggtatgtgcagacatcacaaacaagagatgaaatgagtttagagagttt
tcttggcagatcaggatgcatacatgaatctaaattagaggttacacttgcaaattataa
caaggagaattttacagtgtgggctattaatctacaagaaatggctcaaattagaaggaa
atttgaattgttcacctatactaggtttgattctgaaataaccctagttccatgcatttc
cgcccttagtcaggacattggacacatcacaatgcaatacatgtatgttccaccaggtgc
accggtgcccaatagtagggacgattatgcatggcagtctggcactaatgcctctgtttt
ctggcaacatggacaggcttatccaagattttccttacctttcctaagtgtggcatctgc
ttattacatgttttatgatgggtatgatgaacaagatcaaaactatggtacagcaaacac
aaataacatggggtcactatgctctaggatagtaacagagaaacacattcataaagtaca
tataatgacaagaatctatcacaaggctaaacatgtcaaggcatggtgtccacgcccacc
cagagcgcttgagtatactcgtgctcatcgcactaattttaaaattgaggataggagtat
tcagacagcaattgtgaccagaccaattatcactacagctggccccagtgacatgtatgt
tcatgtaggtaaccttatttatagaaatcttcatcttttcaactctgagatgcatgaatc
tattttggtatcttattcatcagatttaatcatttaccgaacaaacactgtaggtgatga
ttacattccctcttgtgattgtacccaagctacttattattgcaaacataaaaatagata
cttcccaattacagttacaagccatgactggtatgaaatacaggaaagtgagtactatcc
caaacacatacagtacaatttgttgattggtgagggcccttgtgaaccaggtgactgtgg
tggaaagttgctatgcaaacatggtgtcataggtatagtaacagctggtggtgataatca
tgtggcttttattgaccttagacacttccattgtgctgaagaacaaggggttacagatta
tatacatatgctaggagaagcatttggaaatggatttgtggatagtgtaaaagaacatat
acatgccataaacccagtaggaaatatcagcaagaaaattattaaatggatgttgagaat
aatatcagcaatggtcataataattagaaactcttctgacccccaaactatattagcaac
actcacactgattgggtgttctggatcaccctggagatttttaaaggaaaaattctgtaa
atggacacagcttaattatatacacaaagaatcagattcatggttaaagaaatttactga
agcatgcaatgcagctagagggcttgaatggatagggaataagatatctaaatttattga
atggatgaagtcgatgctcccgcaagctcaattgaaggttaagtacttaaacgagcttaa
aaaactcaacctatacgaaaagcaagttgagagcttgcgggtggctgacatgaaaacaca
agaaaaaattaaaatggaaatagacactttacatgatttgtcacgtaaatttctaccttt
gtatgcaagtgaggcaaaaaggataaaaaccctatacattaaatgtgataatatcatcaa
gcagaagaaaagatgtgaaccagtagctatagttattcatggaccacctggtgctggcaa
atctataacaacaaatttcctggccaaaatgataactaatgatagtgacatatactctct
acctcctgatccaaaatattttgatggttatgaccaacagagtgtagtaataatggatga
cattatgcagaatccagccggggatgacatgacactgttctgccaaatggtttctagtgt
tacatttataccaccaatggctgatctaccagataaaggcaaggcttttgattctaggtt
tgtattatgcagcacaaatcattcccttctaacacccccgacaataacttcactacctgc
aatgaatagaagatttttcctagatttagatataatagtacatgataacttcaaagatcc
acagggcaaacttaatgtggcagcagcgtttcgaccatgtgatgtagataatagaatagg
aaatgcacgttgttgtccatttgtgtgtggaaaagcagtttctttcaaagatcgtaactc
ttgcaacaaatacagccttgcgcaggtgtacaacataatgattgaagaagacagacggag
aagacaagtggttgatgtcatgacagctatattccaagggccaattgatatgaaaaaccc
accaccacctgctattactgacttgctccagtctgttagaacccctgaagttattaagta
ttgtgagggtaatagatggataattccagcagaatgcaagatagaaaaggagttgaactt
ggctaacacaatcataacaatcattgcaaatgttattggtatggcgagaataatatatgt
tatttacaaacttttttgcacattacagggaccatattcaggagaaccaaagcccaagac
taaaatcccagaaaggcgtgtagtaacacagggaccagaggaggaatttgggatgtcttt
aattaaacataactcatgtgttattacaacagaaaatgggaaattcacaggtcttggagt
atacgacagatttgtggtcgtaccaacacatgcagatcctggaaaggaaattcaggttga
tggtataactacaaaagtcattgactcatatgacctatacaacaagaatgggataaagct
agaaataacagtacttaaattagatagaaatgaaaaatttagagatatcaggagatatat
acctaacaatgaagatgattaccccaattgcaacttagcactgctagcaaaccagcctga
accaactataatcaatgttggagatgttgtatcctatggcaatatactgctcagtggcaa
ccaaacggctagaatgcttaaatacagttacccaactaaatctggttactgtggaggtgt
cttatacaaaattgggcaagtgcttggaatacatgttgggggcaatggtagggatggttt
ctcagctatgttactcagatcctatttcactgatgttcagggccaaataacgttatcaaa
gaagaccagtgaatgtaacctacccagtatacacaccccatgcaaaaccaaattgcagcc
tagtgttttctatgatgtattccctggttcaaaagaaccagctgtgttgtctgaaaaaga
tgcccggttacaagttgatttcaatgaagcactattttctaaatacaaagggaatacaga
ttgctccattaatgaccacataagaattgcatcatcacattatgcagcacaactcattac
cttagatattgacccaaaacctattacacttgaggacagtgtctttggcactgatggatt
agaggctcttgatttgaacactagcgcaggatttccatatattgcaatgggagttaaaaa
gagagatttaataaacaacaagaccaaggatataagcaaacttaaagaagcaattgacaa
atacggagttgacttacctatggtcaccttcttgaaagatgaactcagaaagcatgaaaa
ggtaattaaaggtaaaactagagttattgaagctagtagtgtgaatgataccctattatt
tagaacaacttttggcaacctcttttcaaagttccacttgaatcctggaattgttactgg
atcagcagttggatgtgatccagaggtgttttggtcaaaaataccagcaatgttggatga
taaatgtattatggcttttgattatacaaattatgatggtagtatacaccctatttggtt
tgaagctcttaaacaggtactggtagatctatcatttaatccaacattaatagatagact
atgcaagtctaaacacatcttcaaaaatacatactatgaagtggagggaggtgtaccatc
tgggtgttcaggtactagtatttttaacactatgatcaataatattatcataaggacctt
agtgttagatgcatacaagaatatagatctagataagcttaagataattgcctatggtga
tgatgtcatattctcatacatacatgaactggacatggaggctatagcaatagagggtgt
taaatatggtttgactataactcctgctgataaatctaacacatttgtaaaattagacta
tagcaatgttacttttttaaaaagagggtttaagcaagatgagaagtataactttctaat
acatccaactttccctgaagatgaaatatttgaatccatcagatggacaaagaaaccatc
acaaatgcatgaacatgtgttgtctctgtgtcacttaatgtggcacaatggacgtgacgc
atacaaaaaatttgtggagaagatacgcagtgtaagcgctggtcgtgcactgtacatccc
tccgtatgatttgcttttgcatgagtggtatgaaaaattttaaagatatagaaatagtaa
actgatagtttattagttttat
"""
        fasta = ">hrv\n{}\n{}".format(seq, struct)
        bg1 = fgb.from_fasta_text(fasta)

        seq = seq.replace('\n', '')
        fasta = ">hrv\n{}\n{}".format(seq, struct)
        bg2 = fgb.from_fasta_text(fasta)

        self.assertEqual(bg1.defines['s56'], [713,717,731,735])
        self.assertEqual(len(bg1.defines), 1167)
        self.assertEqual(bg1.seq, seq.replace("t", "u"))

        self.assertEqual(bg1.defines, bg2.defines)
        self.assertEqual(bg1.seq, bg2.seq)

        #print >>sys.stderr, "bg.defines:", bg.defines


    def test_from_fasta_double_zero_length_ml(self):
        a = """
>test
AAAGGGUUUCCC
((([[[)))]]]
"""
        bg = fgb.from_fasta_text(a)
        self.assertEqual(len(bg.defines), 5)
        self.assertEqual(bg.defines["m0"], [])
        self.assertEqual(bg.defines["m1"], [])
        self.assertEqual(bg.defines["m2"], [])
        self.assertIn("m0", bg.edges["s0"])
        self.assertIn("m1", bg.edges["s0"])
        self.assertIn("m2", bg.edges["s0"])
        self.assertIn("m0", bg.edges["s1"])
        self.assertIn("m1", bg.edges["s1"])
        self.assertIn("m2", bg.edges["s1"])

    def test_from_bpseq_file(self):
        with open('test/forgi/data/1gid.bpseq', 'r') as f:
            lines = f.readlines()

        bpseq_str = "".join(lines)
        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(bpseq_str, dissolve_length_one_stems=True)

        for d in bg.defines:
            self.assertFalse(d[0] == 'x')

        with open('test/forgi/data/1ymo.bpseq', 'r') as f:
            lines = f.readlines()

        bpseq_str = "".join(lines)
        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(bpseq_str, dissolve_length_one_stems=True)

        node = bg.get_node_from_residue_num(25)
        self.assertFalse(node[0] == 'h')

    def test_from_bpseq_error(self):
        bpstr1 = """
1 G 2
2 G 0
"""
        bpstr2 = """
1 G 2
2 G 1
3 A 2
"""
        bg = fgb.BulgeGraph()
        with self.assertRaises(GraphConstructionError):
            bg.from_bpseq_str(bpstr1)
        with self.assertRaises(GraphConstructionError):
            bg.from_bpseq_str(bpstr2)

    def test_from_bpseq(self):
        bg = fgb.BulgeGraph()

        bpstr = """1 G 8
2 G 7
3 C 6
4 A 5
5 U 4
6 G 3
7 C 2
8 C 1
"""
        bg.from_bpseq_str(bpstr, dissolve_length_one_stems=False)
        bg.from_bpseq_str(bpstr, dissolve_length_one_stems=True)

        #seq = 'AAAAAAAAAAAAAAAAAAAAA'
        #db = '.((((..)).))..((..)).'
        #n = '12345678901234567890.'

        bg = fgb.BulgeGraph()
        bpstr = """1 A 0
2 A 12
3 A 11
4 A 9
5 A 8
6 A 0
7 A 0
8 A 5
9 A 4
10 A 0
11 A 3
12 A 2
13 A 0
14 A 0
15 A 20
16 A 19
17 A 0
18 A 0
19 A 16
20 A 15
21 A 0
"""

        bg.from_bpseq_str(bpstr)

        self.assertEqual(bg.defines['i0'], [10, 10])
        self.assertEqual(bg.defines['h0'], [6, 7])
        self.assertEqual(bg.defines['h1'], [17, 18])
        self.assertEqual(bg.defines['s0'], [2, 3, 11, 12])
        self.assertEqual(bg.defines['s1'], [4, 5, 8, 9])
        self.assertEqual(bg.defines['s2'], [15, 16, 19, 20])
        self.assertEqual(bg.defines['t0'], [21, 21])

        bg.get_node_from_residue_num(21)

        bpstr = """1 G 26
2 A 25
3 G 24
4 C 23
5 U 22
6 G 21
7 C 0
8 A 0
9 G 0
10 C 19
11 A 18
12 C 17
13 G 0
14 A 0
15 A 0
16 A 0
17 G 12
18 U 11
19 G 10
20 A 0
21 C 6
22 G 5
23 G 4
24 C 3
25 U 2
26 C 1
"""

        bg.from_bpseq_str(bpstr)
        bg.get_node_from_residue_num(1)

        bpstr = """1 G 8
2 G 7
3 C 6
4 A 5
5 U 4
6 G 3
7 C 2
8 C 1
"""
        bg.from_bpseq_str(bpstr)
        bg.get_node_from_residue_num(1)
        bg.get_node_from_residue_num(2)

        #db = '(.(.(.).).)'
        #nm = '12345678901'

        bpstr = """1 A 11
2 A 0
3 A 9
4 A 0
5 A 7
6 A 0
7 A 5
8 A 0
9 A 3
10 A 0
11 A 1
"""
        bg.from_bpseq_str(bpstr)

        #db = '[.(.].)'
        #nm = '1234567'

        bpstr = """1 A 5
2 A 0
3 A 7
4 A 0
5 A 1
6 A 0
7 A 3
"""

        bg.from_bpseq_str(bpstr)

        db = '[[.((..]]...))'
        nm = '12345678901234'

        bpstr="""1 A 9
2 A 8
3 A 0
4 A 14
5 A 13
6 A 0
7 A 0
8 A 2
9 A 1
10 A 0
11 A 0
12 A 0
13 A 5
14 A 4
"""
        bg.from_bpseq_str(bpstr)

        db='[[.((..]]...)).((..((..)).))'
        nm='1234567890123456789012345678'
        bpstr="""1 A 9
2 A 8
3 A 0
4 A 14
5 A 13
6 A 0
7 A 0
8 A 2
9 A 1
10 A 0
11 A 0
12 A 0
13 A 5
14 A 4
15 A 0
16 A 28
17 A 27
18 A 0
19 A 0
20 A 25
21 A 24
22 A 0
23 A 0
24 A 21
25 A 20
26 A 0
27 A 17
28 A 16
"""
        bg.from_bpseq_str(bpstr)

        db='[[.((..]]..((..)).))'
        nm='12345678901234567890'

        bpstr="""1 A 9
2 A 8
3 A 0
4 A 20
5 A 19
6 A 0
7 A 0
8 A 2
9 A 1
10 A 0
11 A 0
12 A 17
13 A 16
14 A 0
15 A 0
16 A 13
17 A 12
18 A 0
19 A 5
20 A 4
"""
        bg.from_bpseq_str(bpstr)



    def test_from_dotplot(self):
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(self.dotbracket)

        self.assertEquals(bg.seq_length, len(self.dotbracket))

        bg = fgb.BulgeGraph()
        bg.from_dotbracket('....')


    def test_define_residue_num_iterator1(self):
        bg = fgb.BulgeGraph(dotbracket_str='(.(...).)')
        ress = list(bg.define_residue_num_iterator('i0', adjacent=True))
        self.assertEqual(ress, [1,2,3,7,8,9])

        bg = fgb.BulgeGraph(dotbracket_str='((.((...))))')
        ress = list(bg.define_residue_num_iterator('i0', adjacent=True))
        self.assertEqual(ress, [2,3,4,10,11])


    def test_define_residue_num_iterator(self):
        bg = fgb.BulgeGraph(dotbracket_str='((.).)')
        self.assertEqual(list(bg.define_residue_num_iterator('s1', adjacent=True)), [1,2,3,4,5])

        bg = fgb.BulgeGraph(dotbracket_str='((..((..))((..))))')
        drni = bg.define_residue_num_iterator('m2', adjacent=True)
        # the second multiloop should have at least two adjacent nucleotides
        self.assertEqual(len(list(drni)), 2)
        drni = bg.define_residue_num_iterator('m1', adjacent=True)
        # the second multiloop should have at least two adjacent nucleotides
        self.assertEqual(len(list(drni)), 2)

        bg.define_residue_num_iterator('m1', adjacent=True)

        bg = fgb.from_fasta_text('..((..((...))..))..((..))..')

        self.assertEqual(list(bg.define_residue_num_iterator('f0')),
                         [1,2])
        self.assertEqual(list(bg.define_residue_num_iterator('t0')),
                         [26, 27])
        self.assertEqual(list(bg.define_residue_num_iterator('s1')),
                         [7, 8, 12, 13])
        self.assertEqual(list(bg.define_residue_num_iterator('i0')),
                         [5,6,14,15])

        fa=""">blah
AAAAAAAAAA
((((.)).))
"""
        bg = fgb.from_fasta_text(fa, dissolve_length_one_stems=True)
        self.assertEqual(list(bg.define_residue_num_iterator('i0', adjacent=True)),
                         [2,3,7,8,9])

        self.assertEqual(list(x.resid for x in bg.define_residue_num_iterator('i0', adjacent=True, seq_ids=True)),
                         [(' ', 2, ' '), (' ', 3, ' '), (' ', 7, ' '), (' ', 8, ' '), (' ', 9, ' ')])

    def test_define_range_iterator(self):
        fa = """>blah
AAAAAAAAAAAAAAAAAAAAAAAAAAA
..((..((...))..))..((..))..
"""
        bg = fgb.from_fasta_text(fa, dissolve_length_one_stems=False)
        self.assertEqual(list(bg.define_range_iterator('i0')),
                         [[5,6],[14,15]])
        self.assertEqual(list(bg.define_range_iterator('i0', adjacent = True)),
                         [[4,7],[13,16]])

        self.assertEqual(list(bg.define_range_iterator('f0')),
                         [[1,2]])
        self.assertEqual(list(bg.define_range_iterator('t0')),
                         [[26,27]])

    def test_from_dotplot4(self):
        dotbracket = '()'
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket)

        # this structure should have a hairpin
        self.assertTrue('h0' not in bg.defines)

        # or should it?

    def test_from_dotplot3(self):
        dotbracket = '(.(.((((((...((((((....((((.((((.(((..(((((((((....)))))))))..((.......))....)))......))))))))...))))))..)).))))).)..((((..((((((((((...))))))))).))))).......'
        bg = fgb.BulgeGraph()
        self.check_graph_integrity(bg)

        bg.from_dotbracket(dotbracket)
        self.check_graph_integrity(bg)

    def test_from_dotplot2(self):
        bg = fgb.BulgeGraph()

        bg.from_dotbracket('(.(..))')
        self.check_for_overlapping_defines(bg)
        self.check_for_all_nucleotides(bg)

        # secondary structure taken from 1y26
        bg.from_dotbracket('((..))')
        elem_str = bg.to_element_string()
        self.assertEquals(elem_str, "sshhss")

        bg.from_dotbracket('((..))..')
        elem_str = bg.to_element_string()
        self.assertEquals(elem_str, "sshhsstt")

        dotbracket = '..((..))..'
        bg.from_dotbracket(dotbracket)
        elem_str = bg.to_element_string()

        self.assertEquals(elem_str, "ffsshhsstt")

        dotbracket = '..((..))..((..))..'
        bg.from_dotbracket(dotbracket)
        elem_str = bg.to_element_string()

        self.assertEquals(elem_str, "ffsshhssmmsshhsstt")

        dotbracket = '((((((((((..(((((((.......)))))))......).((((((.......))))))..)))))))))'
        bg.from_dotbracket(dotbracket)
        self.check_for_overlapping_defines(bg)
        self.check_for_all_nucleotides(bg)

        dotbracket = '((((((((((..(((((((.......)))))))......).((((((.......))))))..)))))))))'
        bg.from_dotbracket(dotbracket, dissolve_length_one_stems=True)
        self.check_for_overlapping_defines(bg)
        self.check_for_all_nucleotides(bg)


    def test_from_bg_string(self):
        bg = fgb.BulgeGraph()
        bg.from_bg_string(self.bg_string)

        self.assertEquals(bg.seq_length, 71)

    def check_from_and_to_dotbracket(self, dotbracket):
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket)
        self.assertEquals(bg.to_dotbracket_string(), dotbracket)

    def test_to_fasta_string(self):
        filename = 'test/forgi/data/2hoj.fa'
        with open(filename, 'r') as f:
            instr = f.read()
            bg = fgb.from_fasta_text(instr)
            outstr = bg.to_fasta_string()

            self.assertEqual(instr.strip(), outstr.strip())

    def test_get_multiloop_side(self):
        # see page 85 in the notebook
        bg = fgb.BulgeGraph(dotbracket_str='(.().().)')

        s = bg.get_multiloop_side('m0')
        self.assertEqual(s, (1, 0))

        s = bg.get_multiloop_side('m1')
        self.assertEquals(s, (3, 0))

        s = bg.get_multiloop_side('m2')
        self.assertEquals(s, (2, 3))

    def test_get_any_sides(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..((..))..)).((..))')

        self.assertEqual(bg.get_any_sides('s0', 'i0'), (1,0))
        self.assertEqual(bg.get_any_sides('i0', 's0'), (0,1))

        bg = fgb.BulgeGraph(dotbracket_str='((..((..))((..))))')
        bg.log(logging.INFO)
        self.assertEqual(bg.get_any_sides('s1', 'm0'), (0, 1))
        self.assertEqual(bg.get_any_sides('m0', 's1'), (1, 0))

    def test_get_sides(self):
        with open('test/forgi/data/1ymo.bpseq', 'r') as f:
            lines = f.readlines()

        bpseq_str = "".join(lines)
        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(bpseq_str, dissolve_length_one_stems=True)

    def test_get_sides_plus(self):
        bg = fgb.BulgeGraph(dotbracket_str='(.().().)')

        p1 = bg.get_sides_plus('s0', 'm0')
        self.assertEquals(p1[0], 1)

        p1 = bg.get_sides_plus('s0', 'm2')
        self.assertEquals(p1[0], 2)

        p1 = bg.get_sides_plus('s1', 'm0')
        self.assertEquals(p1[0], 0)

        bg = fgb.BulgeGraph(dotbracket_str='(((((((((...(((((((.......)))))))........((((((.......))))))..)))))))))')

        for d in bg.mloop_iterator():
            connections = bg.connections(d)

            (s1c, d1c) = bg.get_sides_plus(connections[0], d)
            (s2c, d2c) = bg.get_sides_plus(connections[1], d)

            self.assertTrue((s1c, s2c) in [(1,0),(3,0),(2,3)])

    def test_to_dotbracket(self):
        self.check_from_and_to_dotbracket('..((..))..')
        self.check_from_and_to_dotbracket('..((..))..((..))')
        self.check_from_and_to_dotbracket('..((..((..))..))')

        pass

    def test_pairing_partner(self):
        # documented
        bg = fgb.BulgeGraph()
        bg.from_dotbracket('((..))')

        self.assertEquals(bg.pairing_partner(1), 6)
        self.assertEquals(bg.pairing_partner(2), 5)
        self.assertEquals(bg.pairing_partner(5), 2)

    def test_big_structure(self):
        bg = fgb.BulgeGraph()
        bg.from_dotbracket('')

    def test_get_bulge_dimensions(self):
        bg = fgb.BulgeGraph(dotbracket_str='((.).)')
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (0,1))

        bg = fgb.BulgeGraph(dotbracket_str='(.(.))')
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (1,0))

        bg = fgb.BulgeGraph(dotbracket_str='().()')
        bd = bg.get_bulge_dimensions('m0')

        dotbracket = '(.(.).(.).(.))'
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket)
        found = col.Counter()
        for loop in bg.defines:
            if loop[0]!="m": continue
            bd = bg.get_bulge_dimensions(loop)
            found[bd]+=1
        self.assertEqual(found[(0,1000)], 1)
        self.assertEqual(found[(1,1000)], 3)

        bg = fgb.BulgeGraph(dotbracket_str='((..((..))....))..((..((..))...))')

        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (2, 4))
        bd = bg.get_bulge_dimensions('i1')
        self.assertEquals(bd, (2, 3))

    def test_get_length(self):
        bg = fgb.BulgeGraph(dotbracket_str='(())')

        bg = fgb.BulgeGraph(dotbracket_str='((..))..(((.)))')

        self.assertEquals(bg.get_length('s0'), 2)
        self.assertEquals(bg.get_length('h0'), 2)
        self.assertEquals(bg.get_length('m0'), 2)
        self.assertEquals(bg.get_length('s1'), 3)

        bg = fgb.BulgeGraph(dotbracket_str='(())(())')
        self.assertEquals(bg.get_length('m0'), 0)

        bg = fgb.BulgeGraph(dotbracket_str='(((((((((..(((..((((.(((((((((.....(((((.(((((....((((....))))....))))).....(((((((((.......)))))))))....))))).((........))...)))))))))))))...)))..))....))))))).')

        self.assertEqual(bg.get_length('i4'), 2)

    def test_get_define_seq_str(self):
        bg = fgb.from_fasta_text(""">1Y26_X
CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
((((((((((..((((((.........))))))......).((((((.......))))))..)))))))))""")
        self.assertEqual(bg.get_define_seq_str("s0"), ['CGCUUCAUA', 'UAUGAAGUG'])
        self.assertEqual(bg.get_define_seq_str("s0", True), ['CGCUUCAUAU', 'UUAUGAAGUG'])

        bg = fgb.BulgeGraph(dotbracket_str="(.(.))")
        bg.seq = 'acgauu'
        self.assertEquals(bg.get_define_seq_str("i0"), ['c', ''])

        self.assertEquals(bg.get_define_seq_str("h0"), ['a'])
        self.assertEquals(bg.get_define_seq_str("h0", True), ['gau'])


        bg = fgb.BulgeGraph(dotbracket_str="(.(.))")
        bg.seq = 'acgauu'
        self.assertEquals(bg.get_define_seq_str("i0", True), ['acg','uu'])

        bg = fgb.BulgeGraph(dotbracket_str='(.(.).(.).)')
        bg.seq = 'acguaaccggu'
        self.assertEquals(bg.get_define_seq_str('m0'), ['c'])
        self.assertEquals(bg.get_define_seq_str('m0', True), ['acg'])

        self.assertEquals(bg.get_define_seq_str('m1'), ['a'])
        self.assertEquals(bg.get_define_seq_str('m1', True), ['aac'])

        self.assertEquals(bg.get_define_seq_str('m2'), ['g'])
        self.assertEquals(bg.get_define_seq_str('m2', True), ['ggu'])

        bg = fgb.BulgeGraph(dotbracket_str=".(.).")
        bg.seq = 'acgau'
        self.assertEquals(bg.get_define_seq_str("f0", True), ['ac'])
        self.assertEquals(bg.get_define_seq_str("f0"), ['a'])

        self.assertEquals(bg.get_define_seq_str("t0", True), ['au'])
        self.assertEquals(bg.get_define_seq_str("t0"), ['u'])

    def check_define_integrity(self, bg):
        """
        Check to make sure that the define regions are always 5' to 3'
        """
        for v in bg.defines.values():
            prev = 0
            i = iter(v)

            # iterate over every other element to make sure that the
            # ones in front involve lower-numbered nucleotides than
            # the ones further on
            for e in i:
                self.assertTrue(e > prev)
                prev = e
                next(i)

    def test_bulge_graph_define_sorting(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..((..))..))..((..((..))...))')

        self.check_define_integrity(bg)

    def test_get_flanking_region(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..))')

        (m1, m2) = bg.get_flanking_region('h0')
        self.assertEqual(m1, 1)
        self.assertEqual(m2, 6)

        bg = fgb.BulgeGraph(dotbracket_str='((.((.)).(.).))')

        (m1, m2) = bg.get_flanking_region('m0')
        self.assertEqual(m1, 1)
        self.assertEqual(m2, 5)

        (m1, m2) = bg.get_flanking_region('m1')
        self.assertEqual(m1, 7)
        self.assertEqual(m2, 10)

        (m1, m2) = bg.get_flanking_region('m2')
        self.assertEqual(m1, 12)
        self.assertEqual(m2, 15)

        bg = fgb.BulgeGraph(dotbracket_str='(.(.).).(.(.))')
        (m1, m2) = bg.get_flanking_region('i1', side=0)
        self.assertEqual(bg.get_flanking_region('i0', side=0),
                         (1,3))
        self.assertEqual(bg.get_flanking_region('i0', side=1),
                         (5,7))
        self.assertEqual(bg.get_flanking_region('i1', side=0),
                         (9,11))
        self.assertEqual(bg.get_flanking_region('i1', side=1),
                         (13,14))

        dotbracket = '...(((((((((((((((((())))))))))))))))))...(((((((((((((((())))))))))))))))'
        seq = fus.gen_random_sequence(len(dotbracket))
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket, seq=seq)
        (m1, m2) = bg.get_flanking_region('m0')

        fa = """>blah
AAAACCGGGCCUUUUACCCCAAAUUGGAA
((((..(((..)))..))))...((..))
"""

    def test_get_flanking_sequence(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..))')
        bg.seq = 'AACCGG'

        self.assertEqual(bg.get_flanking_sequence('h0'),
                         'AACCGG')

        bg = fgb.BulgeGraph(dotbracket_str='((.((.)).(.).))')
        bg.seq = 'AUGCaugcAUGCaug'
        self.assertEqual(bg.get_flanking_sequence('m0'),
                         'AUGCa')
        self.assertEqual(bg.get_flanking_sequence('m1'),
                         'gcAU')
        self.assertEqual(bg.get_flanking_sequence('m2'),
                         'Caug')

        dotbracket = '...(((((((((((((((((())))))))))))))))))...(((((((((((((((())))))))))))))))'
        seq = fus.gen_random_sequence(len(dotbracket))
        print(seq)
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket, seq=seq)
        s = bg.get_flanking_sequence('m0')

    def test_get_flanking_handles(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..))')
        h = bg.get_flanking_handles('h0')

        self.assertEqual(h, (2, 5, 1, 4))

        bg = fgb.BulgeGraph(dotbracket_str='((.((.)).(.).))')

        self.assertEqual(bg.get_flanking_handles('m0'),
                         (2,4,1,3))
        self.assertEqual(bg.get_flanking_handles('m1'),
                         (8,10,1,3))
        self.assertEqual(bg.get_flanking_handles('m2'),
                         (12,14,0,2))

        bg = fgb.BulgeGraph(dotbracket_str='(.(.).).(.(.))')
        self.assertEqual(bg.get_flanking_handles('i0', side=0),
                         (1,3,0,2))
        self.assertEqual(bg.get_flanking_handles('i0', side=1),
                         (5,7,0,2))
        self.assertEqual(bg.get_flanking_handles('i1', side=0),
                         (9,11,0,2))
        self.assertEqual(bg.get_flanking_handles('i1', side=1),
                         (13,14,0,1))

        bg = fgb.BulgeGraph(dotbracket_str='((.((.)).)).((.((.))))')
        #                                   1234567890123456789012
        self.assertEqual(bg.get_flanking_handles('i0', side=0),
                         (2,4,1,3))
        self.assertEqual(bg.get_flanking_handles('i0', side=1),
                         (8,10,1,3))
        self.assertEqual(bg.get_flanking_handles('i1', side=0),
                         (14,16,1,3))
        self.assertEqual(bg.get_flanking_handles('i1', side=1),
                         (20,21,1,2))

    def test_are_adjacent_stems(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..((..))..))..((..))')

        self.assertTrue(bg.are_adjacent_stems('s0', 's1'))
        self.assertTrue(bg.are_adjacent_stems('s0', 's2'))
        self.assertFalse(bg.are_adjacent_stems('s1', 's2'))

        self.assertFalse(bg.are_adjacent_stems('s0', 's2',
                                               multiloops_count=False))

    def test_element_length(self):
        bg = fgb.BulgeGraph(dotbracket_str='.((..(((..))).))((..))')

        self.assertEqual(bg.element_length('s0'), 4)
        self.assertEqual(bg.element_length('i0'), 3)

    def test_stem_length(self):
        bg = fgb.BulgeGraph(dotbracket_str='.((..(((..))).))((..))')

        self.assertEqual(bg.stem_length('s0'), 2)
        self.assertEqual(bg.stem_length('s1'), 3)
        self.assertEqual(bg.stem_length('m0'), 0)
        self.assertEqual(bg.stem_length('i0'), 1)
        self.assertEqual(bg.stem_length('f0'), 1)

    def test_connection_type(self):
        bg = fgb.BulgeGraph(dotbracket_str='(.(.).).(.(.))')

        self.assertEqual(bg.connection_type('m0', ['s0', 's2']), 3)
        self.assertEqual(bg.connection_type('m0', ['s2', 's0']), -3)

        self.assertEqual(bg.connection_type('i0', ['s0', 's1']), 1)
        self.assertEqual(bg.connection_type('i0', ['s1', 's0']), -1)

    def test_random_subgraph(self):
        bg = fgb.BulgeGraph(dotbracket_str='(.(.).).(.(.))..((..((..((..))..))..))')

        sg = bg.random_subgraph()

        # check to make sure there are no duplicate elements
        self.assertEquals(len(sg), len(set(sg)))

    def test_random_subgraph2(self):
        bg = fgb.BulgeGraph(dotbracket_str = "...(((...)))...(((...)))...(((...(((...)))...)))",
                                 seq="AAAGGGAAACCCAAAGGGAAACCCAAAGGGUUUGGGAAACCCUUUCCC")
        for rep in range(10):
            for l in range(3, 8):
                sg = bg.random_subgraph(l)
                self.assertGreaterEqual(len(sg), l)
    def test_has_connection(self):
        bg = fgb.BulgeGraph(dotbracket_str='(())..(())..(())..')

        self.assertTrue(bg.has_connection('m0', 'm1'))
        self.assertTrue(bg.has_connection('m1', 't0'))
        self.assertFalse(bg.has_connection('m0', 't0'))

    def test_compare_hairpins(self):
        bg = fgb.BulgeGraph(dotbracket_str='(())(())')

    def test_create_mst(self):
        """
        Test the creation of a minimum spanning tree from the graph.
        """
        db = '....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).)))).....((((...(((((((((...)))))))))..)))).......'
        bg = fgb.BulgeGraph(dotbracket_str=db)
        mst = bg.get_mst()
        self.assertTrue("m0" in mst)
        build_order = bg.traverse_graph()

        db = '..((.(())..(())...)).'
        bg = fgb.BulgeGraph(dotbracket_str=db)
        mst = bg.get_mst()

        self.assertTrue('m0' in mst)
        self.assertTrue('m1' in mst)

        build_order = bg.traverse_graph()

    def test_create_mst_telomerase(self):
        """
        Test the creation of a minimum spanning tree from the telomerase
        secondary structure.
        """
        bg = fgb.BulgeGraph('test/forgi/data/telomerase.cg')

        mst = bg.get_mst()
        self.assertTrue('m0' not in mst)
        self.assertTrue('m3' not in mst)

    def test_traverse_graph(self):
        # the dotbracket for 1gid
        db = '....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).)))).....((((...(((((((((...)))))))))..)))).......'
        bg = fgb.BulgeGraph(dotbracket_str=db)

        build_order = bg.traverse_graph()
        all_stems = set(bg.stem_iterator())

        for (f, c, t) in build_order:
            if f in all_stems:
                all_stems.remove(f)
            if t in all_stems:
                all_stems.remove(t)

        self.assertTrue(('s0', 'i4', 's1') in build_order)
        self.assertEqual(len(all_stems), 0)

    def test_get_node_from_residue_num(self):
        bg = fgb.BulgeGraph('test/forgi/data/telomerase.cg')

    def test_get_connected_nucleotides(self):
        db = '((..((..))..))'
        bg = fgb.BulgeGraph(dotbracket_str=db)

        cr = bg.get_connected_residues('s0', 's1')

        self.assertEqual(len(cr), 2)
        self.assertEqual(cr[0], [2,5])
        self.assertEqual(cr[1], [13, 10])

        # The order depends on the order of the input stems.
        cr = bg.get_connected_residues('s1', 's0')

        self.assertEqual(len(cr), 2)
        self.assertEqual(cr[0], [5,2])
        self.assertEqual(cr[1], [10, 13])


        db = '((..))..((..))'
        bg = fgb.BulgeGraph(dotbracket_str=db)

        cr = bg.get_connected_residues('s0', 's1')

        self.assertEqual(len(cr), 1)
        self.assertEqual(cr[0], [6,9])

    def test_nucleotides_to_elements(self):
        db = '((..))..((..))'
        bg = fgb.BulgeGraph(dotbracket_str=db)
        self.assertEqual(bg.nucleotides_to_elements([1,5,7,11]), {"s0", "m0", "h1"})

    def test_to_bpseq_str(self):
        bpstr = self.bpseq['1y26']

        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(bpstr)
        bg.to_bpseq_string()

    def test_to_pair_tuples(self):
        bpstr = self.bpseq['1y26']

        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(bpstr)
        pair_tuples = bg.to_pair_tuples()

        self.assertTrue((1,26) in pair_tuples)
        self.assertTrue((26,1) in pair_tuples)
        self.assertTrue((7,0) in pair_tuples)

    def test_to_pairtable(self):
        bpstr = self.bpseq['1y26']

        bg = fgb.BulgeGraph()
        bg.from_bpseq_str(bpstr)

        pt = bg.to_pair_table()

        self.assertEqual(pt[0], 26)
        self.assertEqual(pt[1], 26)
        self.assertEqual(pt[26], 1)
        self.assertEqual(pt[7], 0)

    def test_to_networkx(self):
        fasta = """>1L2X_A
GCGCGGCACCGUCCGCGGAACAAACGG
.(((((..[[[.))))).......]]]
"""
#2345678
        bg = fgb.from_fasta_text(fasta)

        # needs networkx, and for what?
        #bg.to_networkx()


    def test_shortest_bg_loop(self):
        fasta = """>1L2X_A
AAAAAAAAAAAAAAAAAAAAAAAAAA
((...[[...{{...))..]]...}}
"""

        bg = fgb.from_fasta_text(fasta)

        sbl = bg.shortest_bg_loop('m0')
        self.assertIn(3, sbl)
        self.assertIn(4, sbl)
        self.assertIn(18, sbl)
        self.assertNotIn(13, sbl)
        self.assertEqual(len(sbl), 11)

        sbl = bg.shortest_bg_loop('m2')
        self.assertEqual(len(sbl), 14)
        self.assertIn(13, sbl)
        self.assertIn(14, sbl)
        self.assertNotIn(1, sbl)

        sbl = bg.shortest_bg_loop('m4')
        self.assertIn(23, sbl)
        self.assertIn(9, sbl)
        self.assertNotIn(3, sbl)
        self.assertEqual(len(sbl), 12)

    def test_is_loop_pseudoknot(self):
        fasta = """>test
AAAAAAAAAAAAAAAA
((...[[...))..]]
"""
        bg = fgb.from_fasta_text(fasta)
        loops = bg.find_mlonly_multiloops()
        self.assertEqual(len(loops), 1)
        for loop in loops:
            self.assertTrue(bg.is_loop_pseudoknot(loop))

    @unittest.expectedFailure
    def test_is_loop_pseudoknot_2(self):
        """
        The current loop-finding algorithm wrongly decomposes this into 2 loops

        """
        fasta = """>1L2X_A
AAAAAAAAAAAAAAAAAAAAAAAAAA
((...[[...{{...))..]]...}}
"""
        bg = fgb.from_fasta_text(fasta)
        loops = bg.find_mlonly_multiloops()
        for loop in loops:
            log.info("%s, %s", loop, bg.describe_multiloop(loop))
            self.assertTrue(bg.is_loop_pseudoknot(loop))

    def test_is_loop_pseudoknot_no_pk(self):
        dotbracket = '....((((.........)))).((((.(((((..(((((((((....(((.(((..(((..((.((((((((((.....)))))))))).))))).\
.....(((......((((((((..((...(((((((.(((((....((((((....)))))).....)))))....((((.(((((....))))).))))...((((...)))).))))\
)))..))))))))))(((....(((..((((((((.......)))))))))))......)))..((((((((....))))...))))))).(((((............))))).(((((\
..)))))...)))))).).....(.(((...(((((....)))).))))).)).))))))..((((......((((....)))).....))))....(((((...(....((((.....\
)))).....)....)))))......(((((......(((((.....((....)).......))))))))))..)))))))))..........(((.....(.((((...(((.((((((\
(.((((((((((......((((((.....))))))....))))))))..)))))))))..(((((((((...((((((((....((((((....((........)).......))))))\
....).......((....)).)))))))..))))).)))..))))...))))....((((((...((...((((.........))))...))))))))..........((((((..(((\
(((((((...))))))))))...((....)).....)))))))))).(((......((((....))))....)))........(((((.(((((((.((..(((((..((((((((((.\
.....((........))........(.(((((((..(...(............((((....))))...........................)).((.(((...((((((.(....(((\
((((((....)))...(((......)))...)))))).....((((.(((((.(..((...(((.....)))).)...).)))))..(..(((((....))))).....)..))))...\
..).).)))...)).)))))....))))))))..)).)))))))).(...(((((((.....(((..((..((((....))))..))....))).....)))))))......(....((\
(((((........)))))))....)..)..))))).....(((((((.(.....)..)))))))......))...)))))))))).))..(.(..((.(.((((.(((..((((((.((\
((((...(.((((....(((....))).)))).)..)))))).))))))..))).))))..).))...)..)..(((((((((....)))))))))......'

        bg = fgb.from_fasta_text(dotbracket)

        loops = bg.find_mlonly_multiloops()
        for loop in loops:
            if bg.is_loop_pseudoknot(loop):
                log.error("%s is PK", loop)
            self.assertFalse(bg.is_loop_pseudoknot(loop))

    def test_remove_pseudoknots(self):
        fasta = """>1L2X_A
GCGCGGCACCGUCCGCGGAACAAACGG
.(((((..[[[.))))).......]]]
"""
#23456789012345678901234567
        bg = fgb.from_fasta_text(fasta)

        pairs = bg.remove_pseudoknots()
        self.assertTrue((9,27) in pairs)
        self.assertTrue((10,26) in pairs)
        self.assertTrue((11,25) in pairs)

        self.assertEqual(bg.pairing_partner(9), None)
        self.assertEqual(bg.pairing_partner(12), None)

        self.assertEqual(bg.pairing_partner(2), 17)
        self.assertEqual(bg.pairing_partner(2), 17)

#2345678901234

    def test_find_external_loops(self):
        db = '..((..))..'
        bg = fgb.BulgeGraph()

        bg.from_dotbracket(db)
        eloops = bg.find_external_loops()

        self.assertEqual(eloops, ['f0', 't0'])

        bg.from_dotbracket('..((.)).((.))..')
        eloops = bg.find_external_loops()
        self.assertEqual(eloops, ['f0', 't0', 'm0'])

        bg.from_dotbracket('..((.))((.))..')
        eloops = bg.find_external_loops()
        self.assertEqual(eloops, ['f0', 't0', 'm0'])

        bg.from_dotbracket('..(((.))((.)))..(..)')
        eloops = bg.find_external_loops()
        self.assertEqual(eloops, ['f0', 'm3'])

    def test_stem_bp_iterator(self):
        fasta = """>1L2X_A
GCGCGGCACCGUCCGCGGAACAAACGG
.(((((..[[[.))))).......]]]
"""
#23456789012345678901234567
        bg = fgb.from_fasta_text(fasta)

        self.assertEqual(list(bg.stem_bp_iterator("s1")), [(9, 27), (10, 26), (11, 25)])

    def test_ss_distance(self):
        db = '((.((..))..))'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        self.assertEqual(bg.ss_distance('s0', 's1'), 2)
        self.assertEqual(bg.ss_distance('i0', 'h0'), 3)
        self.assertEqual(bg.ss_distance('s0', 's0'), 0)
        self.assertEqual(bg.ss_distance('s0', 'i0'), 1)
        self.assertEqual(bg.ss_distance('s0', 'h0'), 4)

        db = '((..))((..))'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        self.assertEqual(bg.ss_distance('s0', 's1'), 1)

        db = '((((..))..))'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        self.assertEqual(bg.ss_distance('s0', 's1'), 1)

    def test_get_position_in_element(self):
        db = '(((((...))....)))'
        #     12345678901234567

        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        p,l = bg.get_position_in_element(1)
        self.assertEqual(p, 0)
        self.assertEqual(l, 2)
        p,l = bg.get_position_in_element(12)

        self.assertEqual(p, 2)
        self.assertEqual(l, 5)

        db = '((..(((...))....)))'
        #     1234567890123445678
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        p,l = bg.get_position_in_element(3)
        self.assertEqual(p, 1)
        self.assertEqual(l, 3)

        p,l = bg.get_position_in_element(8)
        self.assertEqual(p, 1)
        self.assertEqual(l, 2)

        p,l = bg.get_position_in_element(9)
        self.assertEqual(p, 2)
        self.assertEqual(l, 2)

        db = '((....))'
        #     12345678
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        p,l = bg.get_position_in_element(3)
        self.assertEqual(p, 1)
        self.assertEqual(l, 2)

        p,l = bg.get_position_in_element(6)
        self.assertEqual(p, 1)
        self.assertEqual(l, 2)

        db = '(.).(.)'
        #     12345678
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        p,l = bg.get_position_in_element(4)
        self.assertEqual(p, 1)
        self.assertEqual(l, 2)

        db = '(((......((((.(........).(((..(((((((((.((((.......))))))((((.....))))...(((.....((((.......))))...)))..)))..))))))))))).......))).....'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        p, l = bg.get_position_in_element(14)

        self.assertEqual(p, 1)
        self.assertEqual(l, 2)

        db = '((.(.....).))...........((((..((((((((((((((.......))))).((((.....))))................................)))))..))))..))))................'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        p, l = bg.get_position_in_element(14)

        self.assertEqual(p, 1)
        #self.assertEqual(l, 2)

        db = '(((((...........)))))............'
        #     123456789012345678901234567890123
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        p1, l1 = bg.get_position_in_element(8)
        p2, l2 = bg.get_position_in_element(14)

    def test_connections(self):
        db = '((.(.....).))...........((((..((((((((((((((.......))))).((((.....))))................................)))))..))))..))))................'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        connections = bg.connections('m0')
        self.assertEqual(connections, ['s0', 's2'])

    def test_connected(self):
        db = '((..((..))..((..))..((..))..))'
        # clockwise from the bottom
        # s0 m0 s1 m1 s2 m2 s3 m3
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        self.assertTrue(bg.connected('s0', 'm0'))
        self.assertTrue(bg.connected('m0', 'm1'))

        db = '((..((..)).((..((..))..((..))..))..((..))..))'
        #     123456789012345678901234567890123456789012345
        #     ssmmsshhssmssmmsshhssmmsshhssmmssmmsshhssmmss
        #     000011001112222331133334422444422555533556600

        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.log(logging.INFO)

        # clockwise from the bottom
        # s0 m0 s1 m1 s2 m2 s3 m3 s4 m4 s2 m5 s5 m6 s0

        self.assertTrue(bg.connected('s0', 'm0'))
        self.assertTrue(bg.connected('m0', 'm1'))
        self.assertTrue(bg.connected('m1', 'm5'))
        self.assertTrue(bg.connected('m5', 'm6'))
        self.assertTrue(bg.connected('m6', 'm0'))
        self.assertTrue(bg.connected('m2', 'm3'))
        self.assertTrue(bg.connected('m2', 'm4'))
        self.assertFalse(bg.connected('m1', 'm2'))
        self.assertFalse(bg.connected('m4', 'm5'))
        self.assertFalse(bg.connected('m1', 'm6'))

    def test_cg_min_max_bp_distance(self):
        db = '((..((..))..((..))..((..))..))'
        # clockwise from the bottom
        # s0 m0 s1 m2 s2 m3 s3 m4
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        #fud.pv('bg.to_bg_string()')

        (mi, mx) = bg.min_max_bp_distance('s1', 's2')
        self.assertEqual(mi, 3)
        self.assertEqual(mx, 7)

        (mi, mx) = bg.min_max_bp_distance('s1', 's1')
        self.assertEqual(mi, 0)
        self.assertEqual(mx, 2)

        db =  '..(((..(((..(((..((((((...)))..)))..)))(((...))).(((...(((((((((...))).(((...)))...))).))).)))....))))))..'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        (mi, mx) = bg.min_max_bp_distance('s1', 's10')
        self.assertEqual(mi, 19)
        self.assertEqual(mx, 24)
        (mi, mx) = bg.min_max_bp_distance('s4', 's7')
        self.assertEqual(mi, 18)
        self.assertEqual(mx, 24)
    def test_global_pos_to_stem_pos(self):
        db = '...((((((((...))))))))...'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        self.assertEqual(bg.stem_side_vres_to_resn("s0", 0, 3), 7)
        self.assertEqual(bg.stem_side_vres_to_resn("s0", 1, 3), 19)
        self.assertEqual(bg.stem_resn_to_stem_vres_side("s0",7),(3,0))
        self.assertEqual(bg.stem_resn_to_stem_vres_side("s0",19),(3,1))

        self.assertEqual(bg.stem_side_vres_to_resn("s0", 0, 0), 4)
        self.assertEqual(bg.stem_side_vres_to_resn("s0", 1, 0), 22)
        self.assertEqual(bg.stem_resn_to_stem_vres_side("s0",4),(0,0))
        self.assertEqual(bg.stem_resn_to_stem_vres_side("s0",22),(0,1))

    def test_get_stem_edge(self):
        db = '........(((((((((.(((.((...)).))).)))))).)))..(((....)))....'
        #     123456789012345678901234567890123456789012345678901234567890
        # from 5' end
        # f0 s0 s1 i1 s2 i0 s3 h0 i2 m0 s4 h1 t0
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)

        edge = bg.get_stem_edge('s0',9)
        self.assertEqual(edge, 0)

        edge = bg.get_stem_edge('s0',42)
        self.assertEqual(edge, 1)

        edge = bg.get_stem_edge('s1',12)
        self.assertEqual(edge, 0)

        edge = bg.get_stem_edge('s1',40)
        self.assertEqual(edge, 1)

        edge = bg.get_stem_edge('s2',19)
        self.assertEqual(edge, 0)

        edge = bg.get_stem_edge('s2',33)
        self.assertEqual(edge, 1)

        edge = bg.get_stem_edge('s3',24)
        self.assertEqual(edge, 0)

        edge = bg.get_stem_edge('s3',28)
        self.assertEqual(edge, 1)

        edge = bg.get_stem_edge('s4',48)
        self.assertEqual(edge, 0)

        edge = bg.get_stem_edge('s4',55)
        self.assertEqual(edge, 1)

    def test_shortest_path(self):
        db =  '........(((((((((.(((.((...)).))).)))))).)))..(((....)))....'
        #      123456789012345678901234567890123456789012345678901234567890
        #
        #      ffffffffsssssssssisssisshhhssisssissssssisssmmssshhhhssstttt
        #      000000000001111111222033000330222111111120000044411114440000
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.log(logging.INFO)
        sp = bg.shortest_path('s1', 'h0') # Path traverses stem from base to loop
        self.assertEqual(sp, ['s1', 'i0', 's2', 'i1', 's3', 'h0'])

        sp = bg.shortest_path('s0', 's1') # Includes assymetric bulge w/ length of 0 on one side (i2)
        self.assertEqual(sp, ['s0', 'i2', 's1'])

        sp = bg.shortest_path('f0', 'h0') # Includes assymetric bulge w/ length of 0 on one side (i2)
        self.assertEqual(sp, ['f0', 's0', 'i2', 's1', 'i0', 's2', 'i1', 's3', 'h0'])

        sp = bg.shortest_path('f0', 't0') # Path traverses a multiloop
        self.assertEqual(sp, ['f0', 's0', 'm0', 's4', 't0'])

        sp = bg.shortest_path('h0','h1') # Path traverses stem from loop to base
        self.assertEqual(sp, ['h0', 's3', 'i1', 's2', 'i0', 's1', 'i2', 's0', 'm0', 's4', 'h1'])

        sp = bg.shortest_path('t0','f0') # Shortest path along graph in reverse
        self.assertEqual(sp, ['t0', 's4', 'm0', 's0', 'f0'])

    def test_get_domains(self):
        db =  '..(((..(((..(((..((((((...)))..)))..)))(((...))).(((...(((((((((...))).(((...)))...))).))).)))....))))))..'
        #      1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
        #      ffsssiisssmmsssiisssssshhhsssiisssiisssssshhhsssmsssiiissssssssshhhsssmssshhhsssmmmsssisssisssmmmmsssssstt
        #                                                                             111   111
        #      0000000111002222233344400044411333222225551115552666333777888999222999400033300055588847773666666611100000

        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        dom = bg.get_domains()
        bg.log(logging.INFO)
        log.info("m4: %s %s %s", bg.connections("m4"), bg.get_angle_type("m4"), bg.get_angle_type("m4", allow_broken=True))
        mls = sorted([
                    sorted(['f0', 't0']),
                    sorted(['m0', 'm1', 'm2', 'm6']),
                    sorted(['m3', 'm4', 'm5'])
                    ])
        rods = sorted([
                        sorted(['s0', 'i0', 's1']),
                        sorted(['s2', 'i2', 's3', 'i1', 's4', 'h0']),
                        sorted(['s5', 'h1']),
                        sorted(['s6', 'i3', 's7', 'i4', 's8']),
                        sorted(['s9', 'h2']),
                        sorted(['s10', 'h3'])
        ])
        self.assertEqual(dom["multiloops"], mls)
        self.assertEqual(dom["rods"], rods)
        self.assertEqual(dom["pseudoknots"], [])

    def test_get_domains2(self):
        db =  '(((...(((...((((((...(((...(((...(((((((((...(((...))))))...)))...))))))...))))))...))))))).))'

        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        dom = bg.get_domains()


        self.assertEqual(len(dom["multiloops"]), 0)
        self.assertEqual(len(dom["rods"]), 1)
        self.assertEqual(len(dom["rods"][0]), len(bg.defines))
        self.assertEqual(len(dom["pseudoknots"]), 0)




class MultiloopFinding(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_next_ml_segment(self):
        db="...(((.(((...))).(((...)))...)))..."
        """
        ...(((.(((...))).(((...)))...)))...
        fffsssmssshhhsssmssshhhsssmmmsssttt
        00000001110001112222111222111000000
        """
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        self.assertEqual(bg.get_next_ml_segment("f0"),"t0")
        self.assertEqual(bg.get_next_ml_segment("m0"),"m1")
        self.assertEqual(bg.get_next_ml_segment("m1"),"m2")
        self.assertEqual(bg.get_next_ml_segment("m2"),"m0")
        self.assertEqual(bg.get_next_ml_segment("t0"),None)
    def test_get_next_ml_segment_pk(self):
        db = "..(((..[[[..)))..]]].."
             #ffsssmmsssmmsssmmssstt
             #0000000111110002211100'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        self.assertEqual(bg.get_next_ml_segment("f0"),"m2")
        self.assertEqual(bg.get_next_ml_segment("m2"),"m1")
        self.assertEqual(bg.get_next_ml_segment("m1"),"m0")
        self.assertEqual(bg.get_next_ml_segment("m0"),"t0")
        self.assertEqual(bg.get_next_ml_segment("t0"),None)
    def test_get_next_ml_segment_cofold(self):
        db = "...(((..&...)))"
             #fff   tt fff
             #000   00 111
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.log(logging.INFO)
        self.assertEqual(bg.get_next_ml_segment("f0"),None)
        self.assertEqual(bg.get_next_ml_segment("t0"),None)
        self.assertEqual(bg.get_next_ml_segment("f1"),"t0")
    def test__get_next_ml_segment_no_stem(self):
        bg = fgb.from_fasta_text(".....")
        self.assertEqual(bg.get_next_ml_segment("f0"),None)


    def test_shortest_mlonly_multiloop(self):
        db = "(((..(((...)))(((...))).(((...)))...)))"
        """
        (((..(((...)))(((...))).(((...)))...)))
        sssmmssshhhsssssshhhsssmssshhhsssmmmsss
        000001110001112221112223333222333111000
        """
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        ml = bg.shortest_mlonly_multiloop("m0")
        print(bg.to_dotbracket_string())
        print(bg.to_element_string(True))
        self.assertEqual(ml, ("m0", "m1", "m2", "m3"))
    #@unittest.skip("shortest_mlonly_multiloop")
    def test_shortest_mlonly_multiloop_first_0_length(self):
        db = "((((((...))).(((...))).(((...)))...)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        ml = bg.shortest_mlonly_multiloop("m0")
        print(bg.to_dotbracket_string())
        print(bg.to_element_string(True))
        self.assertEqual(ml, ("m0", "m1", "m2", "m3"))
    #@unittest.skip("shortest_mlonly_multiloop")
    def test_shortest_mlonly_multiloop_first_longest(self):
        db = "(((......(((...))).(((...))).(((...)))...)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        ml = bg.shortest_mlonly_multiloop("m0")
        print(bg.to_dotbracket_string())
        print(bg.to_element_string(True))
        self.assertEqual(ml, ("m0", "m1", "m2", "m3"))
    #@unittest.skip("shortest_mlonly_multiloop")
    def test_shortest_mlonly_multiloop_first_0_length_different_start(self):
        db = "((((((...))).(((...))).(((...)))...)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        ml = bg.shortest_mlonly_multiloop("m2")
        print(bg.to_dotbracket_string())
        print(bg.to_element_string(True))
        self.assertEqual(ml, ("m0", "m1", "m2", "m3"))
    #@unittest.skip("shortest_mlonly_multiloop")
    def test_shortest_mlonly_multiloop_first_longest_different_start(self):
        db = "(((......(((...))).(((...))).(((...)))...)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        ml = bg.shortest_mlonly_multiloop("m2")
        print(bg.to_dotbracket_string())
        print(bg.to_element_string(True))
        self.assertEqual(ml, ("m0", "m1", "m2", "m3"))

    def test_find_mlonly_multiloops_pseudoknotfree(self):
        db = "(((...(((...(((...)))...(((...)))...)))...(((...)))...(((...)))...)))"
            #"sssmmmsssmmmssshhhsssmmmssshhhsssmmmsssmmmssshhhsssmmmssshhhsssmmmsss"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(sorted(mls),sorted([( "m1", "m2", "m3"),("m0", "m4", "m5", "m6" )]))
    def test_find_mlonly_multiloops_pseudoknotfree_0lengthML(self):
        db = "...(((((((((...)))(((...))))))(((...)))(((...))))))..."
             #fffssssssssshhhsssssshhhssssssssshhhsssssshhhssssssttt
             #000000111222000222333111333111444222444555333555000000
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.log(logging.INFO)
        #print("\n".join(map(str,sorted(bg.edges.items()))))
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(sorted(mls),sorted([("f0", "t0"),( "m0", "m4", "m5", "m6"),("m1", "m2", "m3" )]))
    def test_find_mlonly_multiloops_pseudoknot_free_someML0length(self):
        db = "...(((...((((((...)))...(((...)))...)))(((...)))...(((...))))))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(sorted(mls),sorted([("f0", "t0"), ("m0", "m4", "m5", "m6" ), ( "m1", "m2", "m3")]))
    def test_find_mlonly_multiloops_254_pk_witht1(self):
        db = "(((..[[[..)))..]]].."
             #sssmmsssmmsssmmssstt
             #00000111110002211111'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls,[("m0", "t0", "m2", "m1")])
    def test_find_mlonly_multiloops_254_pk_witht1f1(self):
        db = "..(((..[[[..)))..]]].."
             #ffsssmmsssmmsssmmssstt
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls,[("f0", "m2", "m1", "m0", "t0")])


    def test_find_mlonly_multiloops_2534_pk(self):
        db = "(((..[[[..)))(((...)))..]]]"
            # sssmmsssmmsssssshhhsssmmsss
            # 000001111100022200022233111
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls,[("m0", "m2", "m3", "m1")])


    def test_find_mlonly_multiloops_combined(self):
        db = "...(((..[[[..)))(((...)))..]]]...(((...(((...(((.(((..[[[..)))..]]]...)))...(((...)))...)))...(((...)))...(((...)))...)))..."
             #fffsssmmsssmmsssssshhhsssmmsssmmmsssmmmsssmmmsssmsssmmsssmmsssmmsssmmmsssmmmssshhhsssmmmsssmmmssshhhsssmmmssshhhsssmmmsssttt
             #0000000011111000222000222331114443335554446665557666887779966600777111555222888111888333444444999222999555000333000666333000
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.log(logging.INFO)
        mls = bg.find_mlonly_multiloops()
        pprint(mls)
        self.assertEqual(sorted(mls),
                         sorted([("f0", "m2", "m3", "m1", "m0", "m4", "t0"),
                                 ("m7", "m10", "m9", "m8", "m11"),
                                 ("m5", "m14", "m15", "m16"),
                                 ("m6", "m12", "m13")]))

    def test_find_mlonly_multiloops_combined_with_il(self):
        db = "...(((..[[[..)))(((...)))..]]]...(((...(((...(((.(((..[[[..)))..].]]...)))...(((...)))...)))...(((...)))...(((...)))...)))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        pprint(mls)
        self.assertEqual(sorted(mls),
                         sorted([("f0", "m2", "m3", "m1", "m0", "m4", "t0"),
                                 ("m7", "m10", "m9", "m8", "m11"),
                                 ("m5", "m14", "m15", "m16"),
                                 ("m6", "m12", "m13")]))

    def test_find_mlonly_multiloops_pk2455(self):
        db = "(((...(((...)))...(((...[[[...)))...(((...]]]...)))...)))..."
            #"sssmmmssshhhsssmmmsssmmmsssmmmsssmmmsssmmmsssmmmsssmmmsssttt"
            #000000111000111111222222333333222444444555333666444777000000
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.log(logging.INFO)
        mls = bg.find_mlonly_multiloops()
        pprint(mls)
        self.assertEqual(sorted(mls),
                         sorted([("t0",),("m0", "m1", "m4", "m7"), ("m2", "m6", "m5", "m3")]))

    def test_find_mlonly_multiloops_cofold_1(self):
        db = "(((((...&...)))))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls, [("t0", "f0")])

    def test_find_mlonly_multiloops_cofold_2(self):
        db = "(((((...&(((...))))))))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls, [("t0", "m0")])

    def test_find_mlonly_multiloops_f_only(self):
        db = "...(((...)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls, [("f0", )])
    def test_find_mlonly_multiloops_cofold_f_only(self):
        db = "...(((&)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(mls, [("f0", )])

    def test_describe_multiloop_pseudoknot_no_angle_type_5(self):
        # A real world example
        db = "...(((..(((...(((...)))...[[[...(((...)))...)))...(((...)))...]]]...)))..."
        bg = fgb.from_fasta_text(db)
        mls = bg.find_mlonly_multiloops()
        self.assertEqual(len(mls), 2)
        mls.sort(key=lambda x: len(x))
        self.assertEqual(bg.describe_multiloop(mls[0]), set(["open"]))
        self.assertEqual(bg.describe_multiloop(mls[1]), set(["pseudoknot"]))
        # This "strange" pseudoknot has no angle type 5
        for m in bg.mloop_iterator():
            self.assertNotEqual(abs(bg.get_angle_type(m, allow_broken=True)), 5)


class WalkBackboneTests(unittest.TestCase):
    def setUp(self):
        pass
    def test_iter_elements_along_backbones(self):
        db = "...(((...(((...(((...)))...(((...)))...)))...(((...)))...(((...)))...)))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        l = list(bg.iter_elements_along_backbone())
        self.assertEqual(l, ["f0", "s0", "m0", "s1", "m1", "s2", "h0", "s2", "m2", "s3", "h1", "s3",
                            "m3", "s1", "m4", "s4", "h2", "s4", "m5", "s5", "h3", "s5", "m6", "s0",
                            "t0"])
    def test_iter_elements_along_backbones_no_t1(self):
        db = "(((...)))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        l = list(bg.iter_elements_along_backbone())
        self.assertEqual(l, ["s0", "h0", "s0"])
    def test_iter_elements_along_backbones_onesided_i(self):
        db = "((.((...))))"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        l = list(bg.iter_elements_along_backbone())
        self.assertEqual(l, ["s0", "i0", "s1", "h0", "s1", "s0"])
    def test_iter_elements_along_backbones_empty_graph(self):
        bg = fgb.BulgeGraph()
        l = list(bg.iter_elements_along_backbone())
        self.assertEqual(l, [])
    '''
    def test_walk_backbone_pseudoknot_free(self):
        db = "...(((...(((...(((...)))...(((...)))...)))...(((...)))...(((...)))...)))..."
            #"fffsssmmmsssmmmssshhhsssmmmssshhhsssmmmsssmmmssshhhsssmmmssshhhsssmmmsssttt"
            #"111000000111222222000222555333111333333111444444222444666555333555111000111"
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.walk_backbone()
        self.assertEqual(bg.multiloops["pseudoknots"],[])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[])
        self.assertEqual(bg.multiloops["multiloops"],[[ "m2", "m5", "m3"],["m0", "m4", "m6", "m1" ]])

    def test_walk_backbone_pseudoknot_free_0lengthML(self):
        db = "...(((((((((...)))(((...))))))(((...)))(((...))))))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        print("\n".join(map(str,sorted(bg.edges.items()))))
        bg.walk_backbone()
        self.assertEqual(bg.multiloops["pseudoknots"],[])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[])
        self.assertEqual(bg.multiloops["multiloops"],[[ "m1", "m2", "m4"],["m5", "m0", "m3", "m6" ]])
    def test_walk_backbone_pseudoknot_free_someML0length(self):
        db = "...(((...((((((...)))...(((...)))...)))(((...)))...(((...))))))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.walk_backbone()
        self.assertEqual(bg.multiloops["pseudoknots"],[])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[])
        self.assertEqual(bg.multiloops["multiloops"],[[ "m2", "m5", "m3"],["m0", "m4", "m6", "m1" ]])
    def test_walk_backbone_254_pk(self):
        db = "(((..[[[..)))..]]].."
             #sssmmsssmmsssmmssstt
             #00000111110002211111'
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.walk_backbone()
        self.assertEqual(bg.multiloops["pseudoknots"],[["m0", "m1", "m2"]])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[])
        self.assertEqual(bg.multiloops["multiloops"],[])
    def test_walk_backbone_2534_pk(self):
        db = "(((..[[[..)))(((...)))..]]].."
            # sssmmsssmmsssssshhhsssmmssstt
            # 00000111110002220002223311111
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.walk_backbone()
        self.assertEqual(bg.multiloops["pseudoknots"],[["m0", "m1", "m2", "m3"]])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[])
        self.assertEqual(bg.multiloops["multiloops"],[])
    def test_walk_backbone_combined(self):
        db = "...(((..[[[..)))(((...)))..]]]...(((...(((...(((.(((..[[[..)))..]]]...)))...(((...)))...)))...(((...)))...(((...)))...)))..."
             #fffsssmmsssmmsssssshhhsssmmsssmmmsssmmmsssmmmsssmsssmmsssmmsssmmsssmmmsssmmmssshhhsssmmmsssmmmssshhhsssmmmssshhhsssmmmsssttt
             #1110000011111000222000222331114443335554447775550666337774466655777111555222888111888888444999999222999666000333000666333111
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.walk_backbone()
        pprint(bg.multiloops)
        self.assertEqual(bg.multiloops["pseudoknots"],[["m0", "m1", "m2", "m3"], ["m13", "m14", "m15"]])
        self.assertEqual(bg.multiloops["pk_context"],[["m10", "m11"]])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[["m4"]])
        self.assertEqual(bg.multiloops["multiloops"],[["m7", "m12", "m8"], ["m5", "m9", "m16", "m6"]])

    def test_walk_backbone_combined_with_il(self):
        db = "...(((..[[[..)))(((...)))..]]]...(((...(((...(((.(((..[[[..)))..].]]...)))...(((...)))...)))...(((...)))...(((...)))...)))..."
        bg = fgb.BulgeGraph()
        bg.from_dotbracket(db)
        bg.walk_backbone()
        pprint(bg.multiloops)
        self.assertEqual(bg.multiloops["pseudoknots"],[["m0", "m1", "m2", "m3"], ["m13", "m14", "m15"]])
        self.assertEqual(bg.multiloops["pk_context"],[["m10", "m11"]])
        self.assertEqual(bg.multiloops["pseudo_multiloop"],[["m4"]])
        self.assertEqual(bg.multiloops["multiloops"],[["m7", "m12", "m8"], ["m5", "m9", "m16", "m6"]])\
    '''


class SequenceTest(unittest.TestCase):
    def test_subseq_with_cutpoints(self):
        seq = fgb.Sequence("123&456&789")
        self.assertEqual(seq.subseq_with_cutpoints(1,5), "123&4")
        self.assertEqual(seq.subseq_with_cutpoints(1,9), "123&456&78")
        self.assertEqual(seq.subseq_with_cutpoints(5,9), "56&78")
        self.assertEqual(seq.subseq_with_cutpoints(5,None), "56&789")
        seq = fgb.Sequence("12&3&4&56&789")
        self.assertEqual(seq.subseq_with_cutpoints(5,8), "56&7")

    def test_backbone_breaks_after(self):
        seq = fgb.Sequence("123&456&789")
        self.assertEqual(seq.backbone_breaks_after, [3,6])
        seq = fgb.Sequence("12&3&4&56&789")
        self.assertEqual(seq.backbone_breaks_after, [2,3,4,6])

    def test_string_format(self):
        seq_str = "123&456&789"
        seq = fgb.Sequence(seq_str)
        self.assertEqual("{}".format(seq), seq_str)


class BulgeGraphElementNucleotideTests(GraphVerification):
    def test_define_a_s(self):
        bg = fgb.BulgeGraph(dotbracket_str=".((((...)))).")
        self.assertEqual(bg.define_a("s0"), [1,6,8,13])
        bg = fgb.BulgeGraph(dotbracket_str=".((((...))))")
        self.assertEqual(bg.define_a("s0"), [1,6,8,12])
        bg = fgb.BulgeGraph(dotbracket_str="((((...)))).")
        self.assertEqual(bg.define_a("s0"), [1,5,7,12])
    def test_define_a_i(self):
        bg = fgb.BulgeGraph(dotbracket_str="(..(((...))).)")
        self.assertEqual(bg.define_a("i0"), [1,4,12,14])
        bg = fgb.BulgeGraph(dotbracket_str="((((...))).)")
        self.assertEqual(bg.define_a("i0"), [1,2,10,12])
        bg = fgb.BulgeGraph(dotbracket_str="(..(((...))))")
        self.assertEqual(bg.define_a("i0"), [1,4,12,13])
    def test_define_a_m(self):
        bg = fgb.BulgeGraph(dotbracket_str="((...))((...))")
        self.assertEqual(bg.define_a("m0"), [7,8])
        bg = fgb.BulgeGraph(dotbracket_str="((...)).((...))")
        self.assertEqual(bg.define_a("m0"), [7,9])
    def test_define_a_m_pk(self):
        bg = fgb.BulgeGraph(dotbracket_str="(([[[))]]]")
        self.assertEqual(bg.define_a("m0"), [2,3])
        self.assertEqual(bg.define_a("m1"), [5,6])
        self.assertEqual(bg.define_a("m2"), [7,8])
    def test_define_a_ft(self):
        bg = fgb.BulgeGraph(dotbracket_str="..((...))..")
        self.assertEqual(bg.define_a("f0"), [1,3])
        self.assertEqual(bg.define_a("t0"), [9,11])
    def test_all_connections_ft(self):
        bg = fgb.BulgeGraph(dotbracket_str="..((...))..")
        self.assertEqual(bg.all_connections("t0"), ["s0", None])
        self.assertEqual(bg.all_connections("f0"), [None, "s0"])
        self.assertEqual(bg.all_connections("s0"), ["f0", "h0", "h0", "t0"])
        self.assertEqual(bg.all_connections("h0"), ["s0", "s0"])
    def test_all_connections_m_pk(self):
        bg = fgb.BulgeGraph(dotbracket_str="(([[[))]]]")
        #self.assertEqual(bg.define_a("m0"), [2,3])
        #self.assertEqual(bg.define_a("m1"), [5,6])
        #self.assertEqual(bg.define_a("m2"), [7,8])
        self.assertEqual(bg.all_connections("s0"), [None, "m0", "m1", "m2"])
        self.assertEqual(bg.all_connections("s1"), ["m0", "m1", "m2", None])
        self.assertEqual(bg.all_connections("m0"), ["s0", "s1"])
        self.assertEqual(bg.all_connections("m1"), ["s1", "s0"])
        self.assertEqual(bg.all_connections("m2"), ["s0", "s1"])
