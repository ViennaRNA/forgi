import os
import sys
import unittest
import textwrap

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.mcannotate as ftum
import forgi.graph.residue as fgr


class TestMCAnnotate(unittest.TestCase):
    '''
    Tests for the pdb loading utility functions.
    '''

    def setUp(self):
        return

    def test_get_interacting_base_pairs(self):
        self.assertEqual(ftum.get_interacting_base_pairs("'0'1161-'7'1191 : A-A Bh/O2P pairing\n"),
                         ("0", "1161", "7", "1191"))
        self.assertEqual(ftum.get_interacting_base_pairs("'0'1161-A1191 : G-C Ww/Ww pairing antiparallel cis\n"),
                         ("0", "1161", "A", "1191"))
        self.assertEqual(ftum.get_interacting_base_pairs("A1161-B1191 : G-U Hh/O2' adjacent_5p inward pairing\n"),
                         ("A", "1161", "B", "1191"))
        self.assertEqual(ftum.get_interacting_base_pairs("B1161-'0'1191 : G-C Ww/Ww pairing antiparallel cis XIX\n"),
                         ("B", "1161", "0", "1191"))

    def test_parse_base_pair_id_normal(self):
        self.assertEqual(ftum.parse_base_pair_id("'0'1161-'0'1191"),
                         ("0", "1161", "0", "1191"))
        self.assertEqual(ftum.parse_base_pair_id("'0'1161-A1191"),
                         ("0", "1161", "A", "1191"))
        self.assertEqual(ftum.parse_base_pair_id("A1161-B1191"),
                         ("A", "1161", "B", "1191"))
        self.assertEqual(ftum.parse_base_pair_id("B1161-'0'1191"),
                         ("B", "1161", "0", "1191"))

    def test_parse_base_pair_id_negative_seqid(self):
        self.assertEqual(ftum.parse_base_pair_id("'0'-123-'0'1191"),
                         ("0", "-123", "0", "1191"))
        self.assertEqual(ftum.parse_base_pair_id("'0'1161-A-1191"),
                         ("0", "1161", "A", "-1191"))
        self.assertEqual(ftum.parse_base_pair_id("A-1161-B-1191"),
                         ("A", "-1161", "B", "-1191"))
        self.assertEqual(ftum.parse_base_pair_id("B-1161-'0'1191"),
                         ("B", "-1161", "0", "1191"))

    def test_parse_base_pair_id_insertion(self):
        self.assertEqual(ftum.parse_base_pair_id("'0'-123.A-'0'1191.B"),
                         ("0", "-123.A", "0", "1191.B"))
        self.assertEqual(ftum.parse_base_pair_id("'0'1161.C-A-1191"),
                         ("0", "1161.C", "A", "-1191"))
        self.assertEqual(ftum.parse_base_pair_id("A-1161-B-1191.D"),
                         ("A", "-1161", "B", "-1191.D"))
        self.assertEqual(ftum.parse_base_pair_id("B1161.E-'0'1191.F"),
                         ("B", "1161.E", "0", "1191.F"))

    def test_parse_chain_base(self):
        self.assertEqual(ftum.parse_chain_base("'9'12"), ("9", "12"))
        self.assertEqual(ftum.parse_chain_base("A46"), ("A", "46"))
        self.assertEqual(ftum.parse_chain_base("F-65"), ("F", "-65"))
        self.assertEqual(ftum.parse_chain_base("'9'12.R"), ("9", "12.R"))
        self.assertEqual(ftum.parse_chain_base("A46.Q"), ("A", "46.Q"))
        self.assertEqual(ftum.parse_chain_base("F-65.S"), ("F", "-65.S"))

    def test_get_dotplot_neg_residuenumbers(self):
        dotplot_negative_residues = textwrap.dedent("""\
                Residue conformations -------------------------------------------
                R-40 : C C3p_endo anti
                R-39 : C C3p_endo anti
                R-38 : C C2p_endo anti
                R-37 : C C3p_endo anti
                R-36 : G C4p_exo anti
                R-35 : A C2p_endo syn
                R-34 : A C2p_endo syn
                R-33 : G C3p_endo anti
                R-32 : G C3p_endo anti
                R-31 : G C2p_endo anti
                Adjacent stackings ----------------------------------------------
                R-40-R-39 : adjacent_5p upward
                R-36-R-35 : adjacent_5p inward
                R-34-R-33 : adjacent_5p outward
                R-33-R-32 : adjacent_5p upward
                Non-Adjacent stackings ------------------------------------------
                Number of stackings = 4
                Number of adjacent stackings = 4
                Number of non adjacent stackings = 0
                Base-pairs ------------------------------------------------------
                R-40-R-30 : C-G Wh/Bh pairing parallel cis one_hbond 124
                R-39-R-31 : C-G Wh/Wh pairing parallel cis one_hbond 124
                R-38-R-37 : C-C O2'/Hh adjacent_5p pairing
                R-38-R-33 : C-G Ww/Ww pairing antiparallel cis one_hbond 130
                R-32-R-31 : G-G Ww/Hw adjacent_5p pairing antiparallel cis one_hbond
                """)
        dp = ftum.get_dotplot(dotplot_negative_residues.splitlines())
        self.assertEqual(dp[1][0], fgr.RESID("R", (" ", -40, " ")))
        self.assertEqual(
            dp[0], "1 C 0\n2 C 0\n3 C 8\n4 C 0\n5 G 0\n6 A 0\n7 A 0\n8 G 3\n9 G 0\n10 G 0\n")
