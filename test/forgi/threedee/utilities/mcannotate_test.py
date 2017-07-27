import os
import sys
import unittest

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.mcannotate as ftum

class TestMCAnnotate(unittest.TestCase):
    '''
    Tests for the pdb loading utility functions.
    '''

    def setUp(self):
        return

    def test_get_interacting_base_pairs(self):
        self.assertEqual(ftum.get_interacting_base_pairs("'0'1161-'7'1191 : A-A Bh/O2P pairing\n"),
                        ("0", "1161", "7", "1191") )
        self.assertEqual(ftum.get_interacting_base_pairs("'0'1161-A1191 : G-C Ww/Ww pairing antiparallel cis\n"),
                        ("0", "1161", "A", "1191") )
        self.assertEqual(ftum.get_interacting_base_pairs("A1161-B1191 : G-U Hh/O2' adjacent_5p inward pairing\n"),
                        ("A", "1161", "B", "1191") )
        self.assertEqual(ftum.get_interacting_base_pairs("B1161-'0'1191 : G-C Ww/Ww pairing antiparallel cis XIX\n"),
                        ("B", "1161", "0", "1191") )
    def test_parse_base_pair_id_normal(self):
        self.assertEqual(ftum.parse_base_pair_id("'0'1161-'0'1191"),
                        ("0", "1161", "0", "1191") )
        self.assertEqual(ftum.parse_base_pair_id("'0'1161-A1191"),
                        ("0", "1161", "A", "1191") )
        self.assertEqual(ftum.parse_base_pair_id("A1161-B1191"),
                        ("A", "1161", "B", "1191") )
        self.assertEqual(ftum.parse_base_pair_id("B1161-'0'1191"),
                        ("B", "1161", "0", "1191") )
    def test_parse_base_pair_id_negative_seqid(self):
        self.assertEqual(ftum.parse_base_pair_id("'0'-123-'0'1191"),
                        ("0", "-123", "0", "1191") )
        self.assertEqual(ftum.parse_base_pair_id("'0'1161-A-1191"),
                        ("0", "1161", "A", "-1191") )
        self.assertEqual(ftum.parse_base_pair_id("A-1161-B-1191"),
                        ("A", "-1161", "B", "-1191") )
        self.assertEqual(ftum.parse_base_pair_id("B-1161-'0'1191"),
                        ("B", "-1161", "0", "1191") )
    def test_parse_base_pair_id_insertion(self):
        self.assertEqual(ftum.parse_base_pair_id("'0'-123.A-'0'1191.B"),
                        ("0", "-123.A", "0", "1191.B") )
        self.assertEqual(ftum.parse_base_pair_id("'0'1161.C-A-1191"),
                        ("0", "1161.C", "A", "-1191") )
        self.assertEqual(ftum.parse_base_pair_id("A-1161-B-1191.D"),
                        ("A", "-1161", "B", "-1191.D") )
        self.assertEqual(ftum.parse_base_pair_id("B1161.E-'0'1191.F"),
                        ("B", "1161.E", "0", "1191.F") )
    def test_parse_chain_base(self):
        self.assertEqual(ftum.parse_chain_base("'9'12"), ("9", "12"))
        self.assertEqual(ftum.parse_chain_base("A46"), ("A", "46"))
        self.assertEqual(ftum.parse_chain_base("F-65"), ("F", "-65"))
        self.assertEqual(ftum.parse_chain_base("'9'12.R"), ("9", "12.R"))
        self.assertEqual(ftum.parse_chain_base("A46.Q"), ("A", "46.Q"))
        self.assertEqual(ftum.parse_chain_base("F-65.S"), ("F", "-65.S"))
