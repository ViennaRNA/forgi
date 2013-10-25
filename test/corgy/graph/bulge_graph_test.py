import unittest, os

import corgy.graph.bulge_graph as cgb
import corgy.utilities.debug as cud

import copy, time

class TestBulgeGraph(unittest.TestCase):
    '''
    Simple tests for the BulgeGraph data structure.

    For now the main objective is to make sure that a graph is created
    and nothing crashes in the process. In the future, test cases for
    bugs should be added here.
    '''

    def setUp(self):
        self.dotbracket = '....((((((....((.......((((.((((.(((...(((((..........)))))...((.......))....)))......))))))))......))...)).))))......(((....((((((((...))))))))...)))........'
        self.bg_string = '''
name temp
length 71
seq CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
define f1 0 1
define h1 47 55
define s3 42 47 55 60
define s2 13 19 27 33
define h0 19 27
define s0 1 9 63 71
define t1 71 72
define m1 9 13
define m2 33 42
define m0 60 63
connect s3 h1 m0 m2
connect s2 h0 m1 m2
connect s0 f1 m1 m0 t1
'''

    def test_from_dotplot(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket(self.dotbracket)
        #print bg.to_bg_string()

        self.assertEquals(bg.seq_length, len(self.dotbracket))

    def test_from_dotplot2(self):
        bg = cgb.BulgeGraph()

        # secondary structure taken from 1y26
        bg.from_dotbracket('((..))')
        elem_str = bg.to_element_string()
        self.assertEquals(elem_str, "sshhss")

        dotbracket = '..((..))..'
        bg.from_dotbracket(dotbracket)
        elem_str = bg.to_element_string()

        self.assertEquals(elem_str, "ffsshhsstt")

        dotbracket = '..((..))..((..))..'
        bg.from_dotbracket(dotbracket)
        elem_str = bg.to_element_string()

        self.assertEquals(elem_str, "ffsshhssmmsshhsstt")

    def test_from_bg_string(self):
        bg = cgb.BulgeGraph()
        bg.from_bg_string(self.bg_string)
        
        self.assertEquals(bg.seq_length, 71)

    def check_from_and_to_dotbracket(self, dotbracket):
        bg = cgb.BulgeGraph(dotbracket)
        self.assertEquals(bg.to_dotbracket(), dotbracket)

    def test_to_dotbracket(self):
        self.check_from_and_to_dotbracket('..((..))..')
        self.check_from_and_to_dotbracket('..((..))..((..))')
        self.check_from_and_to_dotbracket('..((..((..))..))')

        pass

    def test_define_residue_num_iterator(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket('((..))')

        cud.pv('bg.to_bg_string()')
        
        unfixed = list(bg.define_residue_num_iterator('f1', adjacent=False))
        self.assertEquals(len(unfixed), 0)
        unfixed = list(bg.define_residue_num_iterator('t1', adjacent=False))
        self.assertEquals(len(unfixed), 0)

        bg.from_dotbracket('.((..))')

        cud.pv('bg.to_bg_string()')
        
        unfixed = list(bg.define_residue_num_iterator('f1', adjacent=False))
        self.assertEquals(len(unfixed), 1)
        unfixed = list(bg.define_residue_num_iterator('f1', adjacent=True))
        self.assertEquals(len(unfixed), 2)
        unfixed = list(bg.define_residue_num_iterator('t1', adjacent=False))
        self.assertEquals(len(unfixed), 0)
        unfixed = list(bg.define_residue_num_iterator('t1', adjacent=True))
        self.assertEquals(len(unfixed), 0)

        bg.from_dotbracket('((..)).')

        cud.pv('bg.to_bg_string()')
        
        unfixed = list(bg.define_residue_num_iterator('f1', adjacent=False))
        self.assertEquals(len(unfixed), 0)
        unfixed = list(bg.define_residue_num_iterator('f1', adjacent=True))
        self.assertEquals(len(unfixed), 0)

        unfixed = list(bg.define_residue_num_iterator('t1', adjacent=False))
        self.assertEquals(len(unfixed), 1)
        unfixed = list(bg.define_residue_num_iterator('t1', adjacent=True))
        cud.pv('unfixed')
        self.assertEquals(len(unfixed), 2)

    def test_pairing_partner(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket('((..))')

        self.assertEquals(bg.pairing_partner(1), 6)
        self.assertEquals(bg.pairing_partner(2), 5)
        self.assertEquals(bg.pairing_partner(5), 2)
