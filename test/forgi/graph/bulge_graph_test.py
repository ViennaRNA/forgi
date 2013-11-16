import unittest, os
import itertools as it

import forgi.graph.bulge_graph as cgb
import forgi.utilities.debug as cud

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

    def check_for_overlapping_defines(self, bg):
        '''
        Check to make sure none of the defines overlap.
        '''
        for d1,d2 in it.combinations(bg.defines.keys(), 2):
            for dx in bg.defines[d1]:
                for dy in bg.defines[d2]:
                    self.assertNotEqual(dx, dy)

    def check_for_all_nucleotides(self, bg):
        '''
        Check to make sure that the bulge_graph covers each nucleotide
        in the structure.
        '''
        nucs = [False for i in range(bg.seq_length)]

        for d in bg.defines.keys():
            for r in bg.define_residue_num_iterator(d):
                nucs[r-1] = True

        for i,n in enumerate(nucs):
            self.assertTrue(n)


    def test_dissolve_stem(self):
        '''
        Test to make sure length one stems can be dissolved.
        '''
        bg = cgb.BulgeGraph()
        bg.from_dotbracket('((.(..((..))..).))', dissolve_length_one_stems = True)
        self.assertEquals(bg.to_dotbracket(), '((....((..))....))')
        self.check_for_overlapping_defines(bg)

    def test_from_dotplot3(self):
        dotbracket = '(.(.((((((...((((((....((((.((((.(((..(((((((((....)))))))))..((.......))....)))......))))))))...))))))..)).))))).)..((((..((((((((((...))))))))).))))).......'
        bg = cgb.BulgeGraph()

        bg.from_dotbracket(dotbracket)

    def test_from_dotplot2(self):
        bg = cgb.BulgeGraph()

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
        bg = cgb.BulgeGraph()
        bg.from_bg_string(self.bg_string)
        
        self.assertEquals(bg.seq_length, 71)

    def check_from_and_to_dotbracket(self, dotbracket):
        bg = cgb.BulgeGraph(dotbracket_str=dotbracket)
        self.assertEquals(bg.to_dotbracket(), dotbracket)

    def test_get_multiloop_side(self):
        bg = cgb.BulgeGraph(dotbracket_str='(.().().)')

        s = bg.get_multiloop_side('m0')
        self.assertEqual(s, (1, 0))

        s = bg.get_multiloop_side('m1')
        self.assertEquals(s, (3, 0))

        s = bg.get_multiloop_side('m2')
        self.assertEquals(s, (2, 3))

        cud.pv('bg.to_bg_string()')
        cud.pv("bg.get_bulge_dimensions('m0')")
        cud.pv("bg.get_bulge_dimensions('m1')")
        cud.pv("bg.get_bulge_dimensions('m2')")

    def test_get_sides_plus(self):
        bg = cgb.BulgeGraph(dotbracket_str='(.().().)')

        p1 = bg.get_sides_plus('s0', 'm0')
        self.assertEquals(p1[0], 1)

        p1 = bg.get_sides_plus('s0', 'm2')
        self.assertEquals(p1[0], 2)

        p1 = bg.get_sides_plus('s1', 'm0')
        self.assertEquals(p1[0], 0)

    def test_to_dotbracket(self):
        self.check_from_and_to_dotbracket('..((..))..')
        self.check_from_and_to_dotbracket('..((..))..((..))')
        self.check_from_and_to_dotbracket('..((..((..))..))')

        pass

    def test_pairing_partner(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket('((..))')

        self.assertEquals(bg.pairing_partner(1), 6)
        self.assertEquals(bg.pairing_partner(2), 5)
        self.assertEquals(bg.pairing_partner(5), 2)

    def test_find_multiloop_loops(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket('((..((..))..((..))..))')
        
        bg.find_multiloop_loops()

        bg.from_dotbracket('((..((..((..))..((..))..))..((..))..))')

    def test_big_structure(self):
        bg = cgb.BulgeGraph()
        bg.from_dotbracket('')

    def test_get_bulge_dimensions(self):
        bg = cgb.BulgeGraph(dotbracket_str='(.(.))')
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (1,0))

        bg = cgb.BulgeGraph(dotbracket_str='((.).)')
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (0,1))

        bg = cgb.BulgeGraph(dotbracket_str='().()')
        cud.pv('bg.to_bg_string()')
        bd = bg.get_bulge_dimensions('m0')

        bg = cgb.BulgeGraph(dotbracket_str='(.(.).(.).(.))')
        bd = bg.get_bulge_dimensions('m0')
        self.assertEquals(bd, (0,1000))
        bd = bg.get_bulge_dimensions('m1')
        self.assertEquals(bd, (1,1000))
        bd = bg.get_bulge_dimensions('m2')
        self.assertEquals(bd, (1,1000))

        bg = cgb.BulgeGraph(dotbracket_str='((..((..))....))..((..((..))...))')

        bd = bg.get_bulge_dimensions('i1')
        self.assertEquals(bd, (2, 4))
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (2, 3))

    def test_get_define_seq_str(self):
        bg = cgb.BulgeGraph(dotbracket_str="(.())") 
        bg.seq = 'acguu'
        self.assertEquals(bg.get_define_seq_str("i0"), ['c'])

    def check_define_integrity(self, bg):
        '''
        Check to make sure that the define regions are always 5' to 3'
        '''
        for v in bg.defines.values():
            prev = 0
            i = iter(v)

            # iterate over every other element to make sure that the 
            # ones in front involve lower-numbered nucleotides than
            # the ones further on
            for e in i:
                self.assertTrue(e > prev)
                prev = e
                i.next()

    def test_bulge_graph_define_sorting(self):
        bg = cgb.BulgeGraph(dotbracket_str='((..((..))..))..((..((..))...))')

        self.check_define_integrity(bg)

    def test_get_flanking_region(self):
        bg = cgb.BulgeGraph(dotbracket_str='((..))')

        (m1, m2) = bg.get_flanking_region('h0')
        self.assertEqual(m1, 1)
        self.assertEqual(m2, 6)

        bg = cgb.BulgeGraph(dotbracket_str='()()')

        (m1, m2) = bg.get_flanking_region('m0')
        self.assertEqual(m1, 2)
        self.assertEqual(m2, 3)
