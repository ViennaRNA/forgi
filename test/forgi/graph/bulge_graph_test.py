import unittest, os
import itertools as it
import random
import sys

import forgi.graph.bulge_graph as fgb
import forgi.utilities.debug as fud
import forgi.utilities.stuff as fus

import copy, time

class GraphVerification(object):
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

    def check_node_labels(self, bg):
        '''
        There should be only six types of nodes in the graph. The internal
        representation sometimes uses nodes that start with 'x' or 'y' as
        intermediates, but these should always be removed.
        '''
        for k in bg.defines:
            self.assertTrue(k[0] in ['s', 'h', 'i', 'm', 't', 'f'])
    
    def check_graph_integrity(self, bg):
        self.check_node_labels(bg)
        self.check_for_all_nucleotides(bg)
        self.check_for_overlapping_defines(bg)

class TestBulgeGraph(unittest.TestCase, GraphVerification):
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

    def test_to_bg_string(self):
        self.fasta = """>1y26
CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
(((((((((...((((((.........))))))........((((((.......))))))..)))))))))
"""
        bg = fgb.BulgeGraph()
        bg.from_fasta(self.fasta, dissolve_length_one_stems=True)

    def test_from_fasta(self):
        bg = fgb.BulgeGraph()

        with open('test/forgi/threedee/data/3V2F.fa', 'r') as f:
            text = f.read()
            bg.from_fasta(text, dissolve_length_one_stems=False)

        for s in bg.stem_iterator():
            bg.stem_length(s)

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

    def test_from_bpseq(self):
        bg = fgb.BulgeGraph()

        bpstr="""1 G 8
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

        seq= 'AAAAAAAAAAAAAAAAAAAAA'
        db = '.((((..)).))..((..)).'
        n  = '12345678901234567890.'

        bg = fgb.BulgeGraph()
        bpstr="""1 A 0
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

        self.assertEqual(bg.defines['i0'], [10,10])
        self.assertEqual(bg.defines['h0'], [6,7])
        self.assertEqual(bg.defines['h1'], [17,18])
        self.assertEqual(bg.defines['s0'], [2,3,11,12])
        self.assertEqual(bg.defines['s1'], [4,5,8,9])
        self.assertEqual(bg.defines['s2'], [15,16,19,20])
        self.assertEqual(bg.defines['t1'], [21, 21])

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

        db='(.(.(.).).)'
        nm='12345678901'

        bpstr="""1 A 11
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

        db='[.(.].)'
        nm='1234567'

        bpstr="""1 A 5
2 A 0
3 A 7
4 A 0
5 A 1
6 A 0
7 A 3
"""

        bg.from_bpseq_str(bpstr)

        db='[[.((..]]...))'
        nm='12345678901234'

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
28 A 20
"""
        bg.from_bpseq_str(bpstr)

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


    def test_define_residue_num_iterator(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..((..))((..))))')
        drni = bg.define_residue_num_iterator('m2', adjacent=True)
        # the second multiloop should have at least two adjacent nucleotides
        self.assertEqual(len(list(drni)), 2)
        drni = bg.define_residue_num_iterator('m1', adjacent=True)
        # the second multiloop should have at least two adjacent nucleotides
        self.assertEqual(len(list(drni)), 2)

        drni = bg.define_residue_num_iterator('m1', adjacent=True)

        bg = fgb.BulgeGraph()
        bg.from_dotbracket('..((..((...))..))..((..))..')

        self.assertEqual(list(bg.define_residue_num_iterator('f1')),
                         [1,2])
        self.assertEqual(list(bg.define_residue_num_iterator('t1')),
                         [26, 27])
        self.assertEqual(list(bg.define_residue_num_iterator('s1')),
                         [7, 8, 12, 13])
        self.assertEqual(list(bg.define_residue_num_iterator('i0')),
                         [5,6,14,15])

        fa=""">blah
AAAAAAAAAA
((((.)).))
"""
        bg.from_fasta(fa, dissolve_length_one_stems=True)
        self.assertEqual(list(bg.define_residue_num_iterator('i0', adjacent=True)),
                         [2,3,7,8,9])

        self.assertEqual(list(bg.define_residue_num_iterator('i0', adjacent=True, seq_ids=True)),
                         [(' ', 2, ' '), (' ', 3, ' '), (' ', 7, ' '), (' ', 8, ' '), (' ', 9, ' ')])

    def test_define_range_iterator(self):
        bg = fgb.BulgeGraph()
        fa = """>blah
AAAAAAAAAAAAAAAAAAAAAAAAAAA
..((..((...))..))..((..))..
"""
        bg.from_fasta(fa, dissolve_length_one_stems=False)
        self.assertEqual(list(bg.define_range_iterator('i0')),
                         [[5,6],[14,15]])

        r1 = list(bg.define_range_iterator('i0', seq_ids=True))[0]
        srange = list(bg.iterate_over_seqid_range(r1[0], r1[1]))
        self.assertEqual(srange[0], (' ', 5, ' '))

        self.assertEqual(list(bg.define_range_iterator('f1')),
                         [[1,2]])
        self.assertEqual(list(bg.define_range_iterator('t1')),
                         [[26,27]])

        

    def test_dissolve_stem(self):
        '''
        Test to make sure length one stems can be dissolved.
        '''
        bg = fgb.BulgeGraph()
        bg.from_dotbracket('((.(..((..))..).))', dissolve_length_one_stems = True)
        self.assertEquals(bg.to_dotbracket_string(), '((....((..))....))')
        self.check_graph_integrity(bg)

        bg = fgb.BulgeGraph(dotbracket_str='((..))..((..))')
        self.assertEquals(bg.to_dotbracket_string(), '((..))..((..))')
        bg.dissolve_stem('s0')
        self.check_graph_integrity(bg)

        self.assertEquals(bg.to_dotbracket_string(), '........((..))')

        bg.dissolve_stem('s0')
        self.check_graph_integrity(bg)

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
            bg = fgb.from_fasta(filename)
            outstr = bg.to_fasta_string()

            self.assertEqual(instr, outstr)

    def test_get_multiloop_side(self):
        # see page 85 in the notebook
        bg = fgb.BulgeGraph(dotbracket_str='(.().().)')

        s = bg.get_multiloop_side('m0')
        self.assertEqual(s, (1, 0))

        s = bg.get_multiloop_side('m2')
        self.assertEquals(s, (3, 0))

        s = bg.get_multiloop_side('m1')
        self.assertEquals(s, (2, 3))

    def test_get_any_sides(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..((..))..)).((..))')

        self.assertEqual(bg.get_any_sides('s0', 'i0'), (1,0))
        self.assertEqual(bg.get_any_sides('i0', 's0'), (0,1))

        bg = fgb.BulgeGraph(dotbracket_str='((..((..))((..))))')

        self.assertEqual(bg.get_any_sides('s1', 'm1'), (0, 1))
        self.assertEqual(bg.get_any_sides('m1', 's1'), (1, 0))

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

        p1 = bg.get_sides_plus('s0', 'm1')
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

    def test_find_multiloop_loops(self):
        bg = fgb.BulgeGraph()
        bg.from_dotbracket('((..((..))..((..))..))')
        
        bg.find_multiloop_loops()
        bg.from_dotbracket('((..((..((..))..((..))..))..((..))..))')
        bg.from_dotbracket('(.(.(.(.).(.).).(.).))')

    def test_big_structure(self):
        bg = fgb.BulgeGraph()
        bg.from_dotbracket('')

    def test_get_bulge_dimensions(self):
        bg = fgb.BulgeGraph(dotbracket_str='(.(.))')
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (1,0))

        bg = fgb.BulgeGraph(dotbracket_str='((.).)')
        bd = bg.get_bulge_dimensions('i0')
        self.assertEquals(bd, (0,1))

        bg = fgb.BulgeGraph(dotbracket_str='().()')
        bd = bg.get_bulge_dimensions('m0')

        dotbracket = '(.(.).(.).(.))'
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket)
        bd = bg.get_bulge_dimensions('m0')
        self.assertEquals(bd, (1,1000))
        bd = bg.get_bulge_dimensions('m1')
        self.assertEquals(bd, (0,1000))
        bd = bg.get_bulge_dimensions('m2')
        self.assertEquals(bd, (1,1000))
        bd = bg.get_bulge_dimensions('m3')
        self.assertEquals(bd, (1,1000))

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
        bg = fgb.BulgeGraph(dotbracket_str="(.(.))") 
        bg.seq = 'acgauu'
        self.assertEquals(bg.get_define_seq_str("i0"), ['c', ''])

        bg = fgb.BulgeGraph(dotbracket_str="(.(.))") 
        bg.seq = 'acgauu'
        self.assertEquals(bg.get_define_seq_str("i0", True), ['acg','uu'])

        bg = fgb.BulgeGraph(dotbracket_str='(.(.).(.).)')
        bg.seq = 'acguaaccggu'
        self.assertEquals(bg.get_define_seq_str('m0'), ['c'])
        self.assertEquals(bg.get_define_seq_str('m0', True), ['acg'])

        self.assertEquals(bg.get_define_seq_str('m1'), ['g'])
        self.assertEquals(bg.get_define_seq_str('m1', True), ['ggu'])

        self.assertEquals(bg.get_define_seq_str('m2'), ['a'])
        self.assertEquals(bg.get_define_seq_str('m2', True), ['aac'])

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

        (m1, m2) = bg.get_flanking_region('m2')
        self.assertEqual(m1, 7)
        self.assertEqual(m2, 10)

        (m1, m2) = bg.get_flanking_region('m1')
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
        bg.seq = '123456789012345'
        self.assertEqual(bg.get_flanking_sequence('m0'),
                         '12345')
        self.assertEqual(bg.get_flanking_sequence('m2'),
                         '7890')
        self.assertEqual(bg.get_flanking_sequence('m1'),
                         '2345')

        dotbracket = '...(((((((((((((((((())))))))))))))))))...(((((((((((((((())))))))))))))))'
        seq = fus.gen_random_sequence(len(dotbracket))
        bg = fgb.BulgeGraph(dotbracket_str=dotbracket, seq=seq)
        s = bg.get_flanking_sequence('m0')

    def test_get_flanking_handles(self):
        bg = fgb.BulgeGraph(dotbracket_str='((..))')
        h = bg.get_flanking_handles('h0')

        self.assertEqual(h, (2, 5, 1, 4))

        bg = fgb.BulgeGraph(dotbracket_str='((.((.)).(.).))')

        self.assertEqual(bg.get_flanking_handles('m0'),
                         (2,4,1,3))
        self.assertEqual(bg.get_flanking_handles('m2'),
                         (8,10,1,3))
        self.assertEqual(bg.get_flanking_handles('m1'),
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
        self.assertEqual(bg.stem_length('f1'), 1)

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

        nbg = fgb.bg_from_subgraph(bg, sg)
        self.assertTrue(set(nbg.defines.keys()) == set(sg))

    def test_has_connection(self):
        bg = fgb.BulgeGraph(dotbracket_str='(())..(())..(())..')

        self.assertTrue(bg.has_connection('m0', 'm1'))
        self.assertTrue(bg.has_connection('m1', 't1'))
        self.assertFalse(bg.has_connection('m0', 't1'))

    def test_compare_hairpins(self):
        bg = fgb.BulgeGraph(dotbracket_str='(())(())')

    def test_create_mst(self):
        '''
        Test the creation of a minimum spanning tree from the graph.
        '''
        db = '....((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).)))).....((((...(((((((((...)))))))))..)))).......'
        bg = fgb.BulgeGraph(dotbracket_str=db)
        mst = bg.get_mst()
        self.assertTrue("m0" in mst)
        build_order = bg.traverse_graph()

        db = '..((.(())..(())...)).'
        bg = fgb.BulgeGraph(dotbracket_str=db)
        mst = bg.get_mst()

        self.assertTrue('m0' in mst)
        self.assertTrue('m2' in mst)

        build_order = bg.traverse_graph()

    def test_create_mst_telomerase(self):
        '''
        Test the creation of a minimum spanning tree from the telomerase
        secondary structure.
        '''
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


