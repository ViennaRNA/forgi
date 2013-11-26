import unittest, os
import sys

import forgi.graph.bulge_graph as cgb
import forgi.threedee.model.coarse_grain as cmc
import forgi.utilities.debug as cud
import tempfile as tf

import copy, time

class TestCoarseGrainRNA(unittest.TestCase):
    '''
    Simple tests for the BulgeGraph data structure.

    For now the main objective is to make sure that a graph is created
    and nothing crashes in the process. In the future, test cases for
    bugs should be added here.
    '''

    def setUp(self):
        self.text = '''
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
coord h1 40.4278627091 -20.8674778121 -8.82330482052 56.2509994507 -31.7719993591 -10.3870000839
coord s3 31.3600393191 -17.1294338011 -0.68791944097 40.4278627091 -20.8674778121 -8.82330482052
coord s2 24.2992093432 -28.493700772 -0.129475995658 37.964536177 -34.3139290144 -3.76403225525
coord h0 37.964536177 -34.3139290144 -3.76403225525 50.0929985046 -23.061000824 -16.1790008545
coord s0 2.38464903301 0.227762971329 3.9355757629 19.2429921013 -11.2312887275 4.2510264888
coord m1 24.2992093432 -28.493700772 -0.129475995658 19.2429921013 -11.2312887275 4.2510264888
coord m0 31.3600393191 -17.1294338011 -0.68791944097 19.2429921013 -11.2312887275 4.2510264888
coord m2 31.3600393191 -17.1294338011 -0.68791944097 24.2992093432 -28.493700772 -0.129475995658
twist s3 -0.5773115352 0.273755751171 -0.769265350855 0.683685010863 0.0797372103421 0.725408011542
twist s2 0.354411598448 0.27377713776 0.894113246589 -0.269145115187 -0.0308076403374 -0.96260677136
twist s0 0.0711019690565 0.0772274674423 -0.994474951051 -0.552638293934 -0.807336040837 0.206880238892
'''
        self.bg = cgb.BulgeGraph()
        self.bg.from_bg_string(self.text)
        self.tempfile = tf.NamedTemporaryFile()

        self.tempfile.write(self.bg.to_bg_string())
        self.tempfile.flush()

    def compare_bg_to_cg(self, bg, cg):
        for d in bg.defines.keys():
            self.assertTrue(d in cg.defines.keys())
            self.assertTrue(bg.defines[d] == cg.defines[d])

        for e in bg.edges.keys():
            self.assertTrue(e in cg.edges.keys())
            self.assertTrue(bg.edges[e] == cg.edges[e])

    def test_from_cg_str(self):
        bg = cgb.BulgeGraph()
        cg = cmc.CoarseGrainRNA()
        bg.from_bg_string(self.text)
        cg.from_cg_string(self.text)

        self.compare_bg_to_cg(bg, cg)

    def test_from_file(self):
        cg = cmc.from_file(self.tempfile.name)

        self.compare_bg_to_cg(self.bg, cg)

    def test_from_pdb(self): 
        #cg = cmc.from_pdb('test/forgi/threedee/data/1y26.pdb')
        cg = cmc.from_pdb('test/forgi/threedee/data/RS_118_S_0.pdb', intermediate_file_dir='tmp')

        self.assertTrue(len(cg.defines) > 1)

    def test_from_cg(self):
        cg = cmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        
        self.assertEqual(len(cg.coords), 8)
        for key in cg.defines.keys():
            self.assertTrue(key in cg.coords)

    def test_get_bulge_angle_stats_core(self):
        cg = cmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')

        cud.pv('cg.to_bg_string()')
        cud.pv('cg.to_dotbracket_string()')
        for d in cg.mloop_iterator():
            cg.get_bulge_angle_stats(d)

    def test_read_longrange_interactions(self):
        cg = cmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')

        self.assertGreater(len(cg.longrange), 0)
