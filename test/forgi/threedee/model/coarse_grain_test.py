import numpy as np
import unittest, os
import sys

import itertools as it
import forgi.graph.bulge_graph as cgb
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.utilities.debug as fud
import tempfile as tf

import test.forgi.graph.bulge_graph_test as tfgb

import copy, time


def cg_from_sg(cg, sg):
    '''
    Create a coarse-grain structure from a subgraph.
    
    @param cg: The original structure
    @param sg: The list of elements that are in the subgraph
    '''
    new_cg = ftmc.CoarseGrainRNA()

    for d in sg:
        new_cg.defines[d] = cg.defines[d]

        if d in cg.coords.keys():
            new_cg.coords[d] = cg.coords[d]
        if d in cg.twists.keys():
            new_cg.twists[d] = cg.twists[d]
        if d in cg.longrange.keys():
            new_cg.longrange[d] = cg.longrange[d]

        for x in cg.edges[d]:
            if x in new_cg.defines.keys():
                new_cg.edges[d].add(x)
                new_cg.edges[x].add(d)

    return new_cg


class TestCoarseGrainRNA(tfgb.GraphVerification):
    '''
    Simple tests for the BulgeGraph data structure.

    For now the main objective is to make sure that a graph is created
    and nothing crashes in the process. In the future, test cases for
    bugs should be added here.
    '''

    def setUp(self):
        pass

    def compare_bg_to_cg(self, bg, cg):
        for d in bg.defines.keys():
            self.assertTrue(d in cg.defines.keys())
            self.assertTrue(bg.defines[d] == cg.defines[d])

        for e in bg.edges.keys():
            self.assertTrue(e in cg.edges.keys())
            self.assertTrue(bg.edges[e] == cg.edges[e])

    def check_cg_integrity(self, cg):
        for s in cg.stem_iterator():
            edges = list(cg.edges[s])
            if len(edges) < 2:
                continue

            multiloops = False
            for e in edges:
                if e[0] != 'i':
                    multiloops = True

            if multiloops:
                continue

            self.assertFalse(np.allclose(cg.coords[edges[0]][0],
                                         cg.coords[edges[1]][0]))
            self.assertFalse(np.allclose(cg.coords[edges[0]][0],
                                         cg.coords[edges[1]][1]))
            self.assertFalse(np.allclose(cg.coords[edges[0]][1],
                                         cg.coords[edges[1]][0]))
            self.assertFalse(np.allclose(cg.coords[edges[0]][1],
                                         cg.coords[edges[1]][1]))

    def test_from_cg_str(self):
        pass

        '''
        bg = cgb.BulgeGraph()
        cg = ftmc.CoarseGrainRNA()
        bg.from_bg_string(self.text)
        cg.from_cg_string(self.text)

        self.compare_bg_to_cg(bg, cg)
        '''

    def test_from_file(self):
        pass

    def test_get_node_from_residue_num(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1X8W.pdb',
                           intermediate_file_dir='tmp', chain_id='A')
        self.check_cg_integrity(cg)

        elem_name = cg.get_node_from_residue_num(247, seq_id=True)

    def test_from_pdb(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/RS_363_S_5.pdb')
        self.check_cg_integrity(cg)

        cg = ftmc.from_pdb('test/forgi/threedee/data/1ymo.pdb',
                           intermediate_file_dir='tmp',
                           remove_pseudoknots=False)
        self.check_cg_integrity(cg)

        node = cg.get_node_from_residue_num(25)
        self.assertFalse(node[0] == 'h')

        cg = ftmc.from_pdb('test/forgi/threedee/data/RS_118_S_0.pdb', intermediate_file_dir='tmp')
        self.check_cg_integrity(cg)

        self.assertTrue(len(cg.defines) > 1)

        cg = ftmc.from_pdb('test/forgi/threedee/data/ideal_1_4_5_8.pdb', intermediate_file_dir='tmp')
        self.check_cg_integrity(cg)
        cg = ftmc.from_pdb('test/forgi/threedee/data/ideal_1_4_5_8.pdb', intermediate_file_dir=None)

        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_missing.pdb', intermediate_file_dir='tmp')
        self.check_cg_integrity(cg)

        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb',
                           intermediate_file_dir='tmp', chain_id='Y')
        self.check_cg_integrity(cg)

        cg = ftmc.from_pdb('test/forgi/threedee/data/1X8W.pdb',
                           intermediate_file_dir='tmp', chain_id='A')
        self.check_cg_integrity(cg)

        cg = ftmc.from_pdb('test/forgi/threedee/data/1FJG_reduced.pdb',
                           intermediate_file_dir='tmp')
        self.check_cg_integrity(cg)

        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb',
                           intermediate_file_dir='tmp')

        for d in cg.defines:
            for r in cg.define_residue_num_iterator(d):
                # make sure all the seq_ids are there
                print cg.seq_ids[r - 1]


    def test_from_cg(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

        # self.assertEqual(len(cg.coords), 8)
        for key in cg.defines.keys():
            self.assertTrue(key in cg.coords)

    def test_get_bulge_angle_stats_core(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        self.check_graph_integrity(cg)

        for d in cg.mloop_iterator():
            cg.get_bulge_angle_stats(d)

    def test_read_longrange_interactions(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        self.check_graph_integrity(cg)

        self.assertGreater(len(cg.longrange), 0)

    def test_radius_of_gyration(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        self.check_graph_integrity(cg)

        cg.radius_of_gyration()

    def test_get_coordinates_list(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        self.check_graph_integrity(cg)

        cl = cg.get_coordinates_list()
        self.assertEqual(len(cl), len(cg.defines) * 2)

    def test_get_sides(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1gid.cg')
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

        (s1b, s1e) = cg.get_sides('s0', 'f1')
        (s1b, s1e) = cg.get_sides('s8', 't1')

    def test_cg_from_sg(self):
        bg = ftmc.CoarseGrainRNA(
            dotbracket_str='.(((((..(((.(((((((.((.((((..((((((....))))))..)))).)).))........(((((.....((((...((((....))))...))))...))))).))))).)))...)))))')
        self.check_graph_integrity(bg)

        # bg = cgb.BulgeGraph(dotbracket_str='.(((((........)))))..((((((((.(((.((...........((((((..(((((.((((((((..(((..)))...((((....)))).....))))))))..)))))................((((((...........))))))..((...(((((((...((((((..)))))).....((......))....)))))))...(((((((((.........))))))))).(((....))).))..........(((((.(((((.......))))))))))..........))))..))............(((.((((((((...((.......))...))))))..))))).........((((((((((((..(((((((((......))))))..))).((((.......)))).....)))))..))..))).))....((...............))....))..)))))))))))...')

        for j in range(40):
            sg = bg.random_subgraph()
            new_cg = cg_from_sg(bg, sg)

            for i in it.chain(new_cg.iloop_iterator(), new_cg.mloop_iterator()):
                c = new_cg.connections(i)

                if len(c) != 2:
                    self.assertEqual(len(c), 2)

    def test_define_residue_num_iterator(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2mis.pdb', intermediate_file_dir='tmp')
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

        self.assertEqual(list(cg.define_range_iterator('i0', adjacent=True, seq_ids=True)),
                         [[(' ', 6, ' '), (' ', 10, ' ')], [(' ', 19, ' '), (' ', 21, ' ')]])

    def test_get_stem_stats(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2mis.pdb', intermediate_file_dir='tmp')

        cg.get_stem_stats("s0")

    def test_get_angle_stats(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2mis.pdb', intermediate_file_dir='tmp')

        cg.get_bulge_angle_stats("i0")

    def test_get_loop_stat(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2mis.pdb', intermediate_file_dir='tmp')
        cg.get_loop_stat("h0")

        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/4GXY_A.cg')
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)
        cg.get_loop_stat('h3')

    def test_length_one_stems(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1byj.pdb',
                           intermediate_file_dir='tmp', remove_pseudoknots=False)
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

        cg = ftmc.from_pdb('test/forgi/threedee/data/2QBZ.pdb',
                           intermediate_file_dir='tmp', remove_pseudoknots=False)
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

    def test_pseudoknot(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1YMO.pdb', intermediate_file_dir='tmp')
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3D0U_A.cg')
        self.check_graph_integrity(cg)
        self.check_cg_integrity(cg)

        cg.traverse_graph()
        self.assertEqual(cg.get_angle_type("i3"), 1)

    def test_from_fasta(self):
        cg = ftmc.CoarseGrainRNA()
        with open('test/forgi/threedee/data/1.fa', 'r') as f:
            text = f.read()
            cg.from_fasta(text)

    def test_small_molecule(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2X1F.pdb')

        self.assertTrue('f1' in cg.coords)
