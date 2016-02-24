from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import unittest, sys
import numpy as np
import networkx as nx
import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as ftmp
import forgi.threedee.utilities.vector as ftuv
import matplotlib.pyplot as plt

#TODO: Add tests for the label of the projected segments!!!

class Projection2DBasicTest(unittest.TestCase):
    '''
    Simple tests for the BulgeGraph data structure.

    For now the main objective is to make sure that a graph is created
    and nothing crashes in the process. In the future, test cases for
    bugs should be added here.
    '''

    def setUp(self):
        pass 
    def test_init_projection(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2X1F.pdb')
        proj=ftmp.Projection2D(cg, [1.,1.,1.])
        self.assertTrue(ftuv.is_almost_colinear(proj.proj_direction,np.array([1.,1.,1.])), 
                        msg="The projection direction was not stored correctly. Should be coliniar"
                            " with {}, got {}".format(np.array([1.,1.,1.]),proj.proj_direction))

class Projection2DTestWithData(unittest.TestCase):
    def setUp(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.proj=ftmp.Projection2D(cg, [1.,1.,1.])
        self.proj2=ftmp.Projection2D(cg, [0.,0.,1.])
    def test_get_bounding_square(self):
        s1=self.proj.get_bounding_square()
        s2=self.proj.get_bounding_square(5)
        #is a square
        self.assertAlmostEqual(s1[1]-s1[0], s1[3]-s1[2])
        self.assertAlmostEqual(s2[1]-s2[0], s2[3]-s2[2])
        #s2 around s1
        self.assertTrue(s2[0]<s1[0] and s1[0]<s2[1])
        self.assertTrue(s2[1]>s1[1] and s1[1]>s2[0])
        self.assertTrue(s2[2]<s1[2] and s1[2]<s2[3])
        self.assertTrue(s2[3]>s1[3] and s1[3]>s2[2])
        #All points inside s1
        for point in self.proj.proj_graph.nodes():
            self.assertLessEqual(round(s1[0],7), round(point[0],7), msg="Point {} outside of "
                                              "bounding square {} at the LEFT".format(point, s1))
            self.assertLessEqual( round(point[0],7), round(s1[1], 7), msg="Point {} outside of "
                                              "bounding square {} at the RIGHT".format(point, s1))
            self.assertLessEqual(round(s1[2], 7), round(point[1], 7), msg="Point {} outside of "
                                              "bounding square {} at the BOTTOM".format(point, s1))
            self.assertLessEqual(round(point[1], 7),round(s1[3], 7), msg="Point {} outside of "
                                              "bounding square {} at the TOP".format(point, s1))

    def test_condense_points(self):
        self.proj.condense_points()
        self.assertEqual(len(self.proj.proj_graph.nodes()), 12)
        self.proj.condense_points(100)
        self.assertEqual(len(self.proj.proj_graph.nodes()), 1)
    def test_condense(self):
        self.proj2.condense_points(1)
        self.assertEqual(self.proj2.get_cyclebasis_len(), 1)
        self.assertAlmostEqual(self.proj2.get_longest_arm_length()[0], 35.95081565305)
        #self.proj.plot(show=True)
        self.proj2.condense(5)
        #self.proj2.plot(show=True)
        self.assertEqual(self.proj2.get_cyclebasis_len(), 2)
        self.assertAlmostEqual(self.proj2.get_longest_arm_length()[0], 20.5730384435)
 
class Projection2DTestOnCondensedProjection(unittest.TestCase):
    def setUp(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.proj=ftmp.Projection2D(cg, [1.,1.,1.])
        self.proj.condense_points(1)
    def test_plot(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        for dire in [ (1.,0.,0.), (0.,0.,1.), (0.,1.,0.),(1.,1.,0.), (1.,0.,1.), (0.,1.,1.), (1.,1.,1.) ]:
            self.proj=ftmp.Projection2D(cg, dire)
            self.proj.condense_points(1)
            self.proj.plot(show=True, add_labels=True)
            #plt.show()

    def test_get_cyclebasis_len(self):
        self.assertEqual(self.proj.get_cyclebasis_len(), 2)
    def test_get_branchpoint_count(self):
        self.assertEqual(self.proj.get_branchpoint_count(), 4)
        self.assertEqual(self.proj.get_branchpoint_count(3), 3)
        self.assertEqual(self.proj.get_branchpoint_count(4), 1)
    def test_get_total_length(self):
        self.assertAlmostEqual(self.proj.get_total_length(), 134.377646879)
    def test_get_maximal_path_length(self):
        self.assertAlmostEqual(self.proj.get_maximal_path_length(), 95.29668468051)
    def test_get_longest_arm_length(self):
        self.assertAlmostEqual(self.proj.get_longest_arm_length()[0], 20.301048998)





