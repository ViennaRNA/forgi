import unittest
import numpy as np
import networkx as nx
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.projection2d as ftmp
import forgi.threedee.utilities.vector as ftuv
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
    def test_plot(self):
        self.proj.plot()
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
    def test_len_cycle_basis(self):
        print(self.proj.proj_graph)
        print(nx.cycle_basis(self.proj.proj_graph))
        self.assertEqual(len([x for x in nx.cycle_basis(self.proj.proj_graph) if len(x)>1]), 4)

