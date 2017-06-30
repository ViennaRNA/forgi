from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import unittest, sys
import numpy as np
import numpy.testing as nptest
import networkx as nx
import forgi.threedee.model.coarse_grain as ftmc
import forgi.projection.projection2d as fpp
import forgi.threedee.utilities.vector as ftuv
import matplotlib.pyplot as plt

#TODO: Add tests for the label of the projected segments!!!

class Projection2DBasicTest(unittest.TestCase):
    def setUp(self):
        self.longMessage = True
    def test_init_projection(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2X1F.pdb')
        proj=fpp.Projection2D(cg, [1.,1.,1.])
        self.assertTrue(ftuv.is_almost_colinear(proj.proj_direction,np.array([1.,1.,1.])), 
                        msg="The projection direction was not stored correctly. Should be coliniar"
                            " with {}, got {}".format(np.array([1.,1.,1.]),proj.proj_direction))
    def test_init_projection2(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2X1F.pdb')
        cg.project_from = [1.,1.,1.]
        proj=fpp.Projection2D(cg)
        self.assertTrue(ftuv.is_almost_colinear(proj.proj_direction,np.array([1.,1.,1.])), 
                        msg="The projection direction was not stored correctly. Should be coliniar"
                            " with {}, got {}".format(np.array([1.,1.,1.]),proj.proj_direction))

    def test_init_with_rotate_for_coords(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.])
        projRot=fpp.Projection2D(cg, [1.,1.,1.], rotation = 90)
        define = list(cg.defines.keys())[0]

        self.assertAlmostEqual(np.dot(proj._coords[define][1]- proj._coords[define][0],
                                      projRot._coords[define][1]- projRot._coords[define][0]),
                               0., msg = "Rotating the projection by 90 degrees does not make "
                                         "corresponding coord-vectors orthogonal!")

    def test_init_with_rotate_for_virtual_atoms(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_atoms=True)
        projRot=fpp.Projection2D(cg, [1.,1.,1.], rotation = 90, project_virtual_atoms=True)

        self.assertAlmostEqual(np.dot(proj._virtual_atoms[100]- proj._virtual_atoms[0],
                                      projRot._virtual_atoms[100]- projRot._virtual_atoms[0]),
                               0., msg = "Rotating the projection by 90 degrees does not make "
                                         "corresponding virtual atoms vectors orthogonal!")

    def test_init_with_rotate_for_virtual_residues(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[1,30])
        projRot=fpp.Projection2D(cg, [1.,1.,1.], rotation = 90, project_virtual_residues=[1,30])

        self.assertAlmostEqual(np.dot(proj.get_vres_by_position(1) - proj.get_vres_by_position(30),
                                      projRot.get_vres_by_position(1) - projRot.get_vres_by_position(30)), 
                               0., msg = "Rotating the projection by 90 degrees does not make "
                                          "corresponding virtual residue vectors orthogonal!")

    def test_init_raises_if_no_proj_dir(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/2X1F.pdb')
        with self.assertRaises(ValueError):
            proj=fpp.Projection2D(cg)

    def test_projection_with_virtual_atoms(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_atoms=True)
        all_vas = sum(len(cg.virtual_atoms(i).keys()) for d in cg.defines.keys() 
                                                      for i in cg.define_residue_num_iterator(d))
        self.assertEqual(all_vas, proj._virtual_atoms.shape[0])
        self.assertGreater(proj._virtual_atoms.shape[0], 5)

    def test_projection_with_some_virtual_atoms(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_atoms="selected")
        proj_all=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_atoms=True)
        all_vas = sum(len(cg.virtual_atoms(i).keys()) for d in cg.defines.keys() 
                                                      for i in cg.define_residue_num_iterator(d))
        self.assertLess(proj._virtual_atoms.shape[0], all_vas)
        self.assertLess(proj._virtual_atoms.shape[0], proj_all._virtual_atoms.shape[0])
        self.assertGreater(proj._virtual_atoms.shape[0], 5) #At least some atoms

    def test_projection_with_virtual_residues(self):

        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[1,len(cg.seq)//2])
        elem1 = cg.get_node_from_residue_num(1)
        elem2 = cg.get_node_from_residue_num(len(cg.seq)//2)
        self.assertLess(ftuv.vec_distance(proj._coords[elem1][0], proj.get_vres_by_position(1)), 
                        ftuv.vec_distance(proj._coords[elem2][0], proj.get_vres_by_position(1)),
                        msg =  "Projected virtual residue is closer to projection of a cg-element "
                               "far away  than to the projection of corresponding cg element.")
        self.assertLess(ftuv.vec_distance(proj._coords[elem2][0], 
                                          proj.get_vres_by_position(len(cg.seq)//2)), 
                        ftuv.vec_distance(proj._coords[elem1][0], 
                                          proj.get_vres_by_position(len(cg.seq)//2)),
                        msg =  "Projected virtual residue is closer to projection of a cg-element "
                               "far away than to the corresponding cg element.")

    def test_projection_with_nonexisting_virtual_residues(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        with self.assertRaises(LookupError):
            proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[len(cg.seq)+2])
        with self.assertRaises(LookupError):
            proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[0])

    def test_get_vres_by_position(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[1,4,23])
        self.assertEqual(proj.get_vres_by_position(1).shape, (2,))
        self.assertEqual(proj.get_vres_by_position(4).shape, (2,))
        self.assertEqual(proj.get_vres_by_position(23).shape, (2,))

    def test_get_vres_by_position_raises_if_this_vres_not_projected(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[1,4,23])
        with self.assertRaises(LookupError):
            proj.get_vres_by_position(7)

    def test_get_vres_by_position_raises_if_no_vres_projected(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.])
        with self.assertRaises(LookupError):
            proj.get_vres_by_position(1)

    def test_vres_iterator(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1y26.cg')
        proj=fpp.Projection2D(cg, [1.,1.,1.], project_virtual_residues=[3,7,26])
        self.assertEqual(sorted(list(x[0] for x in proj.vres_iterator)), [3,7,26])

    def test_projection_longest_axis(self):
        cg = ftmc.CoarseGrainRNA()
        cg.from_dotbracket("((()))")
        cg.coords["s0"]=np.array([0.,0,-1]), np.array([0.,0,1])
        proj = fpp.Projection2D(cg, [0.,1.,0.])
        self.assertAlmostEqual(proj.longest_axis, 2)
        proj = fpp.Projection2D(cg, [0.,0.,1.])
        self.assertAlmostEqual(proj.longest_axis, 0)
class Projection2DTestWithData(unittest.TestCase):
    def setUp(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.proj=fpp.Projection2D(cg, [1.,1.,1.])
        self.proj2=fpp.Projection2D(cg, [0.,0.,1.])
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
        self.assertAlmostEqual(self.proj2.get_longest_arm_length()[0], 35.9, places = 1)
        #self.proj.plot(show=True)
        self.proj2.condense(5)
        #self.proj2.plot(show=True)
        self.assertEqual(self.proj2.get_cyclebasis_len(), 2)
        self.assertAlmostEqual(self.proj2.get_longest_arm_length()[0], 20.393064351368327)
 
class Projection2DTestOnCondensedProjection(unittest.TestCase):
    def setUp(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.proj=fpp.Projection2D(cg, [1.,1.,1.])
        self.proj.condense_points(1)
    @unittest.skip("This is a Manual Test")
    def test_plot(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        for dire in [ (1.,0.,0.), (0.,0.,1.), (0.,1.,0.),(1.,1.,0.), (1.,0.,1.), (0.,1.,1.), (1.,1.,1.) ]:
            self.proj=fpp.Projection2D(cg, dire)
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
        self.assertAlmostEqual(self.proj.get_total_length(), 134.45, places = 1)
    def test_get_maximal_path_length(self):
        self.assertAlmostEqual(self.proj.get_maximal_path_length(), 95.4, places = 1)
    def test_get_longest_arm_length(self):
        self.assertAlmostEqual(self.proj.get_longest_arm_length()[0], 20.123454285158793)


class Projection2DImageHelperfunctionTests(unittest.TestCase):
    def setUp(self):
        self.img_bw = np.array([[ 0., 0., 0., 1. ],
                      [ 1., 0., 0., 1. ],
                      [ 1., 0., 0., 0. ],
                      [ 0.5, 0.5, 0., 0. ]])
        self.img_rgb = np.array([[ [0,0,0], [0,0,0], [0,0,0], [255,255,255] ],
                      [ [255,255,255], [0,0,0], [0,0,0], [255,255,255] ],
                      [ [255,255,255], [0,0,0], [0,0,0], [0,0,0] ],
                      [ [127,127,127], [127,127,127], [0,0,0], [0,0,0] ]])
    def test_to_rgb(self):
        np.testing.assert_array_equal(fpp.to_rgb(self.img_bw), self.img_rgb)
    def test_to_grayscale(self):
        np.testing.assert_almost_equal(fpp.to_grayscale(self.img_rgb), self.img_bw, decimal=2)
    def test_rgb_to_rgb(self):
        np.testing.assert_array_equal(fpp.to_rgb(self.img_rgb), self.img_rgb)
    def test_grayscale_to_grayscale(self):
        np.testing.assert_array_equal(fpp.to_grayscale(self.img_bw), self.img_bw)
class Projection2DHelperfunctionsTests(unittest.TestCase):
  def test_rasterize_2d_coordinates_single_point_X1(self):
    a = np.array([[100.,0.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[10,0]]))
  def test_rasterize_2d_coordinates_single_point_X2(self):
    a = np.array([[109.99,0.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[10,0]]))
  def test_rasterize_2d_coordinates_single_point_X3(self):
    a = np.array([[110.,0.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[11,0]]))
  def test_rasterize_2d_coordinates_single_point_Y1(self):
    a = np.array([[0., 100.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[0,10]]))
  def test_rasterize_2d_coordinates_single_point_Y2(self):
    a = np.array([[0., 109.99]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[0,10]]))
  def test_rasterize_2d_coordinates_single_point_Y3(self):
    a = np.array([[0., 110.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[0,11]]))
  def test_rasterize_2d_coordinates_many_points(self):
    a = np.array([[100.,0.],[109.99,0.],[110., 0.], [0., 100.],[0.,109.99], [0.,110.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10)
    nptest.assert_array_equal(a_rast, np.array([[10,0],[10,0], [11,0], [0,10], [0,10], [0,11]]))
  def test_rasterize_2d_coordinates_many_points_with_origin(self):
    a = np.array([[100.,0.],[109.99,0.],[110., 0.], [0., 100.],[0.,109.99], [0.,110.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10, np.array([-0.5,-25]))
    nptest.assert_array_equal(a_rast, np.array([[10,2],[11,2], [11,2], [0,12], [0,13], [0,13]]))
  def test_rasterize_2d_coordinates_single_point_with_origin(self):
    a = np.array([[109.99,0.]])
    b = np.array([[0.,109.99]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10, np.array([-0.5,-25]))
    b_rast = fpp.rasterized_2d_coordinates(b, 10, np.array([-0.5,-25]))
    nptest.assert_array_equal(a_rast, np.array([[11,2]]))
    nptest.assert_array_equal(b_rast, np.array([[0,13]]))
  def test_rasterize_2d_coordinates_many_points_with_rotation(self):
    a = np.array([[100.,0.],[109.99,0.],[110., 0.], [0., 100.],[0.,109.99], [0.,110.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10, rotate = 90)
    print(a_rast)
    nptest.assert_array_equal(a_rast, np.array([[0,-10],[0,-11], [0,-11], [10,0], [10,0], [11,0]]))
  def test_rasterize_2d_coordinates_many_points_with_rotation2(self):
    a = np.array([[100.,0.],[109.99,0.],[110., 0.], [0., 100.],[0.,109.99], [0.,110.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10, rotate = -90)
    print(a_rast)
    nptest.assert_array_equal(a_rast, np.array([[0,10],[0,10], [0,11], [-10,0], [-11,0], [-11,0]]))
  def test_rasterize_2d_coordinates_many_points_with_rotation_and_origin(self):
    a = np.array([[100.,0.],[109.99,0.],[110., 0.], [0., 100.],[0.,109.99], [0.,110.]])
    a_rast = fpp.rasterized_2d_coordinates(a, 10, np.array([-0.5, -25]), rotate = -90)
    print(a_rast)
    nptest.assert_array_equal(a_rast, np.array([[0,12],[0,13], [0,13], [-10,2], [-11,2], [-11,2]]))

