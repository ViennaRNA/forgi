from __future__ import absolute_import
from __future__ import print_function

from builtins import range
import numpy as np
import numpy.testing as nptest
import random

from forgi.threedee.model.linecloud import CoordinateStorage, LineSegmentStorage
import unittest
from math import sin, cos
import copy
import forgi.threedee.model.coarse_grain as ftmc
import itertools as it
import logging
log= logging.getLogger(__name__)

RANDOM_REPETITIONS = 10

def rand(m = 0, d = 10):
    """A random number with a uniform distribution from m-d to m+d"""
    return random.random()*2*d-d+m

class CoordinateStorageTest(unittest.TestCase):
    def setUp(self):
        self.cs = CoordinateStorage(["s1","s2"])
        self.cs2 = CoordinateStorage(['m2', 'h1', 's3', 'i0', 'h0', 's0', 'm1', 'm0', 's2', 's1'])
    def test_set_and_get(self):
        self.cs["s1"]=[1,2,3],[4,5,6]
        self.cs["s2"]=[np.array([5,6,7]), np.array([5.3,61,7.1])]
        nptest.assert_array_equal(self.cs["s1"][0],[1.,2.,3.])
        nptest.assert_array_equal(self.cs["s1"][1],[4,5,6])
        nptest.assert_array_equal(self.cs["s2"][0],[5,6,7])
        nptest.assert_array_equal(self.cs["s2"][1],[5.3,61,7.1])
    def test_set_and_get_with_nparray(self):
        self.cs["s1"]=np.array([[1,2,3],[4,5,6]])
        self.cs["s2"]=np.array([[5,6,7],[5.3,61,7.1]])
        nptest.assert_array_equal(self.cs["s1"][0],[1.,2.,3.])
        nptest.assert_array_equal(self.cs["s1"][1],[4,5,6])
        nptest.assert_array_equal(self.cs["s2"][0],[5,6,7])
        nptest.assert_array_equal(self.cs["s2"][1],[5.3,61,7.1])
    def test_indices_for(self):
        self.assertEqual(self.cs._indices_for("s1"), [0,1])
        self.assertEqual(self.cs._indices_for("s2"), [2,3])
    def test_indices_for2(self):
        self.assertEqual(self.cs2._indices_for("m2"), [0,1])
        self.assertEqual(self.cs2._indices_for("h1"), [2,3])
        self.assertEqual(self.cs2._indices_for("s3"), [4,5])
        self.assertEqual(self.cs2._indices_for("i0"), [6,7])
        self.assertEqual(self.cs2._indices_for("h0"), [8,9])
        self.assertEqual(self.cs2._indices_for("s0"), [10,11])
        self.assertEqual(self.cs2._indices_for("m1"), [12,13])
        self.assertEqual(self.cs2._indices_for("m0"), [14,15])
        self.assertEqual(self.cs2._indices_for("s2"), [16,17])
        self.assertEqual(self.cs2._indices_for("s1"), [18,19])
        #And again in permuted order
        self.assertEqual(self.cs2._indices_for("h1"), [2,3])
        self.assertEqual(self.cs2._indices_for("i0"), [6,7])
        self.assertEqual(self.cs2._indices_for("m1"), [12,13])
        self.assertEqual(self.cs2._indices_for("s1"), [18,19])
        self.assertEqual(self.cs2._indices_for("i0"), [6,7])
        self.assertEqual(self.cs2._indices_for("m2"), [0,1])
        self.assertEqual(self.cs2._indices_for("h0"), [8,9])
        self.assertEqual(self.cs2._indices_for("s0"), [10,11])
        self.assertEqual(self.cs2._indices_for("m1"), [12,13])
        self.assertEqual(self.cs2._indices_for("m0"), [14,15])
        self.assertEqual(self.cs2._indices_for("m0"), [14,15])
        self.assertEqual(self.cs2._indices_for("m2"), [0,1])

    def test_rotate1(self):
        self.cs["s1"]=[0,0,0],[0,0,10]
        rotMat=np.array([[1,0,0],[0,0,-1],[0,1,0]])
        self.cs.rotate(rotMat)
        nptest.assert_almost_equal(self.cs["s1"][1], [0,-10,0])

class CoordinateStorageTest2(unittest.TestCase):
    def setUp(self):
        self.cs = CoordinateStorage(["s1","s2"])
        self.cs2 = CoordinateStorage(['m2', 'h1', 's3', 'i0', 'h0', 's0', 'm1', 'm0', 's2', 's1'])
        self.cs["s1"]=[0,0,0], [0,0,1]
        self.cs["s2"]=[0,0,1], [0,4,1]
        self.cs2["m2"]=[0,0,0],[0,0,1]
        self.cs2["h1"]=[0,1,0],[2,0,1]
        self.cs2["s3"]=[1,2,3],[3,4,1]
        self.cs2["i0"]=[5,6,7.5],[7.8,8,1]
        self.cs2["h0"]=[10,0.5,0],[0,0,17.5]
        self.cs2["s0"]=[0,3,0],[1,0,1]
        self.cs2["m1"]=[0,1,0],[1,1,1]
        self.cs2["m0"]=[9,9,0],[0,9,1]
        self.cs2["s2"]=[0,0,8],[0,0,9]
        self.cs2["s1"]=[6,8,7],[7,6,1]
    def test_get_element_with_list_as_index(self):
        nptest.assert_array_equal(self.cs[["s1"]], np.array([[0.,0.,0.],[0.,0.,1.]]))
        nptest.assert_array_equal(self.cs[["s1", "s2"]], np.array([[0.,0.,0.],[0.,0.,1.],[0.,0.,1.],[0.,4.,1.]]))
        nptest.assert_array_equal(self.cs2[["h1", "i0", "s2", "m0"]], np.array([[0.,1.,0.],[2.,0.,1.],[5.,6.,7.5],[7.8,8.,1.], [0,0,8],[0,0,9], [9,9,0],[0,9,1]]))
    def test_get_element_raises_key_error(self):
        with self.assertRaises(KeyError):
            self.cs["s3"]
        with self.assertRaises(KeyError):
            self.cs2["s2,s1"]
        with self.assertRaises(KeyError):
            self.cs2[None]

class CoordinateStorageEqualityTests(unittest.TestCase):
    def test_eq_with_nans(self):
        cs1 = CoordinateStorage(["s1", "s2", "s3", "s4", "s5"])
        cs2 = CoordinateStorage(["s2", "s3", "s1", "s5", "s4"])
        self.assertEqual(cs1, cs2)
    def test_eq_with_values(self):
        cs1 = CoordinateStorage(["s1", "s2", "s3", "s4", "s5"])
        cs2 = CoordinateStorage(["s2", "s3", "s1", "s5", "s4"])
        cs1["s1"]=[0,0,0], [0,0,1]
        cs1["s2"]=[0,0,1], [0,4,1]
        cs1["s3"]=[0,0,0],[0,0,1]
        cs1["s4"]=[0,1,0],[2,0,1]
        cs1["s5"]=[1,2,3],[3,4,1]
        cs2["s1"]=[0,0,0], [0,0,1]
        cs2["s2"]=[0,0,1], [0,4,1]
        cs2["s3"]=[0,0,0],[0,0,1]
        cs2["s4"]=[0,1,0],[2,0,1]
        cs2["s5"]=[1,2,3],[3,4,1]
        self.assertEqual(cs1, cs2)
    def test_ne_different_keys(self):
        cs1 = CoordinateStorage(["s1", "s2", "s3", "s4", "s5"])
        cs2 = CoordinateStorage(["s2", "s3", "s1", "s5", "i4"])
        self.assertNotEqual(cs1, cs2)
    def test_ne_different_values(self):
        cs1 = CoordinateStorage(["s1", "s2", "s3", "s4", "s5"])
        cs2 = CoordinateStorage(["s2", "s3", "s1", "s5", "s4"])
        cs1["s1"]=[0,0,0], [0,0,1]
        cs1["s2"]=[0,0,1], [0,4,1]
        cs1["s3"]=[0,0,0],[0,0,1]
        cs1["s4"]=[0,1,0],[2,0,1]
        cs1["s5"]=[1,2,3],[3,4,1]
        cs2["s1"]=[0,0,0], [0,0,1]
        cs2["s2"]=[0,0,1], [0,4,1]
        cs2["s3"]=[0,4,5],[0,0,1] #<== DIFFERENT
        cs2["s4"]=[0,1,0],[2,0,1]
        cs2["s5"]=[1,2,3],[3,4,1]
        self.assertNotEqual(cs1, cs2)
    def test_ne_some_nans(self):
        cs1 = CoordinateStorage(["s1", "s2", "s3", "s4", "s5"])
        cs2 = CoordinateStorage(["s2", "s3", "s1", "s5", "s4"])
        cs1["s1"]=[0,0,0], [0,0,1]
        cs1["s2"]=[0,0,1], [0,4,1]
        cs1["s3"]=[0,0,0], [0,0,1]
        cs1["s4"]=[0,1,0], [2,0,1]
        cs1["s5"]=[1,2,3], [3,4,1]
        cs2["s1"]=[0,0,0], [0,0,1]
        cs2["s2"]=[0,0,1], [0,4,1]
        # "s3" is not initialized
        cs2["s4"]=[0,1,0],[2,0,1]
        cs2["s5"]=[1,2,3],[3,4,1]
        self.assertNotEqual(cs1, cs2)
    def test_eq_for_cgs(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg1_a = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')
        cg2_a = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')
        cg3 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1J1U.cg')
        cg3_a = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1J1U.cg')
        self.assertEqual(cg1.coords, cg1.coords)
        self.assertEqual(cg1.coords, cg1_a.coords)
        self.assertEqual(cg2.coords, cg2_a.coords)
        self.assertEqual(cg3.coords, cg3_a.coords)
        self.assertNotEqual(cg1.coords, cg2.coords)
        self.assertNotEqual(cg1.coords, cg3.coords)
class LineSegmentStorageTests(unittest.TestCase):
    def test_elements_closer_than_on_same_line(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,0.],[0,0,1.]
        cs["s1"]=[0,0,2.],[0,0,3.]
        self.assertEqual(cs.elements_closer_than(2),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(1.01),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(1.),[])
        self.assertEqual(cs.elements_closer_than(0.5),[])
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,0.],[1,0,0.]
        cs["s1"]=[2,0,0.],[3,0,0.]
        self.assertEqual(cs.elements_closer_than(2),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.5),[])

    def test_elements_closer_than_overlapping(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,0.],[0,0,1.]
        cs["s1"]=[0,0,0.5],[0,0,3.]
        self.assertEqual(cs.elements_closer_than(2),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.5),[("s0","s1")])

    def test_elements_closer_than_parallel_ol(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,0.],[0,0,1.]
        cs["s1"]=[0,1,0.5],[0,1,3.]
        self.assertEqual(cs.elements_closer_than(2),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.5),[])
    def test_elements_closer_than_parallel_far(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,0.],[0,0,1.]
        cs["s1"]=[0,1,2],[0,1,3.]
        self.assertEqual(cs.elements_closer_than(1.5),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(1),[])
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s1"]=[0,0,0.],[0,0,1.]
        cs["s0"]=[0,1,2],[0,1,3.]
        self.assertEqual(cs.elements_closer_than(1.5),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(1),[])

    def test_elements_closer_than_intersecting(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[-1,-1,-1.],[1,1,1.]
        cs["s1"]=[-1,-1,1.],[1,1,-1.]
        self.assertEqual(cs.elements_closer_than(1.5),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.001),[("s0","s1")])

    def test_elements_closer_than_in_plane(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,3.],[0,1,3.]
        cs["s1"]=[1,0.5,3.],[2,0.5,-3.]
        self.assertEqual(cs.elements_closer_than(1.1),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.9),[])

    def test_elements_closer_than_windschief(self):
        cs = LineSegmentStorage(["s0", "s1"])
        cs["s0"]=[0,0,0.],[0,0,3.]
        cs["s1"]=[1,1,3.],[1,-1,-3.]
        self.assertEqual(cs.elements_closer_than(1.5),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.4),[])

    def test_elements_closer_than_ignore(self):
        cs = LineSegmentStorage(["s0", "s1", "s2"])
        cs["s0"]=[0,0,0.],[0,0,3.]
        cs["s1"]=[1,1,3.],[1,-1,-3.]
        cs["s2"]=[0.,0.,0.],[1,-1,-3.]

        self.assertEqual(cs.elements_closer_than(1.5, ignore=set([("s0","s2"),("s2", "s1")])),[("s0","s1")])
        self.assertEqual(cs.elements_closer_than(0.4),[("s0","s2"),("s1", "s2")])
        self.assertEqual(cs.elements_closer_than(0.4, ignore=set([("s0","s2"),("s2", "s1")])),[])
        self.assertEqual(cs.elements_closer_than(0.4, ignore=[("s2","s0")]),[("s1", "s2")])

    def test_elements_closer_real_world_example(self):
        cs = LineSegmentStorage(["m0", "m1"])
        cs["m0"] = [0.,0.,1.],[-2.76245752, -6.86976093,  7.54094508]
        cs["m1"] = [-27.57744115,   6.96488989, -22.47619655], [-16.93424799,  -4.0631445 , -16.19822301]
        self.assertEqual(cs.elements_closer_than(25),[("m0","m1")])


    def test_elements_closer_than_like_old_confusion_matrix(self):
        BP_DIST = 16
        CUTOFF_DIST = 25.
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')
        ignore = set()
        for n1, n2 in it.combinations(cg2.defines.keys(), r=2):
            if cg2.connected(n1, n2):
                ignore.add((n1,n2))
                continue
            bp_dist = cg2.min_max_bp_distance(n1, n2)[0]
            if bp_dist < BP_DIST:
                ignore.add((n1,n2))
        old = interactions_old(cg2, CUTOFF_DIST, BP_DIST)
        new = set(cg2.coords.elements_closer_than(CUTOFF_DIST, ignore = ignore))
        print("m0-m1: {} A phys, dist, "
              "{} bp dist".format(cg2.element_physical_distance("m0", "m1"),
                                  cg2.min_max_bp_distance("m0", "m1")[0]))
        print("m0: {}".format(cg2.coords["m0"]))
        print("m1: {}".format(cg2.coords["m1"]))
        self.assertEqual(old, new,
                        msg = "ONLY old: {}\n, ONLY new {},\n {} both".format(old-new, new-old, len(old&new)))



    def test_rmsd_to_self(self):
        cs1 = LineSegmentStorage(["s0", "s1", "s2"])
        for r in range(RANDOM_REPETITIONS):
            cs1["s0"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs1["s1"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs1["s2"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            self.assertAlmostEqual(cs1.rmsd_to(cs1), 0)

    def test_rmsd_to_is_a_symmetric_relation(self):
        cs1 = LineSegmentStorage(["s0", "m3", "s1", "i0", "s2"])
        cs2 = LineSegmentStorage(["s0", "s1", "s2", "i0", "m3"])

        for r in range(RANDOM_REPETITIONS):
            cs1["s0"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs1["s1"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs1["s2"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs1["i0"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs1["m3"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs2["s0"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs2["s1"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs2["s2"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs2["i0"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            cs2["m3"]=[rand(),rand(),rand()],[rand(),rand(),rand()]
            self.assertAlmostEqual(cs1.rmsd_to(cs2), cs2.rmsd_to(cs1))


    def test_rmsd_to_offset_ordered(self):
        cs1 = LineSegmentStorage(["s0", "s1", "s2"])
        cs2 = LineSegmentStorage(["s0", "s1", "s2"])
        cs1["s0"]=[0,0,0.],[0,0,3.]
        cs1["s1"]=[1,1,3.],[1,-1,-3.]
        cs1["s2"]=[0.,0.,0.],[1,-1,-3.]
        cs2["s0"]=[0,0,1.],[0,0,4.]
        cs2["s1"]=[1,1,4.],[1,-1,-2.]
        cs2["s2"]=[0.,0.,1.],[1,-1,-2.]
        self.assertAlmostEqual(cs1.rmsd_to(cs2), 0)
        self.assertAlmostEqual(cs2.rmsd_to(cs1), 0)

    def test_rmsd_to_offset_unordered(self):
        cs1 = LineSegmentStorage(["s0", "s1", "s2"])
        cs2 = LineSegmentStorage(["s2", "s0", "s1"])
        cs1["s0"]=[0,0,0.],[0,0,3.]
        cs1["s1"]=[1,1,3.],[1,-1,-3.]
        cs1["s2"]=[0.,0.,0.],[1,-1,-3.]
        cs2["s0"]=[0,0,1.],[0,0,4.]
        cs2["s1"]=[1,1,4.],[1,-1,-2.]
        cs2["s2"]=[0.,0.,1.],[1,-1,-2.]
        self.assertAlmostEqual(cs1.rmsd_to(cs2), 0)

    def test_rmsd_to_rotated_ordered_unordered(self):
        cs1 = LineSegmentStorage(["s0", "s1", "s2"])
        cs2 = LineSegmentStorage(["s2", "s0", "s1"])
        cs1["s0"]=[0,0,0.],[0,0,3.]
        cs1["s1"]=[1,1,3.],[1,-1,-3.]
        cs1["s2"]=[0.,0.,0.],[1,-1,-3.]
        cs_temp = copy.deepcopy(cs1)
        cs_temp.rotate(np.array([[1,0,0],[0, cos(1.2), -sin(1.2)],[0, sin(1.2), cos(1.2)]]))

        self.assertAlmostEqual(cs1.rmsd_to(cs_temp), 0)

        cs2["s0"] = cs_temp["s0"]
        cs2["s1"] = cs_temp["s1"]
        cs2["s2"] = cs_temp["s2"]


        self.assertAlmostEqual(cs1.rmsd_to(cs2), 0)

    def test_rmsd_to_deviating(self):
        cs1 = LineSegmentStorage(["s0", "s1", "s2", "s3"])
        cs2 = LineSegmentStorage(["s2", "s0", "s3", "s1"])
        cs1["s0"]=[0,0,0.],[0,0,3.]
        cs1["s1"]=[1,1,3.],[1,-1,-3.]
        cs1["s2"]=[0.,0.,0.],[1,-1,-3.]
        cs1["s3"]=[0.,10.,0.],[1,-1,-3.]

        cs2["s0"]=[0,0,3.],[0,0,0.]
        cs2["s1"]=[1,-1,-3.],[1,1,3.]
        cs2["s2"]=[-1.,-1.,-3.],[0,0,0.]
        cs2["s3"]=[0.,5.,0.],[1,-1,-3.]
        self.assertGreater(cs1.rmsd_to(cs2), 1)
        self.assertLess(cs1.rmsd_to(cs2), 5)

    def test_get_direction(self):
        cs1 = LineSegmentStorage(["s0", "s1", "s2", "s3"])
        cs1["s1"]=[-1,-1,-2],[1,2,4]
        nptest.assert_almost_equal(cs1.get_direction("s1"), [2,3,6])


def interactions_old(cg, distance, bp_distance=16):
    """Code from Peter's confusion matrix, for verification"""
    '''
    Calculate the true_positive, false_positive,
    true_negative and false_negative rate for the tertiary
    distances of the elements of two structures, cg1 and cg2.


    :param cg: A coarse grain model
    :param distance: The distance to consider for interactions
    :param bp_distance: Only consider elements separated by this many more pair
                        and backbone bonds
    :return: A set of node-pairs.
    '''

    nodes2 = set(cg.defines.keys())

    tp = 0 #true positive
    tn = 0 #true negative

    fp = 0 #false positive
    fn = 0 #false negative

    #fud.pv('nodes1')
    #fud.pv('nodes2')

    hits_cg2 = set()
    for n1, n2 in it.combinations(nodes2, r=2):


        if cg.connected(n1, n2):
            continue


        bp_dist = cg.min_max_bp_distance(n1, n2)[0]
        if bp_dist < bp_distance:
            continue


        dist2 = cg.element_physical_distance(n1, n2)


        if dist2 < distance:
            if n1<n2:
                hits_cg2.add((n1,n2))
            else:
                hits_cg2.add((n2,n1))



    return hits_cg2
