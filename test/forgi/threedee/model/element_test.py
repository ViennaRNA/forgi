from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import numpy.testing as nptest

from forgi.threedee.model.Element import CoordinateStorage
import unittest

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
    def test_get_direction(self):
        self.cs["s1"]=[-1,-1,-2],[1,2,4]
        nptest.assert_almost_equal(self.cs.get_direction("s1"), [2,3,6])
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

