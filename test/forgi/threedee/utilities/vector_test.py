from __future__ import absolute_import

import unittest
import numpy as np
import forgi.threedee.utilities.vector as ftuv

REPEAT_TESTS_CONTAINING_RANDOM = 10

class TestVector(unittest.TestCase):
    """Tests for the threedee.utilities.vector module"""
    def setUp(self):
        return

    def test_get_inter_distances(self):
        vecs=[np.array([1., 0., 0.]), np.array([0., 0., 0.]), np.array([0., 0., 0.]), np.array([-1., 0., 0.])]
        distances=ftuv.get_inter_distances(vecs)
        self.assertEqual(sorted(distances), [0,1,1,1,1,2])
    def test_get_random_vector(self):
        for _ in range(REPEAT_TESTS_CONTAINING_RANDOM):
            vec=ftuv.get_random_vector()
            self.assertLessEqual(abs(vec[0]),1)        
            self.assertLessEqual(abs(vec[1]),1)          
            self.assertLessEqual(abs(vec[2]),1)
            vec1=ftuv.get_random_vector()
            vec2=ftuv.get_random_vector()
            self.assertTrue(all( vec1[j]!=vec2[j] for j in [0,1,2]),
                              msg="Repeated calls should generate different results."
                                  "This tests depends on random values. if it fails, try running it again.")
    def test_get_orthogonal_unit_vector(self):
        vecs=[np.array([1., 0., 0.]), np.array([2.7, 5.6, 8.2]), np.array([11., -40., 0.]), np.array([-1., 0., 0.])]
        for vec in vecs:
            ortVec=ftuv.get_orthogonal_unit_vector(vec)
            self.assertAlmostEqual(np.dot(ortVec, vec), 0, places=10)
            self.assertAlmostEqual(np.linalg.norm(ortVec), 1, places=10)
        #Every vector is orthogonal to the zero-vector:
        vec=np.array([0., 0., 0.])
        #ortVec=ftuv.get_orthogonal_unit_vector(vec)
        #Currently, ortVec==nan, so the assertion fails.
        #self.assertAlmostEqual(np.dot(ortVec, vec), 0, places=10) 
        #self.assertAlmostEqual(np.linalg.norm(ortVec), 1, places=10)
    def test_seg_intersect(self):
        #normal case
        isec=ftuv.seg_intersect(([0.,1.], [0., -1.]), ([-1.,0.], [1.,0.]))
        self.assertEqual(len(isec), 1)
        np.testing.assert_allclose(isec[0], [0., 0.])
        #parallel, no intersection
        isec=ftuv.seg_intersect(([0., 3.],[1., 3.]),([2.,3.], [3.,3.]))
        self.assertEqual(isec, [])
        #one inside other
        isec=ftuv.seg_intersect(([0.,0.],[4.,4.]), ([1.,1.], [2.,2.]))
        self.assertEqual(len(isec), 2)
        isec=sorted(isec, key=lambda x: (x[0], x[1]))
        np.testing.assert_allclose(isec[0], [1., 1.])
        np.testing.assert_allclose(isec[1], [2., 2.])
        isec=ftuv.seg_intersect(([1.,1.], [2.,2.]), ([0.,0.],[4.,4.]))
        self.assertEqual(len(isec), 2)
        isec=sorted(isec, key=lambda x: (x[0], x[1]))
        np.testing.assert_allclose(isec[0], [1., 1.])
        np.testing.assert_allclose(isec[1], [2., 2.])
        #overlapping
        isec=ftuv.seg_intersect(([0.,2.], [2.,4.]), ([1.,3.],[3.,5.]))
        self.assertEqual(len(isec), 2)
        isec=sorted(isec, key=lambda x: (x[0], x[1]))
        np.testing.assert_allclose(isec[0], [1., 3.])
        np.testing.assert_allclose(isec[1], [2., 4.])
        #non-parallel, no intersection
        isec=ftuv.seg_intersect(([0.,2.], [2.,4.]), ([5.,3.],[10,5.]))
        self.assertEqual(isec, [])
        #shared endpoint
        isec=ftuv.seg_intersect(([0.,1.], [0., 4.]), ([0.,4.], [5.,7.]))
        self.assertEqual(len(isec), 1)
        np.testing.assert_allclose(isec[0], [0., 4.])
        isec=ftuv.seg_intersect(([0.,1.], [0., 4.]), ([0.,1.], [-5.,7.]))
        self.assertEqual(len(isec), 1)
        np.testing.assert_allclose(isec[0], [0., 1.])
        #Invalid inputs
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.,1.], [0., 4.]), ([0.,1.], [-5.,7., 5.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.,1., 3.], [0., 4.]), ([0.,1.], [-5.,7.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.,1.], [0., 4., 5.]), ([0.,1.], [-5.,7.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.,1.], [0., 4.]), ([0.,1., 7.], [-5.,7.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.,1.], [0., 4., 6.]), ([0.,1., 7.], [-5.,7.,8.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.], [0., 4.]), ([0.,1.], [-5.,7.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0., 5.], [4.34]), ([0.,1.], [-5.,7.]))
        with self.assertRaises(ValueError):
            ftuv.seg_intersect(([0.3, 5.2], [0.3, 5.2]), ([0.,1.], [-5.,7.]))
    def test_is_almost_colinear(self):
        self.assertTrue(ftuv.is_almost_colinear(np.array([3,6,7]),np.array([9.,18.,21.])))
        self.assertFalse(ftuv.is_almost_colinear(np.array([1,2,3]),np.array([2.,4.,-6.])))
        self.assertFalse(ftuv.is_almost_colinear(np.array([1,2,3]),np.array([3.,4.,6.])))
        self.assertFalse(ftuv.is_almost_colinear(np.array([1,2,3]),np.array([2.,5.,6.])))
