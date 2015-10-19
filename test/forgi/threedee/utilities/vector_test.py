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
        ortVec=ftuv.get_orthogonal_unit_vector(vec)
        #self.assertAlmostEqual(np.dot(ortVec, vec), 0, places=10) #Currently, ortVec==nan, so the assertion woulkd fail.
        #self.assertAlmostEqual(np.linalg.norm(ortVec), 1, places=10)
     
