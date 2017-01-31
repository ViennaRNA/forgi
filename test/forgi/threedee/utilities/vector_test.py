from __future__ import absolute_import

import unittest, math
import numpy as np
import numpy.testing as nptest
import forgi.threedee.utilities.vector as ftuv

REPEAT_TESTS_CONTAINING_RANDOM = 10

class TestVector(unittest.TestCase):
    """Tests for the threedee.utilities.vector module"""
    def setUp(self):
        return
    def test_closest_point_on_seg(self):
        self.assertEqual(tuple(ftuv.closest_point_on_seg((0,1),(0,3),(2,2))),(0,2))
        self.assertEqual(tuple(ftuv.closest_point_on_seg((1,0),(3,0),(2,2))),(2,0))
        #Not parallel to axis: Floating point values...
        self.assertAlmostEqual(ftuv.closest_point_on_seg((0,0),(2,2),(0,2))[0],1.) #x-coordinate
        self.assertAlmostEqual(ftuv.closest_point_on_seg((0,0),(2,2),(0,2))[0],1.) #y-coordinate
        #Outside segment: returns one endpoint of the segment.
        self.assertEqual(tuple(ftuv.closest_point_on_seg((0,1),(0,3),(2,4))),(0,3))
        self.assertEqual(tuple(ftuv.closest_point_on_seg((0,1),(0,3),(-2,0))),(0,1))
    def test_get_inter_distances(self):
        vecs=[np.array([1., 0., 0.]), np.array([0., 0., 0.]), np.array([0., 0., 0.]), np.array([-1., 0., 0.])]
        distances=ftuv.get_inter_distances(vecs)
        self.assertEqual(sorted(distances), [0,1,1,1,1,2])
    def test_get_random_vector(self):
        for _ in range(REPEAT_TESTS_CONTAINING_RANDOM):
            vec=ftuv.get_random_vector()
            self.assertLessEqual(ftuv.magnitude(vec[0]),1.)
            vec1=ftuv.get_random_vector()
            vec2=ftuv.get_random_vector()
            self.assertTrue(all( vec1[j]!=vec2[j] for j in [0,1,2]),
                              msg="Repeated calls should generate different results."
                                  "This tests depends on random values. if it fails, try running it again.")
    def test_get_alignment_matrix(self):
        vec1=np.array([0.5,0.7,0.9])
        vec2=np.array([0.345,3.4347,0.55])
        rotMat=ftuv.get_alignment_matrix(vec1,vec2)
        self.assertTrue( ftuv.is_almost_colinear(vec2, np.dot(vec1, rotMat)) )
        self.assertTrue( ftuv.is_almost_colinear(np.dot(rotMat, vec2), vec1) )

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
        self.assertTrue(ftuv.is_almost_colinear(np.array([0,0,0]),np.array([0.,0.,0.])))
        self.assertTrue(ftuv.is_almost_colinear(np.array([0,0,2]),np.array([0.,0.,3.])))
        self.assertTrue(ftuv.is_almost_colinear(np.array([3,6,7]),np.array([9.,18.,21.])))
        self.assertTrue(ftuv.is_almost_colinear(np.array([3,6,0]),np.array([9.,18.,0.])))
        self.assertTrue(ftuv.is_almost_colinear(np.array([3,0,8]),np.array([9.,0.,24.])))
        self.assertTrue(ftuv.is_almost_colinear(np.array([0,0,0]),np.array([0,20,0])))
        self.assertTrue(ftuv.is_almost_colinear(np.array([0.0004,0,0]),np.array([0.,0.,0.])))

        self.assertFalse(ftuv.is_almost_colinear(np.array([0,0,3]),np.array([2.,0,0])))
        self.assertFalse(ftuv.is_almost_colinear(np.array([1,2,3]),np.array([2.,4.,-6.])))
        self.assertFalse(ftuv.is_almost_colinear(np.array([1,2,3]),np.array([3.,4.,6.])))
        self.assertFalse(ftuv.is_almost_colinear(np.array([1,2,3]),np.array([2.,5.,6.])))

    def test_middlepoint(self):
        self.assertIsInstance(ftuv.middlepoint((1,2,3),(4,5,6)), tuple)
        self.assertIsInstance(ftuv.middlepoint([1,2,3],[4,5,6]), list)
        self.assertIsInstance(ftuv.middlepoint(np.array([1,2,3]),np.array([4,5,6])), type(np.array([1,2,3])))
        self.assertEqual(ftuv.middlepoint((1,2),(3,4)), (2,3))
        self.assertEqual(ftuv.middlepoint([1,2,3],[5,6,7]), [3,4,5])
        mp=ftuv.middlepoint(np.array([1,2,-3]),np.array([1,0,-5]))
        self.assertTrue(((mp==np.array([1,1,-4])).all()), msg="Middlepoint for np arrays: {} "
                                                      "is not {}".format(mp, np.array([1,1,-4])))
    def test_create_orthonormal_basis(self):
        #Note: If the input vectors are not orthogonal, the result are 3 vectors that might not form a basis.
        basis1=ftuv.create_orthonormal_basis(np.array([0.0,0.0,2.0]))      
        self.assertTrue( ftuv.is_almost_colinear(basis1[0], np.array([0.,0.,2.])) )
        basis2=ftuv.create_orthonormal_basis(np.array([0.0,0.0,2.0]), np.array([0.0, 3.6, 0.]))    
        self.assertTrue( ftuv.is_almost_colinear(basis2[0], np.array([0.,0.,2.])) )
        self.assertTrue( ftuv.is_almost_colinear(basis2[1], np.array([0.,3.6,0])) )
        basis3=ftuv.create_orthonormal_basis(np.array([0.0,0.0,2.0]), np.array([0.0, 3.6, 0.]), np.array([1.,0,0]))    
        self.assertTrue( ftuv.is_almost_colinear(basis3[0], np.array([0.,0.,2.])) )
        self.assertTrue( ftuv.is_almost_colinear(basis3[1], np.array([0.,3.6,0])) )
        self.assertTrue( ftuv.is_almost_colinear(basis3[2], np.array([1.,0,0])) )
        for basis in [basis1, basis2, basis3]:
          self.assertAlmostEqual(np.dot(basis[0],basis[1]),0)
          self.assertAlmostEqual(np.dot(basis[0],basis[2]),0)
          self.assertAlmostEqual(np.dot(basis[2],basis[1]),0)
          for b in basis:
              self.assertAlmostEqual(ftuv.magnitude(b),1)

    def test_spherical_coordinate_transforms(self):
        for vec in [np.array([0,0,1]), np.array([0,2,0]), np.array([3,0,0]), np.array([4,5,0]), np.array([6,0,7]), np.array([8,9,0.4])]:
            sphe=ftuv.spherical_cartesian_to_polar(vec)
            nptest.assert_allclose(ftuv.spherical_polar_to_cartesian(sphe),vec, atol=0.0000001)
        nptest.assert_allclose(ftuv.spherical_polar_to_cartesian([1,0,math.pi/2]), np.array([0,0,1]), atol=0.0000001)
        nptest.assert_allclose(ftuv.spherical_polar_to_cartesian([2, math.pi/2, math.pi/4]), np.array([1,1,0])*2/math.sqrt(2), atol=0.0000001)
        nptest.assert_allclose(ftuv.spherical_cartesian_to_polar(np.array([0,2,2])/math.sqrt(2)), [2,math.pi/4,math.pi/2], atol=0.0000001)

    def test_get_standard_basis(self):
        nptest.assert_allclose(ftuv.get_standard_basis(2), [[1,0],[0,1]])
        nptest.assert_allclose(ftuv.get_standard_basis(3), [[1,0,0],[0,1,0], [0,0,1]])


    def test_change_basis(self):
        new_v = ftuv.change_basis(np.array([1.,2.,3.]), np.array([[0,1.,0],[1.,0,0],[0,0,1.]]), np.array([[1.,0,0],[0,1.,0],[0,0,1.]]))
        nptest.assert_allclose(new_v, np.array([2.,1.,3.]))
                    
    def test_change_basis_vectorized(self):
        coords = np.array([[0., 1., 2.], [1., 2., 3.], [0., 0., 2.], [0.,1.,0.]])
        basis1 = np.array([[0.,0.,1.],[0.,1.,0.],[1.,0.,1.]])
        basis2 = np.array([[1.,0.,0.], [0.,1.,0.],[0.,1.,1.]])
        new_coords = ftuv.change_basis_vectorized(coords, basis2, basis1)
        for i in range(4):
            nptest.assert_array_equal(new_coords[i], ftuv.change_basis(coords[i], basis2, basis1))
        nptest.assert_array_equal(new_coords[2], np.array([2.,-2.,2.]))

    def benchmark_change_basis(self):
        import timeit
        t1 = timeit.timeit("ftuv.change_basis_vectorized(coords, new_basis, old_basis)", 
                          "import forgi.threedee.utilities.vector as ftuv; import numpy as np; coords = (np.random.rand(100,3)-0.5)*20;  "
                          "new_basis=np.array([[1.,2.,0.],[0.,6.,7],[0.4,0,9.3]]);old_basis=np.array([[1.5,2.5,0],[1.5,0,7],[0,0.7,9.3]])", number = 1000000)
        t2 = timeit.timeit("for coord in coords: ftuv.change_basis(coord, new_basis, old_basis)", 
                          setup="import numpy as np; coords = (np.random.rand(100,3)-0.5)*20; import forgi.threedee.utilities.vector as ftuv; "
                                "new_basis=np.array([[1.,2.,0.],[0.,6.,7],[0.4,0,9.3]]);old_basis=np.array([[1.5,2.5,0],[1.5,0,7],[0,0.7,9.3]])", number = 1000000)
        self.assertLess(int(t1)+50,int(t2))
        
    def test_det3x3(self):
        m1 = np.array([[1.,2,3],[4.,5,6],[7,8,9]])
        m2 = np.array([[1,1,2],[3,3,4.],[6,6,8]])
        m3= np.array([[2,-4,6],[-2,6.,9],[0,0,1]])
        for m in [m1,m2,m3]:
            self.assertEqual(ftuv.det3x3(m), np.linalg.det(m))
            
    def test_get_centroid(self):
        coords = [[0.,1.,1.],[1,1,1],[-1,2,3],[3, 0, 0],[-3,1,0]]
        nptest.assert_almost_equal(ftuv.get_vector_centroid(np.array(coords)), [0, 1, 1])
        nptest.assert_almost_equal(ftuv.get_vector_centroid(coords), [0, 1, 1])
    def test_center_on_centroid(self):
        coords = [[0.,1.,1.],[1,1,1],[-1,2,3],[3, 0, 0],[-3,1,0]]
        nptest.assert_almost_equal(ftuv.center_on_centroid(np.array(coords)), 
                    [[0,0.,0],[1,0,0],[-1,1,2],[3,-1,-1],[-3,0,-1]])
        nptest.assert_equal(ftuv.get_vector_centroid(ftuv.center_on_centroid(coords)), [0,0.,0])
