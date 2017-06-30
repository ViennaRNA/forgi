from __future__ import print_function
import unittest, os, math
import numpy as np
import forgi.threedee.model.descriptors as ftmd
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

class TestGyration(unittest.TestCase):
    def setUp(self):
        pass

    def test_radius_of_gyration(self):
        a = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        r = ftmd.radius_of_gyration(a)
        self.assertGreater(r, 0)
        self.assertAlmostEqual(r, math.sqrt(2))

    def test_gyration_tensor_vs_rog(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.], [-2,-2,-2]])
        rog = ftmd.radius_of_gyration(a1)
        g_tensor = ftmd.gyration_tensor(a1, diagonalize=True)
        print(rog, g_tensor)
        #the first invariant is the rog
        self.assertAlmostEqual(g_tensor[0,0]+g_tensor[1,1]+g_tensor[2,2], rog**2)

    def test_gyration_tensor_diagonalize(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.], [-2,-2,-2]])
        a2 = np.array([[1., 2., 1.], [0., 0., 0.], [-1., -1., -2.], [-12,-2,-25]])
        for a in [a1,a2]:
            g_tensor = ftmd.gyration_tensor(a1, diagonalize=True)
            self.assertEqual(g_tensor[0,1],0)
            self.assertEqual(g_tensor[0,2],0)
            self.assertEqual(g_tensor[1,2],0)
            self.assertEqual(g_tensor[1,0],0)
            self.assertEqual(g_tensor[2,0],0)
            self.assertEqual(g_tensor[2,1],0)
            self.assertGreaterEqual(g_tensor[0,0],g_tensor[1,1])
            self.assertGreaterEqual(g_tensor[1,1],g_tensor[2,2])

    def test_anisotropy_linear(self):
        linear = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.], [-2,-2,-2]])
        self.assertAlmostEqual(ftmd.anisotropy(linear), 1)       
    def test_anisotropy_star(self): 
        star = np.array([[0.,0,1], [0,1,0], [1,0,0], [-1,0,0], [0,-1,0], [0,0,1], [0,0,0]])
        self.assertLessEqual(ftmd.anisotropy(star), 0.2)
    def test_anisotropy_planar(self):
        # The anisotropy for planar symmetric objects converges to 1/4.
        # See references 36-40 in doi:10.1063/1.4788616
        planar = np.array([[0.,0,1], [0,1,0], [0,-1,0], [0,0,1], [0,0,0], 
                           [0,0.6,0.6], [0,-0.6,0.6],[0,0.6,-0.6],[0,-0.6,-0.6]])
        self.assertAlmostEqual(ftmd.anisotropy(planar), 0.25, places=2)
        

