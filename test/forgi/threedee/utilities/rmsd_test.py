import unittest, os, math
import numpy as np
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.utilities.vector as ftuv

import forgi.utilities.debug as fud

class TestRMSD(unittest.TestCase):
    '''
    Test some of the rmsd-type functions.
    '''
    def setUp(self):
        pass

    def test_rmsd(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., 1., 0.], [-1., -1., -1.]])
        a3 = np.array([[2., 2., 2.], [0., 0., 0.], [-2., -2., -2.]])
        c1 = ftuv.center_on_centroid(a1)
        c2 = ftuv.center_on_centroid(a2)

        r1 = ftur.centered_rmsd(a1, a2)
        r2 = ftur.rmsd(c1,c2)

        self.assertAlmostEqual(r1,r2)

        self.assertAlmostEqual(ftur.rmsd(a1,a3), math.sqrt(2))
    """
    #Currently fails:
    def test_rmsd_doesnt_center(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., -1., 0.], [-1., -1., -1.]])
        r1 = ftur.centered_rmsd(a1, a2)
        self.assertNotEqual(r1, ftur.rmsd(a1,a2))
    """
    def test_drmsd(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., -1., 0.], [-1., -1., -1.]])


        r = ftur.drmsd(a1, a2)

        self.assertGreater(r, 0)


    def test_radius_of_gyration(self):
        a = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        r = ftur.radius_of_gyration(a)

        self.assertGreater(r, 0)

    def test_rmsd_to_same(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., -1., 0.], [-1., -1., -1.]])
        self.assertAlmostEqual(ftur.drmsd(a1, a1), 0)
        self.assertAlmostEqual(ftur.drmsd(a2, a2), 0)
        self.assertAlmostEqual(ftur.rmsd(a1, a1), 0)
        self.assertAlmostEqual(ftur.rmsd(a2, a2), 0)

    def test_rmsd_in_2D(self):
        a1 = np.array([[1., 1.], [0., 0.], [-1., -1.]])
        a2 = np.array([[2., 2.], [0., 0.], [-2., -2.]])
        a3 = np.array([[2., 4.], [0., 6.], [-1., 0.]])
        self.assertAlmostEqual(ftur.centered_rmsd(a1, a1), 0)
        self.assertAlmostEqual(ftur.centered_rmsd(a2, a2), 0)
        self.assertAlmostEqual(ftur.centered_rmsd(a3, a3), 0)
        self.assertAlmostEqual(ftur.centered_rmsd(a1, a2), math.sqrt(4./3.))
        self.assertGreater(ftur.centered_rmsd(a1, a3), ftur.centered_rmsd(a1, a2))

    def test_gyration_tensor_vs_rog(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.], [-2,-2,-2]])
        rog = ftur.radius_of_gyration(a1)
        g_tensor = ftur.gyration_tensor(a1, diagonalize=True)
        print(rog, g_tensor)
        #the first invariant is the rog
        self.assertAlmostEqual(g_tensor[0,0]+g_tensor[1,1]+g_tensor[2,2], rog**2)

    def test_gyration_tensor_diagonalize(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.], [-2,-2,-2]])
        a2 = np.array([[1., 2., 1.], [0., 0., 0.], [-1., -1., -2.], [-12,-2,-25]])
        for a in [a1,a2]:
            g_tensor = ftur.gyration_tensor(a1, diagonalize=True)
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
        self.assertAlmostEqual(ftur.anisotropy(linear), 1)       
    def test_anisotropy_star(self): 
        star = np.array([[0.,0,1], [0,1,0], [1,0,0], [-1,0,0], [0,-1,0], [0,0,1], [0,0,0]])
        self.assertLessEqual(ftur.anisotropy(star), 0.2)
    def test_anisotropy_planar(self):
        # The anisotropy for planar symmetric objects converges to 1/4.
        # See references 36-40 in doi:10.1063/1.4788616
        planar = np.array([[0.,0,1], [0,1,0], [0,-1,0], [0,0,1], [0,0,0], 
                           [0,0.6,0.6], [0,-0.6,0.6],[0,0.6,-0.6],[0,-0.6,-0.6]])
        self.assertAlmostEqual(ftur.anisotropy(planar), 0.25, places=2)
        

