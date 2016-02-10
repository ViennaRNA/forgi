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

    
