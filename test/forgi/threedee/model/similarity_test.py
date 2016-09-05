import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.similarity as ftme
import forgi.utilities.debug as fud
import forgi.threedee.utilities.graph_pdb as ftug 
import unittest
import unittest, os, math
import numpy as np
import forgi.threedee.utilities.vector as ftuv

import forgi.utilities.debug as fud

class CompareTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_confusion_matrix(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = ftme.confusion_matrix(cg1, cg2)
        cm = ftme.confusion_matrix(cg2, cg2)

        self.assertEqual(cm['fp'], 0)
        self.assertEqual(cm['fn'], 0)
        #fud.pv('cm')
        #TODO assert something
        pass

    def test_mcc(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = ftme.confusion_matrix(cg1, cg2)
        mcc = ftme.mcc(cm)

        self.assertTrue(mcc < 1.0)

        cm = ftme.confusion_matrix(cg2, cg2)
        mcc = ftme.mcc(cm)

        self.assertLess(abs(mcc - 1.0), 0.01)

        pass

    def test_new_confusionmatrix_is_like_old(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = ftme.confusion_matrix(cg1, cg2)
        mcc = ftme.mcc(cm)
        cm_new = ftme.AdjacencyCorrelation(cg1) #previousely named ConfusionMatrix
        mcc_n = ftme.mcc(cm_new.evaluate(cg2))
        self.assertAlmostEqual(mcc, mcc_n)
        self.assertLess(abs(mcc_n-0.0761), 0.0001) #0.086 for distance=30, 0.0761 for 25


        cm = ftme.confusion_matrix(cg1, cg1)        
        mcc = ftme.mcc(cm)
        mcc_n = ftme.mcc(cm_new.evaluate(cg1))
        self.assertAlmostEqual(mcc, mcc_n)
        self.assertAlmostEqual(mcc_n, 1.0)

    def test_cg_rmsd(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')
        residues1 = ftug.bg_virtual_residues(cg1) 
        residues2 = ftug.bg_virtual_residues(cg2)
        self.assertAlmostEqual(ftme.rmsd(residues1, residues2), ftme.cg_rmsd(cg1,cg2))

    def test_cg_rmsd2(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')
        self.assertAlmostEqual(ftme.cg_rmsd(cg1,cg2), 25.170088934277373)


class TestRMSD(unittest.TestCase):
    '''
    Test some of the rmsd-type functions.
    '''
    def setUp(self):
        pass

    def test_rmsd(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a3 = np.array([[2., 2., 2.], [0., 0., 0.], [-2., -2., -2.]])
        self.assertAlmostEqual(ftme.rmsd(a1,a3), math.sqrt(2))

    @unittest.expectedFailure
    def test_rmsd_doesnt_center(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., -1., 0.], [-1., -1., -1.]])
        r1 = ftur.centered_rmsd(a1, a2)
        self.assertNotEqual(r1, ftur.rmsd(a1,a2))

    def test_drmsd(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., -1., 0.], [-1., -1., -1.]])
        r = ftme.drmsd(a1, a2)
        self.assertGreater(r, 0)

    def test_rmsd_to_same(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[1., 2., 1.], [0., -1., 0.], [-1., -1., -1.]])
        self.assertAlmostEqual(ftme.drmsd(a1, a1), 0)
        self.assertAlmostEqual(ftme.drmsd(a2, a2), 0)
        self.assertAlmostEqual(ftme.rmsd(a1, a1), 0)
        self.assertAlmostEqual(ftme.rmsd(a2, a2), 0)

    def test_rmsd_to_same_shifted(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[2., 2., 1.], [1., 1., 0.], [0., 0., -1.]])
        self.assertAlmostEqual(ftme.drmsd(a1, a2), 0)
        self.assertAlmostEqual(ftme.rmsd(a1, a2), 0)

    def test_rmsd_to_same_rotated(self):
        a1 = np.array([[1., 1., 1.], [0., 0., 0.], [-1., -1., -1.]])
        a2 = np.array([[-1., -1., -1.], [0., 0., 0.], [1., 1., 1.]])
        self.assertAlmostEqual(ftme.drmsd(a1, a2), 0)
        self.assertAlmostEqual(ftme.rmsd(a1, a2), 0)

    def test_rmsd_in_2D(self):
        a1 = np.array([[1., 1.], [0., 0.], [-1., -1.]])
        a2 = np.array([[2., 2.], [0., 0.], [-2., -2.]])
        a3 = np.array([[2., 4.], [0., 6.], [-1., 0.]])
        self.assertAlmostEqual(ftme.rmsd(a1, a1), 0)
        self.assertAlmostEqual(ftme.rmsd(a2, a2), 0)
        self.assertAlmostEqual(ftme.rmsd(a3, a3), 0)
        self.assertAlmostEqual(ftme.rmsd(a1, a2), math.sqrt(4./3.))
        self.assertGreater(ftme.rmsd(a1, a3), ftme.rmsd(a1, a2))
 
