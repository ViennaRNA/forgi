import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.comparison as ftme
import forgi.threedee.utilities.rmsd as ftur
import forgi.utilities.debug as fud
import forgi.threedee.utilities.graph_pdb as ftug 
import unittest

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
        self.assertAlmostEqual(ftur.centered_rmsd(residues1, residues2), ftme.cg_rmsd(cg1,cg2))

    def test_cg_rmsd2(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')
        self.assertAlmostEqual(ftme.cg_rmsd(cg1,cg2), 25.170088934277373)
