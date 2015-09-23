import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.comparison as ftme
import forgi.utilities.debug as fud
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
        fud.pv('cm')

        pass

    def test_mcc(self):
        cg1 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = ftme.confusion_matrix(cg1, cg2)
        mcc = ftme.mcc(cm)

        self.assertTrue(mcc < 1.0)

        cm = ftme.confusion_matrix(cg2, cg2)
        mcc = ftme.mcc(cm)

        self.assertTrue(abs(mcc - 1.0) < 0.01)

        pass
