from __future__ import division
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.similarity as ftme
import forgi.utilities.debug as fud
import forgi.threedee.utilities.graph_pdb as ftug
import unittest
import unittest, os, math
import numpy as np
import forgi.threedee.utilities.vector as ftuv

import forgi.utilities.debug as fud

import itertools as it


class OldConfisionMatrixTest(unittest.TestCase):
    def test_confusion_matrix(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = confusion_matrix(cg1, cg2)
        cm = confusion_matrix(cg2, cg2)

        self.assertEqual(cm['fp'], 0)
        self.assertEqual(cm['fn'], 0)


class CompareTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_ppv(self):
        self.assertAlmostEqual(ftme.ppv(1,1), 0.5)
        self.assertAlmostEqual(ftme.ppv(1,0), 1)
        self.assertAlmostEqual(ftme.ppv(0,1), 0)
        self.assertTrue(np.isnan(ftme.ppv(0,0)))

    def test_mcc(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = confusion_matrix(cg1, cg2)
        mcc = ftme.mcc(cm)

        self.assertTrue(mcc < 1.0)

        cm = confusion_matrix(cg2, cg2)
        mcc = ftme.mcc(cm)

        self.assertLess(abs(mcc - 1.0), 0.01)

        pass

    def test_new_confusionmatrix_is_like_old(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')

        cm = confusion_matrix(cg1, cg2)
        mcc = ftme.mcc(cm)
        cm_new = ftme.AdjacencyCorrelation(cg1) #previousely named ConfusionMatrix
        mcc_n = ftme.mcc(cm_new.evaluate(cg2))
        self.assertAlmostEqual(mcc, mcc_n)
        self.assertAlmostEqual(mcc_n, 0.6756639246921762)


        cm = confusion_matrix(cg1, cg1)
        mcc = ftme.mcc(cm)
        mcc_n = ftme.mcc(cm_new.evaluate(cg1))
        self.assertAlmostEqual(mcc, mcc_n)
        self.assertAlmostEqual(mcc_n, 1.0)

    def test_cg_rmsd(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')
        residues1 = ftug.bg_virtual_residues(cg1)
        residues2 = ftug.bg_virtual_residues(cg2)
        self.assertAlmostEqual(ftme.rmsd(residues1, residues2), ftme.cg_rmsd(cg1,cg2))

    def test_cg_rmsd2(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        self.assertAlmostEqual(ftme.cg_rmsd(cg1,cg1), 0)
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A_sampled.cg')
        self.assertAlmostEqual(ftme.cg_rmsd(cg1,cg2), 7.684377397812648)

    def test_cg_rmsd3(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2 = ftmc.CoarseGrainRNA.from_bg_file('test/forgi/threedee/data/1GID_A.cg')
        cg2.coords.rotate(ftuv.rotation_matrix([1.,2.,3.], 22))
        cg2.twists.rotate(ftuv.rotation_matrix([1.,2.,3.], 22))

        self.assertAlmostEqual(ftme.cg_rmsd(cg1,cg2), 0)

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

    @unittest.expectedFailure # Currently it DOES center!
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




##############################################################################
##  The old confusion_matrix code by pkerpedjiev. Used as reference in tests
##############################################################################

def confusion_matrix(cg1, cg2, distance=25, bp_distance=16):
    '''
    Calculate the true_positive, false_positive,
    true_negative and false_negative rate for the tertiary
    distances of the elements of two structures, cg1 and cg2.

    :param cg1: The first coarse grain model
    :param cg2: The second coarse grain model
    :param distance: The distance to consider for interactions
    :param bp_distance: Only consider elements separated by this many more pair
                        and backbone bonds
    :return: A dictionary like this: `{"tp": tp, "tn": tn, "fp": fp, "fn": fn}`
    '''
    nodes1 = set(cg1.defines.keys())
    nodes2 = set(cg2.defines.keys())

    tp = 0 #true positive
    tn = 0 #true negative

    fp = 0 #false positive
    fn = 0 #false negative

    #fud.pv('nodes1')
    #fud.pv('nodes2')

    assert(nodes1 == nodes2)

    for n1, n2 in it.combinations(nodes1, r=2):
        if cg1.connected(n1, n2):
            continue

        if cg2.connected(n1, n2):
            raise Exception("{} {} connected in cg2 but not cg1".format(n1, n2))


        bp_dist = cg2.min_max_bp_distance(n1, n2)[0]
        if bp_dist < bp_distance:
            continue

        dist1 = cg1.element_physical_distance(n1, n2)
        dist2 = cg2.element_physical_distance(n1, n2)

        if dist1 < distance:
            # positive
            if dist2 < distance:
                #true positive
                tp += 1
            else:
                # false negative
                fn += 1
        else:
            #negative
            if dist2 < distance:
                # false positive
                fp += 1
            else:
                # true negative
                tn += 1


    return {"tp": tp, "tn": tn, "fp": fp, "fn": fn}
