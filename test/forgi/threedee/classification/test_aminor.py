from __future__ import division

import unittest
import logging
import math
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix


import forgi.threedee.classification.aminor as ftca
import forgi.threedee.model.coarse_grain as ftmc


log=logging.getLogger(__name__)

class TestRelativeOrientation(unittest.TestCase):
    def test_get_relative_orientation(self):
        cg = ftmc.CoarseGrainRNA.from_dotbracket("...(((...)))...", seq="AAAGGGAAACCCAAA")
        cg.coords["s0"] = [0,0,0.], [0,0,10.]
        cg.twists["s0"] = [0,1,0.], [0, -1., 0]
        cg.coords["f0"] = [5.,0,0], [10.,0,0]
        cg.coords["t0"] = [5.,0,5], [10., 0, 5]
        cg.coords["h0"] = [5.,0,10],[10.,0,10]
        d, a1, a2 = ftca.get_relative_orientation(cg, "f0", "s0")
        self.assertAlmostEqual(d, 5)
        self.assertAlmostEqual(a1, math.pi/2)
        self.assertAlmostEqual(a2, -math.pi/2)
        d, a1, a2 = ftca.get_relative_orientation(cg, "h0", "s0")
        self.assertAlmostEqual(d, 5)
        self.assertAlmostEqual(a1, math.pi/2)
        self.assertAlmostEqual(a2, math.pi/2)
        d, a1, a2 = ftca.get_relative_orientation(cg, "t0", "s0")
        self.assertAlmostEqual(d, 5)
        self.assertAlmostEqual(a1, math.pi/2)
        self.assertIn(a2, [math.pi, -math.pi])
        # Now take an example with a non 90 degrees angle a1
        cg.coords["h0"] = [5.,0,15],[10.,0.,15]
        d, a1, a2 = ftca.get_relative_orientation(cg, "h0", "s0")
        self.assertAlmostEqual(d, math.sqrt(50))
        self.assertAlmostEqual(a1, math.pi/4)
        self.assertAlmostEqual(a2, math.pi/2)



class TestJustClassifyFunctions(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/3pdr_X.cg")
    def test_all_interactions(self):
        interactions=ftca.all_interactions(self.cg)
        print(interactions)
        self.assertGreaterEqual(len(interactions), 4)
        self.assertLessEqual(len(interactions), 25) # We have 23 annotated Fr3d interactions, but some may be directly connected.
        self.assertEqual(len(interactions[0]), 2)
        self.assertEqual(interactions[0][1][0], "s") # second element is stem
        self.assertEqual(len(interactions), len(set([e[0] for e in interactions]))) # Loops are unique

class TestClassifier(unittest.TestCase):
    def test_via_crossval_slowtest(self):
        specs=[]
        sens=[]
        loops="ih"
        for loop in loops:
            data, labels = ftca.get_trainings_data(loop)
            X_train, X_test, y_train, y_test = train_test_split(
                data, labels, test_size=0.5, train_size=0.5, stratify=labels)
            clf = ftca._get_default_clf(loop)
            clf.fit(X_train, y_train)
            #X_test[:,0]*=ftca.ANGLEWEIGHT
            y_true, y_pred = y_test, clf.predict(X_test)
            tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
            specificity =  tn/(tn+fp)
            sensitivity = tp/(tp+fn)
            specs.append(specificity)
            sens.append(sensitivity)
        print("Specificities", specs)
        print("Sensitivities", sens)
        for i,s in enumerate(specs):
            assert s>0.5, "Specificity({})={} <0.7".format(loops[i], specificity)
        for i,s in enumerate(sens):
            assert sensitivity>0.5, "Sensitivity({})={} <0.7".format(loops[i], sensitivity)
