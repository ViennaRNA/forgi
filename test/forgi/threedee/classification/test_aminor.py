import unittest
import logging

from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix


import forgi.threedee.classification.aminor as ftca
import forgi.threedee.model.coarse_grain as ftmc


log=logging.getLogger(__name__)

class TestJustClassifyFunctions(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/3pdr_X.cg")
    def test_all_interactions(self):
        interactions=ftca.all_interactions(self.cg)
        self.assertGreaterEqual(len(interactions), 5)
        self.assertLessEqual(len(interactions), 25) # We have 23 annotated Fr3d interactions, but some may be directly connected.
        self.assertEqual(interactions.shape[1], 2)
        self.assertEqual(interactions[0,1][0], "s") # second element is stem
    def test_classify_interaction(self):
        pass

class TestClassifier(unittest.TestCase):
    def test_via_crossval(self):
        specs=[]
        sens=[]
        loops="imh"
        for loop in loops:
            data, labels = ftca.get_trainings_data(loop)
            X_train, X_test, y_train, y_test = train_test_split(
                data, labels, test_size=0.5, train_size=0.5, stratify=labels)
            clf = ftca.AMinorClassifier()
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
