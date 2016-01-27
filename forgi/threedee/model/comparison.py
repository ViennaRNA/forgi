import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.rmsd as ftur

import itertools as it 
import math

def ppv(tp, fp):
    '''
    Calculate the positive predictive value.

    @param tp: The true positives
    @param fp: The false positives
    @return: The ppv
    '''
    return tp / float(tp + fp)

def sty(tp, fn):
    '''
    Calculate the sensitivity.

    @param tp: True positives
    @param fn: False negatives
    @return: The sensitivity
    '''
    return tp / float(tp + fn)

def mcc(confusion_matrix):
    '''
    Calculate the Matthews correlation coefficient for
    the given confusion matrix.

    MCC = sqrt(ppv * sty)

    @param confustion_matrix: A dictionary like this: {"tp": tp, "tn": tn, "fp": fp, "fn": fn}
    @return: The MCC
    '''

    return math.sqrt(ppv(confusion_matrix['tp'],
                         confusion_matrix['fp']) *
                     sty(confusion_matrix['tp'],
                         confusion_matrix['fn']))


class ConfusionMatrix(object):
    """
    A class used for calculating the confusion_matrix.

    It is initialized with a reference structure and a distance for interactions.
    The evaluate() method is used for calculating this correlation matrix.
    """
    def __init__(self, reference_cg, distance=30.0):
        self._distance=distance
        self._reference_interactions=self.get_interactions(reference_cg)

    def get_interactions(self, cg):
        """
        :return: A set of 2-tuples containing elements that pair.
        """
        interactions=set()
        nodes=set(cg.defines.keys())
        for n1, n2 in it.combinations(nodes, r=2):
            n1,n2=sorted((n1,n2)) #TODO: Read itertools documentation, if this is necessary.
            if cg.connected(n1, n2):
                continue
            if ftuv.elements_closer_then(cg.coords[n1][0],
                                       cg.coords[n1][1],
                                       cg.coords[n2][0],
                                       cg.coords[n2][1], self._distance):
            #dist = cg.element_physical_distance(n1, n2)
            #if dist < self._distance:
                interactions.add((n1,n2))
        return interactions
    
    def evaluate(self, cg):
        '''
        Calculate the true_positive, false_positive,
        true_negative and false_negative rate for the tertiary
        distances of the elements of cg and the reference structure stored in this class.

        @param cg: The first coarse grain model
        @return: A dictionary like this: {"tp": tp, "tn": tn, "fp": fp, "fn": fn}
        '''
        interactions=self.get_interactions(cg)
        nodes=set(cg.defines.keys())
        allIA=set()
        for n1, n2 in it.combinations(nodes, r=2):
            allIA.add(tuple(sorted((n1,n2))))

        d={ "tp":0, "tn":0, "fp":0, "fn":0 }
        d["tp"]=len(self._reference_interactions & interactions)
        d["fp"]=len(interactions - self._reference_interactions)
        d["fn"]=len(self._reference_interactions - interactions)
        d["tn"]=len(allIA - (self._reference_interactions | interactions) )
        return d

def confusion_matrix(cg1, cg2, distance=30, bp_distance=16):
    '''
    Calculate the true_positive, false_positive,
    true_negative and false_negative rate for the tertiary
    distances of the elements of two structures, cg1 and cg2.

    @param cg1: The first coarse grain model
    @param cg2: The second coarse grain model
    @param distance: The distance to consider for interactions
    @param bp_distance: Only consider elements separated by this many more pair
        and backbone bonds
    @return: A dictionary like this: {"tp": tp, "tn": tn, "fp": fp, "fn": fn}
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

def mcc_between_cgs(cg_query, cg_native, distance=25, bp_distance=16):
    '''
    Calculate the MCC of the distance between two elements.

    @param cg_query: The second cg structure
    @param cg_native: The native cg structure
    @param distance: The distance between which we consider interactions.
    @param bp_distance: Only consider pairs of elements that are separated by no
        more than this many base pairs
    @return: The MCC for interactions within a certain distance
    '''
    cm = confusion_matrix(cg_query, cg_native, distance, bp_distance=16)
    cm['tp'] += 1
    cm['fp'] += 1
    cm['fn'] += 1
    cm['tn'] += 1
    if cm['tp'] + cm['fp'] == 0:
        return None
    if cm['tp'] + cm['fn'] == 0:
        return None

    my_mcc = mcc(cm)
    return my_mcc

def cg_rmsd(cg1, cg2):
    '''
    Calculate the RMSD between two Coarse Grain models using their
    set of virtual residues.

    @param cg1: The first coarse grain model.
    @param cg2: The second coarse-grain model.
    @return: The RMSD
    '''

    residues1 = ftug.bg_virtual_residues(cg1)
    residues2 = ftug.bg_virtual_residues(cg2)

    return ftur.centered_rmsd(residues1, residues2)

