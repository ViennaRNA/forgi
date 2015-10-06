import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud

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

def confusion_matrix(cg1, cg2, distance=30):
    '''
    Calculate the true_positive, false_positive,
    true_negative and false_negative rate for the tertiary
    distances of the elements of two structures, cg1 and cg2.

    @param cg1: The first coarse grain model
    @param cg2: The second coarse grain model
    @param distance: The distance to consider for interactions
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

def mcc_between_cgs(cg_query, cg_native, distance=30):
    '''
    Calculate the MCC of the distance between two elements.

    @param cg_query: The second cg structure
    @param cg_native: The native cg structure
    @return: The MCC for interactions within a certain distance
    '''
    cm = confusion_matrix(cg_query, cg_native, distance)
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

