import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug

import itertools as it 
import math
import numpy as np

__all__ = ['AdjacencyCorrelation', 'cg_rmsd', 'rmsd', 'drmsd']

"""
This module contains functions for the comparison of two cg objects or two ordered point-clouds.
"""

try:
    profile
except:
    def profile(f): 
        return f


def ppv(tp, fp):
    '''
    Calculate the positive predictive value.

    :param tp: The true positives
    :param fp: The false positives
    :return: The ppv
    '''
    return tp / float(tp + fp)

def sty(tp, fn):
    '''
    Calculate the sensitivity.

    :param tp: True positives
    :param fn: False negatives
    :return: The sensitivity
    '''
    return tp / float(tp + fn)

def mcc(confusion_matrix):
    '''
    Calculate the Matthews correlation coefficient for
    the given confusion matrix.

    `MCC = sqrt(ppv * sty)`

    :param confustion_matrix: A dictionary like this: `{"tp": tp, "tn": tn, "fp": fp, "fn": fn}`
    :return: The MCC
    '''

    return math.sqrt(ppv(confusion_matrix['tp'],
                         confusion_matrix['fp']) *
                     sty(confusion_matrix['tp'],
                         confusion_matrix['fn']))


class AdjacencyCorrelation(object):
    """
    A class used for calculating the ACC.
    
    The adjacency correlation coefficient is calculated as the Matthews correlation 
    coefficient of potential interactions, defined as nucleotides within 25 Angstrom
    from each other. See chapter 9.3 of Peter's thesis.

    This object is initialized with a reference structure and a distance for interactions.
    The evaluate() method is used for calculating this correlation matrix.

    This is significantly faster than the confusion_matrix function, if 
    many structures will be compared to the same reference structure.
    """
    def __init__(self, reference_cg, distance=25.0, bp_distance=16):
        self._distance=distance        
        self._bp_distance=bp_distance
        self._reference_interactions=self.get_interactions(reference_cg)

    @profile
    def get_interactions(self, cg):
        """
        :return: A set of 2-tuples containing elements that pair.
        """
        ignore = set()
        for n1, n2 in it.combinations(cg.defines.keys(), r=2):
            if cg.connected(n1, n2):
                ignore.add((n1,n2))
                continue
            bp_dist = cg.min_max_bp_distance(n1, n2)[0]
            if bp_dist < self._bp_distance:
                ignore.add((n1,n2))

        interactions=set(cg.coords.elements_closer_than(self._distance, ignore))
        return interactions
    
    def evaluate(self, cg):
        '''
        Calculate the true_positive, false_positive,
        true_negative and false_negative rate for the tertiary
        distances of the elements of cg and the reference structure stored in this class.

        :param cg: The first coarse grain model
        :return: A dictionary like this: `{"tp": tp, "tn": tn, "fp": fp, "fn": fn}`
        '''
        interactions=self.get_interactions(cg)
        nodes=set(cg.defines.keys())
        allIA=set()
        for n1, n2 in it.combinations(nodes, r=2):
            if cg.connected(n1, n2):
                continue
            bp_dist = cg.min_max_bp_distance(n1, n2)[0]
            if bp_dist < self._bp_distance:
                continue
            allIA.add(tuple(sorted((n1,n2))))

        d={ "tp":0, "tn":0, "fp":0, "fn":0 }
        d["tp"]=len(self._reference_interactions & interactions)
        d["fp"]=len(interactions - self._reference_interactions)
        d["fn"]=len(self._reference_interactions - interactions)
        d["tn"]=len(allIA - (self._reference_interactions | interactions) )
        return d

#NOTE: could be deprecated in the future. Use AdjacencyCorrelation.
#NOTE: could be moved to tests as reference for Adjacency correlation.
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

# COVERAGE NOTE: Never used in forgi or ernwin.
# This is left in the code as an example for the usage of confusiom_matrix,
# but might be removed in the future.
def mcc_between_cgs(cg_query, cg_native, distance=25, bp_distance=16):
    '''
    Calculate the MCC of the distance between two elements.

    :param cg_query: The second cg structure
    :param cg_native: The native cg structure
    :param distance: The distance between which we consider interactions.
    :param bp_distance: Only consider pairs of elements that are separated by 
                        MORE than this many base pairs
    :return: The MCC for interactions within a certain distance
    '''
    cm = confusion_matrix(cg_query, cg_native, distance, bp_distance=16)
    #cm['tp'] += 1
    #cm['fp'] += 1
    #cm['fn'] += 1
    #cm['tn'] += 1
    if cm['tp'] + cm['fp'] == 0:
        return None
    if cm['tp'] + cm['fn'] == 0:
        return None

    my_mcc = mcc(cm)
    return my_mcc


def optimal_superposition(crds1, crds2):
    """
    Returns best-fit rotation matrix as [3x3] numpy matrix for aligning crds1 onto crds2
    using the Kabsch algorithm
    """        
    assert(crds1.shape == crds2.shape)
    if crds1.shape[1] == 3 or crds1.shape[1] == 2:
        correlation_matrix = np.dot(np.transpose(crds1), crds2)
        v, s, w_tr = np.linalg.svd(correlation_matrix)
        is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
        if is_reflection:
            v[:, -1] = -v[:, -1]
        return np.dot(v, w_tr)
    else:
        raise ValueError("Wrong dimension of crds1. Needs to be an array of "
                         "Points in 2D or 3D space. Found {}D".format(crds1.shape[1]))
def cg_rmsd(cg1, cg2):
    '''
    Calculate the RMSD between two Coarse Grain models using their
    set of virtual residues.

    :param cg1: The first coarse grain model.
    :param cg2: The second coarse-grain model.
    :return: The RMSD
    '''

    residues1 = cg1.get_ordered_virtual_residue_poss()
    residues2 = cg2.get_ordered_virtual_residue_poss()

    return rmsd(residues1, residues2)

def rmsd_kabsch(crds1, crds2, is_centered=False):
    '''
    Center the coordinate vectors on their centroid
    and then calculate the rmsd.
    '''
    if not is_centered:
        crds1 = ftuv.center_on_centroid(crds1)
        crds2 = ftuv.center_on_centroid(crds2)

    os = optimal_superposition(crds1, crds2)
    crds_aligned = np.dot(crds1, os)

    diff_vecs = (crds2 - crds_aligned)
    vec_lengths = np.sum(diff_vecs * diff_vecs, axis=1)

    return math.sqrt(sum(vec_lengths) / len(vec_lengths))

def drmsd(coords1, coords2):
    '''
    Calculate the dRMSD measure.

    This should be the RMSD between all of the inter-atom distances
    in two structures.

    :param coords1: The vectors of the 'atoms' in the first structure.
    :param coords2: The vectors of the 'atoms' in the second structure.
    :return: The dRMSD measure.
    '''
    ds1 = np.array([ftuv.vec_distance(c1, c2) for c1,c2 in it.combinations(coords1, r=2)])
    ds2 = np.array([ftuv.vec_distance(c1, c2) for c1,c2 in it.combinations(coords2, r=2)])

    rmsd = math.sqrt(np.mean((ds1 - ds2) * (ds1 - ds2)))
    #rmsd = math.sqrt(np.mean(ftuv.vec_distance(ds1, ds2)))

    return rmsd

#The function rmsd points to the faster QC version if available, else to our kabsch implementation.
#The name rmsd_kabsch is always available to refer to our kabsch implementation.
try:
    from py_qcprot import rmsd as rmsd_qc
    rmsd = rmsd_qc
except:
    rmsd = rmsd_kabsch

