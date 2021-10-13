from __future__ import division
from builtins import object

import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug

import itertools as it
import logging
import math
import numpy as np
from collections import defaultdict


log = logging.getLogger(__name__)
__all__ = ['AdjacencyCorrelation', 'cg_rmsd', 'rmsd', 'drmsd']

"""
This module contains functions for the comparison of two cg objects or two ordered point-clouds.
"""

try:
    profile
except:
    def profile(f):
        return f


class Incompareable(ValueError):
    """
    Raised if two objects are compared, which are incompareable.

    E.g. Coordinate sets of different length are incompareable with respect to the RMSD.
    """


def ppv(tp, fp):
    '''
    Calculate the positive predictive value.

    :param tp: The true positives
    :param fp: The false positives
    :return: The ppv
    '''
    try:
        return tp / float(tp + fp)
    except ZeroDivisionError:
        return float("nan")


def sty(tp, fn):
    '''
    Calculate the sensitivity.

    :param tp: True positives
    :param fn: False negatives
    :return: The sensitivity
    '''
    try:
        return tp / float(tp + fn)
    except ZeroDivisionError:
        return float("nan")


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
        self._distance = distance
        self._bp_distance = bp_distance
        self._reference_interactions = self._get_interactions(reference_cg)

    @profile
    def _get_interactions(self, cg):
        """
        :return: A set of 2-tuples containing elements that pair.
        """
        ignore = set()
        for n1, n2 in it.combinations(cg.defines.keys(), r=2):
            if cg.connected(n1, n2):
                ignore.add((n1, n2))
                continue
            bp_dist = cg.min_max_bp_distance(n1, n2)[0]
            if bp_dist < self._bp_distance:
                ignore.add((n1, n2))

        interactions = set(
            cg.coords.elements_closer_than(self._distance, ignore))
        return interactions

    def evaluate(self, cg):
        '''
        Calculate the true_positive, false_positive,
        true_negative and false_negative rate for the tertiary
        distances of the elements of cg and the reference structure stored in this class.

        :param cg: The coarse grain model with which self.reference_cg should be compared.
                   Note: The result is only meaningful, if both coarse grained models
                   correspond to the same RNA.
        :return: A dictionary like this: `{"tp": tp, "tn": tn, "fp": fp, "fn": fn}`
        '''
        interactions = self._get_interactions(cg)
        nodes = set(cg.defines.keys())
        allIA = set()
        for n1, n2 in it.combinations(nodes, r=2):
            if cg.connected(n1, n2):
                continue
            bp_dist = cg.min_max_bp_distance(n1, n2)[0]
            if bp_dist < self._bp_distance:
                continue
            allIA.add(tuple(sorted((n1, n2))))

        d = {"tp": 0, "tn": 0, "fp": 0, "fn": 0}
        d["tp"] = len(self._reference_interactions & interactions)
        d["fp"] = len(interactions - self._reference_interactions)
        d["fn"] = len(self._reference_interactions - interactions)
        d["tn"] = len(allIA - (self._reference_interactions | interactions))
        return d


def optimal_superposition(crds1, crds2):
    """
    Returns best-fit rotation matrix as [3x3] numpy matrix for aligning crds1 onto crds2
    using the Kabsch algorithm
    """
    if crds1.shape != crds2.shape:
        raise Incompareable(
            "Cannot superimpose coordinate lists of different length.")
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

def cg_stem_rmsd(cg1, cg2):
    coords1 = cg1.get_ordered_stem_poss()
    coords2 = cg2.get_ordered_stem_poss()
    return rmsd(coords1, coords2)


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
    try:
        # First, try all coordinates, no matter if the seq-ids match
        return rmsd(residues1, residues2)
    except Incompareable:
        # If the number of nts is not equal, issue a warning and try to match based on seq_ids.
        res_ids = set(cg1.seq._seqids) | set(cg2.seq._seqids)
        common_resids = set(cg1.seq._seqids) & set(cg2.seq._seqids)
        residues1 = []
        residues2 = []
        for res in common_resids:
            residues1.append(cg1.get_virtual_residue(res, allow_single_stranded = True))
            residues2.append(cg2.get_virtual_residue(res, allow_single_stranded = True))
        if len(residues1)>len(res_ids)*0.8:
            log.warning("Using only %s common residues for RMSD comparison based on seq_ids.", len(residues1))
            if set(cg1.seq._seqids)-common_resids :
                log.warning("Ignoring from cg1:  %s.", set(cg1.seq._seqids)-common_resids )
            if set(cg2.seq._seqids)-common_resids:
                log.warning("Ignoring from cg2:  %s.", set(cg2.seq._seqids)-common_resids )
            return rmsd(residues1, residues2)
        else:
            log.warning("Cannot compare based on seqids: Intersection of"
                        " resids is too small: %s", common_resids)
        raise Incompareable("Cgs {} and {} cannot be compared according to the RMSD, "
                            "because they do not have the same number of "
                            "virtual residues.".format(cg1.name, cg2.name))


def rmsd_contrib_per_element(cg1, cg2):
    residues1, elems1 = cg1.get_ordered_virtual_residue_poss(
        return_elements=True)
    residues2, elems2 = cg2.get_ordered_virtual_residue_poss(
        return_elements=True)
    if elems1 != elems2:
        raise ValueError(
            "RNAs with different structure are not compareable by RMSD.")
    diff_vecs = _pointwise_deviation(residues1, residues2)
    elem_devs = defaultdict(list)
    for i, elem in enumerate(elems1):
        elem_devs[elem].append(diff_vecs[i])
    return elem_devs


def _pointwise_deviation(crds1, crds2, is_centered=False):
    """
    Helper function for Kabsch RMSD
    """
    if not is_centered:
        crds1 = ftuv.center_on_centroid(crds1)
        crds2 = ftuv.center_on_centroid(crds2)

    os = optimal_superposition(crds1, crds2)
    crds_aligned = np.dot(crds1, os)

    diff_vecs = (crds2 - crds_aligned)

    return diff_vecs


def rmsd_kabsch(crds1, crds2, is_centered=False):
    '''
    Center the coordinate vectors on their centroid
    and then calculate the rmsd.
    '''
    diff_vecs = _pointwise_deviation(crds1, crds2, is_centered)

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
    ds1 = np.array([ftuv.vec_distance(c1, c2)
                    for c1, c2 in it.combinations(coords1, r=2)])
    ds2 = np.array([ftuv.vec_distance(c1, c2)
                    for c1, c2 in it.combinations(coords2, r=2)])

    rmsd = math.sqrt(np.mean((ds1 - ds2) * (ds1 - ds2)))
    #rmsd = math.sqrt(np.mean(ftuv.vec_distance(ds1, ds2)))

    return rmsd


def rmsd_qc_wrap(coords1, coords2, is_centered=False):
    r = rmsd_qc(coords1, coords2, is_centered)
    if np.isnan(r):
        return rmsd_kabsch(coords1, coords2, is_centered)
    return r


try:
    from py_qcprot import rmsd as rmsd_qc  # Faster C version, if available
    rmsd = rmsd_qc_wrap
except:
    rmsd = rmsd_kabsch


def basepair_distance(cg1, cg2):
    # QUESTION Move to forgi.graph? or forgi.utilities
    # Note: An implementation in c (with python bindings) is available in the Vienna RNA package
    dist = 0
    if str(cg1.seq) != str(cg2.seq):  # Compare as strings, to ignore missing and mofified residues
        raise Incompareable(
            "We do not support a basepair distance between rnas with different sequences.")
    for stem in cg1.stem_iterator():
        for bp in cg1.stem_bp_iterator(stem):
            if cg2.pairing_partner(bp[0]) != bp[1]:
                dist += 1
    for stem in cg2.stem_iterator():
        for bp in cg2.stem_bp_iterator(stem):
            if cg1.pairing_partner(bp[0]) != bp[1]:
                dist += 1
    return dist
