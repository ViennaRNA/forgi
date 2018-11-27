from __future__ import division
from builtins import range

import math
import warnings
import sys
import forgi.threedee.utilities.vector as ftuv
import itertools as it
import logging

import numpy as np

log = logging.getLogger(__name__)
"""
This module contains functions for characterising a single point-cloud.
"""


def radius_of_gyration(coords):
    '''
    Calculate the radius of gyration, given a set of coordinates.
    '''
    centroid = sum(coords) / float(len(coords))
    diff_vecs = coords - centroid
    # cud.pv('diff_vecs')
    sums = np.sum(diff_vecs * diff_vecs, axis=1)
    # cud.pv('sums')
    total = sum(sums)
    total /= len(coords)
    rmsd = math.sqrt(total)

    return rmsd


def gyration_tensor(coords, diagonalize=True):
    '''
    Calculate the gyration tensor, given a list of 3D coordinates.

    The gyration tensor is defiend as in doi:10.1063/1.4788616, eq. 4

    :param diagonalize: Diagonalize the tensor to diag(lambda1, lambda2, lambda3)
    '''
    if len(coords) == 0:
        log.warning("Cannot calculate gyration tensor: coords are empty, returning 'nan'")
        return np.zeros((3, 3)) * float("nan")
    if len(coords[0]) != 3:
        raise ValueError("Coordinates for Gyration Tensor must be in 3D space")
    centroid = sum(coords) / float(len(coords))
    diff_vecs = coords - centroid
    tensor = np.zeros((3, 3))
    tensor[0, 0] = sum(diff_vecs[:, 0]**2)
    tensor[1, 1] = sum(diff_vecs[:, 1]**2)
    tensor[2, 2] = sum(diff_vecs[:, 2]**2)
    tensor[0, 1] = tensor[1, 0] = sum(diff_vecs[:, 0] * diff_vecs[:, 1])
    tensor[0, 2] = tensor[2, 0] = sum(diff_vecs[:, 0] * diff_vecs[:, 2])
    tensor[1, 2] = tensor[2, 1] = sum(diff_vecs[:, 1] * diff_vecs[:, 2])
    if not diagonalize:
        return tensor
    tensor /= len(coords)
    eigenvalues = np.linalg.eigvals(tensor)
    assert len(eigenvalues) == 3
    tensor = np.zeros((3, 3))
    for i, v in enumerate(sorted(eigenvalues, reverse=True)):
        tensor[i, i] = v
    return tensor


def anisotropy(coords):
    """
    Calculate the anisotropy of a list of points in 3D space

    See for example doi:10.1063/1.4788616
    """
    g_tensor = gyration_tensor(coords)
    eigVs = [g_tensor[i, i] for i in range(3)]
    pp = 0
    for eV1, eV2 in it.combinations(eigVs, 2):
        pp += eV1 * eV2
    return 1 - 3 * (pp) / (sum(eigVs)**2)


def asphericity(coords):
    """
    Calculate the asphericity of a list of points in 3D space

    See for example doi:10.1063/1.4788616
    """
    g_tensor = gyration_tensor(coords)
    return g_tensor[0, 0] - (g_tensor[1, 1] + g_tensor[2, 2]) / 2.
