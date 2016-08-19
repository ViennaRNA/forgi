#!/usr/bin/python

import itertools as it
import numpy, sys
import numpy as np
import math, warnings
import forgi.utilities.debug as cud
#from forgi.graph.bulge_graph import BulgeGraph
import forgi.threedee.utilities.vector as ftuv

# Shamelessly stolen from:
# http://boscoh.com/protein/rmsd-root-mean-square-deviation

#def rmsd(crds1, crds2):
#    """Returns RMSD between 2 sets of [nx3] numpy array"""
#    #assert(crds1.shape[1] == 3)
#    assert(crds1.shape == crds2.shape)
#    n_vec = numpy.shape(crds1)[0]
#    correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
#    v, s, w_tr = numpy.linalg.svd(correlation_matrix)
#    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
#    if is_reflection:
#        s[-1] = - s[-1]
#    E0 = sum(sum(crds1 * crds1)) + \
#       sum(sum(crds2 * crds2))
#
#    rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
#    rmsd_sq = max([rmsd_sq, 0.0])
#    return numpy.sqrt(rmsd_sq)

def centered_rmsd(c1, c2):
    warnings.warn("forgi.threedee.utilities.rmsd.centered_rmsd is deprecated and will be "  
                  "removed in the future. Use forgi.threedee.utilities.rmsd.rmsd instead.",
                  DeprecationWarning, stacklevel=2)
    return rmsd(c1,c2)

def rmsd(crds1, crds2):
    '''
    Center the coordinate vectors on their centroid
    and then calculate the rmsd.
    '''
    crds1 = ftuv.center_on_centroid(crds1)
    crds2 = ftuv.center_on_centroid(crds2)

    os = optimal_superposition(crds1, crds2)
    crds_aligned = np.dot(crds1, os)

    diff_vecs = (crds2 - crds_aligned)
    vec_lengths = np.sum(diff_vecs * diff_vecs, axis=1)

    return math.sqrt(sum(vec_lengths) / len(vec_lengths))   #rmsd(crds1, crds2)

#def centered_drmsd(crds1, crds2): -> Replaced by drmsd(...)

def optimal_superposition(crds1, crds2):
    """Returns best-fit rotation matrix as [3x3] numpy matrix for aligning crds1 onto crds2"""        
    assert(crds1.shape == crds2.shape)
    if crds1.shape[1] == 3 or crds1.shape[1] == 2:
        correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
        v, s, w_tr = numpy.linalg.svd(correlation_matrix)
        is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
        if is_reflection:
            v[:, -1] = -v[:, -1]
        return numpy.dot(v, w_tr)
    else:
        raise ValueError("Wrong dimension of crds1. Needs to be an array of "
                         "Points in 2D or 3D space. Found {}D".format(crds1.shape[1]))
  

def radius_of_gyration(coords):
    '''
    Calculate the radius of gyration, given a set of coordinates.
    '''
    centroid = sum(coords) / float(len(coords))
    diff_vecs = coords - centroid
    #cud.pv('diff_vecs')
    sums = np.sum(diff_vecs * diff_vecs, axis=1)
    #cud.pv('sums')
    total = sum(sums)
    total /= len(coords)
    rmsd = math.sqrt(total)

    return rmsd

def gyration_tensor(coords, diagonalize = True):
    '''
    Calculate the gyration tensor, given a list of 3D coordinates.

    The gyration tensor is defiens as in doi:10.1063/1.4788616, eq. 4

    :param diagonalize: Diagonalize the tensor to diag(lambda1, lambda2, lambda3)
    '''
    if len(coords[0])!=3:
        raise ValueError("Coordinates for Gyration Tensor must be in 3D space")
    centroid = sum(coords) / float(len(coords))
    diff_vecs = coords - centroid
    tensor = np.zeros((3,3))
    tensor[0,0]=sum(diff_vecs[:,0]**2)
    tensor[1,1]=sum(diff_vecs[:,1]**2)
    tensor[2,2]=sum(diff_vecs[:,2]**2)
    tensor[0,1]=tensor[1,0]=sum(diff_vecs[:,0]*diff_vecs[:,1])
    tensor[0,2]=tensor[2,0]=sum(diff_vecs[:,0]*diff_vecs[:,2])
    tensor[1,2]=tensor[2,1]=sum(diff_vecs[:,1]*diff_vecs[:,2])
    if not diagonalize:
        return tensor
    tensor/=len(coords)
    eigenvalues = np.linalg.eigvals(tensor)
    assert len(eigenvalues)==3
    tensor = np.zeros((3,3))
    for i, v in enumerate(sorted(eigenvalues, reverse=True)):
        tensor[i,i]=v
    return tensor

def anisotropy(coords):
    """
    Calculate the anisotropy of a list of points in 3D space

    See for example doi:10.1063/1.4788616
    """
    g_tensor = gyration_tensor(coords)
    eigVs = [ g_tensor[i,i] for i in range(3) ]
    pp = 0
    for eV1, eV2 in it.combinations(eigVs, 2):
        pp+=eV1*eV2
    return 1-3*(pp)/(sum(eigVs)**2)

def asphericity(coords):
    """
    Calculate the asphericity of a list of points in 3D space

    See for example doi:10.1063/1.4788616
    """
    g_tensor = gyration_tensor(coords)
    return g_tensor[0,0]-(g_tensor[1,1]+g_tensor[2,2])/2
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

