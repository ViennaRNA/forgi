#!/usr/bin/python

import numpy, sys
import numpy as np
#from corgy.graph.bulge_graph import BulgeGraph
from tess.utilities.vector import center_on_centroid

# Shamelessly stolen from:
# http://boscoh.com/protein/rmsd-root-mean-square-deviation

def rmsd(crds1, crds2):
    """Returns RMSD between 2 sets of [nx3] numpy array"""
    assert(crds1.shape[1] == 3)
    assert(crds1.shape == crds2.shape)
    n_vec = numpy.shape(crds1)[0]
    correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
    v, s, w_tr = numpy.linalg.svd(correlation_matrix)
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
    if is_reflection:
        s[-1] = - s[-1]
    E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))

    rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
    rmsd_sq = max([rmsd_sq, 0.0])
    return numpy.sqrt(rmsd_sq)

def centered_rmsd(crds1, crds2):
    '''
    Center the coordinate vectors on their centroid
    and then calculate the rmsd.
    '''
    crds1 = center_on_centroid(crds1)
    crds2 = center_on_centroid(crds2)

    os = optimal_superposition(crds1, crds2)
    crds_aligned = np.dot(crds1, os)

    s2 = sum(sum((crds2 - crds_aligned) * (crds2 - crds_aligned)))
    diff_vecs = (crds2 - crds_aligned)
    sums = np.sum(diff_vecs * diff_vecs, axis=1)
    sqrts = np.sqrt(sums)

    return rmsd(crds1, crds2)

def optimal_superposition(crds1, crds2):
    """Returns best-fit rotation matrix as [3x3] numpy matrix"""
    assert(crds1.shape[1] == 3)
    assert(crds1.shape == crds2.shape)
    correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
    v, s, w_tr = numpy.linalg.svd(correlation_matrix)
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
    if is_reflection:
        v[:, -1] = -v[:, -1]
    return numpy.dot(v, w_tr)
  

'''
def main():
    if len(sys.argv) < 3:
        print "usage: ./rmsd.py coords1 coords2"
        sys.exit(1)

    bg1 = BulgeGraph(sys.argv[1])
    bg2 = BulgeGraph(sys.argv[2])

    coords1 = bg1.get_centers()
    coords2 = bg2.get_centers()

    mat1 = numpy.array(coords1)
    mat2 = numpy.array(coords2)

    centroid1 = numpy.array([0., 0., 0.])
    centroid2 = numpy.array([0., 0., 0.])
    
    for i in range(len(mat2)):
        centroid1 += mat1[i]
        centroid2 += mat2[i]

    centroid1 /= float(len(mat2))
    centroid2 /= float(len(mat2))

    print "centroid1:", centroid1
    print "centroid2:", centroid2

    for i in range(len(mat2)):
        mat1[i] -= centroid1
        mat2[i] -= centroid2

    print "mat1:", mat1
    print "mat2:", mat2

    print "rmsd:", rmsd(mat1, mat2)

if __name__ == '__main__':
    main()
'''
