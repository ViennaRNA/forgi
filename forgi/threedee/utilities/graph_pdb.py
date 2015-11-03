#!/usr/bin/python

import Bio.PDB as bp
import Bio.PDB as bpdb
import operator
import itertools as it
import collections as col

import os.path as op
import forgi.config as cc
import warnings

import math as m

import numpy as np
import numpy.linalg as nl
import random
import sys

#import forgi.threedee.utilities.average_stem_vres_atom_positions as ftus #Depricated
import forgi.utilities.debug as fud
import forgi.threedee.utilities.my_math as ftum
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.rmsd as cur
import forgi.threedee.utilities.vector as cuv
import forgi.threedee.utilities.vector as ftuv
import forgi

import scipy.optimize as so

catom_name = "C1'"

def stem_stem_orientation(cg, s1, s2):
    '''
    Calculate the orientation of stem s2 in relation to stem s1
    as described by 3 parameters:

    1. The distance between the closest points of the two stems.
    2. The angle between s1 and s2 in the plane formed by the axis of
       the first stem and the vector between the two points closest
       to each on both stems.
    3. The angle of s2 out of the plane formed by their axes.

    @param bg: The BulgeGraph containing the stems.
    @param s1: The name of the first stem
    @param s2: The name of the second stem
    @return: (x,y,z) where x,y and z are the parameters described in
        the description above.
    '''
    # shorten the names a little bit
    s1_p0 = cg.coords[s1][0]
    s1_p1 = cg.coords[s1][1]

    s2_p0 = cg.coords[s2][0]
    s2_p1 = cg.coords[s2][1]

    # The vectors of the axes of the cylinders
    s1_vec = cg.coords[s1][1] - cg.coords[s1][0]
    s2_vec = cg.coords[s2][1] - cg.coords[s2][0]

    # the minimum distance between the two stems, which are represented
    # as line segments
    (i1, i2) = cuv.line_segment_distance(s1_p0, s1_p1, s2_p0, s2_p1)
    i_vec = i2 - i1

    i_rej = cuv.vector_rejection(i_vec, s1_vec)
    plane_vec = np.cross(i_rej, s1_vec)
    # s2_proj is in the intersection plane
    s2_proj_in = cuv.vector_rejection(s2_vec, plane_vec)
    # s2 proj_out is out of the intersection plane
    s2_proj_out = cuv.vector_rejection(s2_vec, i_rej)
    # the normal of the plane defined by the two stem vectors

    #ang1 = cuv.vec_angle(s1_vec, s2_proj_out)
    #ang2 = cuv.vec_angle(s1_vec, s2_proj_in)
    ang1 = cuv.vec_angle(s2_proj_in, s1_vec)
    ang2 = cuv.vec_angle(s2_proj_out, s1_vec)

    #ang3 = cuv.vec_angle(s1_vec, s2_proj_in)
    #ang4 = cuv.vec_angle(s1_vec, s2_proj_out)

    # ever so slightly increased to prevent domain errors
    # in the lateral_offset calculation below
    dist = cuv.magnitude(i_vec) + 0.0001

    ortho_offset = cuv.magnitude(i_rej)
    lateral_offset = m.sqrt(dist * dist - ortho_offset * ortho_offset)

    return (cuv.magnitude(i_vec), ang1,
            ang2, cuv.vec_angle(s1_vec, s2_vec), lateral_offset, ortho_offset)


def get_stem_phys_length(coords):
    '''
    Return the physical length of a stem.

    @param coords: The coordinates of the ends of the helix axis.
    '''

    return cuv.magnitude(coords[1] - coords[0])


def base_normals(pdb_filename):
    '''
    Return a list of the normals for each base in the structure.

    As defined by the average of the cross products between the C2-C5
    and C2-C6 vectors and the N3-C6 and N3-C5 vectors. The origin of
    the vector will be the centroid of these four atoms.

    @param pdb_filename: The name of the pdb file containing the structure
    @return: A list of pairs containing the origin the normal as well as the
        normal itself.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        struct = bp.PDBParser().get_structure('t', pdb_filename)
    chain = list(struct.get_chains())[0]
    origin_norms = []

    for res in chain:
        c2 = res['C2'].get_vector().get_array()
        c5 = res['C5'].get_vector().get_array()
        c6 = res['C6'].get_vector().get_array()
        n3 = res['N3'].get_vector().get_array()

        v1 = cuv.normalize(np.cross(c6 - c2, c5 - c2))
        v2 = cuv.normalize(np.cross(c6 - n3, c5 - n3))

        # take the average of the two, for accuracy or something
        v_norm = (v1 + v2) / 2

        origin = (c2 + c5 + c6 + n3) / 4
        origin_norms += [(origin, v_norm)]

    return origin_norms


def get_twist_angle(coords, twists):
    '''
    Get the angle of the twists with respect to each other.

    @param coords: The coordinates of the ends of the stem.
    @param twists: The two twist vectors.
    @return angle: The angle between the two twist vectors.
    '''

    stem_vec = coords[1] - coords[0]
    basis = cuv.create_orthonormal_basis(stem_vec, twists[0])

    twist2 = cuv.change_basis(twists[1], basis, cuv.standard_basis)
    #assert_allclose(twist2[0], 0., rtol=1e-7, atol=1e-7)

    angle = m.atan2(twist2[2], twist2[1])
    return angle


def twist2_from_twist1(stem_vec, twist1, angle):
    '''
    Get an orientation for the second twist which will place it an
    angle of angle from the first twist.

    @param stem_vec: The vector of the stem.
    @param twist1: The vector of the first twist.
    @param angle: The angular difference between the two twists.
    '''
    basis = cuv.create_orthonormal_basis(stem_vec, twist1)

    twist2_new = np.array([0., m.cos(angle), m.sin(angle)])
    twist2 = np.dot(basis.transpose(), twist2_new)
    #twist2 = cuv.change_basis(twist2_new, cuv.standard_basis, basis)

    return twist2


def get_twist_parameter(twist1, twist2, (u, v)):
    '''
    Calculate how much stem1 must be twisted for its twist vector
    to coincide with that of stem2.

    @param twist1: The twist notator of stem1
    @param twist2: The twist notator of stem2
    @param (u,v): The parameters for rotating stem2 onto stem1
    '''

    rot_mat1 = cuv.rotation_matrix(cuv.standard_basis[2], v)
    rot_mat2 = cuv.rotation_matrix(cuv.standard_basis[1], u - m.pi / 2)

    twist2_new = np.dot(rot_mat1, twist2)
    twist2_new = np.dot(rot_mat2, twist2_new)

    #print "get_twist_parameter twist2:", twist2_new

    return m.atan2(twist2_new[2], twist2_new[1])


def get_stem_orientation_parameters(stem1_vec, twist1, stem2_vec, twist2):
    '''
    Return a parameterization of the orientation of stem2 with respect to
    stem1.

    stem1 -> bulge -> stem2

    @param stem1_vec: The vector representing the axis of stem1
    @param twist1: The twist of stem1 closest to the bulge
    @param stem2_vec: The vector representing teh axis of stem2
    '''

    # Since we will denote the orientation of stem2 with respect to stem1
    # We first need to define a new coordinate system based on stem1

    stem1_basis = cuv.create_orthonormal_basis(stem1_vec, twist1)

    # Transform the vector of stem2 to the new coordinate system
    stem2_new_basis = cuv.change_basis(stem2_vec, stem1_basis,
                                       cuv.standard_basis)
    twist2_new_basis = cuv.change_basis(twist2, stem1_basis,
                                        cuv.standard_basis)

    # Convert the cartesian coordinates to polar coordinates
    (r, u, v) = cuv.spherical_cartesian_to_polar(stem2_new_basis)
    t = get_twist_parameter(twist1, twist2_new_basis, (u, v))

    return (r, u, v, t)


def get_stem_separation_parameters(stem, twist, bulge):
    '''
    Parameterize the location of the bulge with respect to the stem.

    @param stem: The stem vector.
    @param bulge: the bulge vector.
    '''

    stem_basis = cuv.create_orthonormal_basis(stem, twist)
    bulge_new_basis = cuv.change_basis(bulge, stem_basis, cuv.standard_basis)

    return cuv.spherical_cartesian_to_polar(bulge_new_basis)


def get_stem_twist_and_bulge_vecs(bg, bulge, connections):
    '''
    Return the vectors of the stems and of the twists between which
    we want to calculate angles.

    The two vectors will be defined as follows:

    s1e -> s1b -> b -> s2b -> s2e

    The twists will be the two closest to the bulge.

    @param bulge: The name of the bulge separating the two helices.
    @param connections: The two stems that are connected to this bulge.
    @return: (stem1, twist1, stem2, twist2, bulge)
    '''

    s1 = connections[0]
    s2 = connections[1]

    #s1d = bg.defines[s1]
    #s2d = bg.defines[s2]

    mids1 = bg.coords[s1]
    twists1 = bg.twists[s1]

    mids2 = bg.coords[s2]
    twists2 = bg.twists[s2]

    # find out which sides of the stems are closest to the bulge
    # the results will be indexes into the mids array
    (s1b, s1e) = bg.get_sides(s1, bulge)
    (s2b, s2e) = bg.get_sides(s2, bulge)

    # Create directional vectors for the stems
    stem1_vec = mids1[s1b] - mids1[s1e]
    bulge_vec = mids2[s2b] - mids1[s1b]
    stem2_vec = mids2[s2e] - mids2[s2b]

    #twists1_vec = [twists1[s1b], twists1[s1e]]
    #twists2_vec = [twists2[s2e], twists2[s2b]]

    return (stem1_vec, twists1[s1b], stem2_vec, twists2[s2b], bulge_vec)

def stem2_pos_from_stem1(stem1, twist1, params):
    '''
    Get the starting point of a second stem, given the parameters
    about where it's located with respect to stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist parameter of stem1
    @param params: The parameters describing the orientaiton of stem2 wrt stem1
    '''
    (r, u, v) = params
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))

    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    stem2_start = np.dot(stem1_basis.transpose(), stem2)

    return stem2_start


def stem2_pos_from_stem1_1(stem1_basis, params):
    '''
    Get the starting point of a second stem, given the parameters
    about where it's located with respect to stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist parameter of stem1
    @param params: The parameters describing the orientaiton of stem2 wrt stem1
    '''
    (r, u, v) = params
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    stem2_start = np.dot(stem1_basis, stem2)

    return stem2_start


def twist2_orient_from_stem1(stem1, twist1, (u, v, t)):
    '''
    Calculate the position of the twist factor of the 2nd stem from its
    parameters and the first stem.

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (u, v, t): The parameters describing how the twist of stem2 is
                      oriented with respect to stem1
    '''

    twist2_new = np.array([0., m.cos(t), m.sin(t)])

    rot_mat1 = cuv.rotation_matrix(cuv.standard_basis[2], v)
    rot_mat2 = cuv.rotation_matrix(cuv.standard_basis[1], u - m.pi / 2)

    rot_mat = np.dot(rot_mat2, rot_mat1)
    twist2_new = np.dot(nl.inv(rot_mat), twist2_new)

    '''
    twist2_new = dot(inv(rot_mat2), twist2_new)
    twist2_new = dot(inv(rot_mat1), twist2_new)
    '''

    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    twist2_new_basis = cuv.change_basis(twist2_new, cuv.standard_basis,
                                        stem1_basis)

    return twist2_new_basis


def twist2_orient_from_stem1_1(stem1_basis, (u, v, t)):
    '''
    Calculate the position of the twist factor of the 2nd stem from its
    parameters and the first stem.

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (u, v, t): The parameters describing how the twist of stem2 is
                      oriented with respect to stem1
    '''

    twist2_new = np.array([0., m.cos(t), m.sin(t)])

    rot_mat1 = cuv.rotation_matrix(cuv.standard_basis[2], v)
    rot_mat2 = cuv.rotation_matrix(cuv.standard_basis[1], u - m.pi / 2)

    rot_mat = np.dot(rot_mat2, rot_mat1)
    twist2_new = np.dot(nl.inv(rot_mat), twist2_new)

    twist2_new_basis = np.dot(stem1_basis, twist2_new)

    return twist2_new_basis


def stem2_orient_from_stem1(stem1, twist1, (r, u, v)):
    '''
    Calculate the orientation of the second stem, given its parameterization
    and the parameterization of stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (r,u,v): The orientation of stem2 wrt stem1
    '''
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    stem2 = cuv.change_basis(stem2, cuv.standard_basis, stem1_basis)

    return stem2


def stem2_orient_from_stem1_1(stem1_basis, (r, u, v)):
    '''
    Calculate the orientation of the second stem, given its parameterization
    and the parameterization of stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (r,u,v): The orientation of stem2 wrt stem1
    '''
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    #stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    #stem2 = cuv.change_basis(stem2, cuv.standard_basis, stem1_basis)
    stem2 = np.dot(stem1_basis, stem2)

    return stem2


def get_centroid(chain, residue_num):
    residue_num = [int(i) for i in residue_num]
    #print >>sys.stderr, "residue_num:", residue_num
    atoms = []
    for i in residue_num:
        try:
            atoms += [chain[i][catom_name]]
        except KeyError:
            # the C1* atom probably doesn't exist
            continue

    vectors = [atom.get_vector().get_array() for atom in atoms]

    return cuv.get_vector_centroid(vectors)


def get_bulge_centroid(chain, define):
    i = 0
    res_nums = []
    while i < len(define):
        res_nums += range(int(define[i]), int(define[i + 1]) + 1)
        i += 2

    #print >>sys.stderr, "res_nums:", res_nums
    return get_centroid(chain, res_nums)


def get_furthest_c_alpha(cg, chain, stem_end, d, seq_ids=True):
    '''
    Get the position of the c-alpha atom furthest from the end of the stem.
    '''
    max_dist = 0
    furthest_pos = None

    res_ids = it.chain(*cg.get_resseqs(d, seq_ids=seq_ids))

    for i in res_ids:
        try:
            c_apos = chain[i][catom_name].get_vector().get_array()
        except KeyError as ke:
            print >>sys.stderr, "Nucleotide %s missing in element %s" % (str(i),d )
            continue

        dist = cuv.magnitude(stem_end - c_apos)

        if dist >= max_dist:
            max_dist = dist
            furthest_pos = c_apos

    return furthest_pos


def estimate_mids_core(chain, start1, start2, end1, end2):
    '''
    Get the start and end points of a helix.

    Presume that start1, start2, end1 and end2 are
    integers. This may not working with the resname mapping
    where there are missing nucleotides and such.

    Otherwise, it should only be called with the ideal stems
    which are numbered in a very orderely manner.
    '''
    start1 = int(start1)
    start2 = int(start2)

    end1 = int(end1)
    end2 = int(end2)
    #assert(abs(end1 - start1) == abs(end2 - start2))

    fragment_length = end1 - start1 + 1

    if fragment_length < 2:
        error_str = "Helix shorter than 1 nucleotide: "
        error_str += "start1: %d start2: %d " % (start1, start2)
        error_str += "end1: %d end2: %d " % (end1, end2)
        error_str += "fragment_length: %d" % (fragment_length)
        raise Exception(error_str)

    # get the vector between the CA atoms of the two starting residues
    # as well as for the two next residues
    start_vec1 = (chain[start1][catom_name].get_vector() -
                  chain[start2][catom_name].get_vector())
    end_vec1 = (chain[end1][catom_name].get_vector() -
                chain[end2][catom_name].get_vector())

    # get the vector between the CA atoms of the two ending residues
    # as  well as for the two previous residues
    start_vec2 = (chain[start1 + 1][catom_name].get_vector() -
                  chain[start2 - 1][catom_name].get_vector())
    end_vec2 = (chain[end1 - 1][catom_name].get_vector() -
                chain[end2 + 1][catom_name].get_vector())

    # a vector kind of pointing in the direction of the helix
    start_norm_vec = bp.Vector(np.cross(start_vec1.get_array(),
                                        start_vec2.get_array()))
    start_norm_vec.normalize()

    # another vector kind of pointing in the directions of the helix
    end_norm_vec = bp.Vector(np.cross(end_vec2.get_array(),
                                      end_vec1.get_array()))
    end_norm_vec.normalize()

    start_vec1 = -start_vec1
    end_vec1 = -end_vec1

    # I guess we're converting them to Vector format in a weird sort of way...
    # otherwise these steps don't really make sense
    start_axis_vec = start_vec1 + bp.Vector([0., 0., 0.])
    start_axis_vec.normalize()

    end_axis_vec = end_vec1 + bp.Vector([0., 0., 0.])
    end_axis_vec.normalize()

    start_origin = chain[start1][catom_name].get_vector()
    end_origin = chain[end1][catom_name].get_vector()

    start_y_vec = bp.Vector(np.cross(start_norm_vec.get_array(),
                                     start_axis_vec.get_array()))
    start_y_vec.normalize()
    start_c_vec = (start_axis_vec + start_y_vec) / 2
    start_c_vec.normalize()
    start_c_norm = start_origin + start_c_vec / (1 / 8.4)

    end_y_vec = bp.Vector(np.cross(end_norm_vec.get_array(),
                                   end_axis_vec.get_array()))
    end_y_vec.normalize()
    end_c_vec = (end_axis_vec + end_y_vec) / 2
    end_c_vec.normalize()
    end_c_norm = end_origin + end_c_vec / (1 / 8.4)

    mid1 = start_c_norm
    mid2 = end_c_norm

    return (mid1, mid2)


def basenormals_mids(chain, start1, start2, end1, end2):
    '''
    Calculate a helix axis based on the base normal vectors.

    See Laederach et al., RNA 2007.
    '''

    # calculate the scatter matrix using the base normal vectors
    scatter = np.zeros((3, 3))
    base_normals = np.array([0., 0., 0.])

    residue_numbers = [i for i in range(start1, end1 + 1)]
    residue_numbers += [i for i in range(end2, start2 + 1)]

    for r in residue_numbers:
        c2 = chain[r]['C2'].get_vector().get_array()
        c4 = chain[r]['C4'].get_vector().get_array()
        c6 = chain[r]['C6'].get_vector().get_array()

        xi = cuv.normalize(np.cross(c2 - c4, c4 - c6))
        base_normals += xi
        scatter += np.dot(xi, xi.T)

    scatter /= float(len(residue_numbers))
    base_normals /= float(len(residue_numbers))

    # compute the eigenvalues and eigenvectors
    w, v = np.linalg.eig(scatter)
    index, value = max(enumerate(w), key=operator.itemgetter(1))

    # estimate the start and end position, which will be converted scaled
    # to the position of the helix axis
    ''''
    start_pos = (chain[start1][catom_name].get_vector().get_array() +
                 chain[start2][catom_name].get_vector().get_array()) / 2.


    start1_catom = chain[start1 + real_stem_length][catom_name]
    start2_catom = chain[start2 + real_stem_length][catom_name]
    end_pos = (chain[start1 + start1_catom.get_vector().get_array() +
                 chain[start2 - start2_catom.get_vector().get_array()) / 2.

    print start_pos, end_pos
    '''


def get_mids_core_a(chain, start1, start2, end1, end2, 
                    use_template=True):
    '''
    Estimate the stem cylinder using the old method and then refine it
    using fitted parameters.
    '''
    if use_template:
        template_stem_length = 30
    else:
        template_stem_length = end1 - start1 + 1

    real_stem_length = end1 - start1

    tstart1 = 1
    tstart2 = template_stem_length * 2
    tend1 = template_stem_length
    tend2 = template_stem_length + 1

    template_filename = 'ideal_1_%d_%d_%d.pdb' % (tend1, tend2, tstart2)
    filename = forgi.threedee.data_file(op.join('data',
                       template_filename))
    ideal_chain = ftup.get_first_chain(filename)

    est_mids = estimate_mids_core(ideal_chain, tstart1, tstart2, tend1, tend2)
    est_mids = [est_mids[0].get_array(), est_mids[1].get_array()]

    start_pos = (chain[start1][catom_name].get_vector().get_array() +
                 chain[start2][catom_name].get_vector().get_array()) / 2.

    start1_catom = chain[start1 + real_stem_length][catom_name]
    start2_catom = chain[start2 - real_stem_length][catom_name]

    end_pos = (start1_catom.get_vector().get_array() +
               start2_catom.get_vector().get_array()) / 2.

    atom_poss = []
    residue_numbers = [i for i in range(tstart1, tend1 + 1)]
    residue_numbers += [i for i in range(tend2, tstart2 + 1)]

    for rn in residue_numbers:
        #atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        pot_atoms = ['P', 'O3*', 'C3*', 'C4*', 'C5*', 'O5*', 'C1*']
        for atom in pot_atoms:
            try:
                atom_poss += [ideal_chain[rn][atom].get_vector().get_array()]
            except KeyError:
                pass

    mids = fit_circle(est_mids, np.array(atom_poss), start_pos, end_pos)
    return mids


def get_mids_core(cg, chain, define,
                  use_template=True, seq_ids=True):
    ######## Debug function
    '''
    vec = mids[1] - mids[0]
    n1 = []
    n2 = []
    for i in range(0, end1-start1+1):
        notch1 = chain[start1+i][catom_name].get_vector().get_array()
        notch2 = chain[start2-i][catom_name].get_vector().get_array()
        basis = cuv.create_orthonormal_basis(vec)
        notch1_n = cuv.change_basis(notch1, basis, cuv.standard_basis)
        notch2_n = cuv.change_basis(notch2, basis, cuv.standard_basis)

        n1 += [notch1_n[0]]
        n2 += [notch2_n[0]]
    dists1 = [j-i for i,j in zip(n1[:-1], n1[1:])]
    dists2 = [j-i for i,j in zip(n2[:-1], n2[1:])]

    for dist in dists1 + dists2:
        print "ladder", dist
    '''
    ######## End debug

    stem_length = cg.stem_length(define)

    #filename =
    template_filename = 'ideal_1_%d_%d_%d.pdb' % (stem_length, stem_length + 1,
                                                  stem_length * 2)
    filename = forgi.threedee.data_file(op.join('data', template_filename))
    ideal_chain = ftup.get_first_chain(filename)

    # extract the coordinates of the stem from the input chain
    # and place them into a new chain
    stem_chain = bpdb.Chain.Chain(' ')
    resnames = cg.get_resseqs(define, seq_ids=seq_ids)

    for strand in resnames:
        for resname in strand:
            stem_chain.add(chain[resname])


    # get the rotation and translation to rotate the ideal chain onto
    # the stem chain
    rotran = ftup.pdb_rmsd(stem_chain, ideal_chain, sidechains=False,
                          superimpose=True, apply_sup=False)[2]

    # get the mids of the ideal chain using the fit method
    '''
    ideal_mids = get_mids_core_a(ideal_chain, 1, stem_length * 2,
                                 stem_length, stem_length + 1,
                                 use_template=use_template)
    '''

    # average length of a base-pair: 2.547
    mult=0.1
    ideal_mids = np.array([[0., 0., mult], 
                  np.array([0., 0., -mult]) + (stem_length - 1) * np.array([0., 0., -2.547])])

    # apply the rotation and translation to get the mids of the
    # target chain
    #chain_new_mids[0] = 

    chain_new_mids = np.dot(ideal_mids, rotran[0]) + rotran[1]
    #chain_new_mids = [chain_new_mids[0] + mult * stem_vec, chain_new_mids[1] - mult * stem_vec]

    # return it as a Bio.PDB.Vector for some strang reason
    return (bpdb.Vector(chain_new_mids[0]), bpdb.Vector(chain_new_mids[1]))


def get_twists_core(cg, chain, define,
                    mids=None, method=cc.Configuration.mids_method,
                   seq_ids=True):
    '''
    Get the vectors indicating the twist of the cylinder. In actuality,
    this will be the projection of the (ca_start1 - mid1) onto the plane
    defined by (mid2 - mid1).
    '''

    if mids is None:
        mids = get_mids(cg, chain, define,
                        method=cc.Configuration.mids_method, 
                        seq_ids=seq_ids)

    resnames = cg.get_resseqs(define, seq_ids)

    # the first nucleotide of the first strand
    # and the last nucleotide of the second strand
    start_vec1 = chain[resnames[0][0]][catom_name].get_vector() - mids[0]
    end_vec1 = chain[resnames[0][-1]][catom_name].get_vector() - mids[1]

    # the last nucleotide of the first strand
    # and the first nucleotide of the second strand
    start_vec1a = chain[resnames[1][-1]][catom_name].get_vector() - mids[0]
    end_vec1a = chain[resnames[1][0]][catom_name].get_vector() - mids[1]

    notch1 = cuv.vector_rejection(start_vec1.get_array(),
                                  (mids[0] - mids[1]).get_array())
    notch2 = cuv.vector_rejection(end_vec1.get_array(),
                                  (mids[1] - mids[0]).get_array())

    notch1a = cuv.vector_rejection(start_vec1a.get_array(),
                                   (mids[0] - mids[1]).get_array())
    notch2a = cuv.vector_rejection(end_vec1a.get_array(),
                                   (mids[1] - mids[0]).get_array())

    #print >>sys.stderr, "twist_angle_1:", cuv.vec_angle(notch1, notch1a)
    #print >>sys.stderr, "twist_angle_2:", cuv.vec_angle(notch2, notch2a)

    return (cuv.normalize(notch1 + notch1a), cuv.normalize(notch2 + notch2a))
    #return (normalize(notch1), normalize(notch2))


def get_mids(cg, chain, define, method=cc.Configuration.mids_method, seq_ids=True):
    '''
    Get the mid points of the abstract cylinder which represents a helix.

    @param chain: The Bio.PDB representation of the 3D structure.
    @param define: The define of the helix, as per the BulgeGraph
                   definition standard.
    @return: An array of two vectors representing the two endpoints of the
             helix.
    '''

    if method == 'template':
        return get_mids_core(cg, chain, define, seq_ids=seq_ids)
    elif method == 'fit':
        return get_mids_fit_method(cg, chain, define, seq_ids=seq_ids)
    elif method == 'superimpose':
        return get_mids_core(cg, chain, define, 
                             use_template=False, seq_ids=seq_ids)
    else:
        print >>sys.stderr, "Unknown mids method:", method
        sys.exit(1)

    '''
    elif method == 'estimate':
        return estimate_mids_core(cg, chain, define,
                                  stem_length=stem_length)
    elif method == 'basenormals':
        return basenormals_mids(cg, chain, define,
                                stem_length=stem_length)
    '''

def get_twists(cg, chain, define, mids=None, method=cc.Configuration.mids_method,
              seq_ids=True):
    '''
    Get the projection of the (ca - mids) vectors onto the helix axis. This,
    in a sense will define how much the helix twists.

    @param cg: The CoarseGrainRNA representation
    @param chain: The Bio.PDB representation of the 3D structure.
    @param define: The name of the define
    @return: Two vectors which represent the twist of the helix.
    '''

    return get_twists_core(cg, chain, define, mids, method, seq_ids)


'''
def get_helix_vector(chain, start1, start2, end1, end2):
    (mid1, mid2) = get_mids(chain, start1, start2, end1, end2)
    return mid2 - mid1
'''


def virtual_res_3d_pos_core(coords, twists, i, stem_len, stem_inv=None):
    '''
    Calculate the virtual position of the i'th nucleotide in the stem.

    The virtual position extrapolates the position of the residues based
    on the twists of the helix.

    @param bg: The BulgeGraph structure
    @param stem: The name of the stem
    @param i: The i'th residue of the stem

    @return: A tuple containing the point located on the axis of the stem
             and a vector away from that point in the direction of the
             residue.
    '''
    #stem_len = bg.defines[stem][1] - bg.defines[stem][0] + 1
    stem_vec = coords[1] - coords[0]

    # the position of the virtual residue along the axis of
    # the stem
    if stem_len == 1:
        vres_stem_pos = coords[0]
    else:
        vres_stem_pos = coords[0] + (i / float(stem_len - 1)) * stem_vec

    # the angle of the second twist with respect to the first
    if stem_inv is None:
        stem_basis = cuv.create_orthonormal_basis(stem_vec, twists[0])
        t2 = cuv.change_basis(twists[1], stem_basis, cuv.standard_basis)
    else:
        t2 = np.dot(stem_inv, twists[1])

    ang = ftum.atan3(t2[2], t2[1])

    # calculated from an ideal length 30 helix
    average_ang_per_nt = 0.636738030735
    expected_ang = (stem_len - 1) * average_ang_per_nt
    expected_dev = expected_ang
    while (expected_dev - (2 * m.pi) > 0):
        expected_dev -= 2 * m.pi

    if ang < expected_dev:
        forward = 2 * m.pi + ang - expected_dev
        backward = expected_dev - ang
    else:
        forward = ang - expected_dev
        backward = 2 * m.pi + expected_dev - ang

    if forward < backward:
        ang = expected_ang + forward
    else:
        ang = expected_ang - backward

    if stem_len == 1:
        ang = 0.
    else:
        ang_per_nt = ang / float(stem_len - 1)
        ang = ang_per_nt * i

    # the basis vectors for the helix along which the
    # virtual residues will residue
    u = twists[0]
    v = cuv.normalize(np.cross(stem_vec, twists[0]))

    ang_offset = 0.9
    # equation for a circle in 3-space
    return (vres_stem_pos,
            u * m.cos(ang) + v * m.sin(ang),
            u * m.cos(ang + ang_offset) + v * m.sin(ang + ang_offset),
            u * m.cos(ang - ang_offset) + v * m.sin(ang - ang_offset))


def virtual_res_3d_pos(bg, stem, i, stem_inv=None, stem_length=None):
    if stem_length is None:
        return virtual_res_3d_pos_core(bg.coords[stem], bg.twists[stem], i,
                                       bg.stem_length(stem), stem_inv)
    else:
        return virtual_res_3d_pos_core(bg.coords[stem], bg.twists[stem], i,
                                       stem_length, stem_inv)

def bg_virtual_residues(bg):
    vress = []

    for s in bg.sorted_stem_iterator():
        for i in range(bg.stem_length(s)):
            vres = virtual_res_3d_pos(bg, s, i)
            vress += [vres[0] + vres[2], vres[0] + vres[3]]

    return np.array(vress)

def numbered_virtual_residues(bg):
    '''
    Return a list of virtual residues, along with their
    nucleotide positions.

    @param bg: A coarse grain RNA
    @return: A list of tuples containing nucleotides numbers and coordinates.
    '''
    vress = []

    for s in bg.sorted_stem_iterator():
        for i in range(bg.stem_length(s)):
            vres = virtual_res_3d_pos(bg, s, i)
            vress += [(bg.defines[s][0] + i, vres[0] + vres[2]), (bg.defines[s][3] - i, vres[0] + vres[3])]

    return vress

def virtual_res_basis_core(coords, twists, i, stem_len, vec=None):
    '''
    Define a basis based on the location of a virtual stem residue.

    The basis will be defined by the direction of the stem, the direction
    of the virtual residue.

    @param bg: The BulgeGraph structure
    @param stem: The name of the stem
    @param i: The i'th residue of the stem

    @return: A 3x3 matrix defining the coordinate system above.
    '''

    if vec is None:
        (pos, vec, vec_l, vec_r) = virtual_res_3d_pos_core(coords, twists,
                                                           i, stem_len)

    stem_vec = coords[1] - coords[0]

    return cuv.create_orthonormal_basis(stem_vec, vec)


def virtual_res_basis(bg, stem, i, vec=None):
    return virtual_res_basis_core(bg.coords[stem], bg.twists[stem], i,
                                  bg.stem_length(stem), vec)


def pos_to_spos(bg, s1, i1, s2, i2):
    '''
    Convert the location of s2, i2 into the coordinate system
    defined by (s1, i1)

    @param bg: The BulgeGraph containing the stems
    @param s1: The basis stem name
    @param i1: The basis res position
    @param s2: The stem containing the nucleotide to be converted
    @param i2: The nucleotide to be converted position
    '''
    sbasis = virtual_res_basis(bg, s1, i1)
    (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = virtual_res_3d_pos(bg, s1, i1)
    (s2_pos, s2_vec, s2_vec_l, s2_vec_r) = virtual_res_3d_pos(bg, s2, i2)

    #rpos = (s2_pos + 7. * s2_vec) - (s1_pos + 7 * s1_vec)
    rpos = (s2_pos + 7. * s2_vec) - (s1_pos)
    #print "sbasis:", sbasis

    spos = cuv.change_basis(rpos, sbasis, cuv.standard_basis)

    '''
    if spos[1] ** 2 + spos[2] ** 2 < 5 and spos[0] > -5 and spos[0] < 5:
        print >>sys.stderr, "spos:", spos, s1, i1, s2, i2
    '''
    return spos


def spos_to_pos(bg, stem, i, spos):
    '''
    Convert the location of spos from the coordinate system
    of (stem, i) into the standard coordinate system.

    @param bg: The BulgeGraph
    @param stem: The name of the stem in the BulgeGraph
    @param i: The i'th residue in 'stem' which will define the coordinate
              system
    @param spos: The position in the alternate coordinate system

    @return: The coordinates in the cartesian coordinate system of the
        rest of the model.
    '''
    sbasis = virtual_res_basis(bg, stem, i)
    (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = virtual_res_3d_pos(bg, stem, i)
    pos = cuv.change_basis(spos, cuv.standard_basis, sbasis)
    return pos + (s1_pos + s1_vec)


def get_residue_type(i, stem_len):
    '''
    Each nucleotide will be classified according to its position
    within the stem. That way, the distribution of surrounding
    nucleotides will be conditioned on the type of nucleotides.

    This is important due to the fact that nucleotides at the end
    of a stem may have other stem nucleotides in the direction
    of the stem vector. Nucleotides, in the middle shoubulge not due
    to the excluded volume of the stem they occupy.

    @param i: The position of the nucleotide.
    @param stem_len: The length of the stem.

    @return: The type of nucleotide position.
    '''
    assert(i < stem_len)

    return 0


def junction_virtual_res_distance(bg, bulge):
    '''
    Compute the distance between the two virtual residues flanking
    a bulge region.

    @param bg: The BulgeGraph containing the bulge.
    @param bulge: The name of the bulge.
    '''
    cs = list(bg.edges[bulge])

    (s1b, s1e) = bg.get_sides(cs[0], bulge)
    (s2b, s2e) = bg.get_sides(cs[1], bulge)

    if s1b == 1:
        res = bg.v3dposs[cs[0]][bg.stem_length(cs[0]) - 1]
    else:
        res = bg.v3dposs[cs[0]][0]
    (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = res

    if s2b == 1:
        res = bg.v3dposs[cs[1]][bg.stem_length(cs[1]) - 1]
    else:
        res = bg.v3dposs[cs[1]][0]

    (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = res

    dist2 = cuv.vec_distance((vr1_p + 7. * vr1_v), (vr2_p + 7. * vr2_v))
    return dist2

"""
def get_strand_atom_vrn(bg, s, i):
    '''
    Return the strand and which atom to use for the adjacent
    nucleotide distance calculation.
    '''
    if i == 0:
        return (0, 'P', 0)

    # this might have to just be bg.stem_length(s)
    if i == 1:
        return (0, 'O3*', bg.stem_length(s) - 1)
    if i == 2:
        return (1, 'P', bg.stem_length(s) - 1)
    if i == 3:
        return (1, 'O3*', 0)
"""

def junction_virtual_atom_distance(bg, bulge):
    '''
    Compute the distance between the O3' atom and P' atom
    of the two residues that flank the junction segment.

    @param bg: The BulgeGraph containing the bulge.
    @param bulge: The name of the bulge

    @return: A single number corresponding to the distance above.
    '''
    connecting_stems = bg.connections(bulge)    
    (i1, k1) = bg.get_sides_plus(connecting_stems[0], bulge)
    (i2, k2) = bg.get_sides_plus(connecting_stems[1], bulge)
    pos1=bg.defines[connecting_stems[0]][i1]
    pos2=bg.defines[connecting_stems[1]][i2]
    coords=virtual_atoms(bg)
    if i1==0 or i1==2:
        a1="P"
    else:
        a1="O3'"
    if i2==0 or i2==2:
        a2="P"
    else:
        a2="O3'"
    assert a1!=a2
    return cuv.magnitude(coords[pos1][a1]-coords[pos2][a2])

def add_virtual_residues(bg, stem):
    '''
    Create all of the virtual residues and the associated
    bases and inverses.

    @param bg: The CoarseGrainRNA bulge graph containing the stem
    @param stem: The name of the stem to be included
    '''
    stem_vec = bg.coords[stem][1] - bg.coords[stem][0]
    stem_basis = cuv.create_orthonormal_basis(stem_vec, bg.get_twists(stem)[0])
    stem_inv = nl.inv(stem_basis.transpose())

    bg.bases[stem] = stem_basis
    bg.stem_invs[stem] = stem_inv

    for i in range(bg.stem_length(stem)):
        vpos = virtual_res_3d_pos(bg, stem, i, stem_inv=stem_inv)
        vbasis = virtual_res_basis(bg, stem, i, vec=vpos[1])
        vinv = nl.inv(vbasis.transpose())

        bg.vposs[stem][i] = vpos[0]
        bg.vvecs[stem][i] = vpos[1]
        bg.v3dposs[stem][i] = vpos
        bg.vbases[stem][i] = vbasis
        bg.vinvs[stem][i] = vinv


def stem_vres_reference_atoms(bg, chain, s, i):
    '''
    Calculate the position of each atom in the reference of the
    stem and virtual residue.

    @param bg: The BulgeGraph
    @param chain: The PDB representation of the chain
    @param s: The stem identifier
    @param i: The i'th base-pair in the stem

    @return (origin, bases, [dict(atoms), dict(atoms)])
        The origin of the coordinate system (vpos)
        The basises (one for each nucleotide)
        Two dictionaries containing the positions of each atom in its
        respective coordinate system.
    '''
    coords = [dict(), dict()]
    (vpos, vvec, vvec_l, vvec_r) = virtual_res_3d_pos(bg, s, i)
    vec1 = cuv.normalize(bg.coords[s][1] - bg.coords[s][0])
    vec2 = cuv.normalize(vvec)
    basis = cuv.create_orthonormal_basis(vec1, vec2)

    for k in range(2):
        if k == 0:
            r = bg.defines[s][0] + i
        else:
            r = bg.defines[s][3] - i

        for atom in ftup.all_rna_atoms:
            try:
                c = chain[r][atom].coord
                new_c = cuv.change_basis(c - vpos, basis, cuv.standard_basis)
                coords[k][atom] = new_c

            except KeyError:
                continue

    return (vpos, basis, coords)


def bounding_boxes(bg, chain, s, i):
    '''
    Return the bounding boxes of the two nucleotides at the
    i'th position on the stem.

    @param bg: The BulgeGraph
    @param chain: The PDB representation of the chain
    @param s: The stem identifier
    @param i: The i'th base-pair in the stem

    @return: (origin, bases, [(c1, c2), (c1, c2)]) The bases
            (one for each nucleotide) and the corners defining the bounding box
            of the two nucleotides
    '''

    (vpos, bases, atoms) = stem_vres_reference_atoms(bg, chain, s, i)
    corners = []

    for k in range(2):
        min_c = [10000., 10000., 10000.]
        max_c = [-10000., -10000., -10000.]

        for atom in atoms[k].values():
            for j in range(3):
                min_c[j] = min(min_c[j], atom[j])
                max_c[j] = max(max_c[j], atom[j])
        n = min_c
        x = max_c
        corners += [(n, x)]
    return (vpos, bases, corners)


def virtual_residue_atoms(bg, s, i, strand=0):
    '''
    Return two sets of atoms for the virtual residue. One for the nucleotide
    on each strand.

    @param bg: The BulgeGraph
    @param s: The stem
    @param i: The virtual residue number
    @param strand: The strand for which to get the virtual atoms
    '''
    '''
    if vpos == None or vvec == None:
        (vpos, vvec, vvec_l, vvec_r) = virtual_res_3d_pos(bg, s, i)
    if basis == None:
        basis = virtual_res_basis(bg, s, i, vvec).transpose()
    '''
    if s[0]!="s":
        raise ValueError("Expected stem (not single-stranded RNA element), got {}".format(s))

    glob_pos = (bg.defines[s][0] + i, bg.defines[s][3] - i)
    glob_pos = glob_pos[strand]


    coords=virtual_atoms(bg)

    return coords[glob_pos]


def calc_R(xc, yc, p):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((p[:, 0] - xc) ** 2 + (p[:, 1] - yc) ** 2)


def f_2(c, p):
    """ calculate the algebraic distance between the data points and the mean
        circle centered at c=(xc, yc) """
    Ri = calc_R(*c, p=p)
    return Ri - Ri.mean()

def circle_fit(p):
    x = p[:, 0]
    y = p[:, 1]
    x_m = np.mean(x)
    y_m = np.mean(y)

    u = x - x_m
    v = y - y_m

        # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv = sum(u * v)
    Suu = sum(u ** 2)
    Svv = sum(v ** 2)
    Suuv = sum(u ** 2 * v)
    Suvv = sum(u * v ** 2)
    Suuu = sum(u ** 3)
    Svvv = sum(v ** 3)

    # Solving the linear system
    A = np.array([[Suu, Suv], [Suv, Svv]])
    B = np.array([Suuu + Suvv, Svvv + Suuv]) / 2.0
    uc, vc = nl.solve(A, B)

    xc_1 = x_m + uc
    yc_1 = y_m + vc

    return (xc_1, yc_1)
    '''
    Ri_1     = sqrt((x-xc_1)**2 + (y-yc_1)**2)
    R_1      = mean(Ri_1)
    residu_1 = sum((Ri_1-R_1)**2)

    return (xc_1, yc_1, R_1)
    '''

def circle_error(c, p):
    errors = f_2(c, p)
    return sum([e ** 2 for e in errors])


def sum_square(nums):
    return sum([n ** 2 for n in nums])


def f_3(vec, points, est):
    """ calculate the optimal circle for the points (p) projected onto
    the plane orthogonal to v """
    basis = cuv.create_orthonormal_basis(vec)
    new_points = cuv.change_basis(points.T, basis, cuv.standard_basis).T
    p = new_points[:, 1:]

    #center_2, ier=so.leastsq(f_2, center_estimate,args=p)
    center_2 = circle_fit(p)

    return f_2(center_2, p)


def fit_circle(mids, points, start_pos, end_pos):
    '''
    Calculate the projection of points on the plane normal to
    vec and fit a circle to them.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v1, ier = so.leastsq(f_3, mids[1] - mids[0],
                             args=(points, mids[0][1:]))

    basis1 = cuv.create_orthonormal_basis(v1)

    points1 = cuv.change_basis(points.T, basis1, cuv.standard_basis).T
    start_pos1 = cuv.change_basis(start_pos, basis1, cuv.standard_basis)
    end_pos1 = cuv.change_basis(end_pos, basis1, cuv.standard_basis)

    center_5 = circle_fit(points1[:, 1:])

    mids_stem_basis = [[start_pos1[0], center_5[0], center_5[1]],
                       [end_pos1[0], center_5[0], center_5[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis, basis1).T
    '''
    # works!
    mids_stem_basis = [[nmids[0][0], center_4[0], center_4[1]],
                       [nmids[1][0], center_4[0], center_4[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis).T
    '''
    return mids_standard_basis


def extract_define_residues(define, chain):
    '''Extract the residues in the define and return them as a new chain.'''
    c = bpdb.Chain.Chain(' ')
    ranges = zip(*[iter(define)] * 2)
    for r in ranges:
        for x in range(r[0], r[1] + 1):
            c.add(chain[x])
    return c


def fit_circle_old(mids, points, start_pos, end_pos, chain, stem_length,
                   define):
    '''
    Calculate the projection of points on the plane normal to
    vec and fit a circle to them.
    '''
    template_filename = 'ideal_1_%d_%d_%d.pdb' % (stem_length, stem_length + 1,
                                                  stem_length * 2)
    filename = op.join(cc.Configuration.stem_fragment_dir,
                       template_filename)
    ideal_chain = ftup.get_first_chain(filename)

    '''
    rotran = ftup.pdb_rmsd(ideal_chain, chain, sidechains=False,
            superimpose=True, apply_sup=False)[2]
    '''
    rotran = ftup.pdb_rmsd(chain, ideal_chain, sidechains=False,
                          superimpose=True, apply_sup=False)[2]

    ideal_mids = get_mids_core(ideal_chain, 1, stem_length * 2,
                               stem_length, stem_length + 1)

    ideal_mids = np.array([ideal_mids[0].get_array(),
                           ideal_mids[1].get_array()])

    chain_new_mids = np.dot(ideal_mids, rotran[0]) + rotran[1]

    return chain_new_mids

    vec = mids[1] - mids[0]
    basis = cuv.create_orthonormal_basis(vec)
    new_points = cuv.change_basis(points.T, basis, cuv.standard_basis).T
    p = new_points[:, 1:]
    f_3(mids[1] - mids[0], points, mids[0][1:])

    v1, ier = so.leastsq(f_3, mids[1] - mids[0], args=(points, mids[0][1:]),
                         maxfev=10000)
    basis1 = cuv.create_orthonormal_basis(v1)

    points1 = cuv.change_basis(points.T, basis1, cuv.standard_basis).T
    start_pos1 = cuv.change_basis(start_pos, basis1, cuv.standard_basis)
    end_pos1 = cuv.change_basis(end_pos, basis1, cuv.standard_basis)

    center_5 = circle_fit(points1[:, 1:])
    r5 = np.mean(calc_R(*center_5, p=points1[:, 1:]))
    center_estimate = mids[0][1:]
    center_2, ier = so.leastsq(f_2, center_estimate, args=(p))
    center_4 = circle_fit(p)
    r4 = np.mean(calc_R(*center_4, p=p))

    center_x = center_4
    rx = r4
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable='box', aspect=1)
    ax.plot(new_points[:, 1], new_points[:, 2], "o")
    ax.plot(points1[:, 1], points1[:, 2], "go")
    ax.plot(center_x[0], center_x[1], 'ro')
    ax.plot(center_5[0], center_5[1], 'co')

    mids = cuv.change_basis(np.array(mids).T, basis, cuv.standard_basis).T
    ax.plot(mids[0][1], mids[0][2], 'yo')

    #ax.plot(mids[0][1], mids[0][2], 'go')
    circle1 = plt.Circle(center_x, rx, color='r', alpha=0.5)
    circle2 = plt.Circle(mids[0][1:], rx, color='y', alpha=0.5)
    circle3 = plt.Circle(center_5, r5, color='g', alpha=0.5)

    fig.gca().add_artist(circle1)
    fig.gca().add_artist(circle2)
    fig.gca().add_artist(circle3)

    mids_stem_basis = [[start_pos1[0], center_5[0], center_5[1]],
                       [end_pos1[0], center_5[0], center_5[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis1).T
    '''
    # works!
    mids_stem_basis = [[nmids[0][0], center_4[0], center_4[1]],
                       [nmids[1][0], center_4[0], center_4[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis).T
    '''

    plt.show()
    return mids_standard_basis


def get_mids_fit_method(cg, chain, define):
    '''
    Estimate the endpoints of the cylinder axis by fitting it and using
    the rmsd of the best fit circle as the function to minimize.
    '''
    atom_poss = []
    residue_numbers = list(it.chain(*cg.get_resseqs(define)))
    ideal_chain = chain

    for rn in residue_numbers:
        #atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        try:
            '''
            for atom in chain[rn].get_list():
                atom_poss += [atom.get_vector().get_array()]
            '''

            atom_poss += [ideal_chain[rn]['P'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['O3*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C3*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C4*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C5*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['O5*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C1*'].get_vector().get_array()]
        except KeyError:
            pass

    points = np.array(atom_poss)
    mids = estimate_mids_core(chain, start1, start2, end1, end2)
    mids = np.array([mids[0].get_array(), mids[1].get_array()])
    vec = mids[1] - mids[0]
    v1, ier = so.leastsq(f_3, vec, args=(points, mids[0][1:]), maxfev=10000)

    start_pos = (chain[start1]['C1*'].get_vector().get_array() +
                 chain[start2]['C1*'].get_vector().get_array()) / 2.
    end_pos = (chain[end1]['C1*'].get_vector().get_array() +
               chain[end2]['C1*'].get_vector().get_array()) / 2.

    basis1 = cuv.create_orthonormal_basis(v1)
    points1 = cuv.change_basis(points.T, basis1, cuv.standard_basis).T
    start_pos1 = cuv.change_basis(start_pos.T, basis1, cuv.standard_basis).T
    end_pos1 = cuv.change_basis(end_pos.T, basis1, cuv.standard_basis).T

    center = circle_fit(points1[:, 1:])
    mids_stem_basis = [[start_pos1[0], center[0], center[1]],
                       [end_pos1[0], center[0], center[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis1).T
    return [bpdb.Vector(mids_standard_basis[0]),
            bpdb.Vector(mids_standard_basis[1])]


def stem_vec_from_circle_fit(bg, chain, stem_name='s0'):
    '''
    Attempt to find the stem direcion vector given a set of atom positions.

    This will be done by solving for the stem_vector, then using that
    to project the atoms onto a plane orthogonal to that vector. On that plane,
    a circle will be fit to the positions of the atoms. The stem vector that
    gives a circle with the least residuals will be considered the ideal
    stem vector.

    @return: stem_vector
    '''
    atom_poss = []
    #stem_name = 's0'
    for rn in bg.stem_res_numbers(stem_name):
        #atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        try:
            atom_poss += [chain[rn]['P'].get_vector().get_array()]
            atom_poss += [chain[rn]['O3*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C3*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C4*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C5*'].get_vector().get_array()]
            atom_poss += [chain[rn]['O5*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        except KeyError:
            pass

    define = bg.defines[stem_name]
    start_pos = (chain[define[0]]['C1*'].get_vector().get_array() +
                 chain[define[3]]['C1*'].get_vector().get_array()) / 2.
    end_pos = (chain[define[1]]['C1*'].get_vector().get_array() +
               chain[define[2]]['C1*'].get_vector().get_array()) / 2.

    mids = get_mids(bg, chain, bg.defines[stem_name])
    # use the original calculation to provide an estimate for the
    # optimized stem position calculation
    resnames = [bg.seq_ids[d-1] for d in bg.defines[stem_name]]
    stem_chain = extract_define_residues(resnames, chain)

    mids = (mids[0].get_array(), mids[1].get_array())
    return fit_circle_old(mids, np.array(atom_poss),
                          start_pos, end_pos, stem_chain,
                          bg.stem_length(stem_name), bg.defines[stem_name])


def receptor_angle(bg, l, s):
    (i1, i2) = cuv.line_segment_distance(bg.coords[l][0],
                                         bg.coords[l][1],
                                         bg.coords[s][0],
                                         bg.coords[s][1])

    stem_len = bg.stem_length(s)
    stem_vec = bg.coords[s][1] - bg.coords[s][0]

    m1 = cuv.magnitude(i2 - bg.coords[s][0])
    m2 = cuv.magnitude(bg.coords[s][1] - bg.coords[s][0])

    res_num = (stem_len - 1.) * m1 / m2
    vres = virtual_res_3d_pos(bg, s, res_num)[1]
    if cuv.magnitude(i2 - i1) == 0.:
        return 0.

    incoming_angle = cuv.vector_rejection(i1 - i2, stem_vec)

    return cuv.vec_angle(vres, incoming_angle)


def add_stem_information_from_pdb_chain(cg, chain, seq_ids=True):
    '''
    Get the 3D information of the stems.

    Output the mid points of the helices as well as the 'twist' vectors
    which describe the projection of the (ca - mids) vectors onto
    the plane perpendicular to the axis of the helix.

    Add all of this information to the BulgeGraph data structure.

    @param bg: The BulgeGraph.
    @param chain: The Bio.PDB chain representation of the 3D structure.
    '''
    chain = ftup.rename_rosetta_atoms(chain)

    for d in cg.defines.keys():
        if d[0] == 's':
            mids = get_mids(cg, chain, d, seq_ids=seq_ids)
            twists = get_twists(cg, chain, d, seq_ids=seq_ids)

            cg.coords[d] = (mids[0].get_array(), mids[1].get_array())
            cg.twists[d] = (twists[0], twists[1])
            cg.sampled[d] = [cg.name] + cg.defines[d]


def add_bulge_information_from_pdb_chain(bg, chain):
    '''
    Add the information about the starts and ends of the bulges. The stems
    have to be created beforehand.

    Modifies the structure bg.

    @param bg: The BulgeGraph.
    @param chain: The Bio.PDB chain representation of the 3D structure.
    '''
    for d in bg.defines.keys():
        if d[0] != 's':
            if len(bg.edges[d]) == 2:
                edges = list(bg.edges[d])

                s1d = bg.defines[edges[0]]
                s2d = bg.defines[edges[1]]

                (s1b, s1e) = bg.get_sides(edges[0], d)
                (s2b, s2e) = bg.get_sides(edges[1], d)

                mids1 = bg.coords[edges[0]] #get_mids(chain, s1d)
                mids2 = bg.coords[edges[1]] #get_mids(chain, s2d)

                bg.coords[d] = (mids1[s1b], mids2[s2b])

def add_loop_information_from_pdb_chain(bg, chain, seq_ids=True):
    for d in it.chain(bg.hloop_iterator(), ['t1', 'f1']):
        if d not in bg.defines:
            continue

        edges = list(bg.edges[d])
        
        if len(edges) == 0:
            # Odd case where there are no stems in the structure
            # We should find the furthest distance from the first
            # nucleotide
            first_res = None
            for res in chain.get_residues():
                if catom_name in res:
                    first_res = res
                    break

            start_point = first_res[catom_name].get_vector().get_array() 
            centroid = get_furthest_c_alpha(bg, chain, 
                                            first_res[catom_name].get_vector().get_array(), 
                                            d, seq_ids=seq_ids)
        else:
            s1 = edges[0]
            s1d = bg.defines[s1]
            bd = bg.defines[d]

            (s1b, s2b) = bg.get_sides(s1, d)

            mids = bg.coords[s1]
            start_point = mids[s1b]
            #centroid = get_bulge_centroid(chain, bd)

            centroid = get_furthest_c_alpha(bg, chain, mids[s1b], d, seq_ids=seq_ids)

            if centroid is None:
                print >>sys.stderr, "No end found for loop %s... using the end of stem %s" % (d, s1)
                centroid = mids[s1b]

        bg.coords[d] = (start_point, centroid)


def bg_rmsd(bg1, bg2, rmsd_function=None):
    '''
    Calculate the rmsd between the virtual residues of bg1 and bg2.
    '''
    if rmsd_function == None:
        rmsd_function = cur.centered_rmsd

    centers1 = bg_virtual_residues(bg1)
    centers2 = bg_virtual_residues(bg2)

    return rmsd_function(centers1, centers2)

def cylinder_works(cg, cylinders_to_stems, tv, c, r= 4.):
    '''
    Check if all of these points are inside the cylinder.

    '''
    points = [cg.coords[tv][0], cg.coords[tv][1]]
    
    for s in cylinders_to_stems[c]:
        points += [cg.coords[s][0], cg.coords[s][1]]
        
    data = np.array(points)
    datamean = data.mean(axis=0)
    
    uu, dd, vv = np.linalg.svd(data - datamean)
    
    n = vv[0]
    p = data
    a = datamean
    
    dist_vec = (a - p) - (np.dot((a-p), n)[:,np.newaxis]) * n
    mags = [ftuv.magnitude(c) for c in dist_vec]
    
    linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]
    linepts += datamean
    
    '''
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as m3d
    
    ax = m3d.Axes3D(plt.figure())
    ax.scatter3D(*data.T)
    ax.plot3D(*linepts.T)
    '''

    if max(mags) > r:
        return False
    return True
                   
def get_encompassing_cylinders(cg, radius=6.):
    visited = set()
    
    # the stems_in_cylinders dictionary will be indexed by stem name and contain
    # the number of the cylinder it contains
    #stems_to_cylinders = {'s0': 0}
    stems_to_cylinders = dict()
    cylinders_to_stems = col.defaultdict(list)
    
    in_cylinder = set()
    #cylinders_to_stems = {0: ['s0']}
    
    # the first cylinder is equal to the first stem
    #cylinders = {0: cg.coords['s0']}
    cylinders = dict()
    to_visit = [random.choice(cg.defines.keys())]
    
    cylinder_counter = 0
    
    while to_visit:
        tv = to_visit.pop(0)
        
        if tv in visited:
            continue
            
        visited.add(tv)
        for e in cg.edges[tv]:
            to_visit.append(e)
        
        # not interested in non- stem, multiloop or interior loop elements
        if tv[0] != 's' and tv[0] != 'm' and tv[0] != 'i':
            continue
            
        #cylinders_to_check = set(cylinders_to_stems.keys())
        cylinders_to_check = set()
        
        # find which cylinders we need to check
        for e in cg.edges[tv]:
            if e in stems_to_cylinders:
                cylinders_to_check.add(stems_to_cylinders[e])  
            
        found = False
        for c in sorted(cylinders_to_check, key=lambda x: -sum([cg.stem_length(k) for k in cylinders_to_stems[x]])):
            # the new node will definitely be at the end of the cylinder
            #print "checking...:", c, tv
            if cylinder_works(cg, cylinders_to_stems, tv, c, radius):
                cylinders_to_stems[c] += [tv]
                stems_to_cylinders[tv] = c
                found = True
                
                break
        
        if not found:
            # no appropriately sized cylinder has been found so we
            # just create new one containing just this stem
            cylinder_counter += 1
            cylinders_to_stems[cylinder_counter] += [tv]
            stems_to_cylinders[tv] = cylinder_counter
                
    return cylinders_to_stems

def element_coord_system(cg, d):
    '''
    Get a coordinate system for a particular coarse grain element.

    If an element has an axis vector, a, twist vectors t1 and t2,
    then the coordinate system will be a normalized version 
    of the axis a, the second, v2,  will be equal to norm((t1 + t2) / 2.)

    And the third will be equal to a x v2.
    '''

    vec_axis = ftuv.normalize(cg.coords[d][1] - cg.coords[d][0])
    twists = cg.get_twists(d)


    mid_twist = ftuv.normalize(twists[0] + twists[1])
    #return (((cg.coords[d][0] + cg.coords[d][1]) / 2.),
    return (((cg.coords[d][0] + cg.coords[d][1]) / 2.),
            ftuv.create_orthonormal_basis(vec_axis, mid_twist))

def virtual_atoms(cg, given_atom_names=None, sidechain=True):
    '''
    Get a list of virtual atoms for this structure.

    @param cg: The coarse grain structure.
    '''
    return new_virtual_atoms(cg, given_atom_names, sidechain)
    import forgi.threedee.utilities.average_atom_positions as ftua
    coords = col.defaultdict(dict)

    for d in cg.defines.keys():
        origin, basis = element_coord_system(cg, d)

        if d[0] == 'i' or d[0] == 'm':
            conn = cg.connections(d)
            conn_type = cg.connection_type(d, conn)
        else:
            conn_type = 0

        for i,r in it.izip(it.count(),
                           cg.define_residue_num_iterator(d)):

            if given_atom_names is None:
                if sidechain:
                    atom_names = ftup.nonsidechain_atoms + [ cg.seq[r-1]+"."+x for x in ftup.side_chain_atoms[cg.seq[r-1]] ]
                else:
                    atom_names = ftup.nonsidechain_atoms
            else:
                atom_names = given_atom_names

            for aname in atom_names:
                identifier = "%s %s %d %d %s" % (d[0],
                                              " ".join(map(str, cg.get_node_dimensions(d))),
                                              conn_type, i, aname)
                try:
                    coords[r][aname] = origin + ftuv.change_basis(np.array(ftua.avg_atom_poss[identifier]), ftuv.standard_basis, basis )
                except KeyError as ke:
                    pass
    return coords

def new_virtual_atoms(cg, given_atom_names=None, sidechain=True):
    '''
    Get a VirtualAtomLookup Object for the virtual atoms for this structure.

    @param cg: The coarse grain structure.
    '''
    return VirtualAtomsLookup(cg, given_atom_names, sidechain)
class VirtualAtomsLookup(object):
    """An object with a dict-like interface that calculated the virtual atom positions on demand (lazy evaluation)"""
    def __init__(self, cg, given_atom_names=None, sidechain=True):
        """
        :param cg: The coarse grain structure, for which the virtual atoms are generated. 
        ..note :: If cg is modified, new virtual atom positions are calculated.
        """
        self.cg=cg  
        self.given_atom_names=given_atom_names
        self.sidechain=sidechain
    def __getitem__(self, position):
        """
        :returns: A dictionary containing all atoms (as keys) and their positions (as values) for the given residue.
        :param position: The position of the residue in the RNA (starting with 1)
        """
        #Find out the stem for which we have to calculate virtual atom positions
        for key, value in self.cg.defines.items():
            if len(value)<2:
                pass #For multiloops of length 0, value is []
            elif position>=value[0] and position<=value[1]:
                return self._getitem_for_element(key, position)
            elif len(value)==4 and position>=value[2] and position<=value[3]:
                return self._getitem_for_element(key, position)
        assert False, "No return for pos {}".format(position)
    def keys(self):
        k=set()
        for value in self.cg.defines.values():
            if len(value)>1:
                for i in range(value[0], value[1]+1):
                    k.add(i)
            if len(value)>3:
                for i in range(value[2], value[3]+1):
                    k.add(i)
        return k     
    def _getitem_for_element(self, d, pos):
        """
        :returns: A dictionary containing all atoms (as keys) and their positions (as values) for the given residue.
        :param d: The coarse grained element (e.g. "s1")
        :param pos: The position of the residue. It has to be in the element d!
        """
        import forgi.threedee.utilities.average_atom_positions as ftua
        e_coords=dict()
        origin, basis = element_coord_system(self.cg, d)
        if d[0] == 'i' or d[0] == 'm':
            conn = self.cg.connections(d)
            conn_type = self.cg.connection_type(d, conn)
        else:
            conn_type = 0
        for i,r in it.izip(it.count(),
                           self.cg.define_residue_num_iterator(d)):
            if r!=pos: continue
            if self.given_atom_names is None:
                if self.sidechain:
                    atom_names = (ftup.nonsidechain_atoms + [ self.cg.seq[r-1]+"."+x for x in ftup.side_chain_atoms[self.cg.seq[r-1]] ])
                else:
                    atom_names = ftup.nonsidechain_atoms
            else:
                atom_names = given_atom_names
            for aname in atom_names:
                identifier = "%s %s %d %d %s" % (d[0],
                                              " ".join(map(str, self.cg.get_node_dimensions(d))),
                                              conn_type, i, aname)
                if "." in aname:
                    _,_,aname=aname.partition(".")
                try:
                    e_coords[aname] = origin + ftuv.change_basis(np.array(ftua.avg_atom_poss[identifier]), ftuv.standard_basis, basis )
                except KeyError as ke:
                    pass
            return e_coords

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

"""def add_atoms(coords, twists, define, side, seq, new_coords):
    stem_len = define[1] - define[0] + 1

    for i in range(stem_len):
        if side == 0:
            resnum = define[0] + i
        else:
            resnum = define[1] - i
        resname = seq[resnum-1]
        
        vbasis = virtual_res_basis_core(coords, twists, i, stem_len)
        vpos = virtual_res_3d_pos_core(coords, twists, i, stem_len)

        for a in ftus.avg_stem_vres_atom_coords[side][resname].items():
            c = a[1]
            new_coords[resnum][a[0]] = np.dot(vbasis.transpose(), c) + vpos[0]
    
    return new_coords

def cg_atoms(cg):
    '''
    Place atoms as if every segment was a helix.
    '''
    new_coords = col.defaultdict(dict)

    for d in cg.defines:
        coords = cg.coords[d]
        twists = cg.get_twists(d)

        if d[0] == 's' or d[0] == 'i':
            for c in chunks(cg.defines[d], 2):
                for side in range(2):
                    new_coords = add_atoms(coords, twists, c, side, cg.seq, new_coords)
        elif d[0] == 'm':
            side = cg.get_strand(d)

            if len(cg.defines[d]) > 0:
                new_coords = add_atoms(coords, twists, cg.defines[d], side, cg.seq, new_coords)

        elif d[0] == 'f':
            new_coords = add_atoms(coords, twists, cg.defines[d], 0, cg.seq, new_coords)
        
        elif d[0] == 't':
            new_coords = add_atoms(coords, twists, cg.defines[d], 1, cg.seq, new_coords)

    return new_coords"""

def element_distance(cg, l1, l2):
    '''
    Calculate the distance between the two closest points of these
    two elements.
    '''
    (i1, i2) = ftuv.line_segment_distance(cg.coords[l1][0],
                                       cg.coords[l1][1],
                                       cg.coords[l2][0],
                                       cg.coords[l2][1])                                                             
    return ftuv.vec_distance(i1, i2)
            

