#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from builtins import map
from builtins import range
from builtins import object
from builtins import zip, str

import itertools as it
import collections as col
import os.path as op
import warnings
import random
import numbers
import sys
import math
import json
import operator
from pprint import pprint
import logging
import uuid


import numpy as np
import numpy.linalg as nl
import numpy.testing as nptest
import scipy.optimize as so
import Bio.PDB as bp
import Bio.PDB as bpdb
import Bio.PDB.PDBExceptions

from logging_exceptions import log_to_exception

import forgi.threedee.utilities.average_stem_vres_atom_positions as ftus
import forgi.utilities.debug as fud
import forgi.threedee.utilities.my_math as ftum
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as cuv
import forgi.threedee.utilities.vector as ftuv
import forgi
from forgi.threedee.utilities.modified_res import change_residue_id
from forgi.utilities.exceptions import CgConstructionError
from forgi.threedee.utilities.pdb import AtomName
try:
  from .cytvec import get_broken_ml_deviation
except ImportError:
  def get_broken_ml_deviation(*args, **kwargs):
    raise ValueError("cython extention not properly installed!")

log = logging.getLogger(__name__)




REFERENCE_CATOM = AtomName("C1'")


try:
    profile  # The @profile decorator from line_profiler (kernprof)
except:
    def profile(x):
        return x


def stem_stem_orientation(cg, s1, s2):
    '''
    Calculate the orientation of stem s2 in relation to stem s1
    as described by 3 parameters:

    1. The distance between the closest points of the two stems.
    2. The angle between s1 and s2 in the plane formed by the axis of
       the first stem and the vector between the two points closest
       to each on both stems.
    3. The angle of s2 out of the plane formed by their axes.

    :param bg: The BulgeGraph containing the stems.
    :param s1: The name of the first stem
    :param s2: The name of the second stem
    :return: (x,y,z) where x,y and z are the parameters described in
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
    lateral_offset = math.sqrt(dist * dist - ortho_offset * ortho_offset)

    return (cuv.magnitude(i_vec), ang1,
            ang2, cuv.vec_angle(s1_vec, s2_vec), lateral_offset, ortho_offset)


def base_normals(pdb_filename):
    '''
    Return a list of the normals for each base in the structure.

    As defined by the average of the cross products between the C2-C5
    and C2-C6 vectors and the N3-C6 and N3-C5 vectors. The origin of
    the vector will be the centroid of these four atoms.

    :param pdb_filename: The name of the pdb file containing the structure
    :return: A list of pairs containing the origin the normal as well as the
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
        v_norm = (v1 + v2) / 2.

        origin = (c2 + c5 + c6 + n3) / 4.
        origin_norms += [(origin, v_norm)]

    return origin_norms


def get_twist_angle(coords, twists):
    '''
    Get the angle of the twists with respect to each other.

    :param coords: The coordinates of the ends of the stem.
    :param twists: The two twist vectors.
    :return angle: The angle between the two twist vectors.
    '''

    stem_vec = coords[1] - coords[0]
    basis = cuv.create_orthonormal_basis(stem_vec, twists[0])

    twist2 = cuv.change_basis(twists[1], basis, cuv.standard_basis)
    #assert_allclose(twist2[0], 0., rtol=1e-7, atol=1e-7)

    angle = math.atan2(twist2[2], twist2[1])
    return angle


def twist2_from_twist1(stem_vec, twist1, angle):
    '''
    Get an orientation for the second twist which will place it an
    angle of angle from the first twist.

    :param stem_vec: The vector of the stem.
    :param twist1: The vector of the first twist.
    :param angle: The angular difference between the two twists.
    '''
    basis = cuv.create_orthonormal_basis(stem_vec, twist1)

    twist2_new = np.array([0., math.cos(angle), math.sin(angle)])
    twist2 = np.dot(basis.transpose(), twist2_new)
    #twist2 = cuv.change_basis(twist2_new, cuv.standard_basis, basis)

    return twist2


def get_twist_parameter(twist1, twist2, u_v):
    '''
    Calculate how much stem1 must be twisted for its twist vector
    to coincide with that of stem2.

    :param twist1: The twist notator of stem1
    :param twist2: The twist notator of stem2
    :param u_v: The parameters u and v for rotating stem2 onto stem1
    '''

    u, v = u_v
    rot_mat1 = cuv.rotation_matrix("z", v)
    rot_mat2 = cuv.rotation_matrix("y", u - math.pi / 2.)

    twist2_new = np.dot(rot_mat1, twist2)
    twist2_new = np.dot(rot_mat2, twist2_new)

    # print "get_twist_parameter twist2:", twist2_new

    return math.atan2(twist2_new[2], twist2_new[1])


def get_stem_orientation_parameters(stem1_vec, twist1, stem2_vec, twist2):
    '''
    Return a parameterization of the orientation of stem2 with respect to
    stem1.

    stem1 -> bulge -> stem2

    :param stem1_vec: The vector representing the axis of stem1
    :param twist1: The twist of stem1 closest to the bulge
    :param stem2_vec: The vector representing the axis of stem2

    :returns: (r,u,v,t) where r,u,v = the stem orientation in polar coordinates
                        and t is the twist parameter.

    '''

    # Since we will denote the orientation of stem2 with respect to stem1
    # We first need to define a new coordinate system based on stem1

    stem1_basis = cuv.create_orthonormal_basis(stem1_vec, twist1)

    log.debug("Stem1 basis \n%s", stem1_basis)
    # Transform the vector of stem2 to the new coordinate system
    stem2_new_basis = cuv.change_basis(stem2_vec, stem1_basis,
                                       cuv.standard_basis)
    log.debug("Stem2 in basis of stem 1 %s", stem2_new_basis)

    twist2_new_basis = cuv.change_basis(twist2, stem1_basis,
                                        cuv.standard_basis)

    # Convert the cartesian coordinates to polar coordinates
    (r, u, v) = cuv.spherical_cartesian_to_polar(stem2_new_basis)
    t = get_twist_parameter(twist1, twist2_new_basis, (u, v))
    log.debug("r %s, u %s, v %s, t %s", r, u, v, t)
    return (r, u, v, t)


def get_stem_separation_parameters(stem, twist, bulge):
    '''
    Parameterize the location of the bulge with respect to the stem.

    :param stem: The stem vector.
    :param bulge: the bulge vector.
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

    :param bulge: The name of the bulge separating the two helices.
    :param connections: The two stems that are connected to this bulge.
    :return: (stem1, twist1, stem2, twist2, bulge)
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
    log.debug("Side of 1st stem %s attached to %s is %s ", s1, bulge, s1b)
    (s2b, s2e) = bg.get_sides(s2, bulge)
    log.debug("Side of 2nd stem %s attached to %s is %s ", s2, bulge, s2b)

    # Create directional vectors for the stems
    #  For ML:           For IL: -> -> ->
    #    |    A
    #    V    |
    #    * -> *
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

    :param stem1: The vector representing the axis of stem1's cylinder
    :param twist1: The twist parameter of stem1
    :param params: The parameters describing the position of stem2 wrt stem1
    '''
    (r, u, v) = params
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))

    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    stem2_start = np.dot(stem1_basis.transpose(), stem2)

    return stem2_start


def stem2_pos_from_stem1_1(transposed_stem1_basis, params):
    '''
    Get the starting point of a second stem, given the parameters
    about where it's located with respect to stem1

    The params of the stat describe the change in the coordinate system of stem1.
    This function converts that to a carthesian vector the standard coordinate system

    :param transposed_stem1_basis: The vtransposed basis of the first stem.
    :param params: The parameters describing the position of stem2 wrt stem1
                   (i.e. the carthesian vector in standard coordinates pointing from stem1 to stem2)
    '''
    (r, u, v) = params
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    stem2_start = np.dot(transposed_stem1_basis, stem2)

    return stem2_start

# Seems to be unused!


def twist2_orient_from_stem1(stem1, twist1, u_v_t):
    '''
    Calculate the position of the twist factor of the 2nd stem from its
    parameters and the first stem.

    :param stem1: The vector representing the axis of stem1's cylinder
    :param twist1: The twist factor of stem1.
    :param u_v_t: The parameters describing how the twist of stem2 is
                      oriented with respect to stem1. A triple `(u, v, t)`
    '''
    u, v, t = u_v_t
    twist2_new = np.array([0., math.cos(t), math.sin(t)])

    rot_mat1 = cuv.rotation_matrix("z", v)
    rot_mat2 = cuv.rotation_matrix("y", u - math.pi / 2.)

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


def twist2_orient_from_stem1_1(stem1_basis, u_v_t):
    '''
    Calculate the position of the twist factor of the 2nd stem from its
    parameters and the first stem.

    :param stem1: The vector representing the axis of stem1's cylinder
    :param twist1: The twist factor of stem1.
    :param u_v_t: The parameters describing how the twist of stem2 is
                      oriented with respect to stem1. A triple `(u, v, t)`
    '''
    u, v, t = u_v_t
    twist2_new = np.array([0., math.cos(t), math.sin(t)])

    rot_mat1 = cuv.rotation_matrix("z", v)
    rot_mat2 = cuv.rotation_matrix("y", u - math.pi / 2.)

    rot_mat = np.dot(rot_mat2, rot_mat1)
    #assert np.allclose(nl.inv(rot_mat), rot_mat.T)
    twist2_new = np.dot(rot_mat.T, twist2_new)

    twist2_new_basis = np.dot(stem1_basis, twist2_new)

    return twist2_new_basis


def stem2_orient_from_stem1(stem1, twist1, r_u_v):
    '''
    Calculate the orientation of the second stem, given its parameterization
    and the parameterization of stem1

    :param stem1: The vector representing the axis of stem1's cylinder
    :param twist1: The twist factor of stem1.
    :param r_u_v: The orientation of stem2 wrt stem1, a triple `(r, u, v)`
    '''
    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    return stem2_orient_from_stem1_1(stem1_basis.transpose(), r_u_v)


def stem2_orient_from_stem1_1(stem1_basis, r_u_v):
    '''
    Calculate the orientation of the second stem, given its parameterization
    and the parameterization of stem1

    :param stem1: The vector representing the axis of stem1's cylinder
    :param twist1: The twist factor of stem1.
    :param r_u_v: The orientation of stem2 wrt stem1, a triple `(r, u, v)`
    '''
    r, u, v = r_u_v
    stem2_in_basis1 = cuv.spherical_polar_to_cartesian((r, u, v))
    stem2 = np.dot(stem1_basis, stem2_in_basis1)

    return stem2


def get_angle_stat_geometry(stem1_vec, twist1, stem2_vec, twist2, bulge_vec):
    """
    :param stem1_vec: The vector of the first stem, pointing TOWARDS the bulge
    :param twist1: The twist vector at the side of stem1 closest to the bulge
    :param stem2_vec: The vector of the second stem, pointing AWAY FROM the bulge
    :param twist2: The twist vector at the side of stem2 closest to the bulge
    :param bulge_vec: The vector from stem1 to stem2

    :returns: T 6-tuple: u,v (the orientation parameters),
                         t (twist parameter) and
                         r1, u1, v1 (the seperation parameters)


        \                A
         \ stem1        /
          \            / stem2
           V          /
            --------->
              bulge

    """
    try:
        # Get the orientations for orienting these two stems
        (r, u, v, t) = get_stem_orientation_parameters(stem1_vec, twist1,
                                                       stem2_vec, twist2)
        (r1, u1, v1) = get_stem_separation_parameters(
            stem1_vec, twist1, bulge_vec)
    except ZeroDivisionError as e:
        with log_to_exception(log, e):
            log.error("Cannot get stat. The 3D coodinates are probably wrong.")
        raise

    return u, v, t, r1, u1, v1


def _get_vstem_coords(cg, broken_ml_name, fixed_stem_name, virtual_stat):
    sides = cg.get_sides(fixed_stem_name, broken_ml_name)
    fixed_s_vec = cg.coords.get_direction(fixed_stem_name)
    if sides[0] == 0:
        fixed_s_vec = -fixed_s_vec
    s_twist = cg.twists[fixed_stem_name][sides[0]]
    fixed_stem_basis = ftuv.create_orthonormal_basis(fixed_s_vec, s_twist)
    vbulge_vec, vstem_vec, vstem_twist = _virtual_stem_from_bulge(
        fixed_stem_basis, virtual_stat)

    vstem_vec *= 5
    vstem_coords0 = cg.coords[fixed_stem_name][sides[0]] + vbulge_vec
    vstem_coords1 = vstem_coords0 + vstem_vec
    return vbulge_vec, vstem_coords0, vstem_vec



def _plot_element(cg, elem, style="o-", name_suffix=""):
    import matplotlib.pyplot as plt
    plt.plot([cg.coords[elem][0][0], cg.coords[elem][1][0]],
             [cg.coords[elem][0][1], cg.coords[elem][1][1]],
             style,
             label=elem + name_suffix)


def _plot_junction_2d(cg, broken_ml):
    """
    TODO: Move this to a proper location
    """
    import matplotlib.pyplot as plt
    plotted = set()
    _plot_element(cg, broken_ml, name_suffix=" broken")
    elem = cg.get_next_ml_segment(broken_ml)
    while elem != broken_ml:
        _plot_element(cg, elem)
        for s in cg.edges[elem]:
            if s not in plotted:
                _plot_element(cg, s)
                plotted.add(s)
        elem = cg.get_next_ml_segment(elem)


def _virtual_stem_from_bulge(prev_stem_basis,  stat):
    """
    Return a virtual stem with length 1 that would be placed
    by stat and prev_stem

    :param prev_stem_basis: The basis of the previous stem.
    :param stat: The angle stat that describes the orientation of the
                 virtual stem from the previous stem.
    """
    transposed_stem1_basis = prev_stem_basis.transpose()
    start_location = stem2_pos_from_stem1_1(
        transposed_stem1_basis, stat.position_params())
    stem_orientation = stem2_orient_from_stem1_1(transposed_stem1_basis,
                                                 [1] + list(stat.orientation_params()))
    twist1 = twist2_orient_from_stem1_1(
        transposed_stem1_basis, stat.twist_params())
    return start_location, stem_orientation, twist1


def get_centroid(chain, residue_num):
    """
    :param residue_num: A list of integers
    """
    residue_num = [int(i) for i in residue_num]
    #print >>sys.stderr, "residue_num:", residue_num
    atoms = []
    for i in residue_num:
        try:
            atoms += [chain[i][REFERENCE_CATOM]]
        except KeyError:
            # the C1* atom probably doesn't exist
            continue

    vectors = [atom.get_vector().get_array() for atom in atoms]

    return cuv.get_vector_centroid(vectors)

# Seems to be unused!


def get_bulge_centroid(chain, define):
    i = 0
    res_nums = []
    while i < len(define):
        res_nums += range(int(define[i]), int(define[i + 1]) + 1)
        i += 2

    #print >>sys.stderr, "res_nums:", res_nums
    return get_centroid(chain, res_nums)


def get_furthest_c_alpha(cg, chain, stem_end, d):
    '''
    Get the position of the c-alpha atom furthest from the end of the stem.
    '''
    seq_ids = True
    max_dist = 0
    furthest_pos = None

    res_ids = it.chain(*cg.get_resseqs(d, seq_ids=seq_ids))

    for chainId, i in res_ids:  # seq_ids now contain chain
        try:
            c_apos = chain[i][REFERENCE_CATOM].get_vector().get_array()
        except KeyError as ke:
            print("Nucleotide %s missing in element %s" %
                  (str(i), d), file=sys.stderr)
            continue

        dist = cuv.magnitude(stem_end - c_apos)

        if dist >= max_dist:
            max_dist = dist
            furthest_pos = c_apos

    return furthest_pos


def stem_from_chains(cg, chains, elem_name):
    """
    This function combines get_mids and get_twists into one more efficient routine.

    :param chains: A dictionary {chain_id: Biopython_PDB_chain}
    :param elem_name: e.g. "s0"
    """
    stem_length = cg.stem_length(elem_name)
    template_filename = 'ideal_1_%d_%d_%d.pdb' % (stem_length, stem_length + 1,
                                                  stem_length * 2)
    filename = forgi.threedee.data_file(op.join('data', template_filename))
    try:
        ideal_chain = ftup.get_first_chain(filename)
    except IOError:
        if stem_length > 40:
            raise CgConstructionError("Cannot create coordinates. "
                                      "Helices with lengths greater than 40 are currently not supported in forgi.")
        else:
            raise
    stem_chain = bpdb.Chain.Chain(' ')
    try:
        residue_ids = cg.get_resseqs(elem_name, seq_ids=True)
    except IndexError as e:
        with log_to_exception(log, e):
            log.error("seq_ids were '%r'", cg.seq_ids)
        raise
    new_residue_ids = []
    for strand in residue_ids:
        for res_id in strand:
            log.debug("Adding residue %s", res_id)
            original_residue = chains[res_id.chain][res_id.resid]
            residue = original_residue.copy()
            try:
                stem_chain.add(residue)
            except Bio.PDB.PDBExceptions.PDBConstructionException as e:
                log.info(
                    "Temporarily changing resid %s to uuid, because this id is present twice (with different chain) in one stem", residue.id)
                change_residue_id(residue, uuid.uuid4())
                stem_chain.add(residue)
            new_residue_ids.append(residue.id)
    rotran = ftup.pdb_rmsd(stem_chain, ideal_chain, sidechains=False,
                           superimpose=True )[2]

    # average length of a base-pair: 2.547
    mult = 0.01  # Stems with 1 bp have a tiny length
    ideal_coords = np.array([[0., 0., mult],
                             np.array([0., 0., -mult]) + (stem_length - 1) * np.array([0., 0., -2.547])])

    coords = np.dot(ideal_coords, rotran[0]) + rotran[1]
    stem_direction = coords[1] - coords[0]

    # the first nucleotide of the first strand
    # and the last nucleotide of the second strand
    first_res = residue_ids[0][0]
    start_vec1 = chains[first_res.chain][first_res.resid][REFERENCE_CATOM].coord - coords[0]
    last_res = residue_ids[0][-1]
    end_vec1 = chains[last_res.chain][last_res.resid][REFERENCE_CATOM].coord - coords[1]

    # the last nucleotide of the first strand
    # and the first nucleotide of the second strand
    first_res_a = residue_ids[1][-1]
    try:
        start_vec1a = chains[first_res_a.chain][first_res_a.resid][REFERENCE_CATOM].coord - coords[0]
    except KeyError as e:
        log.error("Atoms are %s", chains[first_res_a.chain][first_res_a.resid].child_dict)
        raise
    last_res_a = residue_ids[1][0]
    end_vec1a = chains[last_res_a.chain][last_res_a.resid][REFERENCE_CATOM].coord - coords[1]

    notch1 = cuv.vector_rejection(start_vec1, stem_direction)
    notch2 = cuv.vector_rejection(end_vec1, stem_direction)

    notch1a = cuv.vector_rejection(start_vec1a, stem_direction)
    notch2a = cuv.vector_rejection(end_vec1a, stem_direction)

    twists = (cuv.normalize(notch1 + notch1a), cuv.normalize(notch2 + notch2a))

    # Perform some verification
    if False:
        verify_vatom_positions(residue_ids, chains, coords,
                               twists, "stem_{}_from_chain".format(elem_name))

    return coords, twists


def verify_vatom_positions(residue_ids, chains, coords, twists, label=""):
    """
    :param coords: The coords of ONE stem
    """
    res1, res2 = residue_ids[0][0], residue_ids[1][-1]
    res3, res4 = residue_ids[0][-1], residue_ids[1][0]

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = Axes3D(fig)
    strand0 = np.array(
        [chains[r.chain][r.resid]["C1'"].coord for r in residue_ids[0]])
    strand1 = np.array(
        [chains[r.chain][r.resid]["C1'"].coord for r in residue_ids[1]])
    ax.plot(strand0[:, 0], strand0[:, 1],
            strand0[:, 2], "o-", label="forward strand")
    ax.plot(strand1[:, 0], strand1[:, 1], strand1[:, 2],
            "o-", label="backwards strand")
    ax.plot([chains[res1.chain][res1.resid]["C1'"].coord[0], chains[res2.chain][res2.resid]["C1'"].coord[0]],
            [chains[res1.chain][res1.resid]["C1'"].coord[1],
                chains[res2.chain][res2.resid]["C1'"].coord[1]],
            [chains[res1.chain][res1.resid]["C1'"].coord[2], chains[res2.chain][res2.resid]["C1'"].coord[2]], "--", label="bp")

    ax.plot([chains[res3.chain][res3.resid]["C1'"].coord[0], chains[res4.chain][res4.resid]["C1'"].coord[0]],
            [chains[res3.chain][res3.resid]["C1'"].coord[1],
                chains[res4.chain][res4.resid]["C1'"].coord[1]],
            [chains[res3.chain][res3.resid]["C1'"].coord[2], chains[res4.chain][res4.resid]["C1'"].coord[2]], "--", label="bp")
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], "o-", label="STEM")
    twist1 = np.array([coords[0], coords[0] + twists[0],
                       coords[0] + twists[0] * 10])
    twist2 = np.array([coords[1], coords[1] + twists[1],
                       coords[1] + twists[1] * 10])
    ax.plot(twist1[:, 0], twist1[:, 1], twist1[:, 2], "o-", label="Twist1")
    ax.plot(twist2[:, 0], twist2[:, 1], twist2[:, 2], "o-", label="Twist1")
    # Virtual atoms

    vres_pos = []
    vres_bases = []
    c1_vecs = []
    for i, res in enumerate(residue_ids[0]):
        pos = virtual_res_3d_pos_core(
            coords, twists, i, len(residue_ids[0]))[0]
        vres_pos.append(pos)
        basis = virtual_res_basis_core(coords, twists, i, len(residue_ids[0]))
        vres_bases.append(basis)
        c1_vecs.append(ftuv.change_basis(
            chains[res.chain][res.resid]["C1'"].coord - vres_pos[-1], vres_bases[-1], ftuv.standard_basis))

    #print("C1 vecs:", c1_vecs)
    av_c1_vec = np.sum(c1_vecs, axis=0) / len(c1_vecs)
    #print("Av C1' vec = ", av_c1_vec)
    virtual_c1s = []
    for pos, basis in zip(vres_pos, vres_bases):
        virtual_c1s.append(ftuv.change_basis(
            av_c1_vec, ftuv.standard_basis, basis) + pos)

    virtual_c1s = np.array(virtual_c1s)
    ax.plot(virtual_c1s[:, 0], virtual_c1s[:, 1],
            virtual_c1s[:, 2], "o", label="virtual C1'")

    ax.set_title(label)
    ax.legend()
    plt.show()
    #assert False

'''
# TODO: Incomplete. Refactor this out ouf virtual_res_3d_pos_core
def total_helix_rotation(coords, twists, stem_len):
    """
    Calculate the total rotation of the helix in radians from the twists.

    When we calculate the angle between the two twists, we only know
    the true rotation modulo 2*pi (i.e. a rotation of 45 degrees could
    mean 45 degrees or 405 degrees). Depending on the number of nucleotides and
    knowledge of the ideal helix (which turns roughly 30 degrees per base-pair),
    this function outputs the correct result.
    """
    stem_vec = coords[1] - coords[0]

    # the angle of the second twist with respect to the first
    stem_basis = cuv.create_orthonormal_basis(stem_vec, twists[0])
    t2 = cuv.change_basis(twists[1], stem_basis, cuv.standard_basis)
    twist_angle = ftum.atan3(t2[2], t2[1])

    # calculated from an ideal length 30 helix
    average_ang_per_nt = 0.636738030735
    expected_total_ang = (stem_len - 1) * average_ang_per_nt
    expected_twist_ang = expected_total_ang % (2 * math.pi)
'''

def virtual_res_3d_pos_core(coords, twists, i, stem_len, stem_inv=None):
    '''
    Calculate the virtual position of the i'th nucleotide in the stem.

    The virtual position extrapolates the position of the residues based
    on the twists of the helix.

    :param i: The nucleotide number or a list of nt numbers.
    :return: A tuple containing the point located on the axis of the stem
             and a vector away from that point in the direction of the
             residue.
             Or, if i is a list, return a list of tuples.
    '''
    log.debug("virtual_res_3d_pos_core called with %s, %s, %s, %s, %s",
                coords, twists, i, stem_len, stem_inv)
    stem_vec = coords[1] - coords[0]

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
    while (expected_dev - (2 * math.pi) > 0):
        expected_dev -= 2 * math.pi
    # expected_dev uis now between 0 and 360 degrees
    if ang < expected_dev:
        forward = 2 * math.pi + ang - expected_dev
        backward = expected_dev - ang
    else:
        forward = ang - expected_dev
        backward = 2 * math.pi + expected_dev - ang

    if forward < backward:
        total_ang = expected_ang + forward
    else:
        total_ang = expected_ang - backward

    # the basis vectors for the helix along which the
    # virtual residues will residue
    # TODO: Use create_orthonormal_basis for speedup
    u = twists[0]
    v = cuv.normalize(np.cross(stem_vec, twists[0]))

    ang_offset = 0.9

    if isinstance(i, numbers.Number):
        ret_list = False
        i=[i]
    else:
        ret_list = True

    outlist=[]
    for j in i:
        if stem_len == 1:
            ang = 0.
        else:
            ang_per_nt = total_ang / float(stem_len - 1)
            ang = ang_per_nt * j

        # the position of the virtual residue along the axis of
        # the stem
        if stem_len == 1:
            vres_stem_pos = coords[0]
        else:
            vres_stem_pos = coords[0] + (j / float(stem_len - 1)) * stem_vec

        vpos_tuple = (vres_stem_pos,
                u * math.cos(ang) + v * math.sin(ang),
                u * math.cos(ang + ang_offset) + v * math.sin(ang + ang_offset),
                u * math.cos(ang - ang_offset) + v * math.sin(ang - ang_offset))

        outlist.append(vpos_tuple)

    if not ret_list:
        assert len(outlist)==1
        return outlist[0]
    return outlist


def virtual_res_3d_pos(bg, stem, i, stem_inv=None, stem_length=None):
    if stem_length is None:
        return virtual_res_3d_pos_core(bg.coords[stem], bg.twists[stem], i,
                                       bg.stem_length(stem), stem_inv)
    else:
        return virtual_res_3d_pos_core(bg.coords[stem], bg.twists[stem], i,
                                       stem_length, stem_inv)


def virtual_res_basis_core(coords, twists, i, stem_len, vec=None):
    '''
    Define a basis based on the location of a virtual stem residue.

    The basis will be defined by the direction of the stem, the direction
    of the virtual residue.

    :param bg: The BulgeGraph structure
    :param stem: The name of the stem
    :param i: The i'th residue of the stem

    :return: A 3x3 matrix defining the coordinate system above.
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

    :param bg: The BulgeGraph containing the stems
    :param s1: The basis stem name
    :param i1: The basis res position
    :param s2: The stem containing the nucleotide to be converted
    :param i2: The nucleotide to be converted position
    '''
    sbasis = virtual_res_basis(bg, s1, i1)
    (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = virtual_res_3d_pos(bg, s1, i1)
    (s2_pos, s2_vec, s2_vec_l, s2_vec_r) = virtual_res_3d_pos(bg, s2, i2)

    #rpos = (s2_pos + 7. * s2_vec) - (s1_pos + 7 * s1_vec)
    rpos = (s2_pos + 7. * s2_vec) - (s1_pos)
    # print "sbasis:", sbasis

    spos = cuv.change_basis(rpos, sbasis, cuv.standard_basis)

    '''
    if spos[1] ** 2 + spos[2] ** 2 < 5 and spos[0] > -5 and spos[0] < 5:
        print >>sys.stderr, "spos:", spos, s1, i1, s2, i2
    '''
    return spos


"""
def spos_to_pos(bg, stem, i, spos):
    '''
    Convert the location of spos from the coordinate system
    of (stem, i) into the standard coordinate system.

    :param bg: The BulgeGraph
    :param stem: The name of the stem in the BulgeGraph
    :param i: The i'th residue in 'stem' which will define the coordinate
              system
    :param spos: The position in the alternate coordinate system

    :return: The coordinates in the cartesian coordinate system of the
        rest of the model.
    '''
    if stem in bg.vbases and i in bg.vbases[stem]:
        sbasis=bg.vbases[stem][i]
    else:
        sbasis = virtual_res_basis(bg, stem, i)
    pos = cuv.change_basis(spos, cuv.standard_basis, sbasis)

    try:
        (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = bg.v3dposs[stem][i]
    except KeyError as e:
        log.info("in spos_to_pos: KeyError {}. Adding virtual residues for stem {}".format(e, stem))
        add_virtual_residues(bg, stem)
        (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = bg.v3dposs[stem][i]

    #return pos + (s1_pos + s1_vec) #TODO BT: THIS SEEMS WRONG
    return pos + s1_pos
"""


def get_residue_type(i, stem_len):
    '''
    Each nucleotide will be classified according to its position
    within the stem. That way, the distribution of surrounding
    nucleotides will be conditioned on the type of nucleotides.

    This is important due to the fact that nucleotides at the end
    of a stem may have other stem nucleotides in the direction
    of the stem vector. Nucleotides, in the middle shoubulge not due
    to the excluded volume of the stem they occupy.

    :param i: The position of the nucleotide.
    :param stem_len: The length of the stem.

    :return: The type of nucleotide position.
    '''
    assert(i < stem_len)

    return 0


def junction_virtual_res_distance(bg, bulge):
    '''
    Compute the distance between the two virtual residues flanking
    a bulge region.

    :param bg: The BulgeGraph containing the bulge.
    :param bulge: The name of the bulge.
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

    :param bg: The BulgeGraph containing the bulge.
    :param bulge: The name of the bulge

    :return: A single number corresponding to the distance above.
    '''
    connecting_stems = bg.connections(bulge)
    (i1, k1) = bg._get_sides_plus(connecting_stems[0], bulge)
    (i2, k2) = bg._get_sides_plus(connecting_stems[1], bulge)
    pos1 = bg.defines[connecting_stems[0]][i1]
    pos2 = bg.defines[connecting_stems[1]][i2]
    if bulge[0] == "m":
        assert list(sorted([pos1, pos2])) == bg.flanking_nucleotides(bulge)
    if i1 == 0 or i1 == 2:
        a1 = "P"
    else:
        a1 = "O3'"
    if i2 == 0 or i2 == 2:
        a2 = "P"
    else:
        a2 = "O3'"
    assert a1 != a2
    dist = cuv.magnitude(bg.virtual_atoms(
        pos1)[a1] - bg.virtual_atoms(pos2)[a2])
    # if bg.element_length(bulge)==0:
    #    partner1 = bg.pairing_partner(pos1)
    #    partner2 = bg.pairing_partner(pos2)
    #    dist2 = cuv.magnitude(bg.virtual_atoms(pos1)[a2]-bg.virtual_atoms(pos2)[a1])
    #    dist3 = cuv.magnitude(bg.virtual_atoms(partner1)[a1]-bg.virtual_atoms(partner2)[a2])
    #    dist4 = cuv.magnitude(bg.virtual_atoms(partner1)[a2]-bg.virtual_atoms(partner2)[a1])
    #    assert dist < dist2, "{} ({} nts): {} !< {}".format(bulge, bg.element_length(bulge), dist, dist2)
    #    assert dist < dist3, "{} ({} nts): {} !< {}".format(bulge, bg.element_length(bulge), dist, dist3)
    #    assert dist < dist4, "{} ({} nts): {} !< {}".format(bulge, bg.element_length(bulge), dist, dist4)
    return dist


@profile
def add_virtual_residues(bg, element):
    '''
    Create all of the virtual residues and the associated
    bases and inverses for the given stem.

    .. note::
       This is a low-level function used if only the virtual residues of a single
       stems should be added. To add the virtual residues for all stems, use
       `cg.add_all_virtual_residues`

    :param bg: The CoarseGrainRNA bulge graph containing the stem
    :param element: The name of the stem to be included
    '''
    if element[0] == "s":
        return _add_stem_virtual_residues(bg, element)
    else:
        return _add_loop_virtual_residues(bg, element)


def _add_loop_virtual_residues(cg, element):
    if not cg.chains:
        log.info(
            "No virtual residues added for %s, because no pdb chain present", element)
        return
    for i, resid in enumerate(cg.define_residue_num_iterator(element, seq_ids=True)):
        try:
            global_coords = cg.chains[resid.chain][resid.resid]["C1'"].coord
        except KeyError:
            log.warning("Added virtual residue position for residue %s will be "
                      "inaccurate, because no C1' is present. Atoms are %s", resid,
                      list(cg.chains[resid.chain][resid.resid].child_dict.keys()))
            p=np.zeros(3)
            i=0
            for atom in cg.chains[resid.chain][resid.resid]:
                p+=atom.coord
                i+=1
            global_coords = p/i

        origin, basis = element_coord_system(cg, element)
        element_coords = ftuv.change_basis(
            global_coords - origin, basis, ftuv.standard_basis)
        cg.vposs[element][i] = element_coords

def _add_three_points_per_element(cg, element):
    if not cg.chains:
        log.info(
            "No 3 ppints added for %s, because no pdb chain present", element)
        return

    for i, resid in enumerate(cg.define_residue_num_iterator(element, seq_ids=True)):
        base_coords = []
        sugar_coords = []
        phosphor_coord = []
        for atom in cg.chains[resid.chain][resid.resid]:
            if atom.name in ["C1'", "C2'", "C3'", "C4'", "O4'"]:
                sugar_coords.append(atom.coord)
            elif atom.name in ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"]:
                base_coords.append(atom.coord)
            elif atom.name in ["P", "O5'", "OP1", "OP2"]:
                phosphor_coord.append(atom.coord)
            else:
                log.debug("Getting vbase coordinates: Ignoring atom: %s", atom.name)
        for attr, coords in zip([cg.vbase, cg.vsugar, cg.vbackbone],
                               [base_coords, sugar_coords, phosphor_coord]):
            if len(coords) == 0:
                log.warning("Cannot get 3 points for %s, %s (%s): "
                            "only %s sugar atoms, %s base atoms and %s phosphor atoms",
                            element, resid, repr(cg.chains[resid.chain][resid.resid]), len(sugar_coords), len(base_coords), len(phosphor_coord))
                log.warning("atoms are %s", 
                            ", ".join(atom.name for atom in cg.chains[resid.chain][resid.resid]))

            else:
                x=sum(c[0] for c in coords)/len(coords)
                y=sum(c[1] for c in coords)/len(coords)
                z=sum(c[2] for c in coords)/len(coords)
                global_coords = [x,y,z]
                origin, basis = element_coord_system(cg, element)
                element_coords = ftuv.change_basis(
                    global_coords - origin, basis, ftuv.standard_basis)
                attr[element][i] = element_coords

def _add_stem_virtual_residues(bg, stem):
    try:
        from . import cytvec
    except ImportError as e:
        def transposed_inverted(basis):
            return nl.inv(basis.transpose())
    else:
        transposed_inverted = cytvec.transposed_inverted

    stem_vec = bg.coords.get_direction(stem)
    twist_vec = bg.get_twists(stem)[0]
    if stem in bg.bases and np.allclose(stem_vec, bg.bases[stem][0]) and np.allclose(twist_vec, bg.bases[stem][1]):
        stem_inv = bg.stem_invs[stem]
    else:
        stem_basis = cuv.create_orthonormal_basis(stem_vec, twist_vec)
        stem_inv = transposed_inverted(stem_basis)
        bg.bases[stem] = stem_basis
        bg.stem_invs[stem] = stem_inv

    list_of_is = list( range(bg.stem_length(stem)))
    vposlist = virtual_res_3d_pos(bg, stem, list_of_is, stem_inv=stem_inv)
    for i in list_of_is:
        vpos = vposlist[i]
        vbasis = virtual_res_basis(bg, stem, i, vec=vpos[1])
        vinv = transposed_inverted(vbasis)
        bg.vposs[stem][i] = vpos[0]
        bg.vvecs[stem][i] = vpos[1]
        bg.v3dposs[stem][i] = vpos
        bg.vbases[stem][i] = vbasis
        bg.vinvs[stem][i] = vinv


def stem_vres_reference_atoms(bg, s, i):
    '''
    Calculate the position of each atom in the reference of the
    stem and virtual residue.

    :param bg: The BulgeGraph
    :param s: The stem identifier
    :param i: The i'th base-pair in the stem

    :return (origin, basis, [dict(atoms), dict(atoms)])
        The origin of the coordinate system (vpos)
        The basis of the virtual residue
        Two dictionaries containing the positions of each atom in the coordinate system of the virtual residue
    '''
    coords = [dict(), dict()]
    (vpos, vvec, vvec_l, vvec_r) = virtual_res_3d_pos(bg, s, i)
    #vec1 = cuv.normalize(bg.coords[s][1] - bg.coords[s][0])
    #vec2 = cuv.normalize(vvec)
    stem_direction = bg.coords[s][1] - bg.coords[s][0]
    twist = vvec

    basis = cuv.create_orthonormal_basis(stem_direction, twist)

    residue_ids = bg.get_resseqs(s, seq_ids=True)
    for strand in [0, 1]:
        if strand == 0:
            res_id = residue_ids[0][i]
        else:
            res_id = residue_ids[1][-(1 + i)]
        for atom in ftup.all_rna_atoms:
            res = bg.chains[res_id.chain][res_id.resid]
            try:
                c = res[atom].coord
            except KeyError:
                continue
            else:
                new_c = cuv.change_basis(c - vpos, basis, cuv.standard_basis)
                log.debug("Atom %s has coords %s", atom, new_c)
                coords[strand][atom] = new_c

    return (vpos, basis, coords)


def bounding_boxes(bg, s, i):
    '''
    Return the bounding boxes of the two nucleotides at the
    i'th position on the stem.

    :param bg: A BulgeGraph where bg.chains is not None
    :param s: The stem identifier
    :param i: The i'th base-pair in the stem

    :return: (origin, basis, [(c1, c2), (c1, c2)]) The bases
             and the corners defining the bounding box
             of the two nucleotides
    '''

    (vpos, basis, atoms) = stem_vres_reference_atoms(bg, s, i)
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
    return (vpos, basis, corners)


def virtual_residue_atoms(bg, s, i, strand=0):
    '''
    Return the atoms for the virtual residue.

    :param bg: The BulgeGraph
    :param s: The stem
    :param i: The virtual residue number
    :param strand: The strand for which to get the virtual atoms
    '''
    '''
    if vpos == None or vvec == None:
        (vpos, vvec, vvec_l, vvec_r) = virtual_res_3d_pos(bg, s, i)
    if basis == None:
        basis = virtual_res_basis(bg, s, i, vvec).transpose()
    '''
    if s[0] != "s":
        raise ValueError(
            "Expected stem (not single-stranded RNA element), got {}".format(s))

    glob_pos = (bg.defines[s][0] + i, bg.defines[s][3] - i)
    glob_pos = glob_pos[strand]

    return bg.virtual_atoms(glob_pos)


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


def add_stem_information_from_pdb_chains(cg):
    '''
    Get the 3D information of the stems.

    Output the mid points of the helices as well as the 'twist' vectors
    which describe the projection of the (ca - mids) vectors onto
    the plane perpendicular to the axis of the helix.

    Add all of this information to the BulgeGraph data structure.

    :param bg: The BulgeGraph.
    :param chain: The Bio.PDB chain representation of the 3D structure.
    '''
    new_chains = {}
    for name, chain in cg.chains.items():
        new_chains[name] = ftup.rename_rosetta_atoms(chain)

    for d in cg.defines.keys():
        if d[0] == 's':
            coords, twists = stem_from_chains(cg, new_chains, d)
            cg.coords[d] = coords
            stem_dir = cg.coords[d][1] - cg.coords[d][0]
            cg.twists[d] = twists
            assert abs(np.dot(stem_dir, twists[0])) < 10**-10
            assert abs(np.dot(stem_dir, twists[1])) < 10**-10
            #cg.sampled[d] = [cg.name] + cg.defines[d]


def get_incomplete_elements(cg):
    """
    Get an estimated list of cg-elements which have missing residues in the PDB.

    One of many problems with PDB data are residues, for which no
    coordinates could be determined experimentally. This function gives
    an estimated list of cg-elements, which are affected by missing residues.
    """
    incomplete = set()
    for elem in cg.defines:
        for r in cg.define_range_iterator(elem, adjacent=elem[0] != "s"):
            if cg.seq[r[0]:r[1]] != cg.seq.with_missing[r[0]:r[1]]:
                incomplete.add(elem)
    return incomplete


def add_loop_information_from_pdb_chains(bg):
    seq_ids = True
    #log.info("add_loop_information_from_pdb_chains called")
    for d in it.chain(bg.hloop_iterator(), bg.floop_iterator(), bg.tloop_iterator()):
        if d not in bg.defines:
            assert False

        edges = list(bg.edges[d])

        if len(edges) == 0:
            # Odd case where there are no stems in the structure
            # We should find the furthest distance from the first
            # nucleotide
            log.info(
                "add_loop_information_from_pdb_chain: {} has no neighbor".format(d))

            chain_ids = set(x.chain for y in bg.get_resseqs(d) for x in y)
            assert len(chain_ids) == 1
            c, = chain_ids
            chain = bg.chains[c]

            first_res = None
            for res in chain.get_residues():
                if REFERENCE_CATOM in res:
                    first_res = res
                    break
            try:
                start_point = first_res[REFERENCE_CATOM].get_vector().get_array()
            except TypeError:
                if first_res is not None:
                    raise
                else:
                    e = CgConstructionError("The PDB chain does not contain any {} atom (despite containing {} residues).".format(
                        REFERENCE_CATOM, len(list(chain.get_residues()))))
                    with log_to_exception(log, e):
                        log.error(
                            "The chain's last residue only has the following atoms: %s", res.child_list)
                        raise e
            centroid = get_furthest_c_alpha(bg, chain,
                                            first_res[REFERENCE_CATOM].get_vector(
                                            ).get_array(),
                                            d)

        else:
            chain_ids = set(x.chain for y in bg.get_resseqs(d) for x in y)
            assert len(chain_ids) == 1
            c, = chain_ids
            chain = bg.chains[c]

            s1 = edges[0]
            s1d = bg.defines[s1]
            bd = bg.defines[d]

            (s1b, s2b) = bg.get_sides(s1, d)

            mids = bg.coords[s1]
            start_point = mids[s1b]
            #centroid = get_bulge_centroid(chain, bd)

            centroid = get_furthest_c_alpha(bg, chain, mids[s1b], d)

            if centroid is None:
                print("No end found for loop %s... using the end of stem %s" %
                      (d, s1), file=sys.stderr)
                centroid = mids[s1b]

        assert start_point is not None
        assert centroid is not None
        bg.coords[d] = (start_point, centroid)


def _add_loop_vres(cg):
    if len(cg.defines) < 2:
        return  # fifeprime only-cgs have no twists
    log.debug("Adding virtual residues")
    for elem in cg.defines:
        try:
            if elem[0] != "s":
                _add_loop_virtual_residues(cg, elem)
        except:
            log.exception("Could not add virtual residues from PDB for %s, elem %s", cg.name, elem)
        try:
          log.debug("Add 3 points for elem %s", elem)
          _add_three_points_per_element(cg, elem)
        except:
            log.exception("Could not add 3 virtual points from PDB for %s, elem %s", cg.name, elem)

def cylinder_works(cg, cylinders_to_stems, tv, c, r=4.):
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

    dist_vec = (a - p) - (np.dot((a - p), n)[:, np.newaxis]) * n # pylint: disable=invalid-sequence-index
    mags = [ftuv.magnitude(c) for c in dist_vec]

    '''
    linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]
    linepts += datamean


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

    #cylinders_to_stems = {0: ['s0']}

    # the first cylinder is equal to the first stem
    #cylinders = {0: cg.coords['s0']}
    to_visit = [random.choice(list(cg.defines.keys()))]

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
            # print "checking...:", c, tv
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

    assert abs(np.dot(vec_axis, twists[0])) < 10**- \
        10, "{}: {}".format(d, abs(np.dot(vec_axis, twists[0])))
    assert abs(np.dot(vec_axis, twists[1])) < 10**-10
    return (((cg.coords[d][0] + cg.coords[d][1]) / 2.),
            ftuv.create_orthonormal_basis(vec_axis, mid_twist))


def virtual_atoms(cg, given_atom_names=None, sidechain=True):
    '''
    Get a list of virtual atoms for this structure.

    :param cg: The coarse grain structure.
    '''
    return VirtualAtomsLookup(cg, given_atom_names, sidechain)


# Module-level var used for caching.
_average_atom_positions = None


class VirtualAtomsLookup(object):
    """
    An object with a dict-like interface that calculates the virtual atom positions on demand.
    """

    def __init__(self, cg, given_atom_names=None, sidechain=True):
        """
        :param cg: The coarse grain structure, for which the virtual atoms are generated.

        ..note ::
            If cg is modified, new virtual atom positions are calculated.
        """
        self.cg = cg
        self.given_atom_names = given_atom_names
        self.sidechain = sidechain
    #@profile

    def __getitem__(self, position):
        """
        :returns: A dictionary containing all atoms (as keys) and their
                  positions (as values) for the given residue.
        :param position: The position of the residue in the RNA (starting with 1)
        """
        # Find out the stem for which we have to calculate virtual atom positions
        for key, value in self.cg.defines.items():
            if len(value) < 2:
                continue  # For multiloops of length 0, value is []
            elif position >= value[0] and position <= value[1]:
                return self._getitem_for_element(key, position)
            elif len(value) == 4 and position >= value[2] and position <= value[3]:
                return self._getitem_for_element(key, position)
        assert False, "No return for pos {}".format(position)
    #@profile

    def keys(self):
        k = set()
        for value in self.cg.defines.values():
            if len(value) > 1:
                for i in range(value[0], value[1] + 1):
                    k.add(i)
            if len(value) > 3:
                for i in range(value[2], value[3] + 1):
                    k.add(i)
        return k

    def _getitem_for_element(self, d, pos):
        """
        :returns:   A dictionary containing all atoms (as keys) and
                    their positions (as values) for the given residue.
        :param d:   The coarse grained element (e.g. "s1")
        :param pos: The position of the residue. It has to be in the element d!
        """
        global _average_atom_positions
        if d[0] == "s":
            # Use virtual residues for stems.
            return self._getitem_for_stem(d, pos)
        if _average_atom_positions is None:
            log.info("LOADING AV_ATOM_POS")
            import pkgutil
            data = pkgutil.get_data(
                'forgi', 'threedee/data/average_atom_positions.json')
            _average_atom_positions = json.loads(data.decode("ascii"))
        log.debug("Using loaded av_atom_pos")

        e_coords = dict()
        try:
            origin, basis = element_coord_system(self.cg, d)
        except ValueError as e:
            # 0-length hairpin.
            if d[0] == "h" and np.array_equal(self.cg.coords[d][0], self.cg.coords[d][1]):
                warnings.warn(
                    "Returning empty set of virtual atoms for 0-length hairpin")
                return e_coords
            else:
                raise

            print(e, "for position {} in element {} with define {}".format(
                pos, d, self.cg.defines[d]))
            raise
        if d[0] == 'i' or d[0] == 'm':
            conn = self.cg.connections(d)
            conn_type = self.cg.connection_type(d, conn)
        else:
            conn_type = 0
        for i, r in zip(it.count(), self.cg.define_residue_num_iterator(d)):
            if r != pos:
                continue
            if self.given_atom_names is None:
                if self.sidechain:
                    # Seq is now 1-based
                    atom_names = (ftup.nonsidechain_atoms + [
                                  self.cg.seq[r] + "." + x for x in ftup.side_chain_atoms[self.cg.seq[r]]])
                else:
                    atom_names = ftup.nonsidechain_atoms
            else:
                atom_names = self.given_atom_names
            for aname in atom_names:
                identifier = "%s %s %d %d %s" % (d[0],
                                                 " ".join(
                                                     map(str, self.cg.get_node_dimensions(d))),
                                                 conn_type, i, aname)

                if "." in aname:
                    _, _, aname = aname.partition(".")
                try:
                    e_coords[aname] = origin + ftuv.change_basis(
                        np.array(_average_atom_positions[identifier]), ftuv.standard_basis, basis)
                except KeyError as ke:
                    #warnings.warn("KeyError in virtual_atoms. No coordinates found for: {}".format(ke))
                    pass
            return e_coords

    def _getitem_for_stem(self, d, pos):
        log.debug("getitem_for_stem %s, pos %s", d, pos)
        pos_in_stem, side = self.cg.stem_resn_to_stem_vres_side(d, pos)
        assert pos >= 1
        try:
            residue = (self.cg.seq[pos])
        except IndexError as e:
            with log_to_exception(log, e):
                log.error("position {} not in sequence {}".format(
                    pos - 1, self.cg.seq))
            raise
        if self.given_atom_names is None:
            if self.sidechain:
                atom_names = (ftup.nonsidechain_atoms +
                              ftup.side_chain_atoms[residue])
            else:
                atom_names = ftup.nonsidechain_atoms
        else:
            atom_names = self.given_atom_names
        atom_keys = []
        atom_coords = []
        for aname in atom_names:
            if aname[-1] == "*":
                aname_dash = aname[:-1] + "'"
            else:
                aname_dash = aname
            spos = ftus.avg_stem_vres_atom_coords[side][residue][aname_dash]
            # TODO: Maybe we can vectorize this and calculate pos from spos for all atoms of the residue at once.
            atom_keys.append(aname)
            atom_coords.append(spos)
        try:
            # virtual_res_basis(self.cg, d, pos_in_stem)
            vres_basis = self.cg.vbases[d][pos_in_stem]
            # virtual_res_3d_pos(self.cg, d, pos_in_stem)[0]
            vres_pos = self.cg.vposs[d][pos_in_stem]
        except KeyError:
            self.cg.add_all_virtual_residues()
            # virtual_res_basis(self.cg, d, pos_in_stem)
            vres_basis = self.cg.vbases[d][pos_in_stem]
            # virtual_res_3d_pos(self.cg, d, pos_in_stem)[0]
            vres_pos = self.cg.vposs[d][pos_in_stem]

        atom_coords = ftuv.change_basis_vectorized(
            np.array(atom_coords), ftuv.standard_basis, vres_basis) + vres_pos

        return {aname: coord for aname, coord in zip(atom_keys, atom_coords)}


def vres_to_global_coordinates(vres_pos, vres_basis, positions):
    newpos = {}
    for key, v_pos in positions.items():
        pos = ftuv.change_basis(v_pos, ftuv.standard_basis, vres_basis)
        newpos[key] = pos + vres_pos
    return newpos


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


def get_basepair_center(cg, pos):
    """
    The center of a basepair, as defined in doi: 10.1261/rna.305307

    :param pos: The number of one of the two pairing bases
    """
    pos2 = cg.pairing_partner(pos)
    seq1 = cg.seq[pos - 1]
    seq2 = cg.seq[pos2 - 1]
    atoms = {"A": ["C1'", "C8"], "G": ["C1'", "C8"],
             "U": ["C1'", "C6"], "C": ["C1'", "C6"]}
    va1 = cg.virtual_atoms(pos)
    va2 = cg.virtual_atoms(pos2)
    avpos = np.zeros(3)
    try:
        for atom in atoms[seq1]:
            avpos += va1[atom]
        for atom in atoms[seq2]:
            avpos += va2[atom]
    except KeyError:
        log.error("%s\n%s", va1.keys(), va2.keys())
        raise
    avpos /= (len(atoms[seq1]) + len(atoms[seq2]))
    return avpos


def get_basepair_plane(cg, pos):
    """
    The plane of the basepair, as defined in figure 13 of doi: 10.1261/rna.305307

    :param pos: The number of one of the two pairing bases
    """
    pos2 = cg.pairing_partner(pos)
    seq1 = cg.seq[pos - 1]
    seq2 = cg.seq[pos2 - 1]
    va1 = cg.virtual_atoms(pos)
    va2 = cg.virtual_atoms(pos2)
    h_bonds = {"U": {"A": [("O4", "N6"), ("N3", "N1")],
                     "G": [("N3", "O6"), ("O2", "N1")]},
               "A": {"U": [("N6", "O4"), ("N1", "N3")]},
               "G": {"U": [("O6", "N3"), ("N1", "O2")],
                     "C": [("O6", "N4"), ("N1", "N3"), ("N2", "O2")]},
               "C": {"G": [("N4", "O6"), ("N3", "N1"), ("O2", "N2")]}
               }
    #print( seq1, seq2 )
    try:
        hb = h_bonds[seq1][seq2]
    except KeyError:
        # Non-canonical basepair
        warnings.warn("Estimating plane from stem vector for "
                      " non-canonical basepair {}-{} at positions"
                      " {},{}".format(seq1, seq2, pos, pos2))
        # ValueError, if cg.pairing_partner is buggy
        stem, = cg.nucleotides_to_elements([pos, pos2])
        return cg.coords[stem][0] - cg.coords[stem][1]
    else:
        plane = np.zeros(3)
        contribs = 0

        for l1, l2 in it.combinations(hb, 2):
            left_1 = va1[l1[0]]
            left_2 = va1[l2[0]]
            right_1 = va2[l1[1]]
            right_2 = va2[l2[1]]
            add = np.cross(right_1 - left_1, right_2 - left_1)
            if np.any(plane != np.zeros(3)):
                assert ftuv.vec_angle(add, plane) < math.radians(15), ("{}-{}: {}, {}: {}"
                                                                       " degrees".format(seq1, seq2, plane, add, math.degrees(ftuv.vec_angle(add, plane))))
            plane += add
            add = np.cross(right_1 - left_1, left_2 - right_1)
            assert ftuv.vec_angle(add, plane) < math.radians(15), ("{}-{}: {}, {}: {}"
                                                                   " degrees".format(seq1, seq2, plane, add, math.degrees(ftuv.vec_angle(add, plane))))
            plane += add
            add = np.cross(right_2 - left_2, right_2 - left_1)
            assert ftuv.vec_angle(add, plane) < math.radians(15), ("{}-{}: {}, {}: {}"
                                                                   " degrees".format(seq1, seq2, plane, add, math.degrees(ftuv.vec_angle(add, plane))))
            plane += add
            add = np.cross(right_2 - left_2, left_2 - right_1)
            assert ftuv.vec_angle(add, plane) < math.radians(15), ("{}-{}: {}, {}: {}"
                                                                   " degrees".format(seq1, seq2, plane, add, math.degrees(ftuv.vec_angle(add, plane))))
            plane += add
        return ftuv.normalize(plane)
