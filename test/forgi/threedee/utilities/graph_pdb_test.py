from __future__ import division

from builtins import map
from builtins import range
import Bio.PDB as bpdb
import unittest, os
import warnings
import math
import numpy as np
import numpy.testing as nptest
import forgi.threedee.model.coarse_grain as ftmc
import itertools as it

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv
import forgi.graph.bulge_graph as fgb

import forgi.utilities.debug as fud

import logging
log=logging.getLogger(__name__)


class TestGraphPDB(unittest.TestCase):
    '''
    Test some of the rmsd-type functions.
    '''
    def setUp(self):
        pass

    def test_get_incomplete_elements(self):
        db = "(((...(((...).)))))"
        cg = ftmc.CoarseGrainRNA(dotbracket_str=db)
        # Residue number 3 is missing. Probably bulged out from stem 0
        seq_ids = map(str, [1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
        cg.seq_ids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["s0"]))
        # Residue number 4 is missing between s0 and i0
        seq_ids = map(str, [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
        cg.seq_ids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["i0"]))
        # Residue number 5 is missing inside i0
        seq_ids = map(str, [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
        cg.seq_ids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["i0"]))
        # Residue number 17 is missing between s0 and s1, ==> i0
        seq_ids = map(str, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20])
        cg.seq_ids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["i0"]))
        # Residue number 10 is missing  ==> h0
        seq_ids = map(str, [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20])
        cg.seq_ids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["h0"]))
        # Multiple residues are missing
        seq_ids = map(str, [1,2,4,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22])
        cg.seq_ids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["h0", "s0", "i0"]))



    def test_add_loop_information_from_pdb_chain(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1A34.pdb')

    """
    def test_add_stem_information_from_pdb_chains(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1gid.cg')
        pdb_filename = 'test/forgi/threedee/data/1gid_rosetta.pdb'

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s = bpdb.PDBParser().get_structure('temp', pdb_filename)
            chain = list(s.get_chains())[0]

        chain = ftup.renumber_chain(chain, cg.seq_ids)
        ftug.add_stem_information_from_pdb_chains(cg)
    """

    def verify_virtual_twist_angles(self, cg, s):
        sl = cg.stem_length(s)

        for i in range(0, sl):
            (pos, vec, vec_l, vec_r) = ftug.virtual_res_3d_pos_core(cg.coords[s],
                                                                    cg.twists[s],i,sl)

            if i==0:
                nptest.assert_array_equal(pos, cg.coords[s][0])

            if i > 1:
                self.assertGreater(ftuv.vec_angle(vec, prev_vec), 0.53)
                self.assertLess(ftuv.vec_angle(vec, prev_vec), 0.73)

            prev_vec = vec

    def test_first_virtual_res_basis(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')
        basis = ftug.virtual_res_basis(cg, "s0", 0)
        nptest.assert_array_equal(basis, ftuv.create_orthonormal_basis(cg.coords["s0"][1]-cg.coords["s0"][0], cg.twists["s0"][0]))

    def test_angle_between_twists(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb', dissolve_length_one_stems=False)

        self.verify_virtual_twist_angles(cg, 's0')
        self.verify_virtual_twist_angles(cg, 's1')
        self.verify_virtual_twist_angles(cg, 's2')
        self.verify_virtual_twist_angles(cg, 's3')

    def test_coordinates_for_add_virtual_residues(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')
        ftug.add_virtual_residues(cg, 's0')
        #XYZ coordinate for first residue are ok:
        self.assertAlmostEqual(cg.vposs["s0"][0][0], 2.3, delta=3,
                  msg="Wrong x-position for virtual residue 0 of stem s0: {}".format(cg.vposs["s0"][0][0]))
        self.assertAlmostEqual(cg.vposs["s0"][0][1], 1.3, delta=3,
                  msg="Wrong y-position for virtual residue 0 of stem s0: {}".format(cg.vposs["s0"][0][1]))
        self.assertAlmostEqual(cg.vposs["s0"][0][2], 1.0, delta=3,
                  msg="Wrong z-position for virtual residue 0 of stem s0: {}".format(cg.vposs["s0"][0][2]))
        last_residue=cg.stem_length("s0") - 1
        self.assertAlmostEqual(cg.vposs["s0"][last_residue][0], 16, delta=4,
                  msg="Wrong x-position for virtual residue {} of stem s0: {}".format(last_residue,cg.vposs["s0"][last_residue][0]))
        self.assertAlmostEqual(cg.vposs["s0"][last_residue][1], -13, delta=4,
                  msg="Wrong y-position for virtual residue {} of stem s0: {}".format(last_residue,cg.vposs["s0"][last_residue][1]))
        self.assertAlmostEqual(cg.vposs["s0"][last_residue][2], 8, delta=4,
                  msg="Wrong z-position for virtual residue {} of stem s0: {}".format(last_residue,cg.vposs["s0"][last_residue][2]))

    def test_basis_transformation_for_virtual_residues(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')
        ftug.add_virtual_residues(cg, 's0')
        offset=cg.vposs["s0"][0]
        vbasis=cg.vbases["s0"][0]
        local_pos=np.array([0,0,1])
        global_pos = np.dot(vbasis.transpose(), local_pos) + offset
        #Checking dimensions of vectors, just in case...
        self.assertEqual(len(global_pos),3)
        self.assertEqual(len(vbasis),3)
        self.assertEqual(len(vbasis[0]),3)
        #Analytically true:
        self.assertTrue(all(global_pos[x]-vbasis[2][x]-offset[x]<0.0000001 for x in [0,1,2]),
                  msg="global pos for (0,0,1) should be {}+{}={}, but is {} instead.".format(
                                                  vbasis[2], offset, vbasis[2]+offset, global_pos))

    def test_virtual_residue_atoms(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        ftug.add_virtual_residues(cg, 's0')
        ftug.add_virtual_residues(cg, 's1')
        bases_to_test=[]
        bases_to_test.append(ftug.virtual_residue_atoms(cg, 's0', 1, 0))
        bases_to_test.append(ftug.virtual_residue_atoms(cg, 's0', 2, 1))
        bases_to_test.append(ftug.virtual_residue_atoms(cg, 's1', 0, 0))

        #Assert that any two atoms of the same base within reasonable distance to each other
        #(https://en.wikipedia.org/wiki/Bond_length says than a CH-bond is >= 1.06A)
        for va in bases_to_test:
            for k1, v1 in va.items():
                for k2, v2 in va.items():
                    dist=ftuv.magnitude(v1-v2)
                    self.assertLess(dist, 30, msg="Nucleotide too big: "
                                    "Distance between {} and {} is {}".format(k1, k2, dist))
                    if k1!=k2:
                        dist=ftuv.magnitude(v1-v2)
                        self.assertGreater(dist, 0.8, msg="Nucleotide too small: "
                                    "Distance between {} and {} is {}".format(k1, k2, dist))

    def test_virtual_residue_atom_exact_match(self):
        #This test serves to detect unwanted changes in the virtual atom calculation algorithm.
        #It is allowed to fail, if the virtual atom calculation changes.
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')
        ftug.add_virtual_residues(cg, 's0')
        vres= ftug.virtual_residue_atoms(cg, 's0', 1, 0)
        nptest.assert_allclose(vres['C8'], np.array([ 5.23455929, -2.9606417 , -2.18156476]))
        nptest.assert_allclose(vres['N2'], np.array([ 6.99285237,  2.32505693, -1.95868568]))


    """def test_numbered_virtual_residues(self): #Not USED (?)
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        nres = ftug.numbered_virtual_residues(cg)
        #fud.pv('nres')
        #TODO assert something"""

class TestOrientation(unittest.TestCase):
    def setUp(self):
        pass
    def test_get_stem_orientation_parameters(self):
        stem1_vec = np.array([0., 0., 1.])
        twist1    = np.array([0., 1., 0.])
        stem2_vec = np.array([0., 0., 1.])
        twist2    = np.array([1., 0., 0.])
        r, u, v, t = ftug.get_stem_orientation_parameters(stem1_vec, twist1, stem2_vec, twist2)
        self.assertEqual(r, 1)
        self.assertEqual(u, math.pi/2)
        self.assertEqual(v, 0)

        stem2_vec = np.array([0., 1., 0.])
        r, u, v, t = ftug.get_stem_orientation_parameters(stem1_vec, twist1, stem2_vec, twist2)
        self.assertEqual(r, 1)
        self.assertEqual(u, math.pi/2)
        self.assertEqual(v, math.pi/2)

        stem2_vec = np.array([0., 1., 1.])
        r, u, v, t = ftug.get_stem_orientation_parameters(stem1_vec, twist1, stem2_vec, twist2)
        self.assertEqual(r, math.sqrt(2))
        self.assertEqual(u, math.pi/2)
        self.assertEqual(v, math.pi/4)

    def test_stem2_orient_from_stem1(self):
        stem1_vec = np.array([0., 0., 1.])
        twist1    = np.array([0., 1., 0.])
        r = 1
        u = math.pi/2
        v = 0
        nptest.assert_allclose(ftug.stem2_orient_from_stem1(stem1_vec, twist1, (r,u,v)), np.array([0.,0,1]), atol=10**-15)
        v = math.pi/2
        nptest.assert_allclose(ftug.stem2_orient_from_stem1(stem1_vec, twist1, (r,u,v)), np.array([0.,1,0]), atol=10**-15)
        r = math.sqrt(2)
        v = math.pi/4
        nptest.assert_allclose(ftug.stem2_orient_from_stem1(stem1_vec, twist1, (r,u,v)), np.array([0.,1,1]), atol=10**-15)

class TestVirtualStats(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3way.cg')

    def assert_stats_equal(self, stat1, stat2):

        nptest.assert_allclose(stat1.position_params()[1], stat2.position_params()[1], atol=10**-12)
        nptest.assert_allclose(stat1.position_params()[2]%(2*math.pi), stat2.position_params()[2]%(2*math.pi), atol=10**-12)

        nptest.assert_allclose(stat1.twist_params(), stat2.twist_params(), atol=10**-12)
        nptest.assert_allclose(stat1.orientation_params(), stat2.orientation_params(), atol=10**-12)

    def assert_all_angles_different(self, stat1, stat2):
        # Every angle has to change, but the length of the seperation does not change.
        if not np.all( np.logical_not( np.isclose(stat1.position_params()[1:], stat2.position_params()[1:]))):
            assert False, "Not all position parameters changed: {}, {}".format(stat1.position_params(), stat2.position_params())
        if not np.all( np.logical_not( np.isclose(stat1.twist_params(), stat2.twist_params()))):
            log.info("Stat1: Position (carthesian))")
            log.info(ftug.stem2_pos_from_stem1_1(ftuv.standard_basis,stat1.position_params()))
            log.info("Stat1: Orientation (carthesian)")
            log.info(ftug.stem2_orient_from_stem1_1(ftuv.standard_basis, [1]+list(stat1.orientation_params())))
            log.info("Stat2: Position (carthesian))")
            log.info(ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, stat2.position_params()))
            log.info("Stat2: Orientation (carthesian)")
            log.info(ftug.stem2_orient_from_stem1_1(ftuv.standard_basis, [1]+list(stat2.orientation_params())))
            assert False, "Not all twist parameters changed: {}, {}".format(stat1.twist_params(), stat2.twist_params())
        if not np.all( np.logical_not( np.isclose(stat1.orientation_params(), stat2.orientation_params()))):
            assert False, "Not all orientation parameters changed: {}, {}".format(stat1.orientation_params(), stat2.orientation_params())


    def verify_virtual_stat(self, cg, ml1, ml2):
        virtual_stat = ftug.get_virtual_stat(cg, ml1, ml2)

        stems1 = cg.edges[ml1]
        stems2 = cg.edges[ml2]
        middle_stem = stems1 & stems2
        # The involved stems:
        stem1, = stems1 - middle_stem
        stem2, = stems2 - middle_stem
        middle_stem, = middle_stem
        # Now we modify the cg
        both_defines = np.array(cg.define_a(ml1) + cg.define_a(ml2))
        mask = np.array([ (d not in cg.defines[middle_stem]) for d in both_defines ], dtype = bool)
        log.error("Both_defines %s, middle_stem %s, mask %s", both_defines, cg.defines[middle_stem], mask)
        both_defines = both_defines + [+1, -1, +1, -1]
        both_defines = both_defines[mask]
        log.debug("Masked: %s", both_defines)
        both_defines = list(sorted(both_defines))
        assert len(both_defines) ==2
        log.error("New define: %s", both_defines)
        cg.defines[ml1] = both_defines
        cg.remove_vertex(middle_stem)
        cg.remove_vertex(ml2)
        cg._add_edge(ml1, stem2)
        new_stat1, _ = cg.get_stats(ml1)
        self.assert_stats_equal(virtual_stat, new_stat1)

    def test_get_virtual_stat(self):
        cg = self.cg
        mst = cg.get_mst()
        ml1, ml2 = [ m for m in mst if m[0]=="m" ]
        if cg.get_next_ml_segment(ml2)==ml1:
            ml1, ml2 = ml2, ml1
        self.verify_virtual_stat(cg, ml1, ml2)
    def test_get_virtual_stat_reverse(self):
        cg = self.cg
        mst = cg.get_mst()
        ml1, ml2 = [ m for m in mst if m[0]=="m" ]
        if cg.get_next_ml_segment(ml2)==ml1:
            ml1, ml2 = ml2, ml1
        self.verify_virtual_stat(cg, ml2, ml1)

    def test_get_virtual_stat_different_mst_1(self):
        cg = self.cg
        mst = cg.get_mst()
        ml1, ml2 = [ m for m in mst if m[0]=="m" ]
        broken_m, = [m for m in cg.defines if m not in mst ]
        if cg.get_next_ml_segment(ml2)==ml1:
            ml1, ml2 = ml2, ml1
        mst.add(broken_m)
        mst.remove(ml1)
        cg.mst = mst
        cg.traverse_graph()
        #We have to patch the cg object (which is now inconsistent) so stat-exception will not break.
        cg.find_mlonly_multiloops = lambda *args: [ml2, broken_m]
        cg.describe_multiloop = lambda *args: set()

        self.verify_virtual_stat(cg, ml2, broken_m)
    def test_get_virtual_stat_different_mst_2(self):
        cg = self.cg
        mst = cg.get_mst()
        ml1, ml2 = [ m for m in mst if m[0]=="m" ]
        broken_m, = [m for m in cg.defines if m not in mst ]
        if cg.get_next_ml_segment(ml2)==ml1:
            ml1, ml2 = ml2, ml1
        mst.add(broken_m)
        mst.remove(ml2)
        cg.mst = mst
        cg.traverse_graph()
        self.verify_virtual_stat(cg, ml1, broken_m)

    def test_get_virtual_stat_not_part_of_mst1(self):
        cg = self.cg
        mst = cg.get_mst()
        ml1, ml2 = [ m for m in mst if m[0]=="m" ]
        broken_m, = [m for m in cg.defines if m not in mst ]
        #Patch the (inconsistent) cg object
        cg.describe_multiloop = lambda *args: set()
        cg.find_mlonly_multiloops = lambda *args: [ml1, ml2]
        self.verify_virtual_stat(cg, ml1, broken_m)

    def test_get_virtual_stat_not_part_of_mst2(self):
        cg = self.cg
        mst = cg.get_mst()
        ml1, ml2 = [ m for m in mst if m[0]=="m" ]
        broken_m, = [m for m in cg.defines if m not in mst ]
        #Patch the (inconsistent) cg object
        cg.describe_multiloop = lambda *args: set()
        cg.find_mlonly_multiloops = lambda *args: [ml1, ml2]
        self.verify_virtual_stat(cg, ml2, broken_m)

    def test_twice_invert_angle_stat_gives_original_stat(self):
        cg = self.cg
        for ml in [ m for m in cg.get_mst() if m[0]=="m" ]:
            stat, = [ stat for stat in cg.get_stats(ml) if stat.ang_type == cg.get_angle_type(ml) ]

            inverse_stat = ftug.invert_angle_stat(stat)
            # Make sure the inversion changes the stat
            self.assert_all_angles_different(stat, inverse_stat)
            # The second inversion changes it back to the original stat
            self.assert_stats_equal(stat, ftug.invert_angle_stat(inverse_stat))

    def visualize_invert_stat(self):
        import matplotlib.pyplot as plt
        cg = self.cg
        for ml in [ m for m in cg.get_mst() if m[0]=="m" ]:
            stat1, stat2 = cg.get_stats(ml)
            inverse1 = ftug.invert_angle_stat(stat1)
            inverse2 = ftug.invert_angle_stat(stat2)
            sep1 = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, stat1.position_params())
            sep1i = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, inverse1.position_params())
            sep2 = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, stat2.position_params())
            sep2i = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, inverse2.position_params())
            fig, ax = plt.subplots()
            vecs = np.array([sep1, sep1i, sep2, sep2i])
            labels = ["stat1", "inverse1", "stat2", "inverse2"]
            sizes = [10,10,5,5]
            for i in range(4):
                ax.plot(vecs[i,0], vecs[i,1], "o", label=labels[i], markersize = sizes[i])
            ax.plot([0], [0], "s", label="origin")
            ax.legend()
            plt.show()

    def test_invert_angle_stat2(self):
        cg = self.cg
        for ml in [ m for m in cg.get_mst() if m[0]=="m" ]:
            stat1, stat2 = cg.get_stats(ml)
            self.assert_all_angles_different(stat1, stat2)
            self.assert_stats_equal(stat2, ftug.invert_angle_stat(stat1))
            self.assert_stats_equal(stat1, ftug.invert_angle_stat(stat2))

    def test_invert_angle_stat3(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3D0U_A.cg')

        for loop in it.chain(cg.mloop_iterator(), cg.iloop_iterator()):
            log.info("Checking loop %s", loop)
            stat1, stat2 = cg.get_stats(loop)
            inverse1 = ftug.invert_angle_stat(stat1)
            inverse2 = ftug.invert_angle_stat(stat2)

            #self.assert_all_angles_different(stat1, stat2)
            self.assert_stats_equal(stat2, inverse1)
            self.assert_stats_equal(stat1, inverse2)
            self.assert_stats_equal(stat1, ftug.invert_angle_stat(inverse1))
            self.assert_stats_equal(stat2, ftug.invert_angle_stat(inverse2))

    def test_sum_of_stats(self):
        for ml1 in ["m0", "m1", "m2"]:
            cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3way.cg')

            ml2 = cg.get_next_ml_segment(ml1)

            # Get the stat according to the angle type
            at1 = cg.get_angle_type(ml1, allow_broken=True)
            at2 = cg.get_angle_type(ml2, allow_broken=True)
            stat1, = [ stat for stat in cg.get_stats(ml1) if stat.ang_type == at1]
            stat2, = [ stat for stat in cg.get_stats(ml2) if stat.ang_type == at2 ]
            log.error(cg.mst)
            log.error("ml1: %s, %s %s", ml1, cg.get_angle_type(ml1), at1)
            log.error("ml2: %s, %s %s", ml2, cg.get_angle_type(ml2), at2)

            # Special case for angle type 4: a positive sign means it points
            # in the direction OPPOSITE to get_next_ml_segment
            if abs(at1)==4 or abs(at2)==4:
                if at1==-4 or (at1>0 and at1!=4):
                    stat1 = ftug.invert_angle_stat(stat1)
                if at2==-4 or (at2>0 and at2!=4):
                    stat2 = ftug.invert_angle_stat(stat2)
                sum_stat = ftug.sum_of_stats(stat2, stat1)
            else:
                if at1<0:
                    stat1 = ftug.invert_angle_stat(stat1)
                if at2<0:
                    stat2 = ftug.invert_angle_stat(stat2)
                sum_stat = ftug.sum_of_stats(stat1, stat2)

            virtual_stat = ftug.get_virtual_stat(cg, ml1, ml2)
            self.assert_stats_equal(sum_stat, virtual_stat)
            self.assert_stats_equal(sum_stat, virtual_stat)

    def test_sum_of_stats_standard_direction(self):
        for ml1 in ["m0", "m1", "m2"]:
            cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3way.cg')

            ml2 = cg.get_next_ml_segment(ml1)
            ml3 = cg.get_next_ml_segment(ml2)



            # Get the stat according to the angle type
            at1 = cg.get_angle_type(ml1, allow_broken=True)
            at2 = cg.get_angle_type(ml2, allow_broken=True)

            log.error(cg.mst)
            log.error("ml1: %s, %s %s", ml1, cg.get_angle_type(ml1), at1)
            log.error("ml2: %s, %s %s", ml2, cg.get_angle_type(ml2), at2)


            stat1, = [ stat for stat in cg.get_stats(ml1) if stat.ang_type == at1]
            stat2, = [ stat for stat in cg.get_stats(ml2) if stat.ang_type == at2 ]

            sum_stat = ftug.sum_of_stat_in_standard_direction(stat1, stat2)

            virtual_stat = ftug.get_virtual_stat(cg, ml1, ml2)
            self.assert_stats_equal(sum_stat, virtual_stat)
            # The sum of the two stats in standard direction is independent of the build order
            stat3, = [ stat for stat in cg.get_stats(ml3) if stat.ang_type >0 ]

            self.assert_stats_equal(sum_stat, stat3)

    def test_sum_of_stats_true_RNA(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3D0U_A.cg')
        cg.print_debug(logging.INFO)
        for ml1 in cg.mloop_iterator():
            ml2 = cg.get_next_ml_segment(ml1)
            stat1 = cg.get_stats(ml1)[0]
            stat2 = cg.get_stats(ml2)[0]
            virtual_stat = ftug.get_virtual_stat(cg, ml1, ml2)
            sum_stat = ftug.sum_of_stat_in_standard_direction(stat1, stat2)
            log.info("Stats for %s and %s with angle types %s and %s", ml1, ml2, stat1.ang_type, stat2.ang_type)
            self.assert_stats_equal(sum_stat, virtual_stat)

    def visualize_sum_of_stats(self):
        cg = self.cg
        import matplotlib.pyplot as plt
        stat0 = cg.get_stats("m0")[0]
        stat1 = cg.get_stats("m2")[0]
        stat2 = cg.get_stats("m1")[0]


        sep0 = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, stat0.position_params())
        stemvec1 = ftug.stem2_orient_from_stem1_1(ftuv.standard_basis, [1]+list(stat0.orientation_params()))
        twist1 = ftug.twist2_orient_from_stem1_1(ftuv.standard_basis, stat0.twist_params())
        sep1 = ftug.stem2_pos_from_stem1(-stemvec1, twist1, stat1.position_params())
        sep2 = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, stat2.position_params())
        sepSum = ftug.stem2_pos_from_stem1_1(ftuv.standard_basis,
                                             ftug.sum_of_stat_in_standard_direction(stat0, stat1).position_params())

        fig, ax = plt.subplots()
        ax.plot([0,sep0[0]], [0,sep0[1]], "-o", label="m0")
        ax.plot([sep0[0], sep0[0]+sep1[0]], [sep0[1], sep0[1]+sep1[1]], "-o", label="m2")
        ax.plot([0,sep2[0]], [0,sep2[1]], "-o", label="m1")
        ax.plot([0,sepSum[0]], [0,sepSum[1]], "--", label="m0+m2")

        ax.legend()
        plt.show()


class TestDistanceCalculation(unittest.TestCase):
    def setUp(self):
        self.rs_random_281=ftmc.from_pdb('test/forgi/threedee/data/RS_random_281_S_0.pdb')
        for key in self.rs_random_281.defines.keys():
          if key[0] =="s":
            ftug.add_virtual_residues(self.rs_random_281, key)
        self.minimal_multiloop = ftmc.CoarseGrainRNA()
        self.minimal_multiloop.from_file('test/forgi/threedee/data/minimal_multiloop.cg')
        for key in self.minimal_multiloop.defines.keys():
          if key[0] =="s":
            ftug.add_virtual_residues(self.minimal_multiloop, key)

    def test_junction_virtual_atom_distance_minimalMultiloop(self):
        distance1=ftug.junction_virtual_atom_distance(self.minimal_multiloop, "m0")
        self.assertLess(distance1, 4., msg="{} is not < {} for {}".format(distance1, 4., "m0"))
        distance2=ftug.junction_virtual_atom_distance(self.minimal_multiloop, "m1")
        self.assertLess(distance2, 4., msg="{} is not < {} for {}".format(distance2, 4., "m1"))
    def test_junction_virtual_atom_distance_realPDB(self):
        distance=ftug.junction_virtual_atom_distance(self.rs_random_281, "m3")
        self.assertLess(distance, 4.)  #This should always hold!
        self.assertAlmostEqual(distance, 3.4294889373610675) #This value might change, if we change the virtual atom calculation


class TestAtomPosition_VirtualAtoms(unittest.TestCase):
    def setUp(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')
        # cg.defines['s0']==[1,9,63,71]
        self.va1=ftug.virtual_atoms(cg, sidechain=False)
        self.va2=ftug.virtual_atoms(cg, sidechain=True)

    def test_virtual_atoms_num_atom_per_nucleotide(self):
        """ Number of atoms for each nucleotide"""
        for i in range(1,65):
            # non-sidechain atoms
            self.assertEqual(len(self.va1[1].keys()), 9)
            #Side-chain atoms depend on nucleotide
            self.assertGreaterEqual(len(self.va2[1].keys()), 15)


    def test_virtual_atoms_intranucleotide_distances_stem_nosidechain(self):
        """ distance between two atoms of same nucleotide IN STEM """
        for i in it.chain(range(1,10),range(63,72)):
            for k1, a1 in self.va1[i].items():
                for k2, a2 in self.va1[i].items():
                    if k1==k2: continue
                    dist=ftuv.magnitude(a1-a2)
                    self.assertLess(dist, 20, msg="Nucleotide {} too big: Distance between "
                                                  "{} and {} is {}".format(i, k1, k2, dist) )
                    self.assertGreater(dist, 0.8, msg="Nucleotide {} too small: Distance between "
                                                  "{} and {} is {}".format(i, k1, k2, dist) )
    def test_virtual_atoms_intranucleotide_distances_stem_withsidechain(self):
        """ distance between two atoms of same nucleotide IN STEM """
        for i in it.chain(range(1,10),range(63,72)):
            for k1, a1 in self.va2[i].items():
                for k2, a2 in self.va2[i].items():
                    if k1==k2: continue
                    dist=ftuv.magnitude(a1-a2)
                    self.assertLess(dist, 30, msg="Nucleotide {} too big: Distance between "
                                                  "{} and {} is {}".format(i, k1, k2, dist) )
                    self.assertGreater(dist, 0.8, msg="Nucleotide {} too small: Distance between "
                                                  "{} and {} is {}".format(i, k1, k2, dist) )
    def test_virtual_atoms_basepairdistance_in_stem(self):
        """ distance between two atoms that pair in stem """
        for i in range(1,10):
            for j in range (71,62):
                mindist=min(ftuv.magnitude(a1-a2) for a1 in self.va2[i].values() for a2 in self.va2[j].values())
                self.assertLess(dist, 20, msg="Distance between nucleotide {} and {} is too big: "
                                              "the minimal distance is {}".format(i, j, mindist))
                self.assertGreater(dist, 1.1, msg="Distance between nucleotide {} and {} is too small: "
                                              "the minimal distance is {}".format(i, j, mindist))
    def test_virtual_atoms_same_strand_nuc_distance(self):
        """ Distance between virtual atoms of nucleotides on same strand in stem"""
        for i in range(1,9):
            for j in range (i+1,10):
                dist=ftuv.magnitude(self.va1[i]["C1'"]-self.va1[j]["C1'"])
                self.assertLess(dist, 2.2+4.5*(j-i))#, msg="Distance between nucleotide {} and {} "
                                                  #   "is too big: {}".format(i, j, dist))
                self.assertGreater(dist, 2*(j-i))#, msg="Distance between nucleotide {} and {} "
                                                  #   "is too small: {}".format(i, j, dist))
        for i in range(63,71):
            for j in range (i+1,72):
                dist=ftuv.magnitude(self.va1[i]["C1'"]-self.va1[j]["C1'"])
                self.assertLess(dist, 2.2+4.5*(j-i))#, msg="Distance between nucleotide {} and {} "
                                                   #  "is too big: {}".format(i, j, dist))
                self.assertGreater(dist, 2*(j-i))#, msg="Distance between nucleotide {} and {} "
                                                    # "is too small: {}".format(i, j, dist))
    def test_virtual_atoms_distance_neighboring_atoms_in_nucleotide(self):
        # C2' is next to C3'
        for i in range(1,9):
            dist=ftuv.magnitude(self.va1[i]["C3'"]-self.va1[i]["C2'"])
            self.assertLess(dist, 3, msg="Distance between C3' and C2' for nucleotide {} "
                                         "is too big: {}".format(i, dist))
            self.assertGreater(dist, 1, msg="Distance between C3' and C2' for nucleotide {} "
                                         "is too small: {}".format(i, dist))
