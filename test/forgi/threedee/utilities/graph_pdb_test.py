from __future__ import division

from builtins import map
from builtins import range
import Bio.PDB as bpdb
import unittest
import os
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
log = logging.getLogger(__name__)


class TestBrokenMlDeviation(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1GID_A.cg")

    def test_zero_deviation(self):
        # m2 is broken
        stats = self.cg.get_stats("m2")
        dev = ftug.get_broken_ml_deviation(self.cg, "m2", "s6", stats[0])
        dev2 = ftug.get_broken_ml_deviation(self.cg, "m2", "s6", stats[1])
        print(dev)
        print(dev2)
        self.assertLess(dev[0], 10**-3)
        self.assertLess(dev[1], 10**-3)
        self.assertLess(dev[2], 10**-3)


class TestGraphPDB(unittest.TestCase):
    '''
    Test some of the rmsd-type functions.
    '''

    def setUp(self):
        pass

    @unittest.skip("Old version of incomplete elements")
    def test_get_incomplete_elements(self):
        db = "(((...(((...).)))))"
        cg = ftmc.CoarseGrainRNA.from_dotbracket(db)
        # Residue number 3 is missing. Probably bulged out from stem 0
        seq_ids = map(str, [1, 2, 4, 5, 6, 7, 8, 9, 10, 11,
                            12, 13, 14, 15, 16, 17, 18, 19, 20])
        cg.seq._seqids = list(map(fgb.resid_from_str, seq_ids))
        cg.seq._set_missing_residues(
            [{"model": None, "ssseq": 3, "res_name": "G", "chain": "A", "insertion": ""}])
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["s0"]))
        # Residue number 4 is missing between s0 and i0
        seq_ids = map(str, [1, 2, 3, 5, 6, 7, 8, 9, 10, 11,
                            12, 13, 14, 15, 16, 17, 18, 19, 20])
        cg.seq._set_missing_residues(
            [{"model": None, "ssseq": 4, "res_name": "G", "chain": "A", "insertion": ""}])
        cg.seq._seqids = list(map(fgb.resid_from_str, seq_ids))
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["i0"]))
        # Residue number 5 is missing inside i0
        seq_ids = map(str, [1, 2, 3, 4, 6, 7, 8, 9, 10, 11,
                            12, 13, 14, 15, 16, 17, 18, 19, 20])
        cg.seq._seqids = list(map(fgb.resid_from_str, seq_ids))
        cg.seq._set_missing_residues(
            [{"model": None, "ssseq": 5, "res_name": "G", "chain": "A", "insertion": ""}])
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["i0"]))
        # Residue number 17 is missing between s0 and s1, ==> i0
        seq_ids = map(str, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                            11, 12, 13, 14, 15, 16, 18, 19, 20])
        cg.seq._seqids = list(map(fgb.resid_from_str, seq_ids))
        cg.seq._set_missing_residues(
            [{"model": None, "ssseq": 17, "res_name": "G", "chain": "A", "insertion": ""}])
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["i0"]))
        # Residue number 10 is missing  ==> h0
        seq_ids = map(str, [1, 2, 3, 4, 5, 6, 7, 8, 9, 11,
                            12, 13, 14, 15, 16, 17, 18, 19, 20])
        cg.seq._seqids = list(map(fgb.resid_from_str, seq_ids))
        cg.seq._set_missing_residues(
            [{"model": None, "ssseq": 10, "res_name": "G", "chain": "A", "insertion": ""}])
        self.assertEqual(ftug.get_incomplete_elements(cg), set(["h0"]))
        # Multiple residues are missing
        seq_ids = map(str, [1, 2, 4, 6, 7, 8, 9, 10, 11,
                            13, 14, 15, 16, 17, 18, 19, 20, 21, 22])
        cg.seq._seqids = list(map(fgb.resid_from_str, seq_ids))
        cg.seq._set_missing_residues([
            {"model": None, "ssseq": 3, "res_name": "G",
                "chain": "A", "insertion": ""},
            {"model": None, "ssseq": 5, "res_name": "G",
                "chain": "A", "insertion": ""},
            {"model": None, "ssseq": 12, "res_name": "G",
                "chain": "A", "insertion": ""}
        ])
        self.assertEqual(ftug.get_incomplete_elements(cg),
                         set(["h0", "s0", "i0"]))

    def test_add_loop_information_from_pdb_chain(self):
        cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1A34.pdb')

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
                                                                    cg.twists[s], i, sl)

            if i == 0:
                nptest.assert_array_equal(pos, cg.coords[s][0])

            if i > 1:
                self.assertGreater(ftuv.vec_angle(vec, prev_vec), 0.53)
                self.assertLess(ftuv.vec_angle(vec, prev_vec), 0.73)

            prev_vec = vec

    def test_first_virtual_res_basis(self):
        cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1y26.pdb')
        basis = ftug.virtual_res_basis(cg, "s0", 0)
        nptest.assert_array_equal(basis, ftuv.create_orthonormal_basis(
            cg.coords["s0"][1] - cg.coords["s0"][0], cg.twists["s0"][0]))

    def test_angle_between_twists(self):
        cg, = ftmc.CoarseGrainRNA.from_pdb(
            'test/forgi/threedee/data/1y26.pdb', dissolve_length_one_stems=False)

        self.verify_virtual_twist_angles(cg, 's0')
        self.verify_virtual_twist_angles(cg, 's1')
        self.verify_virtual_twist_angles(cg, 's2')
        self.verify_virtual_twist_angles(cg, 's3')

    def test_coordinates_for_add_virtual_residues(self):
        try:
            cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1y26.pdb', annotation_tool="MC-Annotate")
        except ftmc.AnnotationToolNotInstalled:
            self.skipTest("This requires MC-Annotate")

        ftug.add_virtual_residues(cg, 's0')
        # XYZ coordinate for first residue are ok:
        self.assertAlmostEqual(cg.vposs["s0"][0][0], 2.3, delta=3,
                               msg="Wrong x-position for virtual residue 0 of stem s0: {}".format(cg.vposs["s0"][0][0]))
        self.assertAlmostEqual(cg.vposs["s0"][0][1], 1.3, delta=3,
                               msg="Wrong y-position for virtual residue 0 of stem s0: {}".format(cg.vposs["s0"][0][1]))
        self.assertAlmostEqual(cg.vposs["s0"][0][2], 1.0, delta=3,
                               msg="Wrong z-position for virtual residue 0 of stem s0: {}".format(cg.vposs["s0"][0][2]))
        last_residue = cg.stem_length("s0") - 1
        self.assertAlmostEqual(cg.vposs["s0"][last_residue][0], 16, delta=4,
                               msg="Wrong x-position for virtual residue {} of stem s0: {}".format(last_residue, cg.vposs["s0"][last_residue][0]))
        self.assertAlmostEqual(cg.vposs["s0"][last_residue][1], -13, delta=4,
                               msg="Wrong y-position for virtual residue {} of stem s0: {}".format(last_residue, cg.vposs["s0"][last_residue][1]))
        self.assertAlmostEqual(cg.vposs["s0"][last_residue][2], 8, delta=4,
                               msg="Wrong z-position for virtual residue {} of stem s0: {}".format(last_residue, cg.vposs["s0"][last_residue][2]))

    def test_basis_transformation_for_virtual_residues(self):
        cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1y26.pdb')
        ftug.add_virtual_residues(cg, 's0')
        offset = cg.vposs["s0"][0]
        vbasis = cg.vbases["s0"][0]
        local_pos = np.array([0, 0, 1])
        global_pos = np.dot(vbasis.transpose(), local_pos) + offset
        # Checking dimensions of vectors, just in case...
        self.assertEqual(len(global_pos), 3)
        self.assertEqual(len(vbasis), 3)
        self.assertEqual(len(vbasis[0]), 3)
        # Analytically true:
        self.assertTrue(all(global_pos[x] - vbasis[2][x] - offset[x] < 0.0000001 for x in [0, 1, 2]),
                        msg="global pos for (0,0,1) should be {}+{}={}, but is {} instead.".format(
            vbasis[2], offset, vbasis[2] + offset, global_pos))

    def test_virtual_residue_atoms(self):
        cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1y26.pdb')

        ftug.add_virtual_residues(cg, 's0')
        ftug.add_virtual_residues(cg, 's1')
        bases_to_test = []
        bases_to_test.append(ftug.virtual_residue_atoms(cg, 's0', 1, 0))
        bases_to_test.append(ftug.virtual_residue_atoms(cg, 's0', 2, 1))
        bases_to_test.append(ftug.virtual_residue_atoms(cg, 's1', 0, 0))

        # Assert that any two atoms of the same base within reasonable distance to each other
        #(https://en.wikipedia.org/wiki/Bond_length says than a CH-bond is >= 1.06A)
        for va in bases_to_test:
            for k1, v1 in va.items():
                for k2, v2 in va.items():
                    dist = ftuv.magnitude(v1 - v2)
                    self.assertLess(dist, 30, msg="Nucleotide too big: "
                                    "Distance between {} and {} is {}".format(k1, k2, dist))
                    if k1 != k2:
                        dist = ftuv.magnitude(v1 - v2)
                        self.assertGreater(dist, 0.8, msg="Nucleotide too small: "
                                           "Distance between {} and {} is {}".format(k1, k2, dist))

    def test_virtual_residue_atom_exact_match(self):
        # This test serves to detect unwanted changes in the virtual atom calculation algorithm.
        # It is allowed to fail, if the virtual atom calculation changes.
        try:
            cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1y26.pdb',
                                               annotation_tool="MC-Annotate")
        except ftmc.AnnotationToolNotInstalled:
            self.skipTest("Requires MC-Annotate")
        ftug.add_virtual_residues(cg, 's0')
        vres = ftug.virtual_residue_atoms(cg, 's0', 1, 0)
        print(vres['C8'])
        print(vres['N2'])
        nptest.assert_allclose(vres['C8'], np.array(
            [5.17130963, -2.93616933, -2.16725865]))
        nptest.assert_allclose(vres['N2'], np.array(
            [6.91170703,  2.35440312, -1.92080436]))


class TestOrientation(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_stem_orientation_parameters(self):
        stem1_vec = np.array([0., 0., 1.])
        twist1 = np.array([0., 1., 0.])
        stem2_vec = np.array([0., 0., 1.])
        twist2 = np.array([1., 0., 0.])
        r, u, v, t = ftug.get_stem_orientation_parameters(
            stem1_vec, twist1, stem2_vec, twist2)
        self.assertEqual(r, 1)
        self.assertEqual(u, math.pi / 2)
        self.assertEqual(v, 0)

        stem2_vec = np.array([0., 1., 0.])
        r, u, v, t = ftug.get_stem_orientation_parameters(
            stem1_vec, twist1, stem2_vec, twist2)
        self.assertEqual(r, 1)
        self.assertEqual(u, math.pi / 2)
        self.assertEqual(v, math.pi / 2)

        stem2_vec = np.array([0., 1., 1.])
        r, u, v, t = ftug.get_stem_orientation_parameters(
            stem1_vec, twist1, stem2_vec, twist2)
        self.assertEqual(r, math.sqrt(2))
        self.assertEqual(u, math.pi / 2)
        self.assertEqual(v, math.pi / 4)

    def test_stem2_orient_from_stem1(self):
        stem1_vec = np.array([0., 0., 1.])
        twist1 = np.array([0., 1., 0.])
        r = 1
        u = math.pi / 2
        v = 0
        nptest.assert_allclose(ftug.stem2_orient_from_stem1(
            stem1_vec, twist1, (r, u, v)), np.array([0., 0, 1]), atol=10**-15)
        v = math.pi / 2
        nptest.assert_allclose(ftug.stem2_orient_from_stem1(
            stem1_vec, twist1, (r, u, v)), np.array([0., 1, 0]), atol=10**-15)
        r = math.sqrt(2)
        v = math.pi / 4
        nptest.assert_allclose(ftug.stem2_orient_from_stem1(
            stem1_vec, twist1, (r, u, v)), np.array([0., 1, 1]), atol=10**-15)


class TestDistanceCalculation(unittest.TestCase):
    def setUp(self):
        try:
            self.rs_random_281, = ftmc.CoarseGrainRNA.from_pdb(
                'test/forgi/threedee/data/RS_random_281_S_0.pdb', annotation_tool="MC-Annotate")
        except ftmc.AnnotationToolNotInstalled:
            self.skipTest("Requires MC-Annotate")

        for key in self.rs_random_281.defines.keys():
            if key[0] == "s":
                ftug.add_virtual_residues(self.rs_random_281, key)
        self.minimal_multiloop = ftmc.CoarseGrainRNA.from_bg_file(
            'test/forgi/threedee/data/minimal_multiloop.cg')
        for key in self.minimal_multiloop.defines.keys():
            if key[0] == "s":
                ftug.add_virtual_residues(self.minimal_multiloop, key)

    def test_junction_virtual_atom_distance_minimalMultiloop(self):
        distance1 = ftug.junction_virtual_atom_distance(
            self.minimal_multiloop, "m0")
        self.assertLess(
            distance1, 4., msg="{} is not < {} for {}".format(distance1, 4., "m0"))
        distance2 = ftug.junction_virtual_atom_distance(
            self.minimal_multiloop, "m1")
        self.assertLess(
            distance2, 4., msg="{} is not < {} for {}".format(distance2, 4., "m1"))

    def test_junction_virtual_atom_distance_realPDB(self):
        distance = ftug.junction_virtual_atom_distance(
            self.rs_random_281, "m3")
        self.assertLess(distance, 4.)  # This should always hold!
        # This value might change, if we change the virtual atom calculation
        self.assertAlmostEqual(distance, 3.486261134885911)


class TestAtomPosition_VirtualAtoms(unittest.TestCase):
    def setUp(self):
        try:
            cg, = ftmc.CoarseGrainRNA.from_pdb('test/forgi/threedee/data/1y26.pdb', annotation_tool="MC-Annotate")
        except ftmc.AnnotationToolNotInstalled:
            self.skipTest("Requires MC-Annotate")

        # cg.defines['s0']==[1,9,63,71]
        self.va1 = ftug.virtual_atoms(cg, sidechain=False)
        self.va2 = ftug.virtual_atoms(cg, sidechain=True)

    def test_virtual_atoms_num_atom_per_nucleotide(self):
        """ Number of atoms for each nucleotide"""
        for i in range(1, 65):
                # non-sidechain atoms
            self.assertEqual(len(self.va1[1].keys()), 9)
            # Side-chain atoms depend on nucleotide
            self.assertGreaterEqual(len(self.va2[1].keys()), 15)

    def test_virtual_atoms_intranucleotide_distances_stem_nosidechain(self):
        """ distance between two atoms of same nucleotide IN STEM """
        for i in it.chain(range(1, 10), range(63, 72)):
            for k1, a1 in self.va1[i].items():
                for k2, a2 in self.va1[i].items():
                    if k1 == k2:
                        continue
                    dist = ftuv.magnitude(a1 - a2)
                    self.assertLess(dist, 20, msg="Nucleotide {} too big: Distance between "
                                                  "{} and {} is {}".format(i, k1, k2, dist))
                    self.assertGreater(dist, 0.8, msg="Nucleotide {} too small: Distance between "
                                       "{} and {} is {}".format(i, k1, k2, dist))

    def test_virtual_atoms_intranucleotide_distances_stem_withsidechain(self):
        """ distance between two atoms of same nucleotide IN STEM """
        for i in it.chain(range(1, 10), range(63, 72)):
            for k1, a1 in self.va2[i].items():
                for k2, a2 in self.va2[i].items():
                    if k1 == k2:
                        continue
                    dist = ftuv.magnitude(a1 - a2)
                    self.assertLess(dist, 30, msg="Nucleotide {} too big: Distance between "
                                                  "{} and {} is {}".format(i, k1, k2, dist))
                    self.assertGreater(dist, 0.8, msg="Nucleotide {} too small: Distance between "
                                       "{} and {} is {}".format(i, k1, k2, dist))

    def test_virtual_atoms_basepairdistance_in_stem(self):
        """ distance between two atoms that pair in stem """
        for i in range(1, 10):
            for j in range(71, 62):
                mindist = min(ftuv.magnitude(a1 - a2)
                              for a1 in self.va2[i].values() for a2 in self.va2[j].values())
                self.assertLess(dist, 20, msg="Distance between nucleotide {} and {} is too big: "
                                              "the minimal distance is {}".format(i, j, mindist))
                self.assertGreater(dist, 1.1, msg="Distance between nucleotide {} and {} is too small: "
                                   "the minimal distance is {}".format(i, j, mindist))

    def test_virtual_atoms_same_strand_nuc_distance(self):
        """ Distance between virtual atoms of nucleotides on same strand in stem"""
        for i in range(1, 9):
            for j in range(i + 1, 10):
                dist = ftuv.magnitude(self.va1[i]["C1'"] - self.va1[j]["C1'"])
                # , msg="Distance between nucleotide {} and {} "
                self.assertLess(dist, 2.2 + 4.5 * (j - i))
                #   "is too big: {}".format(i, j, dist))
                # , msg="Distance between nucleotide {} and {} "
                self.assertGreater(dist, 2 * (j - i))
                #   "is too small: {}".format(i, j, dist))
        for i in range(63, 71):
            for j in range(i + 1, 72):
                dist = ftuv.magnitude(self.va1[i]["C1'"] - self.va1[j]["C1'"])
                # , msg="Distance between nucleotide {} and {} "
                self.assertLess(dist, 2.2 + 4.5 * (j - i))
                #  "is too big: {}".format(i, j, dist))
                # , msg="Distance between nucleotide {} and {} "
                self.assertGreater(dist, 2 * (j - i))
                # "is too small: {}".format(i, j, dist))

    def test_virtual_atoms_distance_neighboring_atoms_in_nucleotide(self):
        # C2' is next to C3'
        for i in range(1, 9):
            dist = ftuv.magnitude(self.va1[i]["C3'"] - self.va1[i]["C2'"])
            self.assertLess(dist, 3, msg="Distance between C3' and C2' for nucleotide {} "
                                         "is too big: {}".format(i, dist))
            self.assertGreater(dist, 1, msg="Distance between C3' and C2' for nucleotide {} "
                               "is too small: {}".format(i, dist))
