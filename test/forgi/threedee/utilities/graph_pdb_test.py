
import Bio.PDB as bpdb
import unittest, os
import warnings
import numpy as np
import forgi.threedee.model.coarse_grain as ftmc

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv
import forgi.graph.bulge_graph as fgb

import forgi.utilities.debug as fud

class TestGraphPDB(unittest.TestCase):
    '''
    Test some of the rmsd-type functions.
    '''
    def setUp(self):
        pass

    def test_add_loop_information_from_pdb_chain(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1A34.pdb')

    def test_add_stem_information_from_pdb_chain(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1gid.cg')
        pdb_filename = 'test/forgi/threedee/data/1gid_rosetta.pdb'

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s = bpdb.PDBParser().get_structure('temp', pdb_filename)
            chain = list(s.get_chains())[0]

        chain = ftup.renumber_chain(chain, cg.seq_ids)
        ftug.add_stem_information_from_pdb_chain(cg, chain)


    def verify_virtual_twist_angles(self, cg, s):
        sl = cg.stem_length(s)

        for i in range(0, sl):
            (pos, vec, vec_l, vec_r) = ftug.virtual_res_3d_pos_core(cg.coords[s],
                                                                    cg.twists[s],i,sl)

            if i > 1:
                self.assertGreater(ftuv.vec_angle(vec, prev_vec), 0.1)
                self.assertLess(ftuv.vec_angle(vec, prev_vec), 0.95)

            prev_vec = vec

    def test_angle_between_twists(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        self.verify_virtual_twist_angles(cg, 's2')
        self.verify_virtual_twist_angles(cg, 's0')


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
    def test_virtual_atoms(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        ftug.virtual_atoms(cg, sidechain=False)

    def test_numbered_virtual_residues(self): 
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        nres = ftug.numbered_virtual_residues(cg)
        #fud.pv('nres')
        #TODO assert something

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
        self.assertLess(distance2, 4., msg="{} is not < {} for {}".format(distance1, 4., "m1"))
    def test_junction_virtual_atom_distance_realPDB(self):
        distance=ftug.junction_virtual_atom_distance(self.rs_random_281, "m4")
        self.assertLess(distance, 4.)
       











