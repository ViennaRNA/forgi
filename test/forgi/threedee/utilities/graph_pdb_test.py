
import Bio.PDB as bpdb
import unittest, os
import warnings
import numpy as np
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv

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

    def test_virtual_residue_atoms(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        ftug.add_virtual_residues(cg, 's0')
        va = ftug.virtual_residue_atoms(cg, 's0', 1, 0)

    def test_virtual_atoms(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        ftug.virtual_atoms(cg, sidechain=False)

    def test_numbered_virtual_residues(self): 
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26.pdb')

        nres = ftug.numbered_virtual_residues(cg)
        fud.pv('nres')


