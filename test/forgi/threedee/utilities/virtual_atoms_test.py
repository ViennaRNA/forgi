import unittest
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.rmsd as ftur
import Bio.PDB as bp
import warnings
import numpy as np

def realatom_vatom_rmsd(cg, chain):
    """
    The RMSD between the real atoms and the virtual atoms of the stems.
    """
    vposs=[]
    rposs=[]
    for stem in cg.stem_iterator():
        for pos in cg.define_residue_num_iterator(stem):
            vas=cg.virtual_atoms(pos)
            for a, coords in vas.items():
                try:
                    rposs.append(chain[pos][a].get_vector().get_array())
                except KeyError: pass
                else:
                    vposs.append(coords)
    return ftur.rmsd(np.array(vposs), np.array(rposs))

class VirtualAtomsTest(unittest.TestCase):
    def setUp(self):
        pdbfile1 = "test/forgi/threedee/data/3FU2.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pdb1 = list(bp.PDBParser().get_structure('test', pdbfile1).get_chains())[0]
        self.cg1  = ftmc.from_pdb(pdbfile1)

        pdbfile2 = "test/forgi/threedee/data/3V2F.pdb" #Takes some time. Big structure
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pdb2 = list(bp.PDBParser().get_structure('test', pdbfile2).get_chains())[0]
        self.cg2  = ftmc.from_pdb(pdbfile1)

    def test_virtualatoms_close_to_pdb(self):
        self.assertLess(realatom_vatom_rmsd(self.cg2, self.pdb2), 1.0) #1.016386 with average_atom_positions
        self.assertLess(realatom_vatom_rmsd(self.cg1, self.pdb1), 2.5) # 2.2 with average_atom_positions.
        # Conclusion: Using virtual residues improves the stem virtual atom positions
        # on average for the bigger structure self.pdb2, as it takes twists into account.
        # ==> Stick with using virtual residues for stems.
