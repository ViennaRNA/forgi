from __future__ import print_function
import unittest
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.rmsd as ftur
import Bio.PDB as bp
import warnings
import numpy as np
import sys

def realatom_vatom_rmsd(cg, chain):
    """
    The RMSD between the real atoms and the virtual atoms of the stems.
    """
    vposs=[]
    rposs=[]
    for stem in cg.stem_iterator():
        stemv=[]
        stemr=[]
        for pos in cg.define_residue_num_iterator(stem):
            vas=cg.virtual_atoms(pos)
            for a, coords in vas.items():
                try:
                    rposs.append(chain[pos][a].get_vector().get_array())
                    stemr.append(chain[pos][a].get_vector().get_array())
                except KeyError: pass
                else:
                    vposs.append(coords)
                    stemv.append(coords)
        #if ftur.rmsd(np.array(stemv), np.array(stemr))>40:
        #    print(stem, np.array(stemv), np.array(stemr))
    assert len(vposs)==len(rposs)
    #print (np.array(vposs).shape, file=sys.stderr)
    return ftur.rmsd(np.array(vposs), np.array(rposs))

class VirtualAtomsTest(unittest.TestCase):
    def setUp(self):
        pdbfile1 = "test/forgi/threedee/data/3FU2.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pdb1 = list(bp.PDBParser().get_structure('test', pdbfile1).get_chains())[0]
        self.cg1  = ftmc.from_pdb(pdbfile1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pdb1b = list(bp.PDBParser().get_structure('test', pdbfile1).get_chains())[0]
        self.cg1b  = ftmc.from_pdb(pdbfile1)
        #pdbfile2 = "test/forgi/threedee/data/3V2F.pdb"  #Takes some time. Big structure
        #with warnings.catch_warnings():
        #    warnings.simplefilter("ignore")
        #    self.pdb2 = list(bp.PDBParser().get_structure('test', pdbfile2).get_chains())[0]
        #self.cg2  = ftmc.from_pdb(pdbfile2)

        pdbfile3 = "test/forgi/threedee/data/1X8W.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pdb3 = list(bp.PDBParser().get_structure('test', pdbfile3).get_chains())[0]
        self.cg3  = ftmc.from_pdb(pdbfile3)

        pdbfile4 = "test/forgi/threedee/data/2QBZ.pdb" 
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pdb4 = list(bp.PDBParser().get_structure('test', pdbfile4).get_chains())[0]
        self.cg4  = ftmc.from_pdb(pdbfile4)

    def test_virtualatoms_close_to_pdb(self):
        self.assertLess(realatom_vatom_rmsd(self.cg4, self.pdb4), 32) #31.97
        self.assertLess(realatom_vatom_rmsd(self.cg3, self.pdb3), 29.3) #22.23
        #self.assertLess(realatom_vatom_rmsd(self.cg2, self.pdb2), 60)
        self.assertLess(realatom_vatom_rmsd(self.cg1, self.pdb1), 2.39) #2.21 without virtual residues
        #The difference between average_virtual_atoms and stem-based virtual atoms is negelctable.
    def test_va_rmsd_always_the_same(self):
        r1=realatom_vatom_rmsd(self.cg1, self.pdb1)
        r2=realatom_vatom_rmsd(self.cg1, self.pdb1)
        self.assertEqual(r1,r2)
        r1=realatom_vatom_rmsd(self.cg1, self.pdb1)
        r2=realatom_vatom_rmsd(self.cg1b, self.pdb1b)
        self.assertEqual(r1,r2)
        r1=realatom_vatom_rmsd(self.cg1, self.pdb1b)
        r2=realatom_vatom_rmsd(self.cg1b, self.pdb1)
        self.assertEqual(r1,r2)
        r1=realatom_vatom_rmsd(self.cg1, self.pdb1)
        r2=realatom_vatom_rmsd(self.cg1b, self.pdb1)
        self.assertEqual(r1,r2)
