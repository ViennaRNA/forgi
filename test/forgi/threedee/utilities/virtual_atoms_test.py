from __future__ import print_function
import unittest
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.similarity as ftme
import forgi.threedee.utilities.pdb as ftup
import Bio.PDB as bp
import warnings
import itertools as it
import numpy as np
import sys

def realatom_vatom_rmsd(cg):
    """
    The RMSD between the real atoms and the virtual atoms of the stems.
    """
    vposs=[]
    rposs=[]
    #print(cg.seq)
    #print("".join(r.get_segid() for r in chain.get_residues()))
    #print(cg.seq[1], chain[2])
    for stem in cg.stem_iterator():
        stemv=[]
        stemr=[]
        resseqs=cg.get_resseqs(stem)
        resseqs=resseqs[0]+resseqs[1]
        for posCg, posPdb in zip(cg.define_residue_num_iterator(stem), resseqs):
            #try:
            #    print( posCg, posPdb, cg.seq[posCg-1], cg.chain[posPdb] )
            #except KeyError as e: print(e)
            vas=cg.virtual_atoms(posCg)
            for a, coords in vas.items():
                try:
                    rposs.append(cg.chain[posPdb][a].get_vector().get_array())
                    stemr.append(cg.chain[posPdb][a].get_vector().get_array())
                except KeyError: 
                    pass
                else:
                    vposs.append(coords)
                    stemv.append(coords)
    assert len(vposs)==len(rposs)
    #print (np.array(vposs).shape, file=sys.stderr)
    return ftme.rmsd(np.array(vposs), np.array(rposs))

class VirtualAtomsTest(unittest.TestCase):
    def setUp(self):
        pdbfile1 = "test/forgi/threedee/data/3FU2.pdb"
        self.cg1  = ftmc.from_pdb(pdbfile1)

        pdbfile2 = "test/forgi/threedee/data/3V2F.pdb"  #Takes some time. Big structure
        self.cg2  = ftmc.from_pdb(pdbfile2)

        pdbfile3 = "test/forgi/threedee/data/1X8W.pdb"
        self.cg3  = ftmc.from_pdb(pdbfile3)

        pdbfile4 = "test/forgi/threedee/data/2QBZ.pdb" 
        self.cg4  = ftmc.from_pdb(pdbfile4)

    def test_virtualatoms_close_to_pdb(self):
        #print(realatom_vatom_rmsd(self.cg1), file=sys.stderr)
        #print(realatom_vatom_rmsd(self.cg2), file=sys.stderr)
        #print(realatom_vatom_rmsd(self.cg3), file=sys.stderr)
        #print(realatom_vatom_rmsd(self.cg4), file=sys.stderr)
        self.assertLess(realatom_vatom_rmsd(self.cg1), 1.16) # 1.32 using average atoms (no vres)
        self.assertLess(realatom_vatom_rmsd(self.cg3), 1.36) # 1.33
        self.assertLess(realatom_vatom_rmsd(self.cg4), 1.32) # 1.16
        self.assertLess(realatom_vatom_rmsd(self.cg2), 1.42) # 1.30



