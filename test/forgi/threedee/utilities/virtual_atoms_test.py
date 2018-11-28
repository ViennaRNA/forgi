from __future__ import print_function
from builtins import zip
import unittest
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.similarity as ftme
import forgi.threedee.utilities.pdb as ftup
import Bio.PDB as bp
import warnings
import itertools as it
import numpy as np
import sys
import logging
log = logging.getLogger(__name__)


def realatom_vatom_rmsd(cg):
    """
    The RMSD between the real atoms and the virtual atoms of the stems.
    """
    vposs = []
    rposs = []

    for stem in cg.stem_iterator():
        stemv = []
        stemr = []
        resseqs = cg.get_resseqs(stem)
        resseqs = resseqs[0] + resseqs[1]
        for posCg, idPdb in zip(cg.define_residue_num_iterator(stem), resseqs):
            chainPdb, posPdb = idPdb
            vas = cg.virtual_atoms(posCg)
            for a, coords in vas.items():
                try:
                    rp = cg.chains[chainPdb][posPdb][a].get_vector(
                    ).get_array()
                except KeyError:
                    pass
                else:
                    rposs.append(rp)
                    vposs.append(coords)
    assert len(vposs) == len(rposs)
    return ftme.rmsd(np.array(vposs), np.array(rposs))


class VirtualAtomsTest(unittest.TestCase):
    def setUp(self):
        pdbfile1 = "test/forgi/threedee/data/3FU2.pdb"
        self.cg1,  = ftmc.CoarseGrainRNA.from_pdb(pdbfile1, load_chains="biggest")

        # pdbfile2 = "test/forgi/threedee/data/3V2F.pdb"  #Takes some time. Big structure
        #self.cg2,  = ftmc.CoarseGrainRNA.from_pdb(pdbfile2, load_chains="biggest")

        pdbfile3 = "test/forgi/threedee/data/1X8W.pdb"
        self.cg3,  = ftmc.CoarseGrainRNA.from_pdb(
            pdbfile3, load_chains="biggest")

        pdbfile4 = "test/forgi/threedee/data/2QBZ.pdb"
        self.cg4,  = ftmc.CoarseGrainRNA.from_pdb(
            pdbfile4, load_chains="biggest")

    def test_virtualatoms_close_to_pdb(self):
        rmsd = realatom_vatom_rmsd(self.cg1)
        print(rmsd)
        self.assertLess(rmsd, 0.94)  # 1.32 using average atoms (no vres)

        rmsd = realatom_vatom_rmsd(self.cg3)
        print(rmsd)
        self.assertLess(rmsd, 0.94)  # 1.33

        rmsd = realatom_vatom_rmsd(self.cg4)
        print(rmsd)
        self.assertLess(rmsd, 0.95)  # 1.16

        #rmsd = realatom_vatom_rmsd(self.cg2)
        # print(rmsd)
        # self.assertLess(rmsd, 1.04) # 1.30
