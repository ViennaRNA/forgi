import Bio.PDB as bpdb

import os
import sys
import unittest

import forgi.threedee.utilities.pdb as ftup
import forgi.utilities.debug as fud

class TestPDBUtilities(unittest.TestCase):
    '''
    Tests for the pdb loading utility functions.
    '''

    def setUp(self):
        return

    def test_get_particular_chain(self):
        '''
        Test loading a specific chain from a pdb file.
        '''
        chain, mr, ir = ftup.get_particular_chain('test/forgi/threedee/data/1y26_two_chains.pdb', 'Y')
        self.assertEqual(chain.id, 'Y')

    def test_get_biggest_chain(self):
        chain, mr, ir = ftup.get_biggest_chain('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.assertEqual(chain.id, 'X')


    def test_trim_chain_between(self):
        chain, mr, ir = ftup.get_biggest_chain('test/forgi/threedee/data/1y26_two_chains.pdb')

        self.assertTrue(15 in chain)
        self.assertTrue(16 in chain)
        self.assertTrue(17 in chain)

        ftup.trim_chain_between(chain, 15, 17)

        self.assertFalse(15 in chain)
        self.assertFalse(16 in chain)
        self.assertFalse(17 in chain)


    def test_is_covalent(self):
        c, mr, ir = ftup.get_biggest_chain('test/forgi/threedee/data/2mis.pdb')

        self.assertTrue(ftup.is_covalent([c[10]["C3'"], c[10]["C4'"]]))
        self.assertFalse(ftup.is_covalent([c[10]["C3'"], c[10]["C5'"]]))

    def test_pdb_file_rmsd(self):
        r = ftup.pdb_file_rmsd('test/forgi/threedee/data/1GID_native.pdb',
                               'test/forgi/threedee/data/1GID_rosetta.pdb')

        self.assertLess(r[1], 28)
        self.assertAlmostEqual( r[1], 27.700276108193787)

    def test_is_protein(self):
        struct = bpdb.PDBParser().get_structure("temp", 'test/forgi/threedee/data/1MFQ.pdb')
        chains = struct.get_chains()

        for c in chains:
            ftup.is_protein(c)
            #fud.pv('ftup.is_protein(c)')


            # code...

    def test_interchain_contacts(self):
        struct = bpdb.PDBParser().get_structure("temp", 'test/forgi/threedee/data/1MFQ.pdb')

        ftup.interchain_contacts(struct)
