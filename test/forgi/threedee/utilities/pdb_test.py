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
        chain = ftup.get_particular_chain('test/forgi/threedee/data/1y26_two_chains.pdb', 'Y')
        self.assertEqual(chain.id, 'Y')

    def test_get_biggest_chain(self):
        chain = ftup.get_biggest_chain('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.assertEqual(chain.id, 'X')


    def test_trim_chain_between(self):
        chain = ftup.get_biggest_chain('test/forgi/threedee/data/1y26_two_chains.pdb')

        self.assertTrue(15 in chain)
        self.assertTrue(16 in chain)
        self.assertTrue(17 in chain)

        ftup.trim_chain_between(chain, 15, 17)

        self.assertFalse(15 in chain)
        self.assertFalse(16 in chain)
        self.assertFalse(17 in chain)

    def test_extract_subchain(self):
        '''
        chain = ftup.get_biggest_chain('test/forgi/threedee/data/1y26_two_chains.pdb')

        self.assertTrue(15 in chain)
        self.assertTrue(16 in chain)
        self.assertTrue(17 in chain)

        new_chain = ftup.extract_subchain(chain, 15, 17)

        self.assertTrue(15 in chain)
        self.assertTrue(16 in chain)
        self.assertTrue(17 in chain)

        self.assertTrue(15 in new_chain)
        self.assertTrue(16 in new_chain)
        self.assertTrue(17 in new_chain)
        '''

        chain = ftup.get_biggest_chain('test/forgi/threedee/data/RS_308_S_7.pdb')

        new_chain = ftup.extract_subchain(chain, (' ', 35, ' '), (' ', 39, ' '))
        self.assertTrue((' ', 35, ' ') in new_chain)
        self.assertTrue((' ', 36, ' ') in new_chain)
        self.assertTrue((' ', 37, ' ') in new_chain)
        self.assertTrue((' ', 38, ' ') in new_chain)
        self.assertTrue((' ', 39, ' ') in new_chain)

    def test_is_covalent(self):
        c = ftup.get_biggest_chain('test/forgi/threedee/data/2mis.pdb')

        self.assertTrue(ftup.is_covalent([c[10]["C3'"], c[10]["C4'"]]))
        self.assertFalse(ftup.is_covalent([c[10]["C3'"], c[10]["C5'"]]))

    def test_pdb_file_rmsd(self):
        r = ftup.pdb_file_rmsd('test/forgi/threedee/data/1GID_native.pdb',
                               'test/forgi/threedee/data/1GID_rosetta.pdb')

        print "rmsd:", r[1]

    def test_is_protein(self):
        struct = bpdb.PDBParser().get_structure("temp", 'test/forgi/threedee/data/1MFQ.pdb')
        chains = struct.get_chains()

        for c in chains:
            ftup.is_protein(c)
            #fud.pv('ftup.is_protein(c)')
            

            # code...
