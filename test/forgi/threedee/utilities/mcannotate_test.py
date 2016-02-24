import os
import sys
import unittest

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup

class TestMCAnnotate(unittest.TestCase):
    '''
    Tests for the pdb loading utility functions.
    '''

    def setUp(self):
        return

    def test_get_dotplot(lines):
        #cg = ftmc.from_pdb('test/forgi/threedee/data/3IVK.pdb')

        #invalid test since the pdb file has negatively numbered residues
        pass
