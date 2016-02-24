import unittest
import forgi.threedee.utilities.average_atom_positions as ftua
import forgi.utilities.debug as fud

class AverageAtomPositionsTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_key1(self):
        key = "h 9 -1 0 7 C4'"

        self.assertTrue(key in ftua.avg_atom_poss.keys())
