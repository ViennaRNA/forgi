import unittest
import pkgutil
import json

import forgi.utilities.debug as fud

class AverageAtomPositionsTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_key1(self):
        key = "h 9 -1 0 7 C4'"
        data = pkgutil.get_data('forgi', 'threedee/data/average_atom_positions.json')
        avg_atom_poss = json.loads(data)

        self.assertTrue(key in avg_atom_poss.keys())
