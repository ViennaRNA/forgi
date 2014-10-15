import unittest, os
import forgi.utilities.stuff as fus


class TestBulgeGraph(unittest.TestCase):
    def setUp(self):
        """
([)]
1234

([{)]}
123456

([)]()
123456

((([[[)))(.(]]])).
123456789012345678
"""
        self.pt_dbs = [([4,3,4,1,2], "([)]"),
                       ([6,4,5,6,1,2,3], "([{)]}"),
                       ([6,3,4,1,2,6,5], "([)]()"),
											 ([18, 9,8,7, 15,14,13, 3,2,1, 17,0,16,6,5,4,12,10,0], "((([[[)))(.(]]])).")]
        pass

    def test_pairtable_to_dotbracket(self):
        """
        Convert a pair table to a dotbracket string.

        """
        for pt,db in self.pt_dbs:
            self.assertEqual(fus.pairtable_to_dotbracket(pt), db)

    def test_dotbracket_to_pairtable(self):
        """
        Convert a dotbracket string to a pair table.
        """
        for pt,db in self.pt_dbs:
            self.assertEqual(fus.dotbracket_to_pairtable(db), pt)
