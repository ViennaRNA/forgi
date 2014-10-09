import unittest, os

class TestBulgeGraph(unittest.TestCase, GraphVerification):
    def setUp(self):
        """
([)]
1234

([{)]}
123456

([)]()
123456
"""
        self.pt_dbs = [([7,3,4,1,2], "([)]"),
                       ([11,4,5,6,1,2,3], "([{)]}"),
                       ([6,3,4,1,2,6,5], "([)]()")]
        pass

    def test_pairtable_to_dotbracket(self):
        """
        Convert a pair table to a dotbracket string.

        """
        for pt,db in self.pt_dbs:
            self.assertEqual(fus.pairtable_to_dotbracket_string(pt), db)

    def test_dotbracket_to_pairtable(self):
        """
        Convert a dotbracket string to a pair table.
        """
        for pt,db in self.pt_dbs:
            self.assertEqual(fus.dotbracket_to_pairtable(db), pt)
