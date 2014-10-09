import unittest, os

class TestBulgeGraph(unittest.TestCase, GraphVerification):
    def setUp(self):
        pass

    def test_pairtable_to_dotbracket(self):
        """
        Convert a pair table to a dotbracket string.

        """
        ot = [7,3,4,1,2]
        db = "([)]"

        self.assertEqual(fus.pairtable_to_dotbracket_string(ot), db)

        db = "([{)]}"
        ot = [11,4,5,6,1,2,3]

        self.assertEqual(fus.pairtable_to_dotbracket_string(ot), db)

    def test_dotbracket_to_pairtable(self):
        ot = [7,3,4,1,2]
        db = "([)]"

        self.assertEqual(fus.dotbracket_to_pairtable(db), ot)
        
        db = "([{)]}"
        ot = [11,4,5,6,1,2,3]

        self.assertEqual(fus.dotbracket_to_pairtable(db), ot)
