import unittest, os
import forgi.utilities.stuff as fus
import sys

#from nose.tools import raises


class StuffTest(unittest.TestCase):
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

    #@raises(ValueError)
    def test_pairtable_to_dotbracket1(self):
        pairtable = [9, 3, 0, 1, 0, 0, 0, 3, 0, 1]
        with self.assertRaises(ValueError):
            db = fus.pairtable_to_dotbracket(pairtable)

    def test_dotbracket_to_pairtable(self):
        """
        Convert a dotbracket string to a pair table.
        """
        for pt,db in self.pt_dbs:
            self.assertEqual(fus.dotbracket_to_pairtable(db), pt)

    def test_dotbracket_to_pairtable1(self):
        db = "....(((((.....((((((((((((.....(..(((....)))..)....))))))))))))......))))).."
        fus.dotbracket_to_pairtable(db)

    #@raises(ValueError)
    def test_dotbracket_to_pairtable2(self):
        db = "((((....))..)"
        with self.assertRaises(ValueError):
            pt = fus.dotbracket_to_pairtable(db)

    #@raises(ValueError)
    def test_dotbracket_to_pairtable3(self):
        db = "((((....))..)).((..)))"
        with self.assertRaises(ValueError):
            pt = fus.dotbracket_to_pairtable(db)

    def test_pairtable_to_tuples(self):
        """
        Convert a pairtable to base pair tuples.
        """
        pt_tuples = [([4,3,4,1,2], [(1,3),(2,4),(3,1),(4,2)]),
                      ([6,4,5,6,1,2,3], [(1,4),(2,5),(3,6),(4,1),(5,2),(6,3)]),
                      ([6,3,4,1,2,6,5], [(1,3),(2,4),(3,1),(4,2),(5,6),(6,5)]),
                      ([6,2,1,0,0,0,0], [(1,2),(2,1),(3,0),(4,0),(5,0),(6,0)])]

        for pt, tup in pt_tuples:
            self.assertEqual(fus.pairtable_to_tuples(pt), tup)

    def test_tuples_to_pairtable(self):
        """
        Convert a pairtable to base pair tuples.
        """
        pt_tuples = [([4,3,4,1,2], [(1,3),(2,4),(3,1),(4,2)]),
                      ([6,4,5,6,1,2,3], [(1,4),(2,5),(3,6),(4,1),(5,2),(6,3)]),
                      ([6,3,4,1,2,6,5], [(1,3),(2,4),(3,1),(4,2),(5,6),(6,5)]),
                      ([6,2,1,0,0,0,0], [(1,2),(2,1)])]

        for pt, tup in pt_tuples:
            self.assertEqual(fus.tuples_to_pairtable(tup, pt[0]), pt)

    def test_pairtable_to_elements(self):
        db = "((((....))..))"
        #     12345678901234
        pt = fus.dotbracket_to_pairtable(db)
        print >>sys.stderr, "pt:", pt
        elements = fus.pairtable_to_elements(pt, 0, 1, len(db))

        self.assertTrue(['s',2,[1,2,13,14]] in elements)
        self.assertTrue(['s',4,[3,4,9,10]] in elements)

        self.assertTrue(['i',3,[2,3,10,11,12,13,]])
        self.assertTrue(['h',5,[4,5,6,7,8,9]])

        print >>sys.stderr, "pairtable_to_elements:", elements


    def test_pairtable_to_elements1(self):
        db = "(.[.).]"
        #     12345678901234
        pt = fus.dotbracket_to_pairtable(db)
        print >>sys.stderr, "pt:", pt
        #elements = fus.pairtable_to_elements(pt, 0, 1, len(db))

        #print >>sys.stderr, "pairtable_to_elements:", elements
