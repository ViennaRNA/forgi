import unittest

import forgi.graph.sequence as fgs
import forgi.graph.residue as fgr

class TestIndexingWithoutMissing(unittest.TestCase):
    def setUp(self):
        self.seq = fgs.Sequence("CAUAAUUUCCG",
                                list(map(fgr.resid_from_str,
                                     "14,15,15.A,16,18,19,20,21,22,23,A:24".split(",")))
                                , [])
    def test_indexing_with_positive_integer(self):
        self.assertEqual(self.seq[1],"C")
        self.assertEqual(self.seq[2],"A")
        self.assertEqual(self.seq[11],"G")
        with self.assertRaises(IndexError):
            self.seq[0]
        with self.assertRaises(IndexError):
            self.seq[12]

    def test_indexing_with_negative_index(self):
        self.assertEqual(self.seq[-1],"G")
        self.assertEqual(self.seq[-5],"U")
        self.assertEqual(self.seq[-11],"C")

    def test_indexing_with_resid(self):
        self.assertEqual(self.seq[fgr.resid_from_str("14")],"C")
        self.assertEqual(self.seq[fgr.resid_from_str("15")],"A")
        self.assertEqual(self.seq[fgr.resid_from_str("15.A")],"U")
        self.assertEqual(self.seq[fgr.resid_from_str("16")],"A")
        self.assertEqual(self.seq[fgr.resid_from_str("18")],"A")
        self.assertEqual(self.seq[fgr.resid_from_str("A:24")],"G")
        with self.assertRaises(LookupError):
            self.seq[fgr.resid_from_str("B:24")]
        with self.assertRaises(LookupError):
            self.seq[fgr.resid_from_str("15.C")]
        with self.assertRaises(LookupError):
            self.seq[fgr.resid_from_str("13")]

    def test_integer_slice_all_positive(self):
        self.assertEqual(self.seq[2:5], "AUAA")
        self.assertEqual(self.seq[:5], "CAUAA")
        self.assertEqual(self.seq[2:], "AUAAUUUCCG")
        with self.assertRaises(IndexError):
            self.seq[1:4:4]

    def test_integer_slice_with_negative_start_stop(self):
        # Negative start
        with self.assertRaises(IndexError):
            self.seq[-2:]
        # Negative stop is allowed for positive steps
        self.assertEqual(self.seq[:-5], "CAUAAU")

    def test_integer_slice_neg_step(self):
        self.assertEqual(self.seq[7:3:-1], "UUAAU")
        self.assertEqual(self.seq[8::-1], "UUUAAUAC")
        self.assertEqual(self.seq[:5:-1], "GCCUUUA")
        with self.assertRaises(IndexError):
            self.seq[:-5:-1]

    def test_resid_slice_forward(self):
        self.assertEqual(self.seq[fgr.resid_from_str("15"):fgr.resid_from_str("18")], "AUAA")
        self.assertEqual(self.seq[:fgr.resid_from_str("18")], "CAUAA")
        self.assertEqual(self.seq[fgr.resid_from_str("15"):], "AUAAUUUCCG")

    def test_resid_slice_backward(self):
        self.assertEqual(self.seq[fgr.resid_from_str("18"):fgr.resid_from_str("15"):-1], "AAUA")
        self.assertEqual(self.seq[fgr.resid_from_str("18")::-1], "AAUAC")
        self.assertEqual(self.seq[:fgr.resid_from_str("15"):-1], "GCCUUUAAUA")

class TestIndexingWithMissing(unittest.TestCase):
    def setUp(self):
        # Full seq:
        # GGCAUACAUUCGUCCGG
        #   **** ***  ****
        self.seq1 = fgs.Sequence("CAUAAUUUCCG",
                                list(map(fgr.resid_from_str,
                                     "14,15,15.A,16,18,19,20,21,22,23,A:24".split(","))),
                                [{"ssseq":8, "res_name":"G", "chain":"A", "insertion":None},
                                 {"ssseq":10, "res_name":"G", "chain":"A", "insertion":"D"},
                                 {"ssseq":17, "res_name":"C", "chain":"A", "insertion":None},
                                 {"ssseq":20, "res_name":"C", "chain":"A", "insertion":"A"},
                                 {"ssseq":20, "res_name":"G", "chain":"A", "insertion":"B"},
                                 {"ssseq":25, "res_name":"G", "chain":"A", "insertion":None}])
        self.seq2 = fgs.Sequence("AAA&GGG",
                                 list(map(fgr.resid_from_str,
                                     "A:14,A:15,A:15.A,B:12,B:13,B:200.A".split(","))),
                                 [{"ssseq":13, "res_name":"G", "chain":"A", "insertion":None},
                                  {"ssseq":16, "res_name":"G", "chain":"A", "insertion":"D"},
                                  {"ssseq":11, "res_name":"C", "chain":"B", "insertion":None},
                                  {"ssseq":202, "res_name":"C", "chain":"B", "insertion":"A"}])

    def test_indexing_with_resid_without_missing(self):
        self.assertEqual(self.seq1[fgr.resid_from_str("14")],"C")
        with self.assertRaises(IndexError):
            self.seq1[fgr.resid_from_str("8")]

    def test_indexing_with_resid(self):
        self.assertEqual(self.seq1.with_missing[fgr.resid_from_str("14")],"C")
        self.assertEqual(self.seq1.with_missing[fgr.resid_from_str("8")],"G")

    def test_integer_slice_all_positive(self):
        self.assertEqual(self.seq1.with_missing[2:5], "AUACA")
        self.assertEqual(self.seq1.with_missing[:5], "GGCAUACA")
        self.assertEqual(self.seq1.with_missing[2:], "AUACAUUCGUCCGG")
        self.assertEqual(self.seq2.with_missing[:], "GAAAG&CGGGC")


    def test_integer_slice_with_negative_stop(self):
        # Negative stop is allowed for positive steps
        self.assertEqual(self.seq2.with_missing[:-3], "GAAAG&CG")

    def test_integer_slice_neg_step(self):
        self.assertEqual(self.seq1.with_missing[7:3:-1], "UUACAU")
        self.assertEqual(self.seq1.with_missing[8::-1], "UGCUUACAUACGG")
        self.assertEqual(self.seq1.with_missing[:5:-1], "GGCCUGCUUA")
        with self.assertRaises(IndexError):
            self.seq1.with_missing[:-5:-1]

    def test_resid_slice_forward_key_in_seq(self):
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("A:15"):fgr.resid_from_str("B:13")],
                         "AAG&CGG")
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("A:15"):],
                                                                    "AAG&CGGGC")
        self.assertEqual(self.seq2.with_missing[:fgr.resid_from_str("B:13")],
                                                                    "GAAAG&CGG")

    def test_resid_slice_backward_key_in_seq(self):
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("B:13"):fgr.resid_from_str("A:15"):-1],
                         "GGC&GAA")
        self.assertEqual(self.seq2.with_missing[:fgr.resid_from_str("A:15"):-1],
                                                                    "CGGGC&GAA")
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("B:13")::-1],
                                                                    "GGC&GAAAG")
    def test_resid_slice_key_not_in_seq(self):
        # Forward
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("A:13"):fgr.resid_from_str("A:16.D")],
                         "GAAAG")
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("B:11"):],
                                                                    "CGGGC")
        self.assertEqual(self.seq2.with_missing[:fgr.resid_from_str("B:11")],
                                                                    "GAAAG&C")
        # Backwards
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("A:16.D"):fgr.resid_from_str("A:13"):-1],
                         "GAAAG")
        self.assertEqual(self.seq2.with_missing[:fgr.resid_from_str("B:11"):-1],
                                                                    "CGGGC")
        self.assertEqual(self.seq2.with_missing[fgr.resid_from_str("B:11")::-1],
                                                                    "C&GAAAG")
