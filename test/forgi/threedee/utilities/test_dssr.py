import unittest
import json
import forgi.threedee.utilities._dssr as ftud
import forgi.threedee.model.coarse_grain as ftmc
import forgi.graph.residue as fgr

class TestHelperFunctions(unittest.TestCase):
    def test_dssr_to_pdb_atom_id_validIds(self):
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "B.C24"), ("B", (" ", 24, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "1:B.C24"), ("B", (" ", 24, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "LYS124"), (None, (" ", 124, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "Z12.U13"), ("Z12", (" ", 13, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "A.5BU36"), ("A", (" ", 36, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "C.C47^M"), ("C", (" ", 47, "M")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            "C.5BU47^M"), ("C", (" ", 47, "M")))
        self.assertEqual(ftud.dssr_to_pdb_resid(u'A.C1'), ("A", (" ", 1, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            u'B.U-1'), ("B", (" ", -1, " ")))
        self.assertEqual(ftud.dssr_to_pdb_resid(
            u'A.A-2'), ("A", (" ", -2, " ")))


class TestCoaxialStacks(unittest.TestCase):
    def setUp(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/1J1U.cg")
        with open("test/forgi/threedee/data/1J1U.json") as f:
            j = json.load(f)
        self.dssr = ftud.DSSRAnnotation(j, cg)

    def test_coaxial_stacks(self):
        self.assertEqual(sorted(self.dssr.coaxial_stacks()),
                         sorted([["s2", "s1"], ["s0", "s3"]]))

    @unittest.skip("Currently not working. TODO")
    def test_compare_coaxial_stacks(self):
        forgi, dssr = self.dssr.compare_coaxial_stack_annotation()
        self.assertEqual(len(dssr), 2)
        self.assertGreaterEqual(len(forgi), 1)
        self.assertGreaterEqual(len(forgi & dssr), 1)
        self.assertIn(("s0", "s5"), (x.stems for x in forgi))
        for x in forgi:
            self.assertEqual(x.forgi, "stacking")
        for x in dssr:
            self.assertEqual(x.dssr, "stacking")

    def test_stacking_nts(self):
        stacks = self.dssr.stacking_nts()
        self.assertIn((fgr.RESID("B:544"), fgr.RESID("B:545")), stacks)
        self.assertNotIn((fgr.RESID("B:549"), fgr.RESID("B:544")), stacks)
