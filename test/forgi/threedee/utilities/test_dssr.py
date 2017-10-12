import unittest
import forgi.threedee.utilities._dssr as ftud
import forgi.threedee.model.coarse_grain as ftmc

class TestHelperFunctions(unittest.TestCase):
    def test_dssr_to_pdb_atom_id_validIds(self):
        self.assertEqual(ftud.dssr_to_pdb_atom_id("B.C24"), ("B",(" ", 24, " ")))
        self.assertEqual(ftud.dssr_to_pdb_atom_id("LYS124"), (None, (" ", 124, " ")))
        self.assertEqual(ftud.dssr_to_pdb_atom_id("Z12.U13"), ("Z12", (" ", 13, " ")))
        self.assertEqual(ftud.dssr_to_pdb_atom_id("A.5BU36"), ("A", (" ", 36, " ")))
        self.assertEqual(ftud.dssr_to_pdb_atom_id("C.C47^M"), ("C", (" ", 47, "M")))
        self.assertEqual(ftud.dssr_to_pdb_atom_id("C.5BU47^M"), ("C",(" ", 47, "M")))
        self.assertEqual(ftud.dssr_to_pdb_atom_id(u'A.C1'), ("A",(" ", 1, " ")))

@unittest.skip("DSSR module not yet ported to multiple chain cg model.")
class TestCoaxialStacks(unittest.TestCase):
    def setUp(self):
        cg = ftmc.CoarseGrainRNA("test/forgi/threedee/data/1J1U.cg")
        self.dssr = ftud.DSSRAnnotation("test/forgi/threedee/data/1J1U.json", cg)
    def test_coaxial_stacks(self):
        self.assertEqual(self.dssr.coaxial_stacks(), [["s0","s5"],["s4","s1"]])

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
