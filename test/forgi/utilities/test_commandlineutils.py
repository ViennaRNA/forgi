from __future__ import print_function, unicode_literals, division
import unittest

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import forgi.utilities.commandline_utils as fuc
import forgi.threedee.model.coarse_grain as ftmc
import forgi.graph.bulge_graph as fgb

class TestSniffFiletype(unittest.TestCase):
    def test_sniff_fasta(self):
        with open("test/forgi/data/2hoj.fa") as f:
            self.assertEqual(fuc.sniff_filetype(f), "fasta")
    def test_sniff_bpseq(self):
        with open("test/forgi/data/1gid.bpseq") as f:
            self.assertEqual(fuc.sniff_filetype(f), "bpseq")
    def test_sniff_bg(self):
        with open("test/forgi/data/telomerase.cg") as f:
            self.assertEqual(fuc.sniff_filetype(f), "forgi")
    def test_sniff_dotbracket(self):
        db = StringIO("(((...(((...))))))")
        self.assertEqual(fuc.sniff_filetype(db), "other")
    def test_sniff_pdb(self):
        with open("test/forgi/threedee/data/1A34.pdb") as f:
            self.assertEqual(fuc.sniff_filetype(f), "pdb")

class TestLoadRNA(unittest.TestCase):
    def test_db_direct(self):
        db = "(((..[[[..)))..(((..]]].)))"
        result = fuc.load_rna(db, "any", allow_many=True)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertIsInstance(result[0], fgb.BulgeGraph)
        with self.assertRaises(IOError): #No such file.
            fuc.load_rna(db, "cg", allow_many=True)
        result = fuc.load_rna(db, "any", allow_many=False)
        self.assertIsInstance(result, fgb.BulgeGraph)
    def test_fasta(self):
        with self.assertRaises(ValueError):
            fuc.load_rna("test/forgi/data/2hoj.fa", "3d")
        result = fuc.load_rna("test/forgi/data/2hoj.fa", "cg", allow_many=False)
        self.assertIsInstance(result, ftmc.CoarseGrainRNA)
        result = fuc.load_rna("test/forgi/data/2hoj.fa", "any", allow_many=False)
        self.assertIsInstance(result, fgb.BulgeGraph)
        result = fuc.load_rna("test/forgi/data/2hoj.fa", "any", allow_many=True)
        self.assertIsInstance(result, list)

class TestCommanldineUtils(unittest.TestCase):
    def test_load_rna_pdb_simple(self):
        cg = fuc.load_rna("test/forgi/threedee/data/1y26.pdb", "pdb",
                          allow_many=False)
        self.assertIsInstance(cg, ftmc.CoarseGrainRNA)
        # With allow_many we get a list of length 1
        cgs = fuc.load_rna("test/forgi/threedee/data/1y26.pdb", "pdb",
                           allow_many=True)
        self.assertEqual(len(cgs), 1)
        self.assertEqual(cg.defines, cgs[0].defines)

    def test_load_rna_pdb_allow_many_and_chain_option(self):
        # More than 1 chain in file
        with self.assertRaises(ValueError):
            fuc.load_rna("test/forgi/threedee/data/1DUQ.pdb", "pdb", False)
        cg = fuc.load_rna("test/forgi/threedee/data/1DUQ.pdb",
                          "pdb", False, pdb_chain="A")
        self.assertIsInstance(cg, ftmc.CoarseGrainRNA)
        # Two chains, forming one stem:
        cg = fuc.load_rna("test/forgi/threedee/data/1DUQ.pdb",
                          "pdb", False, pdb_chain=["A", "B"])
        self.assertIsInstance(cg, ftmc.CoarseGrainRNA)
        cgs = fuc.load_rna("test/forgi/threedee/data/1DUQ.pdb", "pdb", True)
        self.assertIsInstance(cgs, list)
        # 4 molecules with 2 chains each
        self.assertEqual(len(cgs), 4)
        self.assertIsInstance(cgs[0], ftmc.CoarseGrainRNA)

    def test_load_rna_pdb_with_secondary_structure(self):
        with self.assertRaises(ValueError):
            # Need chain for pdb_dotbracket
            fuc.load_rna("test/forgi/threedee/data/1FUF.pdb", "pdb", False,
                         dissolve_length_one_stems=False,
                         pdb_dotbracket="(((((.(.(((((&)))))))))))..")
        cg1 = fuc.load_rna("test/forgi/threedee/data/1FUF.pdb", "pdb", False,
                           dissolve_length_one_stems=False, pdb_chain=[
                               "A", "B"],
                           pdb_dotbracket="(((((.(.(((((&)))))))))))..")
        cg2 = fuc.load_rna("test/forgi/threedee/data/1FUF.pdb", "pdb", False,
                           dissolve_length_one_stems=True, pdb_chain=[
                               "A", "B"],
                           pdb_dotbracket="(((((.(.(((((&)))))))))))..")
        self.assertEqual(len(cg1.defines), 6)
        self.assertEqual(cg1.defines, {"s0": [1, 5, 20, 24], "s1": [7, 7, 19, 19],
                                       "s2": [9, 13, 14, 18], "i0": [6, 6],
                                       "i1": [8, 8], "t0": [25, 26]})
        self.assertEqual(len(cg2.defines), 4)

    def test_load_rna_pdb_dissolve_length_one_stems(self):
        cg1 = fuc.load_rna("test/forgi/threedee/data/1FUF.pdb", "pdb", False,
                           dissolve_length_one_stems=False)
        cg2 = fuc.load_rna("test/forgi/threedee/data/1FUF.pdb", "pdb", False,
                           dissolve_length_one_stems=True)
        self.assertGreater(len(list(cg1.stem_iterator())),
                           len(list(cg2.stem_iterator())))

    def test_load_fasta_types(self):
        with self.assertRaises(fuc.WrongFileFormat):
            fuc.load_rna("test/forgi/data/pk.fa", "3d", False)
        cg = fuc.load_rna("test/forgi/data/pk.fa", "cg", False)
        self.assertIsInstance(cg, ftmc.CoarseGrainRNA)

    def test_load_fasta_cleaning_sec_stru(self):
        bg = fuc.load_rna("test/forgi/data/pk.fa", "bg",
                          False, dissolve_length_one_stems=False)
        bg2 = fuc.load_rna("test/forgi/data/pk.fa", "bg",
                           False, dissolve_length_one_stems=True)
        self.assertLess(len(bg2.defines), len(bg.defines))

    def test_sniff_malformed_file(self):
        file = StringIO("\n>fasta header\nAAAGGGCCC\n.........")
        self.assertEqual(fuc.sniff_filetype(file), "fasta")

    @unittest.skip("Requires ViennaRNApackage")
    def test_with_missing_refolded(self):
        bg, = fuc.load_rna("test/forgi/threedee/data/1FJG_reduced.pdb")
        db = bg.to_dotbracket_string(include_missing=True)
        self.assertIn("-", db)
        cg2 = fuc.with_missing_refolded(bg)
        db2 = cg2.to_dotbracket_string()
        self.assertEqual(len(db2), 208)
        self.assertNotIn("-", db2)
        self.assertEqual(db2, ".........((((.........))))...................................(((..((.((((((((((.....)))))))))).)))))......................................(((((....((((((....)))))).....))))).........(((((....)))))............")

    @unittest.skip("Requires ViennaRNApackage")
    def test_with_missing_refolded_pk(self):
        bg, = fuc.load_rna("test/forgi/threedee/data/3DHS.cif", pdb_remove_pk=False)
        db = bg.to_dotbracket_string(include_missing=True)
        self.assertIn("-", db)
        cg2 = fuc.with_missing_refolded(bg)
        db2 = cg2.to_dotbracket_string()
        self.assertEqual(len(db2), len(db))
        self.assertNotIn("-", db2)
        self.assertEqual(db2, ".((((((((((((.(((((((.(((((((((....))))))))).....[[.[[[[[(((((((((.((......)))))).....(((((....))))))))))...(((...........)))...(((((((................)))))))((((((..........))))))........)))))))(((........(((((......)))))..........)))......]]]]]]]....))))))))))))....")



    def test_insert_pk_into_stru(self):
        self.assertEqual(fuc.insert_pk_into_stru('((..(...)))]]', '(([[-----))]]'),
                            '(([[(...)))]]')
