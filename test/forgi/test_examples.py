
from builtins import str

import sys
import unittest
import subprocess as sp
import tempfile
import os.path
import shutil
import glob
try:
    import io
except ImportError:
    import StringIO as io


class TestRnaConvert(unittest.TestCase):
    def test_normalize_dotbracket(self):
        """
        Normalize a dotbracket string.
        """
        new_db = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "dotbracket", "({[...]})"],
                                 universal_newlines=True)
        self.assertEqual(new_db.strip(), "(((...)))")

    def test_db_to_fasta_raises(self):
        with self.assertRaises(sp.CalledProcessError):
            sp.check_call([sys.executable, "examples/rnaConvert.py",
                           "-T", "fasta", "({[...]})"],
                          universal_newlines=True)
    def test_to_fasta(self):
        fasta = sp.check_output([sys.executable, "examples/rnaConvert.py",
                     "-T", "fasta", "test/forgi/data/temp.bpseq"],
                    universal_newlines=True)
        self.assertEqual(fasta[0], ">")

    def test_to_bpseq(self):
        bpseq = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "bpseq", "test/forgi/data/2hoj.fa"],
                                 universal_newlines=True)
        self.assertEqual(bpseq[0], "1")

    def test_all_pdb_components(self):
        cg_out = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "forgi", "test/forgi/threedee/data/1DUQ.pdb"],
                                 universal_newlines=True)
        self.assertIn("========================================", cg_out)
        self.assertIn("coord", cg_out)

    def test_to_directory(self):
        out_dir = tempfile.mkdtemp()
        cg_out = sp.check_call([sys.executable, "examples/rnaConvert.py",
                                  "-T", "forgi", "test/forgi/threedee/data/1DUQ.pdb",
                                  "--filename", out_dir],
                                 universal_newlines=True)
        #Filenames from RNA name
        self.assertEqual(len(glob.glob(os.path.join(out_dir,"1DUQ_*.cg"))), 4)
        cg_out = sp.check_call([sys.executable, "examples/rnaConvert.py",
                                          "-T", "dotbracket", "((..))",
                                          "--filename", out_dir],
                                         universal_newlines=True)
        # Unnamed RNA. Use "rna" as filename.
        print(glob.glob(out_dir+"/*"))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"rna001.dotbracket")))
        shutil.rmtree(out_dir)

    def test_to_file(self):
        out_dir = tempfile.mkdtemp()
        filename = os.path.join(out_dir, "hallo")
        cg_out = sp.check_call([sys.executable, "examples/rnaConvert.py",
                                  "-T", "fasta", "test/forgi/threedee/data/1DUQ.pdb",
                                  "--filename", filename],
                                 universal_newlines=True)
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo001.fa")))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo002.fa")))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo003.fa")))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo004.fa")))
        self.assertFalse(os.path.isfile(os.path.join(out_dir,"hallo005.fa")))
        shutil.rmtree(out_dir)

    def test_to_neato(self):
        neato = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "neato", "test/forgi/data/2hoj.fa", ],
                                 universal_newlines=True)
        p = sp.Popen(["neato", "-Tsvg" ], stdin=sp.PIPE, stdout=sp.PIPE,
                                stderr=sp.PIPE, universal_newlines=True)
        svg, err = p.communicate(neato)
        self.assertIn("<title>i1</title>", svg)
