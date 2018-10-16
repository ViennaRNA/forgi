
from builtins import str

import sys
import unittest
import subprocess as sp
import tempfile
import os
import shutil
import glob
import math
try:
    import io
except ImportError:
    import StringIO as io

import pandas as pd

import forgi.threedee.model.coarse_grain as ftmc

FORGI_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
assert FORGI_DIR.endswith("forgi")

subprocess_env = os.environ.copy()
subprocess_env["PYTHONPATH"] = FORGI_DIR+os.pathsep+subprocess_env["PATH"]

class TestRnaConvert(unittest.TestCase):
    def test_normalize_dotbracket(self):
        """
        Normalize a dotbracket string.
        """
        new_db = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "dotbracket", "({[...]})"],
                                 universal_newlines=True, env = subprocess_env)
        self.assertEqual(new_db.strip(), "(((...)))")

    def test_db_to_fasta_raises(self):
        with self.assertRaises(sp.CalledProcessError):
            sp.check_call([sys.executable, "examples/rnaConvert.py",
                           "-T", "fasta", "({[...]})"],
                          universal_newlines=True, env = subprocess_env)
    def test_to_fasta(self):
        fasta = sp.check_output([sys.executable, "examples/rnaConvert.py",
                     "-T", "fasta", "test/forgi/data/temp.bpseq"],
                    universal_newlines=True, env = subprocess_env)
        self.assertEqual(fasta[0], ">")

    def test_to_bpseq(self):
        bpseq = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "bpseq", "test/forgi/data/2hoj.fa"],
                                 universal_newlines=True, env = subprocess_env)
        self.assertEqual(bpseq[0], "1")

    def test_all_pdb_components(self):
        cg_out = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "forgi", "test/forgi/threedee/data/1DUQ.pdb"],
                                 universal_newlines=True, env = subprocess_env)
        self.assertIn("========================================", cg_out)
        self.assertIn("coord", cg_out)

    def test_to_directory(self):
        out_dir = tempfile.mkdtemp()
        cg_out = sp.check_call([sys.executable, "examples/rnaConvert.py",
                                  "-T", "forgi", "test/forgi/threedee/data/1DUQ.pdb",
                                  "--filename", out_dir],
                                 universal_newlines=True, env = subprocess_env)
        #Filenames from RNA name
        self.assertEqual(len(glob.glob(os.path.join(out_dir,"1DUQ_*.cg"))), 4)
        cg_out = sp.check_call([sys.executable, "examples/rnaConvert.py",
                                          "-T", "dotbracket", "((..))",
                                          "--filename", out_dir],
                                         universal_newlines=True, env = subprocess_env)
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
                                 universal_newlines=True, env = subprocess_env)
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo001.fa")))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo002.fa")))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo003.fa")))
        self.assertTrue(os.path.isfile(os.path.join(out_dir,"hallo004.fa")))
        self.assertFalse(os.path.isfile(os.path.join(out_dir,"hallo005.fa")))
        shutil.rmtree(out_dir)

    def test_to_neato(self):
        neato = sp.check_output([sys.executable, "examples/rnaConvert.py",
                                  "-T", "neato", "test/forgi/data/2hoj.fa", ],
                                 universal_newlines=True, env = subprocess_env)
        p = sp.Popen(["neato", "-Tsvg" ], stdin=sp.PIPE, stdout=sp.PIPE,
                                stderr=sp.PIPE, universal_newlines=True,
                                env = subprocess_env)
        svg, err = p.communicate(neato)
        self.assertIn("<title>i1</title>", svg)

class TestCompareRNA(unittest.TestCase):
    def test_compare_identity(self):
        cg_out = sp.check_output([sys.executable, "examples/compare_RNA.py",
                                  "test/forgi/threedee/data/1y26.pdb",
                                  "test/forgi/threedee/data/1y26.pdb"],
                                 universal_newlines=True, env = subprocess_env)
        lines = cg_out.split("\n")
        self.assertEqual(lines[0], "ACC:\t1.000")
        self.assertEqual(lines[1], "RMSD:\t0.000")
        self.assertEqual(lines[2], "PDB-RMSD (chain X):\t0.000")
    def test_compare_different(self):
        with self.assertRaises(sp.CalledProcessError):
            sp.check_call([sys.executable, "examples/compare_RNA.py",
                          "test/forgi/threedee/data/1y26.pdb",
                          "test/forgi/threedee/data/1gid.cg"],
                          universal_newlines=True, env = subprocess_env)

class TestOtherScripts(unittest.TestCase):
    def test_burial(self):
        """Just test that burial.py does not crash."""
        sp.check_call([sys.executable, "examples/burial.py",
                       "test/forgi/threedee/data/1y26.pdb"],
                      universal_newlines=True, env = subprocess_env)

    @unittest.skip("2D RMSD no longer supported")
    def test_projection_rmsd(self):
        rmsd1 = sp.check_output([sys.executable, "examples/projection_rmsd.py",
                               "test/forgi/threedee/data/1y26.cg",
                               "test/forgi/threedee/data/1y26.cg", "--directions",
                               "1.0,1.0,0","1.0, 1.5, 0"],
                               universal_newlines=True, env = subprocess_env)
        rmsd2 = sp.check_output([sys.executable, "examples/projection_rmsd.py",
                             "test/forgi/threedee/data/1y26.cg",
                             "test/forgi/threedee/data/1y26.cg", "--directions",
                             "1.0,1.0,0","1.0, 1.1, 0"],
                             universal_newlines=True, env = subprocess_env)
        self.assertGreater(float(rmsd1), float(rmsd2))
        self.assertLess(float(rmsd2), 1.)

    def test_describe_cg(self):
        df_out = sp.check_output([sys.executable, "examples/describe_cg.py",
                               "test/forgi/threedee/data/1y26.cg",
                               "test/forgi/data/2hoj.fa"],
                               universal_newlines=True, env = subprocess_env)
        self.assertIn("num_f", df_out)
        self.assertNotIn("...", df_out) #pandas should not truncate the tabel

    def test_describe_cg_to_file(self):
        out_dir = tempfile.mkdtemp()
        filename = os.path.join(out_dir, "test.csv")
        sp.check_call([sys.executable, "examples/describe_cg.py",
                                   "test/forgi/threedee/data/1y26.cg",
                                   "test/forgi/data/2hoj.fa",
                                    "--csv", filename],
                                   universal_newlines=True, env = subprocess_env)
        df = pd.read_csv(filename)
        jun3, = df[df["name"]=="1Y26"]["3-way-junctions"]
        self.assertEqual(jun3, 1)
        rog, = df[df["name"]=="2HOJ"]["rog_vres"]
        assert math.isnan(rog)
        rog, = df[df["name"]=="1Y26"]["rog_vres"]
        self.assertGreater(rog, 15.)
        shutil.rmtree(out_dir)

    def test_fix_twists(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/unfixed_twists.cg")
        try:
            cg1.add_all_virtual_residues()
        except AssertionError:
            pass
        else:
            raise ValueError("This Test assumes that the input file has wrong twists which "
                             "will cause ftuv to raise an assertion error.")
        fixed_cg_str = sp.check_output([sys.executable, "examples/fix_twists.py",
                                   "test/forgi/threedee/data/unfixed_twists.cg"],
                                   universal_newlines=True, env = subprocess_env)
        cg2 = ftmc.CoarseGrainRNA.from_bg_string(fixed_cg_str)
        cg2.add_all_virtual_residues() #Raises AssertionError, if twists were not fixed
