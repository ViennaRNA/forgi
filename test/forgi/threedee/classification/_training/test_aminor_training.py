from __future__ import print_function, division, absolute_import, unicode_literals

import unittest
import itertools as it
import math
import warnings
import contextlib
import logging
try:
    from io import StringIO  # py3K
except ImportError:
    from StringIO import StringIO  # py2k
try:
    from unittest.mock import mock_open, patch  # Python 3
except:
    from mock import mock_open, patch  # Python 2

import numpy as np
import numpy.testing as nptest

import forgi.threedee.classification._training.aminor_training as ftcta
import forgi.threedee.classification.aminor as ftca

import forgi.threedee.model.coarse_grain as ftmc

log = logging.getLogger(__name__)

@contextlib.contextmanager
def ignore_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield None

A_MULTICHAIN_FR3D_OUTPUT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.0340  A  104  A  957  U 1009 900 ----  ----   cWW  AAA  1909  1885    24                                       -     -     -
""")
A_FR3D_OUTPUT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4119  A  187  C  154  G  182 000 ----  ----   cWW  AAA    33     5    28                                                                             -     -     -
""")

ANOTHER_FR3D_OUTPUT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4119  A  187  C  154  G  182 000 ----  ----   cWW  AAA    33     5    28                                                                             -     -     -
    1S72      0.4247  A 1476  G 1867  C 1862 000  ---- ----   cWW  AAA   390   385     5                                     n2BR                                    -     -     -
    1S72      0.4416  A  395  A  216  U  224 000 ----  ----   cWW  AAA   179   171     8                                                                             -     -     -
    1S72      0.4487  A 2840  C 2047  G 1728 000 ----  ----   cWW  AAA   682   ./5AN9_N_relax/simulation_01989   307                                                                             -     -     -
""")

FR3D_ADJACENT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4292  A   49  G  149  C   42 000 ----  ----   cWW  AAA    98     7   105                                                                             -     -     -
    1S72      0.4358  A   99  C   82  G   92 000 ----  ----   cWW  AAA    17     7    10                                                                             -     -     -
""")

FR3D_NO_STEM = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4590  A 2019  C 1830  G  820 000 ----  ----   ncWW AAS   177  1158   981                                      2BR                                    -     -     -
""")
# The following problematic FR3D annotations were encountered
FR3D_PROBLEMATIC_STRING = StringIO("""
      # Adjacent elements
      4IOA      0.3195  A 1603  G 1430  C 1598 XXX ----  ----   cWW  AAA   173     5   168                                                                             -     -     -
      5DM6      0.3420  A   53  G  158  C   46 XXX ----  ----   cWW  AAA   104     7   111                                                                             -     -     -
      # Illegal residue id (not an integer)
      5DM6      0.3618  A  680  C  180  G 165B XXX ----  ----   cWW  AAA   430   450    20                                                                             -     -     -
      # Stem-Stem interaction
      4IOA      0.3836  A  693  C  185  G  165 XXX ----  ----   cWW  AAA   430   450    20                                                                             -     -     -
      # No file with this pdb-id found
wrongPdbId      0.3989  A    4  A    6  U   44 BDC ----  ----   cWW  AAA    25     4    21                                                                             -     -     -
      # No canonical stem
      2F4V      0.4685  A 1252  G 1355  C 1367 AAA  ntsS ----   cWW  AAA   103   116    13                                      2BR                                    -     -     -
      # Illegal residue number. CAUTION with future: newint converts "24L" to 24!
      1HMH      0.1476  A  24L  G  113  C  103 ACC  tSs   cSs   cWW  AAA    36    29     7                                     n2BR        n3BR                        -     -     -

""")
# Cif chain-ids are multiple characters
FR3D_CIFCHAIN_STRING = StringIO("""
      # Interaction between chains
      4TUE      0.4987  A 1912  C 1407  G 1494 RAQAQA  tSs   ncSs  cWW  AAA  2129  2047    82                                                                             -     -     -
      # legal interaction in bundle1
      4TUE      0.2085  A  767  G 1511  C 1524 QAQAQA  tSs   cSs   cWW  AAA   738   751    13                                                 n3BR                        -     -     -
""")

unittest.TestCase.longMessage = True


class TestFr3dParsing(unittest.TestCase):
    def setUp(self):
        A_MULTICHAIN_FR3D_OUTPUT.seek(0)
        A_FR3D_OUTPUT.seek(0)
        FR3D_ADJACENT.seek(0)
        FR3D_NO_STEM.seek(0)
        ANOTHER_FR3D_OUTPUT.seek(0)

    def test_fr3d_problematic(self):
        cgs = {
            "4IOA": [ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/4IOA_X.cg")],
            "5DM6": [ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/5DM6_X.cg")],
            "2F4V": [ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/2F4V_A.cg")],
            "1HMH": [ftmc.CoarseGrainRNA.from_bg_file("test/forgi/threedee/data/1HMH_A.cg")]
        }
        i = 0
        for line in FR3D_PROBLEMATIC_STRING:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            log.info("Testing line %s", line)
            with ignore_warnings():
                self.assertEqual(ftcta._parse_fred_line(
                    line, cgs, "?", "test/forgi/threedee/data/chain_id_mappings"), None)

    def test_fr3d_orientation(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1S72_0.cg")
        a, skipped = ftcta.parse_fred(
            30, {"1S72": [cg]}, ANOTHER_FR3D_OUTPUT, "test/forgi/threedee/data/chain_id_mappings")
        log.info(a)
        for geo in a:
            self.assertGreater(
                geo.dist, 5, msg="AME-interactions with distance below 5 are unlikely.")
            self.assertLess(geo.angle1, np.pi / 2,
                            msg="Angle 1 is defined between 0 and 90 degrees.")
            self.assertLess(geo.angle2, np.pi,
                            msg="Angle 2 is defined between 0 and 180 degrees.")
        for geo1, geo2 in it.combinations(a, 2):
            self.assertLess(abs(geo1.dist - geo2.dist), 20,
                            msg="Different FR3D hits should have similar coarse-grained geometry: dist. {}, {}".format(geo1, geo2))
            self.assertLess(abs(geo1.angle1 - geo2.angle1), 2.,
                            msg="Different FR3D hits should have similar coarse-grained geometry: angle1. {}, {}".format(geo1, geo2))
            self.assertLess(abs(geo1.angle2 - geo2.angle2), 3.,
                            msg="Different FR3D hits should have similar coarse-grained geometry: angle2. {}, {}".format(geo1, geo2))

    def test_parse_fred_missing_chain(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1S72_0.cg")
        print("Before filter", warnings.filters)
        with self.assertWarns(UserWarning) if hasattr(self, "assertWarns") else ignore_warnings():
            print("With filter", warnings.filters)
            a, skipped = ftcta.parse_fred(30, {"1S72": [
cg]}, A_MULTICHAIN_FR3D_OUTPUT, "test/forgi/threedee/data/chain_id_mappings")
        # Chain 9 is not connected to chain 0, thus not present in the cg-file.
        self.assertEqual(len(a), 0)
        self.assertEqual(skipped, 1)

    def test_parse_fred_adjacent(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1S72_0.cg")
        a, skipped = ftcta.parse_fred(
            30, {"1S72": [cg]}, FR3D_ADJACENT, "test/forgi/threedee/data/chain_id_mappings")
        # The loop is adjacent to the stem.
        self.assertEqual(len(a), 0)
        self.assertEqual(skipped, 2)

    def test_parse_fred_no_stem(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1S72_0.cg")
        a, skipped = ftcta.parse_fred(
            30, {"1S72": [cg]}, FR3D_NO_STEM, "test/forgi/threedee/data/chain_id_mappings")
        # The receptor is no canonical stem.
        self.assertEqual(len(a), 0)
        self.assertEqual(skipped, 1)

    def test_parse_fred1(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1S72_0.cg")
        a, skipped = ftcta.parse_fred(
            30, {"1S72": [cg]}, A_FR3D_OUTPUT, "test/forgi/threedee/data/chain_id_mappings")
        self.assertEqual(len(a), 1)
        self.assertEqual(skipped, 0)
        geometry, = a
        self.assertEqual(geometry.pdb_id, "1S72_0")

    def test_parse_fred2(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/1S72_0.cg")
        a, skipped = ftcta.parse_fred(
            30, {"1S72": [cg]}, ANOTHER_FR3D_OUTPUT, "test/forgi/threedee/data/chain_id_mappings")
        self.assertEqual(len(a), 4)
        self.assertEqual(skipped, 0)

    def test_parse_fred_cif(self):
        cg1 = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/4tue-pdb-bundle1_A.cg")
        cg2 = ftmc.CoarseGrainRNA.from_bg_file(
            "test/forgi/threedee/data/4tue-pdb-bundle2_A.cg")
        a, skipped = ftcta.parse_fred(30, {"4TUE": [
cg2, cg1]}, FR3D_CIFCHAIN_STRING, "test/forgi/threedee/data/chain_id_mappings")
        self.assertEqual(len(a), 1)
        # Interaction between not-connected chains
        self.assertEqual(skipped, 1)
        a, = a
        self.assertEqual(a.pdb_id, "4tue-pdb-bundle1_A")
