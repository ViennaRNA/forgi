from __future__ import unicode_literals
import unittest
import logging

import forgi.graph.bulge_graph as fgb
import forgi.graph.transform_graphs as fgt
import forgi.graph.sequence as fgs
import forgi.graph.residue as fgr
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

fasta = """
>graph1
GCAGGCUAUGCCGGGUCCCGC
((.(((...)))(((.)))))
>graph1_condensed
GAGUCGUCC
(.(.)(.))
>graph2
GCAGGCUAUGCC&GGGUCCCGC
((.(((...)))&(((.)))))
>graph3
AAAAAAAAAAAAAAAA&AAAAAA
((((([[[[[[)))))&]]]]]]
"""


class TestBulgeGraphCondensed(unittest.TestCase):
    def setUp(self):
        bgs = fgb.BulgeGraph.from_fasta_text(fasta)
        self.trafo1 = fgt.BGTransformer(bgs[0])
        self.trafo2 = fgt.BGTransformer(bgs[2])
        self.trafo3 = fgt.BGTransformer(bgs[3])

        self.bg1_condensed = bgs[1]

        missing_residues = list(map(
            lambda x: fgs.MissingResidue(fgr.resid_from_str(x[:-2]), x[-1]),
            ["A:2_C", "A:5_G", "A:6_C", "A:8_A", "A:9_U", "A:10_G", "A:11_C",
             "A:14_G", "A:15_G", "A:17_C", "A:18_C", "A:20_G"]
        ))
        self.bg1_condensed.seq._set_missing_residues(missing_residues)

    def test_seqids(self):
        self.assertEqual(self.trafo1.condensed().seq._seqids,
                         list(map(fgr.resid_from_str, ["A:1", "A:3", "A:4", "A:7", "A:12", "A:13",
                                                       "A:16", "A:19", "A:21"])))

    def test_condensed_missing_correct(self):
        self.assertEqual(self.trafo1.condensed().seq.with_missing[:],
                         str(self.trafo1.bg.seq)
                         )

    def test_condensed_defines_edges_ok(self):
        self.assertNotEqual(self.trafo1.condensed().defines,
                            self.trafo1.bg.defines)
        self.assertEqual(self.trafo1.condensed().edges,
                         self.trafo1.bg.edges)
        self.assertEqual(self.trafo1.condensed().defines,
                         self.bg1_condensed.defines)
        self.assertEqual(self.trafo1.condensed().edges,
                         self.bg1_condensed.edges)

    def test_condense_cofold(self):
        s = self.trafo2.condensed().seq
        log.error("Comparing now")
        self.assertEqual(s, "GAGUC&GUCC")

    def test_condense_cofold2(self):
        db = self.trafo3.condensed().to_dotbracket_string()
        self.assertEqual(db, "([)&]")
