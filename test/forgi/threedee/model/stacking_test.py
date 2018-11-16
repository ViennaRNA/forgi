from __future__ import print_function, absolute_import, unicode_literals, division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip, object)
import unittest
import sys
import forgi.threedee.model.coarse_grain as ftmc
import warnings
from ddt import ddt, data, file_data, unpack

"""
class Test_1B23(unittest.TestCase):
    def setUp(self): 
        cg = ftmc.CoarseGrainRNA("test/forgi/threedee/data/stacking/1B23.cg")
        self.cg = cg
        try:
            #Tyagi has stacking between the basepairs (10-25)&(27-43) and (49-65)&(66-7)
            #The asserts in the setup should make sure we use the same secondary structure in the relevant parts.


            #First stacking is (10-25)&(27-43)
            stem1 = cg.nucleotides_to_elements([cg.seq_id_to_pos(11), cg.seq_id_to_pos(24)])
            self.assertEqual(cg.pairing_partner(cg.seq_id_to_pos(11)),cg.seq_id_to_pos(24))
            self.assertEqual(len(stem1), 1)
            stem1, =stem1
            stem2 = cg.nucleotides_to_elements([cg.seq_id_to_pos(27), cg.seq_id_to_pos(43)])
            self.assertEqual(cg.pairing_partner(cg.seq_id_to_pos(27)),cg.seq_id_to_pos(43))
            self.assertEqual(len(stem2), 1)
            stem2, =stem2
            self.assertNotEqual(stem1, stem2)
            bulges = (cg.edges[stem1] & cg.edges[stem2])
            self.bulge1 = bulges.pop()

            #second_stacking is (49-65)&(66-7)
            stem1 = cg.nucleotides_to_elements([cg.seq_id_to_pos(49), cg.seq_id_to_pos(65)])
            self.assertEqual(cg.pairing_partner(cg.seq_id_to_pos(49)),cg.seq_id_to_pos(65))
            self.assertEqual(len(stem1), 1)
            stem1, =stem1
            stem2 = cg.nucleotides_to_elements([cg.seq_id_to_pos(66), cg.seq_id_to_pos(7)])
            self.assertEqual(cg.pairing_partner(cg.seq_id_to_pos(66)),cg.seq_id_to_pos(7))
            self.assertEqual(len(stem2), 1)
            stem2, =stem2
            self.assertNotEqual(stem1, stem2)
            bulges = (cg.edges[stem1] & cg.edges[stem2])
            self.bulge2 = bulges.pop()
        except AssertionError as e:
            print(cg.seq)
            print(cg.to_dotbracket_string())
            resnums = [ cg.seq_ids.index((" ", 10, " "))+1, cg.seq_ids.index((" ", 25, " "))+1,cg.seq_ids.index((" ", 27, " "))+1,cg.seq_ids.index((" ", 43, " "))+1]
            stri=""
            for i in range(len(cg.seq)):
                pos=i+1
                if pos in resnums:
                    stri+="^"
                else:
                    stri+=" "
            print(stri)
            raise
    def test_stack1(self):
        #warnings.warn("First stack in 1B23 not tested (different secondary structure)")
        self.assertTrue(self.cg.is_stacking(self.bulge1, "CG"))
    def test_stack2(self):
        self.assertTrue(self.cg.is_stacking(self.bulge2, "CG"))
    def test_no_other_stacks(self):
        for d in self.cg.defines:
            if d[0] in "mi" and d!=self.bulge1 and d!=self.bulge2 :
                flk = self.cg.get_flanking_region(d)
                seqstri = "\n"+self.cg.seq+"\n"+self.cg.to_dotbracket_string()+"\n"
                resnums = [ flk[0], self.cg.pairing_partner(flk[0]), flk[1], self.cg.pairing_partner(flk[1]) ]
                for i in range(len(self.cg.seq)):
                    pos=i+1
                    if pos in resnums:
                       seqstri+="^"
                    else:
                       seqstri+=" "
                self.assertFalse(self.cg.is_stacking(d), 
                    msg="Stacking along {} ({}-{})&{}-{}".format(
                                                    d, 
                                                    self.cg.seq_ids[flk[0]-1], 
                                                    self.cg.seq_ids[self.cg.pairing_partner(flk[0])-1],
                                                    self.cg.seq_ids[flk[1]-1], 
                                                    self.cg.seq_ids[self.cg.pairing_partner(flk[1])-1])+
                        seqstri)
"""


class TyagiData(object):
    def __init__(self, filename, bp1, bp2):
        self.bp1 = bp1
        self.bp2 = bp2
        self.filename = filename

    @property
    def __name__(self):
        na = self.filename.split(".")[0]
        if self.bp1 is None:
            return "{}_{}_".format(na, None)
        else:
            return "{}_{}_{}_".format(na, self.bp1[0], self.bp1[1])


tyagiPredictions = [
    ("1B23.cg", [[(49, 65), (66, 7)], [(10, 25), (27, 43)]]),
    ("1QF6.cg", [[(49, 65), (66, 7)]]),
    ("1QU2.cg", [[(10, 25), (27, 43)], [(49, 65), (66, 7)]]),
    ("1C0A.cg", [[(610, 625), (627, 643)], [(650, 664), (666, 607)]]),
    ("1EIY.cg", [[(10, 25), (27, 43)], [(49, 65), (66, 7)]]),
    ("1GAX.cg", [[(910, 924), (926, 942)]]),
    ("1F7U.cg", [[(949, 965), (966, 907)]]),
    ("1H4S.cg", [[(10, 25), (27, 43)], [(49, 65), (66, 7)]]),
    ("1N78.cg", [[(510, 525), (527, 543)], [(549, 565), (566, 507)]]),
    ("1J1U.cg", [[(510, 526), (528, 544)], [(550, 566), (567, 507)]]),
    ("1J2B.cg", [[(950, 966), (967, 907)]]),
    ("1QTQ.cg", [[(910, 925), (927, 943)], [(949, 965), (966, 907)]]),
    ("1EHZ.cg", [[(10, 25), (27, 43)], [(49, 65), (66, 7)]]),
    ("1FIR.cg", [[(49, 65), (66, 7)]]),
    ("1VTQ.cg", [[(10, 25), (27, 43)], [(48, 64), (65, 7)]]),
    ("1YFG.cg", [[(10, 25), (27, 43)], [(49, 65), (66, 7)]]),
    ("1NBS.cg", [[(91, 112), (113, 131)], [(132, 234), (235, 90)]]),
    ("1GID.cg", []),
    ("1QA6.cg", [[(137, 152), (154, 105)]]),
    ("1E8O.cg", [[(128, 143), (144, 103)]]),
    ("1LNG.cg", [[(187, 231), (233, 145)]]),
    ("1MFQ.cg", [[(175, 221), (222, 128)]]),
    ("1U9S.cg", [[(80, 94), (95, 109)], [
     (110, 230), (231, 79)], [(176, 194), (195, 219)]]),
    ("1FFK.cg", [[(116, 124), (125, 129)], [(109, 52), (53, 66)], [(431, 239), (240, 379)],
                 [(636, 1365), (1366, 2058)], [(905, 1300),
                                               (1301, 1354)], [(699, 727), (728, 743)],
                 [(780, 866), (868, 884)], [(1089, 1267),
                                            (1268, 1290)], [(1382, 1400), (1401, 1721)],
                 [(1725, 2050), (2051, 1374)], [(1535, 1650),
                                                (1651, 1655)], [(1634, 1551), (1552, 1569)],
                 [(1570, 1627), (1628, 1633)], [(1820, 2029),
                                                (2031, 1752)], [(1785, 1807), (1808, 1812)],
                 [(2381, 2407), (2409, 2418)], [(2682, 2712), (2713, 2767)], [
        (28, 480), (516, 27)], [(32, 451), (479, 29)],
     [(153, 184), (440, 41)], [(621, 634), (2059, 537)], [
        (1703, 1715), (1718, 1404)], [(1415, 1680), (1697, 1412)],
     [(1891, 1946), (1974, 2008)], [(2293, 2315),
                                    (2316, 2464)], [(2134, 2241), (2245, 2256)],
     [(2832, 2848), (2909, 2831)], [(780, 866),
                                    (887, 774)], [(1045, 1069), (1295, 910)],
     [(1494, 1511), (1514, 1672)]]),
    ("1J5E.cg", [[(367, 393), (395, 46)], [(316, 337), (339, 350)], [(289, 311), (312, 115)],
                 [(821, 879), (880, 575)], [(946, 1235),
                                            (1237, 1337)], [(993, 1045), (1047, 1210)],
                 [(1102, 1073), (1074, 1083)], [
        (198, 219), (221, 142)], [(588, 651), (654, 754)],
        [(673, 717), (734, 672)], [(406, 436), (442, 492)]]),
    ("1NKW.cg", [[(18, 70), (72, 109)], [(520, 30), (31, 485)], [(825, 1209), (1210, 1263)], [(1289, 1307), (1308, 1662)],
                 [(1665, 1992), (1993, 1283)], [(1443, 1580),
                                                (1581, 1584)], [(1563, 1458), (1460, 1481)],
                 [(1755, 1971), (1973, 1691)], [(2326, 2349),
                                                (2351, 2360)], [(2625, 2653), (2655, 2711)],
                 [(35, 457), (484, 32)], [(160, 191), (446, 44)], [
        (56, 69), (112, 55)], [(588, 1274), (1275, 2000)],
        [(835, 848), (849, 954)], [(1643, 1656),
                                   (1659, 1311)], [(1322, 1621), (1638, 1319)],
        [(1724, 1742), (1743, 1747)], [(1827, 1888),
                                       (1916, 1950)], [(2238, 2260), (2261, 2406)],
        [(312, 327), (330, 334)], [(2076, 2179), (2183, 2202)]]),
    ("1KH6.cg", [[(8, 23), (24, 37)], [(39, 48), (49, 5)]])
    # Note: Entry 31 corresponds to an obsolete PDB id.
]

TYAGI_HITS = []

for name, bps in tyagiPredictions:
    for bp in bps:
        TYAGI_HITS.append(TyagiData(name, bp[0], bp[1]))


KNOWN_EXCEPTIONS = {
    "1GAX.cg": [{"bp1": (906, 966), "bp2": (948, 964),
                 "comment": "Different secondary structure. Tyagi has a basepair at 965-907. "
                 "They list this stack as a false positive. They use the secondary structure"
                 " from NDB. (Probably it was generated with FR3D ???)"
                 }],
    "1B23.cg": [{"bp1": (11, 24), "bp2": (29, 41),
                 "comment": "Different secondary structure. Tyagi has a basepair at 10-25, "
                 "that is also stacking. They use the secondary structure"
                 " from NDB. (Probably it was generated with FR3D ???)"
                 }]

}


@unittest.skip("Skipping Stacking Tests. Currently our stacking predictions differ from literature.")
@ddt
class TestAllFilesTyagi(unittest.TestCase):
    def bulge_from_seqids(self, cg, bp1, bp2, asserts=True):
        try:
            stem1 = cg.nucleotides_to_elements(
                [cg.seq_id_to_pos(bp1[0]), cg.seq_id_to_pos(bp1[1])])
            if asserts:
                self.assertEqual(cg.pairing_partner(
                    cg.seq_id_to_pos(bp1[0])), cg.seq_id_to_pos(bp1[1]))
            if asserts:
                self.assertEqual(len(stem1), 1)
            stem1, = stem1
            stem2 = cg.nucleotides_to_elements(
                [cg.seq_id_to_pos(bp2[0]), cg.seq_id_to_pos(bp2[1])])
            if asserts:
                self.assertEqual(cg.pairing_partner(
                    cg.seq_id_to_pos(bp2[0])), cg.seq_id_to_pos(bp2[1]))
            if asserts:
                self.assertEqual(len(stem2), 1)
            stem2, = stem2
            if asserts:
                self.assertNotEqual(stem1, stem2)
            bulges = (cg.edges[stem1] & cg.edges[stem2])
            if asserts:
                self.assertGreaterEqual(len(bulges), 1)
            if len(bulges) == 0:
                return None
            bulge = bulges.pop()
            return bulge
        except ValueError as e:
            if asserts:
                raise
            return None

    @data(*TYAGI_HITS)
    def test_stack_found(self, data):
        cg = ftmc.CoarseGrainRNA(
            "test/forgi/threedee/data/stacking/" + data.filename)
        if data.bp1 is None:
            return
        try:
            bulge = self.bulge_from_seqids(cg, data.bp1, data.bp2)
        except (AssertionError, ValueError) as e:
            warnings.warn("Different secondary structure for {}".format(
                data.filename) + str(e))
            return
        self.assertTrue(cg.is_stacking(bulge, "CG", verbose=True))


"""
    @data(*tyagiPredictions)
    def test_no_other_stack(self, data):
        cg = ftmc.CoarseGrainRNA("test/forgi/threedee/data/stacking/"+data[0])
        for d in cg.defines:
            if d[0] not in "mi":
                continue
            flk = cg.get_flanking_region(d)
            for other in data[1]:
                if d == self.bulge_from_seqids(cg, other[0], other[1], asserts=False):
                    to_cont=True
                    break
            else:
                for exception in KNOWN_EXCEPTIONS.get(data[0], []):
                    if d == self.bulge_from_seqids(cg, exception["bp1"], exception["bp2"], asserts=False):
                        warnings.warn("\nKnown false positive in {}: {}({}-{}).\n{}\n".format(
                                                   data[0], d, exception["bp1"],
                                                   exception["bp2"], exception["comment"]))
                        break
                else:
                    seqstri = "\n"+cg.seq+"\n"+cg.to_dotbracket_string()+"\n"
                    resnums = [ flk[0], cg.pairing_partner(flk[0]), flk[1], cg.pairing_partner(flk[1]) ]
                    definerange = list(cg.define_residue_num_iterator(d))
                    for i in range(len(cg.seq)):
                        pos=i+1
                        if pos in resnums:
                           seqstri+="^"
                        elif pos in definerange:
                            seqstri+="*"
                        else:
                           seqstri+=" "
                    
                    self.assertFalse(cg.is_stacking(d, "CG"), 
                        msg="Stacking along {} ({}-{})&({}-{}) [({}-{})&({}-{})]".format(
                                                            d, 
                                                            cg.seq_ids[flk[0]-1], 
                                                            cg.seq_ids[cg.pairing_partner(flk[0])-1],
                                                            cg.seq_ids[flk[1]-1], 
                                                            cg.seq_ids[cg.pairing_partner(flk[1])-1],
                                                            flk[0],
                                                        cg.pairing_partner(flk[0]),
                                                        flk[1],
                                                        cg.pairing_partner(flk[1]))+
                              seqstri)"""
