from __future__ import division
from __future__ import print_function
import unittest

import math
import logging
import itertools as it

import numpy as np
import numpy.testing as nptest

import forgi.threedee.model.coarse_grain as ftmc
import forgi.graph.bulge_graph as fgb
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.utilities.debug as fud

log = logging.getLogger(__name__)


class TestStats(unittest.TestCase):
    '''
    '''

    def setUp(self):
        self.angle_stats = ftms.get_angle_stats(
            'test/forgi/threedee/data/real.stats')
        pass

    def testRandom(self):
        ftms.RandomAngleStats(self.angle_stats)

    @unittest.skip("Sampling of stats moved to fess.builder.stat_storage")
    def test_filtered_stats(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3pdr_X.cg')
        cs = ftms.ConformationStats('test/forgi/threedee/data/real.stats')
        ftms.FilteredConformationStats(
            'test/forgi/threedee/data/real.stats', 'test/forgi/threedee/data/filtered_3pdr_X.stats')

        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        for d in cg.defines:
            if d[0] == 'm' and cg.get_angle_type(d) is None:
                continue
            else:
                self.assertGreater(len(cs.sample_stats(cg, d)), 0)

        '''
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1GID_A.cg')
        cs = ftms.ConformationStats('test/forgi/threedee/data/real.stats')
        fcs = ftms.FilteredConformationStats('test/forgi/threedee/data/real.stats', 'test/forgi/threedee/data/filtered_stats_1gid.csv')

        '''

    def test_conformation_stats(self):
        '''
        cs = ftms.ConformationStats()
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/1gid.cg')

        h2_stats = cs.sample_stats(cg, 'h2')
        self.assertGreater(len(h2_stats), 0)

        s1_stats = cs.sample_stats(cg, 's1')
        self.assertGreater(len(s1_stats), 0)

        i1_stats = cs.sample_stats(cg, 'i1')
        self.assertGreater(len(i1_stats), 0)

        m1_stats = cs.sample_stats(cg, 'm1')
        self.assertGreater(len(m1_stats), 0)

        f1_stats = cs.sample_stats(cg, 'f1')
        self.assertGreater(len(f1_stats), 0)

        t1_stats = cs.sample_stats(cg, 't1')
        self.assertGreater(len(t1_stats), 0)
        '''

    def test_angle_stat_deviation_from(self):
        as1 = ftms.AngleStat(u=1.57, v=0., r1=1, u1=1.57, v1=0)
        as2 = ftms.AngleStat(u=1.57, v=0., r1=1, u1=1.57, v1=0)

        self.assertTrue(np.allclose([as1.deviation_from(as2)], (0, 0, 0, 0)))
        as2 = ftms.AngleStat(u=0, v=0., r1=1, u1=1.57, v1=0)

        self.assertTrue(np.allclose(
            [as1.deviation_from(as2)], (0, 1.57, 0, 0), atol=0.01))

    def test_angle_stat_deviation_from_around_circle(self):
        as1 = ftms.AngleStat(u=0.17, v=0.17, r1=1, u1=1.57,
                             v1=0)  # 0.17 rad ~ 10 degrees
        # 2.96 rad ~ 170 deg, 3.32 rad ~ 190 deg
        as2 = ftms.AngleStat(u=2.96, v=3.32, r1=2, u1=1.57, v1=0)
        self.assertTrue(np.allclose(
            [as1.deviation_from(as2)], (1, 0.35, 0, 0), atol=0.05))  # 0.35 rad ~ 20 deg

    def test_angle_stat_deviation_from_zero_r(self):
        as1 = ftms.AngleStat(u=0.17, v=0.17, r1=0, u1=1.57, v1=0)
        as2 = ftms.AngleStat(u=0.17, v=0.17, r1=0, u1=2.7, v1=0.2)
        self.assertTrue(np.allclose([as1.deviation_from(as2)], (0, 0, 0, 0)))

    def test_3gx5_stats(self):
        cs = ftms.ConformationStats('test/forgi/threedee/data/3gx5.stats')

        self.assertTrue((8, 8) in cs.stem_stats)

    def test_angle_stat_get_angle(self):
        as1 = ftms.AngleStat(u=math.pi / 2, v=math.pi /
                             4, r1=4, u1=1.27, v1=1.5)
        self.assertAlmostEqual(as1.get_angle(), 3 * math.pi / 4)
        # If u=0 or 180 deg, the angle is always 90 deg
        as1 = ftms.AngleStat(u=0, v=1.23456, r1=3, u1=1.27, v1=1.5)
        self.assertAlmostEqual(as1.get_angle(), math.pi / 2)
        as1 = ftms.AngleStat(u=math.pi, v=0.987654, r1=3, u1=1.27, v1=1.5)
        self.assertAlmostEqual(as1.get_angle(), math.pi / 2)

    def test_angle_stat_get_angle_from_cg(self):
        fa_text = """>1
        AAACCGGGCCCCCCAAUUU
        (((..(((...)))..)))
        """

        cg, = ftmc.CoarseGrainRNA.from_fasta_text(fa_text)

        cg.coords["s0"] = np.array([0., 0., 0.]), np.array([0., 0., 1.])
        cg.twists["s0"] = np.array([0., -1., 0]), np.array([0., 1., 0.])

        cg.coords["s1"] = np.array([0., 0., 2.]), np.array([0., 1., 3.])
        cg.twists["s1"] = np.array([-1., 0., 0.]), np.array([1., 0., 0.])

        cg.coords["h0"] = np.array([0, 1, 3]), np.array([0, 2, 4])
        cg.add_bulge_coords_from_stems()

        print (cg.coords["i0"])
        print (cg.twists)

        as1, as2 = cg.get_bulge_angle_stats("i0")

        self.assertAlmostEqual(as1.get_angle(),
                               ftuv.vec_angle(cg.coords["s0"][0] - cg.coords["s0"][1],
                                              cg.coords["s1"][1] - cg.coords["s1"][0])
                               )
        self.assertAlmostEqual(as2.get_angle(),
                               ftuv.vec_angle(cg.coords["s1"][1] - cg.coords["s1"][0],
                                              cg.coords["s0"][0] - cg.coords["s0"][1])
                               )
        self.assertAlmostEqual(as1.get_angle(), math.radians(135))
        self.assertAlmostEqual(as2.get_angle(), math.radians(135))

    def test_angle_stat_get_angle_from_cg_1(self):
        fa_text = """>1
        AAACCGGGCCCCCCAAUUU
        (((..(((...)))..)))
        """

        cg, = ftmc.CoarseGrainRNA.from_fasta_text(fa_text)

        cg.coords["s0"] = np.array([0., 0., 0.]), np.array([0., 0., 1.])
        cg.twists["s0"] = np.array([0., -1., 0]), np.array([0., 1., 0.])

        cg.coords["s1"] = np.array([0., 0., 2.]), np.array([0., 0., 3.])
        cg.twists["s1"] = np.array([-1., 0., 0.]), np.array([1., 0., 0.])

        cg.coords["h0"] = np.array([0, 1, 3]), np.array([0, 2, 4])
        cg.add_bulge_coords_from_stems()

        print (cg.coords["i0"])
        print (cg.twists)

        as1, as2 = cg.get_bulge_angle_stats("i0")

        self.assertAlmostEqual(as1.get_angle(), math.radians(180))
        self.assertAlmostEqual(as2.get_angle(), math.radians(180))


class StatComparisonMixin:
    def assert_stats_equal(self, stat1, stat2):
        log.info("Asserting equality of %s and %s", stat1, stat2)
        nptest.assert_allclose(stat1.position_params()[
                               0], stat2.position_params()[0], atol=10**-12)
        # In polar coordinates, if the length is 0, the angles do not matter.
        if stat1.position_params()[0] != 0:
            nptest.assert_allclose(stat1.position_params()[
                                   1], stat2.position_params()[1], atol=10**-12)
            nptest.assert_allclose(stat1.position_params()[2] % (
                2 * math.pi), stat2.position_params()[2] % (2 * math.pi), atol=10**-12)
        nptest.assert_allclose(stat1.twist_params(),
                               stat2.twist_params(), atol=10**-12)
        nptest.assert_allclose(stat1.orientation_params(
        ), stat2.orientation_params(), atol=10**-12)

    def assert_all_angles_different(self, stat1, stat2):
        # Every angle has to change, but the length of the seperation does not change.
        if not np.all(np.logical_not(np.isclose(stat1.position_params()[1:], stat2.position_params()[1:]))):
            assert False, "Not all position parameters changed: {}, {}".format(
                stat1.position_params(), stat2.position_params())
        if not np.all(np.logical_not(np.isclose(stat1.twist_params(), stat2.twist_params()))):
            log.info("Stat1: Position (carthesian))")
            log.info(ftug.stem2_pos_from_stem1_1(
                ftuv.standard_basis, stat1.position_params()))
            log.info("Stat1: Orientation (carthesian)")
            log.info(ftug.stem2_orient_from_stem1_1(
                ftuv.standard_basis, [1] + list(stat1.orientation_params())))
            log.info("Stat2: Position (carthesian))")
            log.info(ftug.stem2_pos_from_stem1_1(
                ftuv.standard_basis, stat2.position_params()))
            log.info("Stat2: Orientation (carthesian)")
            log.info(ftug.stem2_orient_from_stem1_1(
                ftuv.standard_basis, [1] + list(stat2.orientation_params())))
            assert False, "Not all twist parameters changed: {}, {}".format(
                stat1.twist_params(), stat2.twist_params())
        if not np.all(np.logical_not(np.isclose(stat1.orientation_params(), stat2.orientation_params()))):
            assert False, "Not all orientation parameters changed: {}, {}".format(
                stat1.orientation_params(), stat2.orientation_params())


class TestStatSimilarTo(unittest.TestCase):
    def test_is_similar_to_v(self):
        stat1 = ftms.AngleStat(u=math.radians(30),
                               v=math.radians(40),
                               t=math.radians(20),
                               r1=15,
                               u1=math.radians(15),
                               v1=math.radians(15))
        stat2 = ftms.AngleStat(u=math.radians(30),
                               v=math.radians(45),
                               t=math.radians(20),
                               r1=15,
                               u1=math.radians(15),
                               v1=math.radians(15))
        self.assertTrue(stat1.is_similar_to(stat1, 10**-5))
        self.assertTrue(stat1.is_similar_to(stat2, 6))
        self.assertFalse(stat1.is_similar_to(stat2, 4))

    def test_is_similar_to_zero_length(self):
        stat1 = ftms.AngleStat(u=math.radians(30),
                               v=math.radians(40),
                               t=math.radians(20),
                               r1=0,
                               u1=math.radians(15),
                               v1=math.radians(15))
        stat2 = ftms.AngleStat(u=math.radians(32),
                               v=math.radians(40),
                               t=math.radians(20),
                               r1=0,
                               u1=math.radians(22),
                               v1=math.radians(65))
        self.assertTrue(stat1.is_similar_to(stat2, 3))
        self.assertFalse(stat1.is_similar_to(stat2, 1))
