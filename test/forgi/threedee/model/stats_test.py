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
        self.angle_stats = ftms.get_angle_stats('test/forgi/threedee/data/real.stats')
        pass

    def testRandom(self):
        ftms.RandomAngleStats(self.angle_stats)

    @unittest.skip("Sampling of stats moved to fess.builder.stat_storage")
    def test_filtered_stats(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3pdr_X.cg')
        cs = ftms.ConformationStats('test/forgi/threedee/data/real.stats')
        ftms.FilteredConformationStats('test/forgi/threedee/data/real.stats', 'test/forgi/threedee/data/filtered_3pdr_X.stats')

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

    def test_angle_stat_difference(self):
        as1 = ftms.AngleStat(u=1.57, v=0., r1=1, u1=1.57, v1=0)
        as2 = ftms.AngleStat(u=1.57, v=0., r1=1, u1=1.57, v1=0)

        self.assertTrue(np.allclose([as1.diff(as2)],[0]))
        as2 = ftms.AngleStat(u=0, v=0., r1=1, u1=1.57, v1=0)

        self.assertTrue(np.allclose([as1.diff(as2)],[math.sqrt(2)],0.01))

        as2 = ftms.AngleStat(u=1.57, v=0., r1=1, u1=0, v1=0)
        fud.pv('as1.diff(as2)')

    def test_3gx5_stats(self):
        cs = ftms.ConformationStats('test/forgi/threedee/data/3gx5.stats')

        self.assertTrue((8,8) in cs.stem_stats)

    @unittest.expectedFailure
    def test_angle_stat_get_angle(self):
        as1 = ftms.AngleStat(u=math.pi/2, v=0., r1=1, u1=1.57, v1=1.1)
        #self.assertAlmostEqual(as1.get_angle(), 0)

        as1 = ftms.AngleStat(u=math.pi/2, v=math.pi/4, r1=4, u1=1.27, v1=1.5)
        self.assertAlmostEqual(as1.get_angle(), math.pi/4)

    def test_angle_stat_get_angle_from_cg(self):
        fa_text = """>1
        AAACCGGGCCCCCCAAUUU
        (((..(((...)))..)))
        """

        bg = fgb.from_fasta_text(fa_text)
        cg = ftmc.from_bulge_graph(bg)

        cg.coords["s0"]=np.array([0.,0.,0.]), np.array([0.,0.,1.])
        cg.twists["s0"]=np.array([0.,-1.,0]), np.array([0.,1.,0.])

        cg.coords["s1"]=np.array([0.,0.,2.]), np.array([0.,1.,3.])
        cg.twists["s1"]=np.array([-1.,0.,0.]), np.array([1.,0.,0.])

        cg.coords["h0"]=np.array([0,1,3]), np.array([0,2,4])
        cg.add_bulge_coords_from_stems()

        print (cg.coords["i0"])
        print (cg.twists)

        as1, as2 = cg.get_bulge_angle_stats("i0")

        self.assertAlmostEqual(as1.get_angle(),
                    ftuv.vec_angle(cg.coords["s0"][0]-cg.coords["s0"][1],
                                   cg.coords["s1"][1]-cg.coords["s1"][0])
                              )
        self.assertAlmostEqual(as2.get_angle(),
                    ftuv.vec_angle(cg.coords["s1"][1]-cg.coords["s1"][0],
                                   cg.coords["s0"][0]-cg.coords["s0"][1])
                              )
        self.assertAlmostEqual(as1.get_angle(), math.radians(135))
        self.assertAlmostEqual(as2.get_angle(), math.radians(135))
    def test_angle_stat_get_angle_from_cg_1(self):
        fa_text = """>1
        AAACCGGGCCCCCCAAUUU
        (((..(((...)))..)))
        """

        bg = fgb.from_fasta_text(fa_text)
        cg = ftmc.from_bulge_graph(bg)

        cg.coords["s0"]=np.array([0.,0.,0.]), np.array([0.,0.,1.])
        cg.twists["s0"]=np.array([0.,-1.,0]), np.array([0.,1.,0.])

        cg.coords["s1"]=np.array([0.,0.,2.]), np.array([0.,0.,3.])
        cg.twists["s1"]=np.array([-1.,0.,0.]), np.array([1.,0.,0.])

        cg.coords["h0"]=np.array([0,1,3]), np.array([0,2,4])
        cg.add_bulge_coords_from_stems()

        print (cg.coords["i0"])
        print (cg.twists)

        as1, as2 = cg.get_bulge_angle_stats("i0")

        self.assertAlmostEqual(as1.get_angle(), math.radians(180))
        self.assertAlmostEqual(as2.get_angle(), math.radians(180))


class StatComparisonMixin:
    def assert_stats_equal(self, stat1, stat2):
        log.info("Asserting equality of %s and %s", stat1, stat2)
        nptest.assert_allclose(stat1.position_params()[0], stat2.position_params()[0], atol=10**-12)
        if stat1.position_params()[0] !=0: # In polar coordinates, if the length is 0, the angles do not matter.
            nptest.assert_allclose(stat1.position_params()[1], stat2.position_params()[1], atol=10**-12)
            nptest.assert_allclose(stat1.position_params()[2]%(2*math.pi), stat2.position_params()[2]%(2*math.pi), atol=10**-12)
        nptest.assert_allclose(stat1.twist_params(), stat2.twist_params(), atol=10**-12)
        nptest.assert_allclose(stat1.orientation_params(), stat2.orientation_params(), atol=10**-12)

    def assert_all_angles_different(self, stat1, stat2):
        # Every angle has to change, but the length of the seperation does not change.
        if not np.all( np.logical_not( np.isclose(stat1.position_params()[1:], stat2.position_params()[1:]))):
            assert False, "Not all position parameters changed: {}, {}".format(stat1.position_params(), stat2.position_params())
        if not np.all( np.logical_not( np.isclose(stat1.twist_params(), stat2.twist_params()))):
            log.info("Stat1: Position (carthesian))")
            log.info(ftug.stem2_pos_from_stem1_1(ftuv.standard_basis,stat1.position_params()))
            log.info("Stat1: Orientation (carthesian)")
            log.info(ftug.stem2_orient_from_stem1_1(ftuv.standard_basis, [1]+list(stat1.orientation_params())))
            log.info("Stat2: Position (carthesian))")
            log.info(ftug.stem2_pos_from_stem1_1(ftuv.standard_basis, stat2.position_params()))
            log.info("Stat2: Orientation (carthesian)")
            log.info(ftug.stem2_orient_from_stem1_1(ftuv.standard_basis, [1]+list(stat2.orientation_params())))
            assert False, "Not all twist parameters changed: {}, {}".format(stat1.twist_params(), stat2.twist_params())
        if not np.all( np.logical_not( np.isclose(stat1.orientation_params(), stat2.orientation_params()))):
            assert False, "Not all orientation parameters changed: {}, {}".format(stat1.orientation_params(), stat2.orientation_params())

class TestIdentityStat(unittest.TestCase, StatComparisonMixin):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3way.cg')
    def test_identity_stat(self):
        loop, = self.cg.find_mlonly_multiloops()
        sum_stat = None
        for elem in loop:
            stat, = [ s for s in self.cg.get_stats(elem)
                      if s.ang_type in [2,3,-4] ]
            if sum_stat is None:
                sum_stat = stat
            else:
                sum_stat = ftug.sum_of_stats(sum_stat, stat)
        self.assert_stats_equal(sum_stat, ftms.IDENTITY_STAT)

class TestStatInversion(unittest.TestCase, StatComparisonMixin):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3way.cg')


    def test_twice_invert_angle_stat_gives_original_stat(self):
        cg = self.cg
        for ml in [ m for m in cg.get_mst() if m[0]=="m" ]:
            stat, = [ stat for stat in cg.get_stats(ml) if stat.ang_type == cg.get_angle_type(ml) ]
            inverse_stat = -stat
            # Make sure the inversion changes the stat
            self.assert_all_angles_different(stat, inverse_stat)
            # The second inversion changes it back to the original stat
            self.assert_stats_equal(stat, -inverse_stat)

    def test_invert_angle_stat2(self):
        cg = self.cg
        for ml in [ m for m in cg.get_mst() if m[0]=="m" ]:
            stat1, stat2 = cg.get_stats(ml)
            self.assert_all_angles_different(stat1, stat2)
            self.assert_stats_equal(stat2, ftug.invert_angle_stat(stat1))
            self.assert_stats_equal(stat1, ftug.invert_angle_stat(stat2))

    def test_invert_angle_stat3(self):
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3D0U_A.cg')
        for loop in it.chain(cg.mloop_iterator(), cg.iloop_iterator()):
            log.info("Checking loop %s", loop)
            stat1, stat2 = cg.get_stats(loop)
            inverse1 = -stat1
            inverse2 = -stat2

            #self.assert_all_angles_different(stat1, stat2)
            self.assert_stats_equal(stat2, inverse1)
            self.assert_stats_equal(stat1, inverse2)
            self.assert_stats_equal(stat1, -inverse1)
            self.assert_stats_equal(stat2, -inverse2)

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


def statsum_around_junction(cg, first_elem):
    elem = first_elem
    stat, = [ s for s in cg.get_stats(elem) if s.ang_type in [2, 3, -4]]
    sum_stat = stat
    elem = cg.get_next_ml_segment(elem)
    i=0
    while elem != first_elem:
        i+=1
        stat, = [ s for s in cg.get_stats( elem) if s.ang_type in [2, 3, -4]]
        sum_stat = ftug.sum_of_stats(sum_stat, stat)
        elem = cg.get_next_ml_segment(elem)
        if i>100:
            assert False
    return sum_stat

class TestStatSum(unittest.TestCase, StatComparisonMixin):
    def test_sum_of_stats_around_a_junction(self):
        """
        Independent of the start, the sum of stats all around a junction aalways gives the same.
        """
        cg = ftmc.CoarseGrainRNA('test/forgi/threedee/data/3D0U_A.cg')
        cg.log(logging.INFO)
        for ml in ["m0", "m1", "m2", "m7", "m8"]:
            assert statsum_around_junction(cg, ml).is_similar_to(ftms.IDENTITY_STAT, 1)

    def test_statsum_doesnot_commutes(self):
        # Rotations in 3D space do not commute.
        stat1 = ftms.AngleStat(pdb_name="A",
                               u=math.radians(30),
                               v=math.radians(40),
                               t=math.radians(20),
                               r1=12,
                               u1=math.radians(15),
                               v1=math.radians(15))
        stat2 = ftms.AngleStat(pdb_name = "B",
                               u=math.radians(32),
                               v=math.radians(5),
                               t=math.radians(12),
                               r1=10,
                               u1=math.radians(22),
                               v1=math.radians(65))
        sum1 = ftug.sum_of_stats(stat1, stat2)
        sum2 = ftug.sum_of_stats(stat2, stat1)
        log.info("%s =?= %s", sum1, sum2)
        with self.assertRaises(AssertionError):
            self.assert_stats_equal(sum1, sum2)

    def test_effect_of_statsum_on_difference(self):
        # The A-B != (C+A)-(C+B) for stats
        stat1 = ftms.AngleStat(pdb_name="A",
                               u=math.radians(30),
                               v=math.radians(40),
                               t=math.radians(20),
                               r1=12,
                               u1=math.radians(15),
                               v1=math.radians(15))
        stat2 = ftms.AngleStat(pdb_name = "B",
                               u=math.radians(32),
                               v=math.radians(5),
                               t=math.radians(12),
                               r1=10,
                               u1=math.radians(22),
                               v1=math.radians(65))
        stat3 = ftms.AngleStat(pdb_name = "C",
                               u=math.radians(12),
                               v=math.radians(8),
                               t=math.radians(7),
                               r1=5,
                               u1=math.radians(13),
                               v1=math.radians(65))
        assert stat1.deviation_from(stat2) != ftug.sum_of_stats(stat3, stat1).deviation_from(ftug.sum_of_stats(stat3, stat2))
