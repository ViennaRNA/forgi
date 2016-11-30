from __future__ import division
import unittest

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.vector as ftuv

import forgi.utilities.debug as fud

import math
import numpy as np

class TestStats(unittest.TestCase):
    '''
    '''

    def setUp(self):
        self.angle_stats = ftms.get_angle_stats('test/forgi/threedee/data/real.stats')
        pass

    def testRandom(self):
        ftms.RandomAngleStats(self.angle_stats)

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

    def test_sample_stats(self):
        fa_text = """>1
        AGAGGUUCUAGCUACACCCUCUAUAAAAAACUAAGG
        (((((............)))))..............
        """
        
        cg = ftmc.CoarseGrainRNA()
        cg.from_fasta(fa_text)

        conf_stats = ftms.get_conformation_stats()
        stats = conf_stats.sample_stats(cg, 't0')

        fud.pv('cg.to_cg_string()')
        fud.pv('stats')

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

        cg = ftmc.CoarseGrainRNA()
        cg.from_fasta(fa_text)

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

        cg = ftmc.CoarseGrainRNA()
        cg.from_fasta(fa_text)

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

