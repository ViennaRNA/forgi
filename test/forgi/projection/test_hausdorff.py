from __future__ import print_function
import unittest, math
import numpy as np
import numpy.testing as nptest
import forgi.projection.hausdorff as fph

import sys

class TestOffsetbasedIteration(unittest.TestCase):
    def setUp(self):
        return

    def test_offsets_increasing_norm(self):
        olddx=0
        olddy=0
        for dd in fph.offsets():
            self.assertGreaterEqual(fph.norm(dd), fph.norm((olddx,olddy)))
            olddx, olddy = dd
            if olddx>120: 
                break

    def test_all_covered_once(self):
        img=np.zeros((60,60))
        for dx, dy in fph.offsets():
            if 30+dx in range(len(img)) and 30+dy in range(len(img[0])):
                img[30+dx, 30+dy]=img[30+dx, 30+dy]+1
            if (dx>30 or dx<-30) and (dy>30 or dy<-30):
                break        
        self.assertEqual(np.min(img), 1, " - ".join(map(str,np.transpose(np.where(img==0)))))
        self.assertEqual(np.max(img), 1)

    def test_two_iterators(self):
      olddd1=(0,0)
      for dd1 in fph.offsets():
          olddd2=(0,0)
          for dd2 in fph.offsets():
              self.assertGreaterEqual(fph.norm(dd2), fph.norm(olddd2))
              olddd2=dd2
              if olddd2[0]>15: 
                  break
          self.assertGreaterEqual(fph.norm(dd1), fph.norm(olddd1))
          olddd1=dd1
          if olddd1[0]>30:
              break

    def test_all_covered_once_two_iterators(self):
        img=np.zeros((80,80))
        img2=np.zeros((80,80))
        for dx, dy in fph.offsets():
            if 40+dx>=0 and 40+dy>=0:
                try:
                    img[40+dx, 40+dy]=img[40+dx, 40+dy]+1
                except IndexError: 
                    pass
            if (dx>40 or dx<-40) and (dy>40 or dy<-40):
                break
        for dx, dy in fph.offsets():
            if 40+dx>=0 and 40+dy>=0:
                try:
                    img2[40+dx, 40+dy]=img2[40+dx, 40+dy]+1
                except IndexError: 
                    pass
            if (dx>40 or dx<-40) and (dy>40 or dy<-40):
                break
        self.assertEqual(np.max(img), 1)
        self.assertEqual(np.min(img), 1, " - ".join(map(str,np.transpose(np.where(img==0)))))
        self.assertEqual(np.max(img2), 1)
        self.assertEqual(np.min(img2), 1," - ".join(map(str,np.transpose(np.where(img2==0)))))

class TestHausdorffDistances(unittest.TestCase):
    def setUp(self):
        img=np.zeros((80,80))
        img[10,10]=1
        img[11,10]=1
        img[12,10]=1
        img[13,10]=1
        img[13,11]=1
        img[13,12]=1
        img[12,14]=1
        img[11,14]=1
        img[10,14]=1
        self.img=img
        img2=np.zeros((80,80))
        img2[20,10]=1
        img2[21,10]=1
        img2[22,10]=1
        img2[23,10]=1
        img2[23,11]=1
        img2[23,12]=1
        img2[22,14]=1
        img2[21,14]=1
        img2[20,14]=1
        self.img2=img2
    def test_hausdorff_helperdistance(self):
        dist=fph.hausdorff_helperdist([10,10],self.img)
        self.assertEqual(dist,0)
        dist=fph.hausdorff_helperdist([10,15],self.img)
        self.assertEqual(dist,1)
        dist=fph.hausdorff_helperdist([8,8],self.img)
        self.assertEqual(dist,math.sqrt(8))
        #Make sure no wrapping is done
        dist=fph.hausdorff_helperdist([75,10],self.img)
        self.assertEqual(dist,62)
    def test_hausdorff_distance(self):
        self.assertEqual(fph.hausdorff_distance(self.img, self.img2),10)
        self.assertEqual(fph.hausdorff_distance(self.img, self.img),0)
        self.assertEqual(fph.hausdorff_distance(self.img2, self.img2),0)
    def test_modified_hausdorff_distance(self):
        self.assertAlmostEqual(fph.modified_hausdorff_distance(
                                              self.img, self.img2),8.7495406219787473)
        self.assertEqual(fph.modified_hausdorff_distance(self.img, self.img),0)
        self.assertEqual(fph.modified_hausdorff_distance(self.img2, self.img2),0)



