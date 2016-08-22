from __future__ import print_function
import unittest, math
import numpy as np
import numpy.testing as nptest
import forgi.projection.hausdorff as fph
import forgi.projection.projection2d as fpp
import forgi.threedee.model.coarse_grain as ftmc
import matplotlib.pyplot as plt
import sys

class TestOffsetbasedIteration(unittest.TestCase):
    def setUp(self):
        return

    def test_offsets_increasing_norm(self):
        olddx=0
        olddy=0
        for dd, norm in fph.offsets():
            self.assertAlmostEqual(fph.norm(dd), norm)
            self.assertGreaterEqual(fph.norm(dd), fph.norm((olddx,olddy)))
            olddx, olddy = dd
            if olddx>120: 
                break

    def test_all_covered_once(self):
        img=np.zeros((60,60))
        for dd, norm in fph.offsets():
            if 30+dd[0] in range(len(img)) and 30+dd[1] in range(len(img[0])):
                img[30+dd[0], 30+dd[1]]=img[30+dd[0], 30+dd[1]]+1
            if (dd[0]>30 or dd[0]<-30) and (dd[1]>30 or dd[1]<-30):
                break        
        self.assertEqual(np.min(img), 1, " - ".join(map(str,np.transpose(np.where(img==0)))))
        self.assertEqual(np.max(img), 1)

    def test_two_iterators(self):
      olddd1=(0,0)
      for dd1, norm1 in fph.offsets():
          olddd2=(0,0)
          for dd2, norm2 in fph.offsets():
              self.assertGreaterEqual(norm2, fph.norm(olddd2))
              olddd2=dd2
              if olddd2[0]>15: 
                  break
          self.assertGreaterEqual(norm1, fph.norm(olddd1))
          olddd1=dd1
          if olddd1[0]>30:
              break

    def test_all_covered_once_two_iterators(self):
        img=np.zeros((80,80))
        img2=np.zeros((80,80))
        for dd, norm in fph.offsets():
            dx,dy=dd
            if 40+dx>=0 and 40+dy>=0:
                try:
                    img[40+dx, 40+dy]=img[40+dx, 40+dy]+1
                except IndexError: 
                    pass
            if (dx>40 or dx<-40) and (dy>40 or dy<-40):
                break
        for dd, norm in fph.offsets():
            dx,dy=dd
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

class TestHelperFunctions(unittest.TestCase):
    def setUp(self):
        self.img = np.zeros((10,10))
        self.img[3,3]=1
        self.img[3,4]=1
        self.img[3,5]=1
        self.img[3,6]=1
        self.img[3,7]=1
    def test_get_longest_img_diameter_straight(self):
        self.assertAlmostEqual(fph.get_longest_img_diameter(self.img, 10), math.sqrt(26))
        self.assertAlmostEqual(fph.get_longest_img_diameter(self.img, 100), 10*math.sqrt(26))
    def test_get_longest_img_diameter_diagonal(self):
        self.img[9,9]=1
        self.assertAlmostEqual(fph.get_longest_img_diameter(self.img, 10), 7*math.sqrt(2))
    def test_get_longest_img_diameter_resolution_invariant(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        ref_proj =  fpp.Projection2D(cg, [1., 1.,   1.   ], project_virtual_atoms=True)
        ref_box=ref_proj.get_bounding_square(margin=30)
        scale=ref_box[1]-ref_box[0]
        img1, _=ref_proj.rasterize(70, bounding_square=ref_box, rotate=0) 
        img2, _=ref_proj.rasterize(40, bounding_square=ref_box, rotate=0) 
        img3, _=ref_proj.rasterize(60, bounding_square=ref_box, rotate=10)
        d1 = fph.get_longest_img_diameter(img1, scale)
        d2 = fph.get_longest_img_diameter(img2, scale)
        d3 = fph.get_longest_img_diameter(img3, scale)
        self.assertAlmostEqual(d1, d2, places=-1 )
        self.assertAlmostEqual(d1, d3, places=-1 )
        self.assertAlmostEqual(d3, d2, places=-1 )
    def test_proj_longest_axis_vs_img_diameter(self):
        cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        ref_proj =  fpp.Projection2D(cg, [1., 1.,   1.   ], project_virtual_atoms=True)
        ref_box=ref_proj.get_bounding_square(margin=30)
        scale=ref_box[1]-ref_box[0]
        ref_img, _=ref_proj.rasterize(70, bounding_square=ref_box, rotate=0)
        self.assertAlmostEqual(ref_proj.longest_axis, fph.get_longest_img_diameter(ref_img, scale))

class TestDistanceCgToImg(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.from_pdb('test/forgi/threedee/data/1y26_two_chains.pdb')
        self.ref_proj =  fpp.Projection2D(self.cg, [1., 1.,   1.   ], project_virtual_atoms=True)
        self.ref_proj_na =  fpp.Projection2D(self.cg, [1., 1.,   1.   ], project_virtual_atoms=False)
    def test_local_search(self):
        ref_box=self.ref_proj.get_bounding_square(margin=30)
        ref_img, _=self.ref_proj.rasterize(70, bounding_square=ref_box, rotate=0)
        scale=ref_box[1]-ref_box[0]
        distance, img, params = fph.locally_minimal_distance(ref_img, scale, self.cg, proj_dir=list(fph.to_polar([1,1.1,1.2]))[1:])
        self.assertLessEqual(distance, 2)
        self.assertLessEqual(abs(params[1]%360), 5)
        nptest.assert_allclose(params[0], list(fph.to_polar([1,1.1,1.2]))[1:], atol=5)
    def test_try_parameters(self):
        ref_box=self.ref_proj.get_bounding_square(margin=30)
        ref_img, _=self.ref_proj.rasterize(70, bounding_square=ref_box, rotate=45)
        scale=ref_box[1]-ref_box[0]
        distance, img, params = fph.try_parameters(ref_img, scale, self.cg, rotations=[0, 30, 60, 90, 180, 270 ], proj_directions=[list(fph.to_polar([1,1.1,1.2]))[1:]])
        self.assertLessEqual(distance, 2)
        self.assertLessEqual(abs(params[1]-45), 5)
        nptest.assert_allclose(params[0], list(fph.to_polar([1,1.1,1.2]))[1:], atol=5)
    def test_global_search(self):
        ref_box=self.ref_proj_na.get_bounding_square(margin=30)
        ref_img, _=self.ref_proj_na.rasterize(70, bounding_square=ref_box, rotate=45)
        scale=ref_box[1]-ref_box[0]
        distance, img, params = fph.globally_minimal_distance(ref_img, scale, self.cg, virtual_atoms=False, verbose = True)
        fig, ax=plt.subplots(2)
        ax[0].imshow(ref_img, interpolation="none", cmap='gray')
        ax[1].imshow(img, interpolation="none", cmap='gray')
        ax[1].set_title("{} distance".format(distance))
        plt.show()
        self.assertLessEqual(distance, 3)
        #self.assertLessEqual(abs(params[1]-45), 5)
        nptest.assert_allclose(params[0], fph.to_polar([2,0,-1.2])[1:], atol=5)

