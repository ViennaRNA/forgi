from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__ = type #New style classes in Python 2.x

import argparse, random, math
import itertools as it
import numpy as np

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.projection2d as ftmp
 
import matplotlib.pyplot as plt
import scipy.signal
import scipy.ndimage

RASTER=50
def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('files', nargs=2, help='Two *.cg/*.coord files to compare.')
    #Options
    #parser.add_argument('-d', '--direction', action='store', default="1.0,0.0,0.0", help='The direction for the projection. A comma-seperated triple.', type=str)
    return parser


def get_box(projection, width=500, offset_x=0, offset_y=0):
    left, right, down, up=projection.get_bounding_square()
    center_hor=left+(right-left)/2
    center_ver=down+(up-down)/2
    box=(center_hor-width/2+offset_x, center_hor+width/2+offset_x, center_ver-width/2+offset_y, center_ver+width/2+offset_y)
    return box

def score_1(ref_img, img):
    ref=(ref_img>np.zeros_like(ref_img))
    test=(img>np.zeros_like(img))
    #TruePositives: White in 1 and in 2
    tp=(np.logical_and(ref, test))
    #True Negatives: Black in 1 and in 2
    tn=(np.logical_and(np.logical_not(ref), np.logical_not(test)))
    #False Positives: White in 2 but not in Ref
    fp=(np.logical_and(np.logical_not(ref), test))
    #False Negatives: Black in 2 but white in 1
    fn=(np.logical_and(ref, np.logical_not(test)))
    return np.count_nonzero(tp)/(np.count_nonzero(fp)+np.count_nonzero(fn)+1)
def score_2(ref_img, img):
    ref=(ref_img>np.zeros_like(ref_img))
    test=(img>np.zeros_like(img))
    #TruePositives: White in 1 and in 2
    tp=(np.logical_and(ref, test))
    #True Negatives: Black in 1 and in 2
    tn=(np.logical_and(np.logical_not(ref), np.logical_not(test)))
    #False Positives: White in 2 but not in Ref
    fp=(np.logical_and(np.logical_not(ref), test))
    #False Negatives: Black in 2 but white in 1
    fn=(np.logical_and(ref, np.logical_not(test)))
    return np.count_nonzero(tp)- np.count_nonzero(fn)

def score_3(ref_img, img):
    ref=(ref_img>np.zeros_like(ref_img))
    test=(img>np.zeros_like(img))
    return hausdorff_distance(img, ref_img)

def hausdorff_helperdist(p, img):
    """
    @param d: a point in matrix coordinates
    @param img: A binary matrix
    """
    d=min((p[0]-x)**2+(p[1]-y)**2 for x,y in it.product(range(len(img)), range(len(img))) if img[x,y] )
    return math.sqrt(d)

def hausdorff_distance(img, ref_img):
    #Source: https://de.wikipedia.org/wiki/Hausdorff-Metrik
    coms=it.product(range(len(img[0])), range(len(img[0]))) 
    h1=max(hausdorff_helperdist([x,y], img) for x,y in coms if ref_img[x,y] )
    h2=max(hausdorff_helperdist([x,y], ref_img) for x,y in it.product(range(len(img)), range(len(img))) if img[x,y] )
    return max( h1, h2)
    
def compare(ref_img, proj, rotation=0, offset_x=0, offset_y=0, show=False):
    proj.rotate(rotation)
    box=get_box(proj, 500, offset_x, offset_y)
    img,_=proj.rasterize(RASTER, bounding_square=box)
    r_z=np.zeros_like(ref_img)
    score=score_3(img, ref_img)
    if show:
        fig, axarr=plt.subplots(2,2)    
        axarr[0,0].imshow(ref_img, cmap='gray', interpolation='none')
        axarr[0,1].imshow(img, cmap='gray', interpolation='none')
        axarr[0,1].set_title(score)
        plt.show()
    """
    score=score_2(ref_img, img)
    #print ("==================\nInitial Score: {}", score)
    for zoom in [ 1/2, 1/4, 1/8 ]:
        ref_zoom=scipy.ndimage.interpolation.zoom(ref_img, zoom)
        i_zoom=scipy.ndimage.interpolation.zoom(img, zoom)
        add_score=score_2(ref_zoom, i_zoom)*zoom
        #print("Zoom {}: {}".format(zoom, add_score))
        score+=add_score
    """
    """
    score=0
    s=scipy.ndimage.generate_binary_structure(2,2)
    labelled_tp, num_f = scipy.ndimage.label(tp, s)
    if show:
        fig, axarr=plt.subplots(2,2)    
        axarr[0,0].imshow(img, cmap='gray', interpolation='none')
        axarr[0,1].imshow(ref_img, cmap='gray', interpolation='none')
        axarr[1,0].imshow(labelled_tp, interpolation='none')
    for feature in range(1, num_f+1):
        coords=np.where(labelled_tp==feature)
        xsize=max(coords[0])-min(coords[0])+1
        ysize=max(coords[1])-min(coords[1])+1
        score+=(xsize+ysize)**2
        if show: print("TP: - ({}+{})**2".format(xsize, ysize))
    labelled_fn, num_f = scipy.ndimage.label(fn, s)
    if show:
        axarr[1,1].imshow(labelled_fn, interpolation='none')
    for feature in range(1, num_f+1):
        coords=np.where(labelled_fn==feature)
        xsize=max(coords[0])-min(coords[0])+1
        ysize=max(coords[1])-min(coords[1])+1
        if show: print("FN: + ({}+{})**2".format(xsize, ysize))
        score-=(xsize+ysize)**2
    if show:
        axarr[1,1].set_title(-score)
        plt.show()
    """
    return score, img


class Optimizer:
    def __init__(self, ref_img):
        self.ref_img=ref_img
        self.start_points=self.get_start_points(60)
    def get_start_points(self, numPoints):
        """
        Return numPoints equally-distributed points on half the unit sphere.

        Implements the 2nd algorithm from https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
        """
        numPoints=2*numPoints
        a=4*math.pi/numPoints
        d=math.sqrt(a)
        Mt=int(math.pi/d)
        dt=math.pi/Mt
        df=a/dt
        points=[]
        for m in range(Mt):
            theta=math.pi*(m+0.5)/Mt
            Mf=int(2*math.pi*math.sin(theta)/df)
            for n in range(Mf):
                phi=2*math.pi*n/Mf
                points.append((math.sin(theta)*math.cos(theta), math.sin(theta)*math.sin(phi), math.cos(theta)))
        return [np.array(p) for p in points if p[0]>=0]  
    def evaluate(self, cg, strategy="M"):
        """
        @param strategy: One of "M". 
                         "M" = mixed: optimize all variables at once.
                         "O" = Other    
        """
        assert strategy in "MO"       
        self.scores=[]
        self.directions=[]
        self.rots=[]
        self.offsets=[]
        self.last_img=None
        if strategy=="M":
            return self._evaluate_M(cg)
        elif strategy=="O":
            return self._evaluate_O(cg)

    def _evaluate_O(self, cg):
        self.cg=cg
        best_score=float("inf")
        for start in self.start_points:
            score, rot, xoff, yoff=self._scoreTR([start[0], start[1], start[2]])
            if score<best_score:
                best_score=score
                best_x=[start[0], start[1], start[2], rot, xoff, yoff]
        proj=ftmp.Projection2D(self.cg, proj_direction=[best_x[0], best_x[1], best_x[2]])
        score, img = compare(self.ref_img, proj, rotation=best_x[3], offset_x=best_x[4], offset_y=best_x[5])
        self.other_img=img
        c1={'type':'eq', 'fun':lambda x: x[0]**2+x[1]**2+x[2]**2-1}
        c2={'type':'ineq', 'fun':lambda x: 180-abs(x[3])}
        c3={'type':'ineq', 'fun':lambda x: 250-abs(x[4])}
        c4={'type':'ineq', 'fun':lambda x: 250-abs(x[5])}
        opt=scipy.optimize.minimize(self._optimizeAll, best_x, constraints=[c1,c2,c3,c4], options={"maxiter":500} )
        if opt.success:
            proj=ftmp.Projection2D(self.cg, proj_direction=[opt.x[0], opt.x[1], opt.x[2]])
            score, img = compare(self.ref_img, proj, rotation=opt.x[3], offset_x=opt.x[4], offset_y=opt.x[5])
            otherscore= self._optimizeAll(opt.x)
            self.last_img=img
            assert otherscore == score == opt.fun, "{}, {}, {}".format(otherscore, score, opt.fun)
            assert abs(score-opt.fun)<1, "{},{}".format(score, opt.fun)
            #print (opt.fun)
            return opt.fun
        else:
            return 100000


    def _evaluate_M(self, cg):
        self.cg=cg
        best_score=float("inf")
        #The longest extention has been rotated in the correct direction.
        for rot in  [ 180, 175, -175, 170, -170, 90, -90, 10, -10, 5, -5, 0 ]:                
            for start in self.start_points:
                #print("X= ", [start[0], start[1], start[2], rot])
                score, xoff, yoff=self._score_withCrossCorr([start[0], start[1], start[2], rot])
                if score<best_score:
                    best_score=score
                    best_x=[start[0], start[1], start[2], rot, xoff, yoff]        
        proj=ftmp.Projection2D(self.cg, proj_direction=[best_x[0], best_x[1], best_x[2]])
        score, img = compare(self.ref_img, proj, rotation=best_x[3], offset_x=best_x[4], offset_y=best_x[5])
        self.other_img=img
        #x0,x1,x2: proj-dir; x3=rot, x4=xoffset, x5=yoffset
        c1={'type':'eq', 'fun':lambda x: x[0]**2+x[1]**2+x[2]**2-1}
        c2={'type':'ineq', 'fun':lambda x: 180-abs(x[3])}
        c3={'type':'ineq', 'fun':lambda x: 250-abs(x[4])}
        c4={'type':'ineq', 'fun':lambda x: 250-abs(x[5])}
        opt=scipy.optimize.minimize(self._optimizeAll, best_x, constraints=[c1,c2,c3,c4], options={"maxiter":500} )
        if opt.success:
            proj=ftmp.Projection2D(self.cg, proj_direction=[opt.x[0], opt.x[1], opt.x[2]])
            score, img = compare(self.ref_img, proj, rotation=opt.x[3], offset_x=opt.x[4], offset_y=opt.x[5])
            otherscore= self._optimizeAll(opt.x)
            self.last_img=img
            assert otherscore == score == opt.fun, "{}, {}, {}".format(otherscore, score, opt.fun)
            assert abs(score-opt.fun)<1, "{},{}".format(score, opt.fun)
            #print (opt.fun)
            return opt.fun
        else:
            return 100000

    def _scoreTR(self,x):
        assert len(x)==3
        best_opt=float("inf")
        for startrot in [0, 180]:
            c2={'type':'ineq', 'fun':lambda x: 180-abs(x[0])}
            c3={'type':'ineq', 'fun':lambda x: 250-abs(x[1])}
            c4={'type':'ineq', 'fun':lambda x: 250-abs(x[2])}
            opt=scipy.optimize.minimize(self._optimizeTR, [startrot, 0, 0], constraints=[c2,c3,c4], options={"maxiter":50}, args=(x[:3]) )
            if opt.fun<best_opt:
                best_opt=opt.fun
                best_x=opt.x
        return best_opt, best_x[0], best_x[1], best_x[2]
    def _optimizeTR(self, x, direct):
        proj=ftmp.Projection2D(self.cg, proj_direction=[direct[0], direct[1], direct[2]])
        score, img = compare(self.ref_img, proj, rotation=x[0], offset_x=x[1], offset_y=x[2])
        return score

    def _score_withCrossCorr(self,x):
        assert len(x)==4
        proj=ftmp.Projection2D(self.cg, proj_direction=[x[0], x[1],x[2]])
        proj.rotate(x[3])
        box=get_box(proj, 500, 0,0)
        img,_=proj.rasterize(RASTER, bounding_square=box)        
        corr=scipy.signal.correlate2d(ref_img, img)
        off_y, off_x = np.unravel_index(np.argmax(corr), corr.shape)
        off_y=off_y-img.shape[0]+1
        off_x=off_x-img.shape[1]+1
        off_x*=500/len(img) #convert to Angstrom
        off_y*=500/len(img)
        proj=ftmp.Projection2D(self.cg, proj_direction=[x[0], x[1],x[2]])
        score, img = compare(self.ref_img, proj, rotation=x[3], offset_x=off_x, offset_y=off_y)
        return score, off_x, off_y

    def _optimizeAll(self, x):
        assert len(x)==6
        proj=ftmp.Projection2D(self.cg, proj_direction=[x[0], x[1],x[2]])
        score, img = compare(self.ref_img, proj, rotation=x[3], offset_x=x[4], offset_y=x[5])
        return score

# Parser is made available if this file is loaded as a module, to allow
# the sphinxcontrib.autoprogram extension to work properly
parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    cg1=ftmc.CoarseGrainRNA(args.files[0])
    cg2=ftmc.CoarseGrainRNA(args.files[1])
    #direction=args.direction.split(",")
    #if len(direction)!=3:
    #    parser.error("The projection direction should have 3 coordinates. e.g.: '1.0,0.0,0.0'")

    #Minimize over rotation and offset_x and offset_y (and proj-dir?)


    scores=[]
    rots=[-180, -90, -45, -22, -10, -5, -2, 0, 2 ,5, 10, 22, 45, 90, 180]
    ref_proj=ftmp.Projection2D(cg1)    
    ref_box=get_box(ref_proj)
    ref_img, _=ref_proj.rasterize(RASTER, bounding_square=ref_box)
    proj=ftmp.Projection2D(cg2, proj_direction=[1,1,1])
    """for rot in rots:
        score, img = compare(ref_img, proj, rotation=rot)
        scores.append(score)

    offsets=[-100, -50, -25, -10, -5, -2, -1, 0, 1, 2, 5, 10, 25, 50, 100]
    o_scores=[]
    for xoff in offsets:
        score, img = compare(ref_img, proj, offset_x=xoff)
        o_scores.append(score)
    """
    optimizer=Optimizer(ref_img)
    print("Final score is:", optimizer.evaluate(cg2, "M"))
    #Plotting 
    #fig, axarr=plt.subplots(2,2)    
    #axarr[0,0].imshow(ref_img, cmap='gray', interpolation='none')
    #axarr[0,1].imshow(optimizer.last_img, cmap='gray', interpolation='none')
    #axarr[1,0].plot(np.arange(len(optimizer.scores)), optimizer.scores )
    #corr=scipy.signal.correlate2d(optimizer.last_img, ref_img)
    #axarr[1,0].imshow(corr)
    #axarr[1,1].imshow(optimizer.other_img, cmap='gray', interpolation='none')
    #plt.show()
