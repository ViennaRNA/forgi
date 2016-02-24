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
import forgi.projection.projection2d as ftmp
 
import matplotlib.pyplot as plt
import scipy.signal
import scipy.ndimage
import scipy.misc

RASTER=55
WIDTH=520
class HausdorffError(ArithmeticError):
    pass

class Hausdorff:
    def __init__(self, length):
        self.offsets=self.get_griddist_iterator(length)
    @staticmethod
    def get_griddist_iterator(length):
        dists={}
        for dx in range(-length-1, length+1):
            for dy in range(-length-1, length+1):
                dists[(dx,dy)]=dx**2+dy**2
        return sorted(dists.keys(), key=lambda x: dists[x])

    def hausdorff_helperdist(self, p, img):
        """
        Returns the shorthest distance from a given point p to any non-zero cell in img.

        @param p: a point in matrix coordinates
        @param img: A binary matrix
        """
        for (dx,dy) in self.offsets:
            try:
                if img[p[0]+dx, p[1]+dy]:
                    return math.sqrt(dx**2+dy**2)
            except IndexError: 
                pass
        else:
            #this will rarely occur, which is why we do not check for it before the for loop.
            if not np.any(img): 
                return float('inf')
            else:
                raise HausdorffError("Cannot calculated distance between point ({},{}) and image {}".format(p[0], p[1], img))

    def hausdorff_distance(self, img, ref_img):
        #Source: https://de.wikipedia.org/wiki/Hausdorff-Metrik
        #print(img)
        #print(ref_img, np.min(ref_img), np.max(ref_img))
        h1=max(self.hausdorff_helperdist([x,y], img) for x,y in np.transpose(np.where(ref_img) ))
        h2=max(self.hausdorff_helperdist([x,y], ref_img) for x,y in np.transpose(np.where(img) ))
        return max( h1, h2)

#GLOBAL VARIABLE FOR NOW
hdCalc=Hausdorff(RASTER)

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


def get_box(projection, width=WIDTH, offset_x=0, offset_y=0):
    left, right, down, up=projection.get_bounding_square()
    #print(up-down)
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
    return hdCalc.hausdorff_distance(ref, test)
    
def compare(ref_img, proj, rotation=0, offset_x=0, offset_y=0, show=False):
    proj.rotate(rotation)
    box=get_box(proj, WIDTH, offset_x, offset_y)
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

def to_polar(x):
  theta=math.atan2(x[1],x[0])
  phi=math.atan2(math.sqrt(x[0]**2+x[1]**2),x[2])
  return theta, phi
def from_polar(theta, phi):
  return [math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)]
class Optimizer:
    def __init__(self, ref_img):
        self.ref_img=ref_img
        self.start_points=self.get_start_points(60)
    def get_start_points(self, numPoints):
        """
        Return numPoints equally-distributed points on half the unit sphere in polar coordinates.

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
                points.append((theta, phi))
        return [np.array(p) for p in points if p[0]<=math.pi/2]

    def evaluate(self, cg, strategy="O"):
        """
        @param strategy: One of "M". 
                         "M" = mixed: optimize all variables at once.
                         "O" = Other    
        """
        return self._evaluate_O(cg)

    def _evaluate_O(self, cg):
        self.cg=cg
        best_score=float("inf")
        start=to_polar(cg.project_from)
        for start in self.start_points+[start]:
            score, rot, xoff, yoff=self._scoreTR([start[0], start[1]]) #Minimize over translation and rotation
            if score<best_score:
                best_score=score
                best_x=[start[0], start[1], rot, xoff, yoff]
        proj=ftmp.Projection2D(self.cg, proj_direction=from_polar(best_x[0], best_x[1]))
        score, img = compare(self.ref_img, proj, rotation=best_x[2], offset_x=best_x[3], offset_y=best_x[4])
        print("Intermediate score", score)
        self.other_img=img
        # x1= phi, x0=theta
        #bounds=[(0,math.pi/2),(0, math.pi*2),(-90,270),(-250,250),(-250,250)]
        opt=scipy.optimize.minimize(self._optimizeAll, best_x, options={"maxiter":500}, method="Powell" ) #or COBYLA
        if opt.success:
            proj=ftmp.Projection2D(self.cg, proj_direction=from_polar(opt.x[0],opt.x[1]))
            score, img = compare(self.ref_img, proj, rotation=opt.x[2], offset_x=opt.x[3], offset_y=opt.x[4])
            self.last_img=img
            print ("After final optimization:", opt.x)
            #print(opt)
            assert score == opt.fun
            return opt.fun
        else:
            print(opt)
            return 100000


    def _scoreTR(self,x):
        assert len(x)==2
        best_opt=float("inf")
        for startrot in [0, 180]: 
            #bounds=[(-90,270),(-250,250),(-250,250)],
            opt=scipy.optimize.minimize(self._optimizeTR, [startrot, 0, 0], options={"maxiter":50}, args=(x), method="Powell" )
            if opt.fun<best_opt:
                best_opt=opt.fun
                best_x=opt.x
                #print(opt)
        #print ("Best translation and rotation: ", best_x)
        return best_opt, best_x[0], best_x[1], best_x[2]
    def _optimizeTR(self, x, direct):
        proj=ftmp.Projection2D(self.cg, proj_direction=from_polar(*direct))
        score, img = compare(self.ref_img, proj, rotation=x[0], offset_x=x[1], offset_y=x[2])
        return score

    """
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
    """

    def _optimizeAll(self, x):
        assert len(x)==5
        proj=ftmp.Projection2D(self.cg, proj_direction=from_polar(x[0], x[1]))
        score, img = compare(self.ref_img, proj, rotation=x[2], offset_x=x[3], offset_y=x[4])
        return score

# Parser is made available if this file is loaded as a module, to allow
# the sphinxcontrib.autoprogram extension to work properly
parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()    
    cg2=ftmc.CoarseGrainRNA(args.files[1])    
    proj=ftmp.Projection2D(cg2)#, proj_direction=[1,1,1])
    if args.files[0][-3:]=="png":
        ref_img=scipy.misc.imread(args.files[0], flatten=True)
        ref_img=scipy.misc.imresize(ref_img,(50,50), "nearest")
        ref_box=get_box(proj)
    else:
        cg1=ftmc.CoarseGrainRNA(args.files[0])    
        ref_proj=ftmp.Projection2D(cg1)    
        ref_img, _=ref_proj.rasterize(RASTER, bounding_square=ref_box)    
        ref_box=get_box(ref_proj)


    scores=[]
    rots=[-180, -90, -45, -22, -10, -5, -2, 0, 2 ,5, 10, 22, 45, 90, 180]



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
    #print("boundingSquare", proj.get_bounding_square())
    score=optimizer.evaluate(cg2, "M")
    print("Final score is:", score)
    #Plotting 
    fig, axarr=plt.subplots(2,2)    
    axarr[0,0].imshow(ref_img, cmap='gray', interpolation='none')
    axarr[0,1].imshow(optimizer.last_img, cmap='gray', interpolation='none')
    axarr[0,1].set_title(score)
    axarr[1,0].imshow(proj.rasterize(RASTER, bounding_square=ref_box)[0], cmap='gray', interpolation='none')
    #corr=scipy.signal.correlate2d(optimizer.last_img, ref_img)
    #axarr[1,0].imshow(corr)
    axarr[1,1].imshow(optimizer.other_img, cmap='gray', interpolation='none')
    plt.show()
