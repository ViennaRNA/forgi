import pkgutil
import logging
import warnings
try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import numpy as np

import pandas as pd

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.neighbors.kde import KernelDensity
from sklearn.metrics import confusion_matrix

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug

import matplotlib.pyplot as plt

log = logging.getLogger(__name__)


CUTOFFDIST=30
ANGLEWEIGHT=10

"""
This module contains code for classifying a coarse-grained geometry as A-Minor
interaction.

:warning: This is intended for low-resolution data , as it does not take the
          orientation of individual bases into account. If you have all-atom
          data, dedicated tools like FR3D will be more accurate.

If you just want to classify interactions in a given structure, you only need
the functions `classify_interaction` or `all_interactions`lassify_interaction`.

To train your own version of the classifier, modify its parameters or perform
cross-validation, use the AMinorClassifier.

Access the default trainings data with the get_trainings_data(loop_type)
function.
"""


P_INTERACTION=0.05 #0.03644949066213922

################# Just classify my structure ############################

def classify_interaction(cg, stem, loop, clf=None):
    """
    :param clf: A fitted AMinorClassifier instance. If not given,
                the default pretrained classifier is used.
    """
    if clf in None:
        clf=_get_default_clf(loop)
    geo = np.array(get_relative_orientation(cg, loop, stem))
    geo.reshape(1, -1)
    y,=clf.predict(geo)
    return y

def all_interactions(cg, clfs=None):
    """
    Get a list of all predicted A-Minor interactions in a cg-object.

    This is more efficient than using classify_interaction iteratively,
    because it uses vectorization.
    :param clfs: A dictionary {loop_type: AMinorClassifier} where
                 loop_type is one of "i", "h", "m".
                 If clfs is None or a key is missing, uses the default
                 pretrained classifier.
    """
    interactions=[]
    for loop_type in ["m", "i", "h"]:
        if clfs is None or loop_type not in clfs:
            clf=_get_default_clf(loop_type)
        else:
            clf=clfs[loop_type]
        labels=[]
        geos=[]
        for stem in cg.stem_iterator():
            for loop in cg.defines:
                if 'A' not in "".join(cg.get_define_seq_str(loop)):
                    continue
                if stem in cg.edges[loop]:
                    continue
                if loop[0]!=loop_type:
                    continue
                # To save computation time
                if not ftuv.elements_closer_than(cg.coords[loop][0],
                                                 cg.coords[loop][1],
                                                 cg.coords[stem][0],
                                                 cg.coords[stem][1],
                                                 CUTOFFDIST):
                    continue
                geos.append(get_relative_orientation(cg, loop, stem))
                labels.append([loop, stem])
        y = clf.predict(geos)
        log.info("Classifying %s", labels)
        log.info("# hits=%s", sum(y))
        if sum(y)==0:
            log.info("geos=%s", geos)
            log.info("probas=%s", clf.predict_proba(geos))
        labels=np.array(labels)
        interactions.extend(labels[y])
    return np.array(interactions)

############# Roll your own classifier ###########################
def get_trainings_data( loop):
    loop_type=loop[0]
    return _DefaultClf._get_data( loop_type)

class AMinorClassifier(BaseEstimator, ClassifierMixin):
    """
    A classifier that predicts A-Minor interactions based on the following formula:

    .. math:: P(I|geo)=f(geo|I)/f(geo)*P(I)

    where :math:`f` is a probability density and :math:`P` is a probability
    and :math:`I` means interaction.

    Since there are no interactions with more than 30 Angstrom distance,
    interactions require an A in the loop sequence and estimating numerator
    and denominator seperately could lead to probabilities greater than 1,
    we use the following formula:

    .. math:: P(I| (geo, d<30, A in seq)) = num/denom

    where
    .. math::

        num=f(geo|I)*P(I|(d<30, A \\in seq))

        denom=num+f(geo|(\\not I, d<30))*(1-P(I|(d<30, A \\in seq))

    We estimate the probability densities as kernel density estimates.
    For :math:`f(geo|I)` we use annotations created by FR3D-searches for
    all 4 types of A-minor motifs.
    We assume that :math:`f(geo|(\\not I, A \\in seq))== f(geo|(\\not I, A \\in seq))`
    and that FR3D might miss some true interactions. Thus we estimate
    :math:`f(geo|(\\not I, A \\in seq, d<30))` as
    :math:`f(geo|(\\not I, A \\notin seq, d<30))`

    Since both densities are normalized, it does not matter that we use a
    different number of datapoints for their estimation.

    We estimate P(I|(d<30, A \\in seq)) as the number of occurrences where
    a loop with A is in an A-minor interaction, over the number of all
    loop-stem pairs with less than 30 angstrom distance and an A in the sequence.
    """
    def __init__(self, kernel="linear", bandwidth=0.09, symmetric=True, p_I=P_INTERACTION):
        self.p_I = p_I
        self.symmetric=symmetric
        self.kernel = kernel
        self.bandwidth = bandwidth

    def fit(self, X, y):
        """
        Train the model.

        :param X: A Nx3 array, where the features are
                  distance(Angstrom)/10, angle1(rad), angle2(rad)
                  The **distance** is the closest distance between the two line
                  segments (i.e. coarse grained elementts)
                  **angle1** is the angle between the line along the stem vector
                  and the line along the shortest connection between the two
                  elements. A an angle between two straight lines, it is
                  defined between 0 and 90 degrees.
                  **angle2** is the angle between the connecting vector
                  (pointing from the stem to the loop), projected onto the
                  plane normal to the stem direction and the twist vector
                  (location of minor groove) at the point closest to the
                  interaction. As an angle between two vectors, it is
                  defined between 0 and 180 degrees.
        :param y: An array of length N. 0 means no interaction,
                  1 means interaction.
        """
        # Check that X and y have correct shape
        X, y = check_X_y(X, y)
        log.info("Trainings-data has shape %s", X.shape)
        log.info("We have %s known interactions ", sum(y))
        if X.shape[1]!=3:
            raise TypeError("Expect exactly 3 features, found {}".format(X.shape[1]))
        if not all(yi in [0,1] for yi in y):
            raise ValueError("y should only contain the values 1 and 0")
        X[:,0]/=ANGLEWEIGHT
        #print("First datapoint", X[0], y[0])
        ame = X[np.where(y)]
        non_ame=X[np.where(y==0)]
        if self.symmetric:
            ame=self._make_symmetric(ame)
            non_ame = self._make_symmetric(non_ame)
        #print(ame)
        log.info("Fitting\n%s", X[np.where(y)])
        self.ame_kde_ = KernelDensity(kernel=self.kernel,
                                      bandwidth=self.bandwidth).fit(ame).score_samples
        self.non_ame_kde_ = KernelDensity(kernel=self.kernel,
                                          bandwidth=self.bandwidth).fit(non_ame).score_samples
        #fig,ax = plt.subplots()
        #plt.plot(np.linspace(0, np.pi/2, 200), self.ame_kde_(np.array([[1.0,a1,1.5] for a1 in np.linspace(0, np.pi/2, 200)])))
        #plt.plot(np.linspace(0, np.pi/2, 200), self.non_ame_kde_(np.array([[1.0,a1,1.5] for a1 in np.linspace(0, np.pi/2, 200)])))
        #plt.plot(np.linspace(0, np.pi/2, 200), self.predict_proba(np.array([[1.0,a1,1.5] for a1 in np.linspace(0, np.pi/2, 200)])))
        #plt.show()
    @staticmethod
    def _make_symmetric(geos):
        """
        Make the trainingsdata symmetric around those bounds of the radial
        variables, where we expect a high density.

        Kernel density estimation does not work too well out of the box on
        bounded support: Samples close to the boundary have create a
        significant density outside the support. Due to symmetry considerations
        of angular data, we can avoid this issue by adding symmetric datapoints
        the following way: Mirror all datapoints along the plane
        angle1=180 degrees (both angles between the two lines are equivalent) and
        around the plane angle2==0 degrees (both rotational
        directions are equivalent).
        By creating 4 times as many datapoints, we avoid bias near the boundary.
        """
        geos = np.concatenate([geos,
                               [(d, np.pi-a1, a2) for d, a1, a2 in geos]])
        geos = np.concatenate([geos,
                               [(d, a1, -a2) for d, a1, a2 in geos]])
        return geos

    def predict(self, X):
        return self.predict_proba(X)>0.5

    def score(self, X, y):
        """
        The average between specificity and sensitivity
        """
        y_pred = self.predict(X)
        tn, fp, fn, tp = confusion_matrix(y, y_pred).ravel()
        specificity =  tn/(tn+fp)
        sensitivity = tp/(tp+fn)
        return 0.5*specificity + 0.5*sensitivity
        #precision = tp/(tp+fp)
        return f1_score(y, y_pred)
    def predict_proba(self, X):
        check_is_fitted(self, ['ame_kde_', 'non_ame_kde_'])
        X = check_array(X)
        #print("PI", self.p_I, "bw", self.bandwidth)
        #print()
        X[:,0] /= ANGLEWEIGHT
        numerator = np.exp(self.ame_kde_(X))*self.p_I
        denom = numerator+np.exp(self.non_ame_kde_(X))*(1-self.p_I)
        #print(np.exp(self.non_ame_kde_(X)))
        log.info("%s/%s=%s", numerator, denom, np.nan_to_num(numerator/denom))
        with warnings.catch_warnings(): # division by 0
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nan_to_num(numerator/denom)


############## get orientation #########################

def get_relative_orientation(cg, loop, stem):
    '''
    Return how loop is related to stem in terms of three parameters.

    The stem is the receptor of a potential A-Minor interaction, whereas the
    loop is the donor.

    The 3 parameters are:

        1. Distance between the closest points of the two elements
        2. The angle between the stem and the vector between the two
        3. The angle between the minor groove of l2 and the projection of
        the vector between stem and loop onto the plane normal to the stem
        direction.
    '''
    point_on_stem, point_on_loop = ftuv.line_segment_distance(cg.coords[stem][0],
                                                              cg.coords[stem][1],
                                                              cg.coords[loop][0],
                                                              cg.coords[loop][1])
    conn_vec = point_on_loop-point_on_stem
    dist = ftuv.magnitude(conn_vec)
    angle1 = ftuv.vec_angle(cg.coords.get_direction(stem),
                            conn_vec)
    # The direction of the stem vector is irrelevant, so
    # choose the smaller of the two angles between two lines
    if angle1>np.pi/2:
        angle1 = np.pi-angle1
    tw = cg.get_twists(stem)
    if dist==0:
        angle2=float("nan")
    else:
        if stem[0] != 's':
            raise ValueError("The receptor needs to be a stem, not {}".format(stem))
        else:
            stem_len = cg.stem_length(stem)
            # Where along the helix our A-residue points to the minor groove.
            # This can be between residues. We express it as floating point nucleotide coordinates.
            # So 0.0 means at the first basepair, while 1.5 means between the second and the third basepair.
            pos = ftuv.magnitude(point_on_stem - cg.coords[stem][0]) / ftuv.magnitude(cg.coords.get_direction(stem)) * (stem_len - 1)
            # The vector pointing to the minor groove, even if we are not at a virtual residue (pos is a float value)
            virt_twist = ftug.virtual_res_3d_pos_core(cg.coords[stem], cg.twists[stem], pos, stem_len)[1]
            # The projection of the connection vector onto the plane normal to the stem
            conn_proj = ftuv.vector_rejection(conn_vec, cg.coords.get_direction(stem))
            try:
                # Note: here the directions of both vectors are well defined,
                # so angles >90 degrees make sense.
                angle2 = ftuv.vec_angle(virt_twist, conn_proj)
            except ValueError:
                if np.all(virt_twist==0):
                    angle2=float("nan")
                else:
                    raise
    return dist, angle1, angle2

############# Private ##########################################
class _DefaultClf(object):
    """Just to put global (cached) variables into a seperate scope."""
    _clfs={} # The pretrained classifiers, lazily loaded.
    _geo_df = None
    @classmethod
    def get_default_clf(cls, loop_type):
        if loop_type not in cls._clfs:
            clf = AMinorClassifier()
            X, y = cls._get_data(loop_type)
            clf.fit(X, y)
            cls._clfs[loop_type]=clf
        return cls._clfs[loop_type]
    @classmethod
    def get_dataframe(cls):
        if cls._geo_df is None:
            rawdata = pkgutil.get_data('forgi', 'threedee/data/aminor_geometries.csv')
            cls._geo_df = pd.read_csv(StringIO(rawdata.decode("ascii")), comment="#", sep=" ")
        return cls._geo_df

    @classmethod
    def _get_data(cls, loop_type):
        df = cls.get_dataframe()
        print(df.columns.values)
        df=df[df.loop_type==loop_type]
        df=df[df.dist<CUTOFFDIST]
        positive = df[df["is_interaction"]]
        negative = df[(df["is_interaction"]==False)&(~df["loop_sequence"].str.contains("A").astype(bool))]
        data = np.concatenate( [
                [[ x.dist, x.angle1, x.angle2 ]
                          for x in positive.itertuples()],
                [[ x.dist, x.angle1, x.angle2 ]
                          for x in negative.itertuples()]])
        return data, np.array([1]*len(positive)+[0]*len(negative))

def _get_default_clf(loop):
    """
    :param loop: A element name, e.g. "i0" or a loop-type, e.g. "i"
    """
    loop_type=loop[0]
    return _DefaultClf.get_default_clf(loop_type)
