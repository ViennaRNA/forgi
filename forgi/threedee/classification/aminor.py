
"""
This module contains code for classifying a coarse-grained geometry as A-Minor
interaction.

.. warning:: This is intended for low-resolution data , as it does not take the
             orientation of individual bases into account. If you have all-atom
             data, dedicated tools like FR3D will be more accurate.

If you just want to classify interactions in a given structure, you only need
the functions `classify_interaction` or `all_interactions`.

To train your own version of the classifier, modify its parameters or perform
cross-validation, use the AMinorClassifier.

Access the default trainings data with the get_trainings_data(loop_type)
function.
"""

from __future__ import print_function, division, absolute_import

import pkgutil
import logging
import warnings
import json

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import numpy as np

import pandas as pd

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
try:
    from sklearn.neighbors import KernelDensity
except:
    from sklearn.neighbors.kde import KernelDensity

from sklearn.metrics import confusion_matrix

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug


log = logging.getLogger(__name__)


CUTOFFDIST = 30
ANGLEWEIGHT = 10



P_INTERACTION = 0.05  # 0.03644949066213922

################# Just classify my structure ############################


def loop_potential_interactions(cg, loop, domain=None):
    """
    Iterate over all stems and return those loop-stem pairs that will be passed
    to the AMinor classification and are not ruled out beforehand.
    """
    geos = []
    labels = []
    if domain is not None:
        stems = (s for s in domain if s[0] == "s")
    else:
        stems = cg.stem_iterator()
    for stem in stems:
        if stem in cg.edges[loop]:
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
    geos = np.array(geos)
    if len(geos) > 0:
        geos[:, 0] /= ANGLEWEIGHT
    return geos, labels


def potential_interactions(cg, loop_type, domain=None):
    """
    :returns: A tuple `geos`, `labels`.
              `geos` is an Nx3 array, where N is the number of
              potential interactions and the inner dimension is dist, angle1, angle2.
              `geos` can be passed to the AMinor classifier.
              `labels` is an Nx2 array, where the inner dimension is loopname, stemname
    """
    labels = []
    geos = []
    for loop in cg.defines:
        if domain is not None and loop not in domain:
            continue
        if loop[0] != loop_type:
            continue
        if 'A' not in "".join(cg.get_define_seq_str(loop)):
            continue
        loop_geos, loop_labels = loop_potential_interactions(cg, loop, domain)
        geos.extend(loop_geos)
        labels.extend(loop_labels)
    return np.array(geos), np.array(labels)


def all_interactions(cg, clfs=None):
    """
    Get a list of all predicted A-Minor interactions in a cg-object.

    This is more efficient than using classify_interaction iteratively,
    because it uses vectorization.

    :param clfs: A dictionary {loop_type: AMinorClassifier} where
                 loop_type is one of "i", "h", "m".
                 If clfs is None or a key is missing, uses the default
                 pretrained classifier.

    :returns: A list of  tuples (loop, stem)
    """
    interactions = []
    for loop_type in ["i", "h"]:
        if clfs is not None and loop_type not in clfs:
            warnings.warn("No classifier specified for loop type %s (only for %s), "
                          "using default classifier.", loop_type, ",".join(clfs.keys()))
        if clfs is None or loop_type not in clfs:
            clf = _get_default_clf(loop_type)
        else:
            clf = clfs[loop_type]
        geos, labels = potential_interactions(cg, loop_type)
        interactions.extend(
            _classify_potential_interactions(clf, geos, labels))
    return interactions


def classify_interaction(cg, loop, stem=None, clf=None):
    """
    Returns the interaction pair loop, stem as a tuple or False if no interaction exists.

    :param cg: The CoarseGrainRNA
    :param loop: The loop name, e.g. "i0"
    :param stem: A stem name, e.g. "s0" to consider interactions of loop with another stem as False
                 and only an interaction with this stem as True.

                ..warning:: Our statistical modelling allows for at most 1 interaction per loop.
                            This means that we have calculate the interaction probability of this
                            loop with all stems, even if stem is given. If another stem has a
                            higher interaction probability than the given stem,
                            this function will return False, regardless of the interaction
                            probability stem-loop.

    :param clf: A trained AMinorClassifier or None (use default classifier for loop type)


    ..note:: `all_interactions` is more efficient, if you are interested in all loop-stem pairs.

    :returns: A tuple (loop, stem) or False
    """
    if clf is None:
        clf = _get_default_clf(loop[0])
    geos, labels = loop_potential_interactions(cg, loop)
    interactions = _classify_potential_interactions(clf, geos, labels)
    if not interactions:
        return False
    if stem is not None and interactions != [(loop, stem)]:
        return False
    interaction, = interactions
    return interaction

############# Roll your own classifier ###########################


def get_trainings_data(loop):
    loop_type = loop[0]
    return _DefaultClf._get_data(loop_type)


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

    def __init__(self, kernel="linear", bandwidth=0.3, symmetric=True, p_I=P_INTERACTION):
        self.p_I = p_I
        self.symmetric = symmetric
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
        if X.shape[1] != 3:
            raise TypeError(
                "Expect exactly 3 features, found {}".format(X.shape[1]))
        if not all(yi in [0, 1] for yi in y):
            raise ValueError("y should only contain the values 1 and 0")
        ame = X[np.where(y)]
        non_ame = X[np.where(y == 0)]
        if self.symmetric:
            ame = self._make_symmetric(ame)
            non_ame = self._make_symmetric(non_ame)
        log.info("Fitting. First positive sample: %s", X[np.where(y)][0])
        self.ame_kde_ = KernelDensity(kernel=self.kernel,
                                      bandwidth=self.bandwidth).fit(ame).score_samples
        self.non_ame_kde_ = KernelDensity(kernel=self.kernel,
                                          bandwidth=self.bandwidth).fit(non_ame).score_samples
        self.X_ = X
        self.y_ = y

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
                               [(d, np.pi - a1, a2) for d, a1, a2 in geos]])
        # We allow negative angles for angle2, so there is no need to make it symmetric.
        # geos = np.concatenate([geos,
        #                       [(d, a1, -a2) for d, a1, a2 in geos]])
        return geos

    def predict(self, X):
        return self.predict_proba(X) > 0.5

    def score(self, X, y):
        """
        The average between specificity and sensitivity
        """
        y_pred = self.predict(X)
        tn, fp, fn, tp = confusion_matrix(y, y_pred).ravel()
        specificity = tn / (tn + fp)
        sensitivity = tp / (tp + fn)
        return (0.8 * sensitivity + 0.2 * specificity)

    def predict_proba(self, X):
        check_is_fitted(self, ['ame_kde_', 'non_ame_kde_'])
        X = check_array(X)
        log.debug("Predicting for %s", X)
        numerator = np.exp(self.ame_kde_(X)) * self.p_I
        denom = numerator + np.exp(self.non_ame_kde_(X)) * (1 - self.p_I)
        with warnings.catch_warnings():  # division by 0
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nan_to_num(numerator / denom)

    def set_params(self, **kwargs):
        """"""
        super(AMinorClassifier, self).set_params(**kwargs)
        # If it was fitted, we must propagate the parameter changes
        # to the child-KDEs by refitting to the same data
        if hasattr(self, "X_"):
            self.fit(self.X_, self.y_)
        return self


############## get orientation #########################

def get_loop_flexibility(cg, loop):
    """
    Unused. We tried to see if the length of the loop vs # bases had an effect on ointeraction probability.
    """
    assert loop[0] == "i"
    d = cg.define_a(loop)
    nt1, nt2 = d[1] - d[0], d[3] - d[2]
    max_nts = max(nt1, nt2)
    loop_length = ftuv.magnitude(cg.coords.get_direction(loop))
    # As number of nucleotide-links (or phosphate groups) per Angstrom
    # 9.2 is the sum of average bond lengths for bonds in the nucleotide linkage.
    # Bond lengths taken from: DOI: 10.1021/ja9528846
    # A value of 1 means, all bonds are stretched.
    # Ideal helices have a value of: 4.41
    # A value below 1 should be rare.
    # Higher values mean higher flexibility.
    return (max_nts) / loop_length * 9.2


def get_relative_orientation(cg, loop, stem):
    '''
    Return how loop is related to stem in terms of three parameters.

    The stem is the receptor of a potential A-Minor interaction, whereas the
    loop is the donor.

    The 3 parameters are:

        1.  Distance between the closest points of the two elements
        2.  The angle between the stem and the vector between the two
        3.  The angle between the minor groove of l2 and the projection of
            the vector between stem and loop onto the plane normal to the stem
            direction.
    '''
    point_on_stem, point_on_loop = ftuv.line_segment_distance(cg.coords[stem][0],
                                                              cg.coords[stem][1],
                                                              cg.coords[loop][0],
                                                              cg.coords[loop][1])
    conn_vec = point_on_loop - point_on_stem
    dist = ftuv.magnitude(conn_vec)
    angle1 = ftuv.vec_angle(cg.coords.get_direction(stem),
                            conn_vec)
    # The direction of the stem vector is irrelevant, so
    # choose the smaller of the two angles between two lines
    if angle1 > np.pi / 2:
        angle1 = np.pi - angle1
    tw = cg.get_twists(stem)
    if dist == 0:
        angle2 = float("nan")
    else:
        if stem[0] != 's':
            raise ValueError(
                "The receptor needs to be a stem, not {}".format(stem))
        else:
            stem_len = cg.stem_length(stem)
            # Where along the helix our A-residue points to the minor groove.
            # This can be between residues. We express it as floating point nucleotide coordinates.
            # So 0.0 means at the first basepair, while 1.5 means between the second and the third basepair.
            pos = ftuv.magnitude(point_on_stem - cg.coords[stem][0]) / ftuv.magnitude(
                cg.coords.get_direction(stem)) * (stem_len - 1)
            # The vector pointing to the minor groove, even if we are not at a virtual residue (pos is a float value)
            virt_twist = ftug.virtual_res_3d_pos_core(
                cg.coords[stem], cg.twists[stem], pos, stem_len)[1]
            # The projection of the connection vector onto the plane normal to the stem
            conn_proj = ftuv.vector_rejection(
                conn_vec, cg.coords.get_direction(stem))
            try:
                # Note: here the directions of both vectors are well defined,
                # so angles >90 degrees make sense.
                angle2 = ftuv.vec_angle(virt_twist, conn_proj)
            except ValueError:
                if np.all(virt_twist == 0):
                    angle2 = float("nan")
                else:
                    raise
            # Furthermore, the direction of the second angle is meaningful.
            # We call use a positive angle, if the cross-product of the two vectors
            # has the same sign as the stem vector and a negative angle otherwise
            cr = np.cross(virt_twist, conn_proj)
            sign = ftuv.is_almost_parallel(cr,  cg.coords.get_direction(stem))
            #assert sign != 0, "{} vs {} not (anti) parallel".format(
            #    cr, cg.coords.get_direction(stem))
            angle2 *= sign

    return dist, angle1, angle2

############# Private ##########################################


def _get_masks(df, loop_type):
    lt_mask = df.loop_type == loop_type
    cutoff_mask = df.dist < CUTOFFDIST
    in_support = lt_mask & cutoff_mask
    has_a = df.loop_sequence.str.contains("A")
    is_i = df.is_interaction
    mask_ame = in_support & is_i & has_a
    mask_non_ame = in_support & np.invert(is_i | has_a)
    mask_non_fred = in_support & np.invert(is_i) & has_a
    return mask_ame, mask_non_ame, mask_non_fred


def _classify_potential_interactions(clf, geos, labels):
    if len(geos) == 0:
        return []
    score = clf.predict_proba(geos)
    y = score > 0.5
    log.info("Classifying %s", labels)
    log.info("score %s", score)
    log.info("# hits (not unique)=%s", sum(y))
    # We modelled the probabilities in a way that each loop can have only 1 interaction.
    # So for multiple interactions, use only the best one.
    best_interactions = {}
    for i, label in enumerate(labels):
        loop = label[0]
        if score[i] < 0.5:
            continue
        if loop not in best_interactions or score[i] > best_interactions[loop][0]:
            best_interactions[loop] = (score[i], label)
    log.info("%s unique interacting loops.", len(best_interactions))
    log.debug(best_interactions)
    return [best_i[1] for best_i in best_interactions.values()]


class _DefaultClf(object):
    """Just to put global (cached) variables into a seperate scope."""
    _clfs = {}  # The pretrained classifiers, lazily loaded.
    _geo_df = None

    @classmethod
    def get_default_clf(cls, loop_type):
        if loop_type not in cls._clfs:
            clf = AMinorClassifier()
            rawdata = pkgutil.get_data(
                'forgi', 'threedee/data/aminor_params.json')
            params = json.load(StringIO(rawdata.decode("ascii")))
            clf.set_params(**params[loop_type])
            X, y = cls._get_data(loop_type)
            clf.fit(X, y)
            cls._clfs[loop_type] = clf
        return cls._clfs[loop_type]

    @classmethod
    def get_dataframe(cls):
        if cls._geo_df is None:
            rawdata = pkgutil.get_data(
                'forgi', 'threedee/data/aminor_geometries.csv')
            cls._geo_df = pd.read_csv(
                StringIO(rawdata.decode("ascii")), comment="#", sep=" ")
        return cls._geo_df

    @classmethod
    def _get_data(cls, loop_type):
        df = cls.get_dataframe()
        return df_to_data_labels(df, loop_type)


def _get_default_clf(loop):
    """
    :param loop: A element name, e.g. "i0" or a loop-type, e.g. "i"
    """
    loop_type = loop[0]
    return _DefaultClf.get_default_clf(loop_type)


def df_to_data_labels(df, loop_type):
    """
    Create the trainings data as two arrays X and y (or data and labels)
    from the initial dataframe

    :returns: X, y
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df["loop_name"] = df.pdb_id + \
            df.interaction.str.split("-").apply(lambda x: x[0])
    mask_ame, mask_non_ame, mask_non_fred = _get_masks(df, loop_type)
    data = pd.concat([df[mask_ame].drop_duplicates(["loop_name"]),
                      df[mask_non_ame]])
    data = np.array([[x.dist / ANGLEWEIGHT, x.angle1, x.angle2, x.is_interaction]
                     for x in data.itertuples()])
    return data[:, :3], data[:, 3]
