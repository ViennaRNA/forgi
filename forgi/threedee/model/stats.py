#!/usr/bin/python

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import zip
from builtins import map
from builtins import object

import csv
import itertools as it
import sys
import warnings

import random as rand
import numpy as np
import numpy.random as nr
import collections as c
import math as m
import math

import logging
log = logging.getLogger(__name__)

import forgi.utilities.debug as fud
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.virtual_residues as ftuvres
import forgi.threedee.utilities.graph_pdb as ftug

# The two constants seem to be unused.
avg_stem_bp_length = 2.24
avg_twist_rotation_per_bp = 360 / 11.


class LoopStat(object):
    '''
    Class for storing the individual statistics about loops.

    phys_length: The length between the start and the centroid of the loop.
    '''

    def __init__(self, line='', s_type="loop"):
        self.stat_type = s_type
        self.pdb_name = ''

        self.bp_length = 0
        self.phys_length = 0.
        self.r = 0.
        self.u = 0.
        self.v = 0.

        self.define = []
        self.seq = ""
        self.vres = {}
        self.vbase = {}
        self.vsugar = {}
        self.vbackbone = {}
        if len(line) > 0:
            try:
                self.parse_line(line)
            except:
                print("Error parsing line:", line, file=sys.stderr)
                raise

    def parse_line(self, line):
        '''
        Parse a line containing statistics about the shape of a stem.

        :param line: The line from the statistics file.
        '''

        parts = line.strip().split()

        self.pdb_name = parts[1]

        self.bp_length = int(parts[2])
        self.phys_length = float(parts[3])
        self.r = float(parts[3])
        self.u = float(parts[4])
        self.v = float(parts[5])
        if len(parts) > 6:
            self.define = list(map(int, [parts[6], parts[7]]))
            self.seq = parts[8]
        if len(parts) > 9:
            self.vres, self.vbase , self.vsugar, self.vbackbone = ftuvres.parse_vres(parts[9:])

    def __str__(self):
        out = ("{stat_type} {pdb_name} {bp_length} {phys_length}"
               " {u} {v} ".format(**self.__dict__))
        out += " ".join(map(str, self.define)) + " " + self.seq
        out += " " + ftuvres.serialize_vres(self.vres)
        out += " vbase " + ftuvres.serialize_vres(self.vbase)
        out += " vsugar " + ftuvres.serialize_vres(self.vsugar)
        out += " vbackbone " + ftuvres.serialize_vres(self.vbackbone)
        return out

    def __eq__(self, other):
        if type(self) == type(other):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other


class StemStat(object):
    '''
    Class for storing the individual statistics about helices.

    Each stem will be defined by its base pair length. Two
    parameters are associated with each base pair length:

    phys_length: The physical length of such a helix
    twist_angle: The angle between its two twist segments
    '''

    def __init__(self, line=''):
        self.pdb_name = ''

        self.bp_length = 0
        self.phys_length = 0.

        self.twist_angle = 0.
        self.define = []
        self.seqs = []
        self.vbase = {}
        self.vsugar = {}
        self.vbackbone = {}
        if len(line) > 0:
            self.parse_line(line)


    def parse_line(self, line):
        '''
        Parse a line containing statistics about the shape of a stem.

        :param line: The line from the statistics file.
        '''

        parts = line.strip().split()

        self.pdb_name = parts[1]

        self.bp_length = int(parts[2])
        self.phys_length = float(parts[3])
        self.twist_angle = float(parts[4])
        if len(parts) > 5:
            self.define = [int(parts[5]), int(parts[6]),
                           int(parts[7]), int(parts[8])]
        try:
            self.seqs = [parts[9], parts[10]]
        except IndexError:
            pass
        else:
            _, self.vbase , self.vsugar, self.vbackbone = ftuvres.parse_vres(parts[11:])

    def __str__(self):
        try:
            return ("stem {pdb_name} {bp_length} {phys_length} {twist_angle} ".format(**self.__dict__)
                    + " ".join(map(str, self.define)) + " " + " ".join(self.seqs)
                    + " vbase {} vsugar {} vbackbone {} ".format(
                        ftuvres.serialize_vres(self.vbase),
                        ftuvres.serialize_vres(self.vsugar),
                        ftuvres.serialize_vres(self.vbackbone)))
        except (KeyError, IndexError):
            warnings.warn(
                "Could not print define '{}' for StemStat".format(self.define))
            return "stem {} {} {} {}".format(self.pdb_name, self.bp_length, self.phys_length, self.twist_angle)


    def __eq__(self, other):
        if type(self) == type(other):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        return not self == other


class AngleStat(object):
    '''
    Class for storing an individual statistic about inter-helical angles.
    '''

    def __init__(self, stat_type="angle", pdb_name='', dim1=0, dim2=0, u=0, v=0, t=0, r1=0, u1=0, v1=0, ang_type='x',
                 define=[], seq="", vres={}, vbase={}, vsugar={}, vbackbone={}):
        #log.debug("Stat init called")
        self.pdb_name = pdb_name
        self.dim1 = dim1
        self.dim2 = dim2
        self.stat_type = stat_type

        self.u = u
        self.v = v
        self.t = t

        self.r1 = r1
        self.u1 = u1
        self.v1 = v1

        self.ang_type = ang_type

        self.define = define
        self.seq = seq
        self.vres = vres
        self.vbase = vbase
        self.vsugar = vsugar
        self.vbackbone = vbackbone

    def __eq__(self, a_s):
        '''
        Is this AngleStat equal to another one.
        '''
        log.debug("Comparing angle stats {} and {}".format(
            self.pdb_name, a_s.pdb_name))
        if self.dim1 != a_s.dim1:
            log.debug("Dim1 {} != {}".format(self.dim1, a_s.dim1))
            return False
        if self.dim2 != a_s.dim2:
            log.debug("Dim2 {} != {}".format(self.dim2, a_s.dim2))
            return False
        if not np.allclose(self.u, a_s.u):
            log.debug("u {} != {}".format(self.u, a_s.u))
            return False
        if not np.allclose(self.v, a_s.v):
            log.debug("v {} != {}".format(self.v, a_s.v))
            return False
        if not np.allclose(self.t, a_s.t):
            log.debug("t {} != {}".format(self.t, a_s.t))
            return False
        if not np.allclose(self.r1, a_s.r1):
            log.debug("r1 {} != {}".format(self.r1, a_s.r1))
            return False
        if not np.allclose(self.u1, a_s.u1):
            log.debug("u1 {} != {}".format(self.u1, a_s.u1))
            return False
        if not np.allclose(self.v1, a_s.v1):
            log.debug("v1 {} != {}".format(self.v1, a_s.v1))
            return False
        log.debug("Angle stats comparing equal")
        return True

    def parse_line(self, line):
        parts = line.strip().split()

        self.stat_type = parts[0]
        self.pdb_name = parts[1]

        self.dim1 = int(parts[2])
        self.dim2 = int(parts[3])

        self.u = float(parts[4])
        self.v = float(parts[5])
        self.t = float(parts[6])

        self.r1 = float(parts[7])
        self.u1 = float(parts[8])
        self.v1 = float(parts[9])

        self.ang_type = int(parts[10])

        if self.dim2 == 1000:
            # multiloop or single stranded region
            if self.dim1 == 0:
                # no unpaired bases means no define
                def_len = 0
            else:
                def_len = 2
        else:
            def_len = 4
            if self.dim1 == 0:
                def_len -= 2
            if self.dim2 == 0:
                def_len -= 2

            # for interior loops, at least one strand has to
            # have an unpaired base
            assert(def_len > 0)

        self.define = list(map(int, parts[11:11 + def_len]))
        self.seq = parts[11 + def_len]
        log.debug("len(parts)=%s, def_len=%s", len(parts), def_len)
        if len(parts) > 12 + def_len:
            self.vres, self.vbase , self.vsugar, self.vbackbone = ftuvres.parse_vres(parts[12 + def_len:])

    def orientation_params(self):
        '''
        Return a tuple containing the parameters which specify the orientation
        of one stem with respect to another.

        :return: (u, v)
        '''

        return (self.u, self.v)

    def twist_params(self):
        '''
        Returns a tuple containing the parameters which specify the difference in
        twist between two stems.

        :return (u, v, t)
        '''
        return (self.u, self.v, self.t)

    def position_params(self):
        '''
        Return a tuple containing the parameters which specify the position
        of one stem with respect to another.

        :return: `(r1, u1, v1)`
        '''
        return (self.r1, self.u1, self.v1)

    def __str__(self):
        out_str = "%s %s %d %d %f %f %f %f %f %f %d %s %s %s vbase %s vsugar %s vbackbone %s" % (self.stat_type,
                                                                 self.pdb_name,
                                                                 self.dim1,
                                                                 self.dim2,
                                                                 self.u,
                                                                 self.v,
                                                                 self.t,
                                                                 self.r1,
                                                                 self.u1,
                                                                 self.v1,
                                                                 self.ang_type,
                                                                 " ".join(
                                                                     map(str, self.define)),
                                                                 self.seq,
                                                                 ftuvres.serialize_vres(self.vres),
                                                                 ftuvres.serialize_vres(self.vbase),
                                                                 ftuvres.serialize_vres(self.vsugar),
                                                                 ftuvres.serialize_vres(self.vbackbone))
        return out_str

    def __hash__(self):
        return hash(str(self))

    def get_angle(self):
        '''
        Return the angle between the two connected stems.
        '''
        return ftuv.vec_angle(np.array([-1., 0., 0.]), ftuv.spherical_polar_to_cartesian([1, self.u, self.v]))

    def deviation_from(self, stat2):
        """
        How much does the other stat differ from this stat?

        :param stat2: Another AngleStat
        :returns: A 4-tuple: The positional deviation in Angstrom, and 3 absolute angular deviations in radians.
                  The  angular deviations are u, v and t
        """
        ret = []
        pos1 = ftuv.spherical_polar_to_cartesian(self.position_params())
        pos2 = ftuv.spherical_polar_to_cartesian(stat2.position_params())
        ret.append(ftuv.magnitude(pos1 - pos2))
        log.debug("Position difference is %f", ret[-1])

        for attr in ["u", "v", "t"]:
            raw_diff = getattr(self, attr) - getattr(stat2, attr)
            # Map the difference to a value between 0 and pi
            raw_diff_on_circle = abs(
                (raw_diff + math.pi / 2) % (math.pi) - math.pi / 2)
            log.debug("Angular difference for %s is %f, mapped to %f",
                      attr, raw_diff, raw_diff_on_circle)
            ret.append(raw_diff_on_circle)
        return tuple(ret)

    def is_similar_to(self, stat2, position_cutoff=4, angular_cutoff=None):
        """
        Returns True, if theAngleStat is similar
        (according to the specified cutoff) to the other Angle Stat.

        :param position_cutoff: in angstrom
        :param angular_cutoff: in radians. If not given, uses the position cutoff as a value in degrees
        """
        if angular_cutoff is None:
            angular_cutoff = math.radians(position_cutoff)
        deviation = self.deviation_from(stat2)
        if deviation[0] > position_cutoff:
            log.debug(
                "Dissimilar, because of position deviation = %f", deviation[0])
            return False
        for dev in deviation[1:]:
            if dev > angular_cutoff:
                log.debug("Dissimilar, because of angular deviation %s (%s) > %f",
                          dev, deviation[1:], angular_cutoff)
                return False
        log.debug("%s < (%f, %f)", deviation, position_cutoff, angular_cutoff)
        return True


class RandomAngleStats(object):
    '''
    Store all of the angle stats.
    '''

    def __init__(self, discrete_angle_stats):
        self.cont_stats = dict()
        self.make_random(discrete_angle_stats)

    def create_random_function(self, data):
        '''
        Create a function that returns a random value for
        each column in the statistics.

        :param stats: A table containing n rows and m columns
        :return: A function returning m random values with a
                 maximum and minimum no greater than the largest
                 and least values in that column, respectively.
        '''
        mins = map(min, data.T)
        maxs = map(max, data.T)

        def bounded_uniform():
            return [nr.uniform(i, x) for i, x in zip(mins, maxs)]

        return bounded_uniform

    def make_random(self, discrete_angle_stats):
        '''
        Create a set of statistsics that is random in each direction.

        The maximum and minimum u and v values will be taken
        from the discrete statistics.
        '''
        import scipy.stats as ss
        for key1, key2, key3 in discrete_angle_stats.keys():
            dims = (key1, key2, key3)
            data = []

            for d in discrete_angle_stats[(key1, key2, key3)]:
                data += [[d.u, d.v, d.t, d.r1, d.u1, d.v1]]

            '''
            if len(data) < 3:
                continue
            '''

            try:
                self.cont_stats[dims] = self.create_random_function(
                    np.array(data))
            except np.linalg.LinAlgError as lae:
                print("Singular matrix, dimensions:", dims, file=sys.stderr)

    def sample_stats(self, dims):
        '''
        Sample a set of statistics.

        :param dims: The dimensions of the bulge for which to sample.
        '''
        new_stats = self.cont_stats[dims]()
        s = AngleStat()
        (s.u, s.v, s.v, s.r1, s.u1, s.v1) = new_stats
        return s


class ConstructionStats(object):
    angle_stats = None
    stem_stats = None
    loop_stats = None
    fiveprime_stats = None
    threeprime_stats = None
    conf_stats = None


def get_angle_stats(filename, refresh=False):
    '''
    Load the statistics about inter the helix-helix orientations from a file.

    The file format should be as follows:

    `angle pdb_name dim1 dim2 r u v t r1 u1 v1 s1b s2b`

    Where the parameters are as follows:

    * `angle`: identifier for a type of statistics... should always just be 'angle'
    * `pdb_name`: the name of the pdb file these statistics came from
    * `dim1`: the smaller dimension of the bulge
    * `dim2`: the larger dimension of the bulge
    * `u`: the polar angle of the orientation of the 2nd stem
    * `v`: the azimuth of the orientation of the 2nd stem
    * `t`: the orientation of the twist of the second stem
    * `r1`: the distance of the start of the 2nd stem helix from the end of the 1st
            stem helix
    * `u1`: the polar angle of the separation vector of the two helices
    * `v1`: the azimuth of the separation vector of the two helices
    * `s1b`: The side of the first stem closest to the bulge
    * `s2b`: The side of the second stem closest to the bulge

    The azimuth is always defined with respect to the coordinate system defined
    by the stem1 helix axis vector and it's twist vector (the one adjacent to the
    bulge element).
    '''
    if ConstructionStats.angle_stats != None and not refresh:
        return ConstructionStats.angle_stats

    '''
    import pickle
    ConstructionStats.angle_stats = pickle.load(open('fess/stats/angle_stats.pickle', 'r'))

    print >>sys.stderr, "done loading stats"
    return ConstructionStats.angle_stats

    '''
    ConstructionStats.angle_stats = c.defaultdict(list)
    # DefaultDict(DefaultDict([]))

    f = open(filename, 'r')

    count = 0
    for line in f:
        line = line.strip()
        if line.startswith('angle'):
            angle_stat = AngleStat()
            angle_stat.parse_line(line)

            if len(angle_stat.define) > 0 and angle_stat.define[0] == 1:
                continue

            ConstructionStats.angle_stats[(
                angle_stat.dim1, angle_stat.dim2, angle_stat.ang_type)] += [angle_stat]
            ConstructionStats.angle_stats[(
                angle_stat.dim2, angle_stat.dim1, -angle_stat.ang_type)] += [angle_stat]
            count += 1

    f.close()

    '''
    pickle.dump(ConstructionStats.angle_stats, open('fess/stats/angle_stats.pickle','w'))
    '''

    #print >>sys.stderr, "done loading stats"
    return ConstructionStats.angle_stats


def get_angle_stat_dims(s1, s2, angle_type, min_entries=1):
    '''
    Return a list of tuples which indicate the dimensions for which angle
    stats are avilable.

    :param s1: The first size
    :param s2: The second size
    :param angle_type: The type of the angle.
    :param min_entries: The minimum number of stats that have to be available
    '''
    available_stats = []
    angle_stats = get_angle_stats()

    for (k1, k2, k3) in angle_stats.keys():
        # and len(angle_stats[(k1,k2,k3)]) >= min_entries: #BT: I think this should't be here.
        if k3 == angle_type:
            dist = m.sqrt((k1 - s1) ** 2 + (k2 - s2) ** 2)
            available_stats += [(dist, (k1, k2, k3))]

    available_stats.sort()
    return available_stats


def get_one_d_stat_dims(d, stats, min_entries=1):
    available_stats = []

    for k in stats.keys():
        available_stats += [(abs(d - k), k)]

    available_stats.sort()
    return available_stats


def get_stem_stats(filename, refresh=False):
    '''
    Load the statistics from the file.

    format:

    stem pdb_name bp_length phys_length twist_angle

    :param filename: The name of the file.
    '''
    if ConstructionStats.stem_stats is not None and not refresh:
        return ConstructionStats.stem_stats

    ConstructionStats.stem_stats = c.defaultdict(list)

    f = open(filename, 'r')
    #print("Opening file", filename)
    for line in f:
        if line.strip().find('stem') == 0:
            stem_stat = StemStat(line)
            '''
            if stem_stat.define[0] == 1:
                continue
            '''

            ConstructionStats.stem_stats[(
                stem_stat.bp_length, stem_stat.bp_length)] += [stem_stat]

    f.close()

    return ConstructionStats.stem_stats


def get_fiveprime_stats(filename, refresh=False):
    '''
    Load the statistics from the file.

    format:

    `fiveprime pdb_name bp_length phys_length`

    :param filename: The name of the file.
    '''
    if ConstructionStats.fiveprime_stats is not None and not refresh:
        return ConstructionStats.fiveprime_stats

    ConstructionStats.fiveprime_stats = c.defaultdict(list)

    f = open(filename, 'r')

    for line in f:
        if line.strip().find('5prime') == 0:
            fiveprime_stat = LoopStat(line)
            ConstructionStats.fiveprime_stats[fiveprime_stat.bp_length] += [
                fiveprime_stat]

    f.close()

    return ConstructionStats.fiveprime_stats


def get_threeprime_stats(filename, refresh=False):
    '''
    Load the statistics from the file.

    format:

    `threeprime pdb_name bp_length phys_length`

    :param filename: The name of the file.
    '''
    if ConstructionStats.threeprime_stats is not None and not refresh:
        return ConstructionStats.threeprime_stats

    ConstructionStats.threeprime_stats = c.defaultdict(list)

    f = open(filename, 'r')

    for line in f:
        if line.strip().find('3prime') == 0:
            threeprime_stat = LoopStat(line)
            ConstructionStats.threeprime_stats[threeprime_stat.bp_length] += [
                threeprime_stat]

    f.close()

    return ConstructionStats.threeprime_stats


def get_loop_stats(filename, refresh=False):
    '''
    Load the statistics from the file.

    format:

    `loop pdb_name bp_length phys_length`

    :param filename: The name of the file.
    '''
    if ConstructionStats.loop_stats != None and not refresh:
        return ConstructionStats.loop_stats

    ConstructionStats.loop_stats = c.defaultdict(list)

    f = open(filename, 'r')

    for line in f:
        if line.strip().find('loop') == 0:
            loop_stat = LoopStat(line)
            ConstructionStats.loop_stats[loop_stat.bp_length] += [loop_stat]

    f.close()

    return ConstructionStats.loop_stats


class ClusteredAngleStats(object):
    def __init__(self, filename):
        """
        A collection of clustered angle_stats. Whenever stats for a certain key are retrieved,
        only one (randomly choosen) representative per cluster is returned.

        Requires a clustered angle stats file (created by fess/scripts/cluster_stats.py).
        The file should look like this:

        '
        # Cluster 0 for (1, 1, -1):
        angle RS_1788_S_000008_A 1 1 1.520877 1.837257 -1.380352 11.237569 0.867246 2.057771 1 19 19 30 30 UCG CAA
        # Cluster 1 for (1, 1, -1):
        angle RS_1140_S_000002_A 1 1 1.563667 1.869580 -0.967966 16.241881 0.896602 2.177000 -1 13 13 24 24 ACG CUU
        angle RS_1183_S_000008_A 1 1 1.533551 1.914482 -1.108396 15.780723 0.829850 2.255556 -1 4 4 46 46 UAG CAG
        '

        :param filename: The filename of the clustered angle stats file.
        """
        #: A dict `(dim1, dim2, ang_type)` : list of lists.
        #: Each value is a list of clusters, and each cluster is a list.
        self._stats_dict = c.defaultdict(list)
        lastkey = None
        with open(filename) as f:
            for line in f:
                if line.startswith('# Cluster'):
                    fields = line.split()
                    dim1 = int(fields[4].lstrip("(").rstrip(","))
                    dim2 = int(fields[5].rstrip(","))
                    ang_type = int(fields[6].rstrip("):"))
                    self._stats_dict[(dim1, dim2, ang_type)].append([])
                    lastkey = (dim1, dim2, ang_type)
                elif line.startswith('#'):
                    continue
                else:
                    angle_stat = AngleStat()
                    angle_stat.parse_line(line)
                    # I don't know what this does, I just copied the following 2 lines from
                    # Peter's code in get_angle_stats
                    if len(angle_stat.define) > 0 and angle_stat.define[0] == 1:
                        continue

                    assert lastkey == (angle_stat.dim1, angle_stat.dim2, angle_stat.ang_type) or lastkey == (
                        angle_stat.dim2, angle_stat.dim1, -angle_stat.ang_type)
                    self._stats_dict[lastkey][-1].append(angle_stat)

    def __getitem__(self, key):
        """
        Returns a list of stats with only one randomly choosen stat for each cluster.
        """
        return [rand.choice(x) for x in self._stats_dict[key]]

    def keys(self):
        return self._stats_dict.keys()

    def lookup_stat(self, stat):
        key = (stat.dim1, stat.dim2, stat.ang_type)
        clusters = self._stats_dict[key]
        total_length = sum(len(cluster) for cluster in clusters)
        cluster_length = -1
        num_clusters = len(clusters)
        for cluster in clusters:
            if stat in cluster:
                cluster_length = len(cluster)
                break
        return cluster_length, total_length, num_clusters

    def cluster_of(self, stat):
        key = (stat.dim1, stat.dim2, stat.ang_type)
        clusters = self._stats_dict[key]
        for i, cluster in enumerate(clusters):
            if stat in cluster:
                return i
        return -1

    def get_angle_stat_dims(self, dim0, dim1, ang_type):
        """
        Returns a list of pairs,`(dist, key)` ordered by increasing distance of the key to the query key.

        :param dim1: dim1 of the query key
        :param dim2: dim2 of the query key
        :param ang_type: angle_type of the query key.
        """
        available_stats = []
        for (k1, k2, k3) in self.keys():
            if k3 == ang_type:
                dist = m.sqrt((k1 - dim0) ** 2 + (k2 - dim1) ** 2)
                available_stats += [(dist, (k1, k2, k3))]

        available_stats.sort()
        return available_stats


class ConformationStats(object):
    def __init__(self, stats_file, clustered_angle_stats_file=None):
        if clustered_angle_stats_file is None:
            self.angle_stats = get_angle_stats(stats_file, refresh=True)
        else:
            self.angle_stats = ClusteredAngleStats(clustered_angle_stats_file)
        self.stem_stats = get_stem_stats(stats_file, refresh=True)
        self.fiveprime_stats = get_fiveprime_stats(stats_file, refresh=True)
        self.threeprime_stats = get_threeprime_stats(stats_file, refresh=True)
        self.loop_stats = get_loop_stats(stats_file, refresh=True)

        self.constrained_stats = c.defaultdict(list)

    def constrain_stats(self, constraint_file):
        '''
        Constrain the statistics for certain regions of the molecule. This
        is created for use with the JAR3D annotations for loop regions.

        :param constraint_file: A file containing the allowed statistics for
                                a particular loop.
        :return: Nothing
        '''
        warnings.warn("Function does nothing!!!")
        raise NotImplementedError("TODO")

    def sample_stats(self, bg, elem, min_entries=10):
        '''
        Return a set of statistics compatible with this element.

        :param bg: The graph representation we're using.
        :param elem: The name of the element
        :return: A list of compatible statistics
        '''
        if elem in self.constrained_stats:
            return self.constrained_stats[elem]

        dims = bg.get_node_dimensions(elem)

        if elem[0] == 's':
            dims = [dims]
            stats = self.stem_stats
            if stats[dims[0]]:
                return stats[dims[0]]
            else:
                raise LookupError("No stats for element {} with dimensions {}. Stats keys are {}".format(
                    elem, dims[0], sorted(stats.keys())))
        elif elem[0] == 'i' or elem[0] == 'm':
            stats = self.angle_stats

            ang_type = bg.get_angle_type(elem)
            try:
                if isinstance(stats, ClusteredAngleStats):
                    dims = stats.get_angle_stat_dims(
                        dims[0], dims[1], ang_type)
                else:  # assert isinstance(stats, defaultdict)
                    dims = get_angle_stat_dims(dims[0], dims[1],
                                               ang_type, min_entries=min_entries)
            except IndexError:
                print("Error in sample_stats:", file=sys.stderr)
                print("elem:", elem, "dims:", dims,
                      "ang_type:", ang_type, file=sys.stderr)
                raise
        elif elem[0] == 'h':
            dims = dims[0]
            stats = self.loop_stats
            dims = get_one_d_stat_dims(dims, stats)
        elif elem[0] == 't':
            dims = dims[0]
            stats = self.threeprime_stats
            dims = get_one_d_stat_dims(dims, stats)

        elif elem[0] == 'f':
            dims = dims[0]
            stats = self.fiveprime_stats
            dims = get_one_d_stat_dims(dims, stats)

        all_stats = []
        for dim in dims:
            if len(all_stats) > min_entries:
                continue

            all_stats += stats[dim[-1]]

        if len(all_stats) == 0:
            msg = "No statistics for bulge {} with dims {}".format(elem, dims)
            if elem[0] in "mi":
                msg += " and ang_type {}".format(ang_type)
            raise LookupError(msg)

        return all_stats


class FilteredConformationStats(ConformationStats):
    def __init__(self, stats_file, filter_filename=None, filter_prob=1):
        """
        :param filter_prob: Return the filtered stats with this probability, else all stats.
                            Default=1 (100%)
        """
        super(FilteredConformationStats, self).__init__(stats_file)

        self.filtered = None
        self.filtered_stats = None
        self.filter_prob = filter_prob
        if filter_filename is not None:
            self.from_file(filter_filename)

    def from_file(self, filename):
        '''
        Read the statistics in from a file, with the following formatting::

            sampled i3 1X8W_A 4 27 31 101 104 ""
            sampled i2 3U5F_6 4 1348 1348 1365 1367 "cWW AG or UU"
            sampled i2 2QBG_B 4 1013 1014 1103 1103 "cWW AG or UU"

        '''
        self.filtered = dict()
        self.filtered_stats = c.defaultdict(list)

        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter=' ', quotechar='"')
            for row in reader:
                if not row:
                    continue  # Empty line
                elem_name = row[0]
                define_len = int(row[2])
                pdb_id = row[1]
                dims = tuple(map(int, row[3:5]))
                define = map(int, row[5:5 + define_len])

                # get filtered stats for each type of angle
                ang_types = [1, -1]
                for at in ang_types:
                    # for interior loops, the mirror loop should induce the same
                    # bend as if the direction of the stem was reversed
                    # TODO: make sure this is true
                    for stat in it.chain(self.angle_stats[(dims[0], dims[1], at)],
                                         self.angle_stats[(dims[1], dims[0], -at)]):
                        if stat.pdb_name == pdb_id:
                            # print stat.define, define
                            pass

                        if stat.pdb_name == pdb_id and stat.define == define:
                            #print >>sys.stderr, "found filtered stats:", stat.pdb_name, define
                            self.filtered_stats[(elem_name, at)] += [stat]

    def sample_stats(self, bg, elem):
        r = rand.random()
        if self.filtered_stats is not None and r <= self.filter_prob:
            ang_type = bg.get_angle_type(elem)
            if (elem, ang_type) in self.filtered_stats:
                if len(self.filtered_stats[(elem, ang_type)]) > 0:
                    return self.filtered_stats[(elem, ang_type)]
        # No filtered stats found. Return normal stats for this element.
        return super(FilteredConformationStats, self).sample_stats(bg, elem)


def get_conformation_stats(stats_file, angle_stats_file=None):
    if ConstructionStats.conf_stats is not None:
        return ConstructionStats.conf_stats
    else:
        return ConformationStats(stats_file, angle_stats_file)


def set_conformation_stats(conf_stats):
    ConstructionStats.conf_stats = conf_stats
