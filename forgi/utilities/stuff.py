from builtins import next
from builtins import range
import itertools as it
import contextlib
import random
import shutil
import tempfile as tf
import collections as col
import sys
import subprocess

import logging
log = logging.getLogger(__name__)

import forgi
import forgi.utilities.debug as cud
from .exceptions import GraphConstructionError


bracket_left = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
bracket_right = ")]}>abcdefghijklmnopqrstuvwxyz"


def is_string_type(stri):
    if sys.version_info < (3,):
        return isinstance(stri, (str, unicode))
    else:
        return isinstance(stri, str)


def get_version_string():
    """
    If installed from github, print the version obtained by `git describe`
    On my local machine, when run within a git directory, get the commit
    hash directly using `git describe`
    """
    try:
        # Installed with setup.py from a gitrepo
        label = "forgi {}".format(forgi.__complete_version__)
    except:
        try:
            # On my local machine, run from git directory.
            repo = subprocess.check_output(
                ["git", "rev-parse", "--show-toplevel"]).decode('ascii')
            if "forgi" not in repo:
                raise OSError("")
            label = subprocess.check_output(
                ["git", "describe", "--dirty"]).decode('ascii')
            label = "forgi {}".format(label)
        except OSError:
            # In production, use the version variable
            label = "forgi {}".format(forgi.__version__)
    return label

# COVERAGE: Not used


def grouped(iterable, n):
    '''
    Return a list of every n elements in iterable.

    http://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list

    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    '''
    return zip(*[iter(iterable)] * n)


def merge_intervals(intervals, diff=0):
    '''
    Take a set of intervals, and combine them whenever the endpoints
    match.

    I.e. [(42,47), (55,60), (60,63), (1,9), (63,71)]

    Should yield

    [(1,9),(42,47), (55,71)]

    There should be no overlapping intervals.

    @param intervals: A set of tuples indicating intervals
    @return: A list of merged intervals
    '''
    intervals.sort()
    iter_intervals = iter(intervals)

    # get the first interval
    curr_interval = list(next(iter_intervals))

    merged_intervals = []

    for i in iter_intervals:
        if abs(i[0] - curr_interval[1]) <= diff:
            # the start of this interval is equal to the end of the
            # current merged interval, so we merge it
            curr_interval[1] = i[1]
        else:
            # start a new interval and add the current merged one
            # to the list of intervals to return
            merged_intervals += [curr_interval]
            curr_interval = list(i)

    merged_intervals += [curr_interval]
    return merged_intervals


def gen_random_sequence(l):
    '''
    Generate a random RNA sequence of length l.
    '''
    return "".join([random.choice(['A', 'C', 'G', 'U']) for i in range(l)])


@contextlib.contextmanager
def make_temp_directory():
    '''
    Yanked from:

    http://stackoverflow.com/questions/13379742/right-way-to-clean-up-a-temporary-folder-in-python-class
    '''
    temp_dir = tf.mkdtemp()
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)


def insert_into_stack(stack, i, j):
    # print "add", i,j
    k = 0
    while len(stack[k]) > 0 and stack[k][len(stack[k]) - 1] < j:
        k += 1
    stack[k].append(j)
    return k


def delete_from_stack(stack, j):
    # print "del", j
    k = 0
    while len(stack[k]) == 0 or stack[k][len(stack[k]) - 1] != j:
        k += 1
    stack[k].pop()
    return k


def pairtable_to_dotbracket(pt):
    """
    Converts arbitrary pair table array (ViennaRNA format) to structure in dot bracket format.
    """
    stack = col.defaultdict(list)
    seen = set()
    res = ""
    for i in range(1, pt[0] + 1):
        if pt[i] != 0 and pt[i] in seen:
            raise ValueError('Invalid pairtable contains duplicate entries')

        seen.add(pt[i])

        if pt[i] == 0:
            res += '.'
        else:
            # '(' check if we can stack it...
            if pt[i] > i:
                res += bracket_left[insert_into_stack(stack, i, pt[i])]
            else:                                    # ')'
                res += bracket_right[delete_from_stack(stack, i)]

    return res


def inverse_brackets(bracket):
    res = col.defaultdict(int)
    for i, a in enumerate(bracket):
        res[a] = i
    return res


def dotbracket_to_pairtable(struct):
    """
    Converts arbitrary structure in dot bracket format to pair table (ViennaRNA format).


    """
    if len(struct) == 0:
        raise ValueError("Cannot convert empty structure to pairtable")
    pt = [0] * ((len(struct) + 1) - struct.count("&"))
    pt[0] = len(struct) - struct.count("&")

    stack = col.defaultdict(list)
    inverse_bracket_left = inverse_brackets(bracket_left)
    inverse_bracket_right = inverse_brackets(bracket_right)

    i = 0
    for a in struct:
        if a == '&':
            continue
        i += 1
        # print i,a, pt
        log.debug("Parsing bracket %r", a)
        if a == ".":
            pt[i] = 0
        else:
            if a in inverse_bracket_left:
                stack[inverse_bracket_left[a]].append(i)
            else:
                assert a in inverse_bracket_right
                if len(stack[inverse_bracket_right[a]]) == 0:
                    raise ValueError('Too many closing brackets!')
                j = stack[inverse_bracket_right[a]].pop()
                pt[i] = j
                pt[j] = i

    if len(stack[inverse_bracket_left[a]]) != 0:
        raise ValueError('Too many opening brackets!')

    return pt


def pairtable_to_tuples(pt):
    '''
    Convert a pairtable to a list of base pair tuples.

    i.e. [4,3,4,1,2] -> [(1,3),(2,4),(3,1),(4,2)]

    :param pt: A pairtable
    :return: A list paired tuples
    '''
    pt = iter(pt)

    # get rid of the first element which contains the length
    # of the sequence. We'll figure it out after the traversal
    next(pt)

    tuples = []
    for i, p in enumerate(pt):
        tuples += [(i + 1, p)]

    return tuples


def tuples_to_pairtable(pair_tuples, seq_length=None):
    '''
    Convert a representation of an RNA consisting of a list of tuples
    to a pair table:

    i.e. [(1,3),(2,4),(3,1),(4,2)] -> [4,3,4,1,2]

    :param tuples: A list of pair tuples
    :param seq_length: How long is the sequence? Only needs to be passed in when
                       the unpaired nucleotides aren't passed in as (x,0) tuples.
    :return: A pair table
    '''
    if seq_length is None:
        max_bp = max([max(x) for x in pair_tuples])
    else:
        max_bp = seq_length

    pt = [0] * (max_bp + 1)
    pt[0] = max_bp

    for tup in pair_tuples:
        pt[tup[0]] = tup[1]

    return pt


def pairtable_to_elements(pt, level, i, j):
    '''
    Convert a pair table to a list of secondary structure
    elements:

     [['s',1,[2,3]]

      The 's' indicates that an element can be a stem. It can also be
      an interior loop ('i'), a hairpin loop ('h') or a multiloop ('m')

      The second number (1 in this case) indicates the depth or
      how many base pairs have to be broken to get to this element.

     Finally, there is the list of nucleotides which are part of
     of this element.
    '''
    elements = []
    u5 = [i - 1]
    u3 = [j + 1]

    if (i > j):
        return []

    # iterate over the unpaired regions on either side
    # this is either 5' and 3' unpaired if level == 0
    # or an interior loop or a multiloop
    while (pt[i] == 0):
        u5.append(i)
        i += 1
    while (pt[j] == 0):
        u3.append(j)
        j -= 1

    if (i > j):
        # hairpin loop or one large unpaired molecule
        u5.append(i)
        if (level == 0):
            return [['e', level, sorted(u5)]]
        else:
            # check to see if we have chain breaks due
            # to multiple strands in the input
            external = False
            left = []
            right = []
            for k in range(0, len(u5)):
                if (external):
                    right.append(u5[k])
                else:
                    left.append(u5[k])

            return [['h', level, sorted(u5)]]

    if (pt[i] != j):
        # multiloop
        m = u5
        k = i

        # the nucleotide before and the starting nucleotide
        m.append(k)
        while (k <= j):
            # recurse into a stem
            elements += pairtable_to_elements(pt, level, k, pt[k])

            # add the nucleotides between stems
            m.append(pt[k])
            k = pt[k] + 1
            while (pt[k] == 0 and k <= j):
                m.append(k)
                k += 1

            m.append(k)

        m.pop()
        m += u3

        if (len(m) > 0):
            if (level == 0):
                elements.append(['e', level, sorted(m)])
            else:
                elements.append(['m', level, sorted(m)])

        return elements

    if (pt[i] == j):
        # interior loop
        u5.append(i)
        u3.append(j)

        combined = u5 + u3
        if len(combined) > 4:
            if (level == 0):
                elements.append(['e', level, sorted(u5 + u3)])
            else:
                elements.append(['i', level, sorted(u5 + u3)])

    s = []
    # go through the stem
    while (pt[i] == j and i < j):
        # one stem
        s.append(i)
        s.append(j)

        i += 1
        j -= 1

        level += 1

    u5 = [i - 1]
    u3 = [j + 1]
    elements.append(['s', level, sorted(s)])

    return elements + pairtable_to_elements(pt, level, i, j)


def bpseq_to_tuples_and_seq(bpseq_str):
    """
    Convert a bpseq string to a list of pair tuples and a sequence
    dictionary. The return value is a tuple of the list of pair tuples
    and a sequence string.

    :param bpseq_str: The bpseq string
    :return: ([(1,5),(2,4),(3,0),(4,2),(5,1)], 'ACCAA')
    """
    lines = bpseq_str.split('\n')
    seq = []
    tuples = []
    pairing_partner = {}
    for line in lines:
        parts = line.split()

        if len(parts) == 0:
            continue

        (t1, s, t2) = (int(parts[0]), parts[1], int(parts[2]))
        if t2 in pairing_partner and t1 != pairing_partner[t2]:
            raise GraphConstructionError("Faulty bpseq string. {} pairs with {}, "
                                         "but {} pairs with {}".format(t2, pairing_partner[t2], t1, t2))
        if t1 in pairing_partner and t2 != pairing_partner[t1]:
            raise GraphConstructionError("Faulty bpseq string. {} pairs with {}, "
                                         "but {} pairs with {}".format(pairing_partner[t1], t1,  t1, t2))

        pairing_partner[t1] = t2
        if t2 != 0:
            pairing_partner[t2] = t1
        tuples += [(t1, t2)]
        seq += [s]

    seq = "".join(seq).upper().replace('T', 'U')

    return (tuples, seq)


def renumber_bpseq(bpseq_triples):
    """
    :param bpseq_triples: A list of triples (from, res, to)
    """
    out = []
    mapping = {}
    for i, triple in enumerate(bpseq_triples):
        mapping[triple[0]] = i + 1
    for triple in bpseq_triples:
        from_, res, to_ = triple
        if to_ not in [0, '0']:
            to_ = mapping[to_]
        from_ = mapping[from_]
        out.append("{} {} {}".format(from_, res, to_))
    return "\n".join(out)
