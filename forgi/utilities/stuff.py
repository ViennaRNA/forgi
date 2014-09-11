import itertools as it
import contextlib
import random
import shutil
import tempfile as tf

import forgi.utilities.debug as cud

def grouped(iterable, n):
    '''
    Return a list of every n elements in iterable.

    http://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list

    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    '''
    return it.izip(*[iter(iterable)]*n)

def merge_intervals(intervals, diff = 0):
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
    return "".join([random.choice(['A','C','G','U']) for i in range(l)])

@contextlib.contextmanager
def make_temp_directory():
    '''
    Yanked from:

    http://stackoverflow.com/questions/13379742/right-way-to-clean-up-a-temporary-folder-in-python-class
    '''
    temp_dir = tf.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)
