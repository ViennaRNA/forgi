import itertools as it
import contextlib
import random
import shutil
import tempfile as tf
import collections as col

import forgi.utilities.debug as cud

bracket_left =  "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
bracket_right = ")]}>abcdefghijklmnopqrstuvwxyz"

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

def insert_into_stack(stack, i, j):
	#print "add", i,j
	k = 0
	while len(stack[k])>0 and stack[k][len(stack[k])-1] < j:
		k+=1
	stack[k].append(j)
	return k

def delete_from_stack(stack, j):
	#print "del", j 
	k = 0
	while len(stack[k])==0 or stack[k][len(stack[k])-1] != j:
		k+=1
	stack[k].pop()
	return k	

def pairtable_to_dotbracket(pt):
	"""
	Converts arbitrary pair table array (ViennaRNa format) to structure in dot bracket format.
	"""
	stack = col.defaultdict(list)
	res = ""
	for i in range(1, pt[0]+1):
		if pt[i]==0: 
			res += '.'
		else:
			if pt[i]>i:						# '(' check if we can stack it...
				res += bracket_left[insert_into_stack(stack, i, pt[i])]
			else:									# ')'
				res += bracket_right[delete_from_stack(stack, i)]

	return res

def inverse_brackets(bracket):
	res = col.defaultdict(int)
	for i,a in enumerate(bracket):
		res[a] = i
	return res

def dotbracket_to_pairtable(struct):
	"""
	Converts arbitrary structure in dot bracket format to pair table (ViennaRNA format).
	"""	
	pt = [0] * (len(struct)+1)
	pt[0] = len(struct)

	stack = col.defaultdict(list)
	inverse_bracket_left = inverse_brackets(bracket_left)
	inverse_bracket_right = inverse_brackets(bracket_right)

	for i,a in enumerate(struct):
		i += 1
		#print i,a, pt
		if a == ".": pt[i] = 0
		else: 
			if a in inverse_bracket_left: stack[inverse_bracket_left[a]].append(i)
			else: 
				j = stack[inverse_bracket_right[a]].pop()
				pt[i] = j
				pt[j] = i

	return pt

