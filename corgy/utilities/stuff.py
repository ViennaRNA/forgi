import corgy.utilities.debug as cud

def merge_intervals(intervals):
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
        if i[0] == curr_interval[1]:
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
