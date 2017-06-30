#!/usr/bin/python

from __future__ import print_function
import collections as c
import itertools as it

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv

import forgi.utilities.debug as fud
import numpy as np
import sys
from optparse import OptionParser

def main():
    usage = """
    python average_atom_positions.py atom_positions1.csv atom_positions2.csv ... etc
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    poss = c.defaultdict(list)
    sources = c.defaultdict(list)

    poss = c.defaultdict(list)

    lines = []

    if args[0] == '-':
        lines = sys.stdin.readlines()
    else:
        for i,arg in enumerate(args):
            with open(arg, 'r') as f:
                lines += f.readlines()

    for line in lines:
        parts = line.strip().split(':')

        identifier = parts[0]
        pos = map(float,parts[1].split(','))
        poss[identifier] += [pos]

    print("import collections as co")

    print("avg_atom_poss = dict()")
    for key in poss.keys():
        pos = np.mean(poss[key], axis=0)
        print('avg_atom_poss["%s"] = [%s] #%d' % (key, ",".join(map(str, pos)), len(poss[key])))

    '''
    print "sources = dict()"
    for key in sources.keys():
        print 'sources["%s"] = [%s]' % (key, ",".join(sources[key]))
    '''

if __name__ == '__main__':
    main()

