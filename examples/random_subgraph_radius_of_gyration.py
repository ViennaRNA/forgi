#!/usr/bin/python

from __future__ import print_function
import sys
import forgi.utilities.debug as cud
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.descriptors as ftmd
from optparse import OptionParser

def main():
    usage = """
    python random_subgraph_radius_of_gyration.py file.cg
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-i', '--iterations', dest='iterations', default=1, help="The number of iterations to perform", type='int')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg = ftmc.CoarseGrainRNA(args[0])
    for i in range(options.iterations):
        sg = cg.random_subgraph()
        stems = [s for s in sg if s[0] == 's']

        if len(stems) == 0:
            continue

        coords = []
        for s in stems:
            coords += [cg.coords[s][0]]
            coords += [cg.coords[s][1]]

        rmsd = ftmd.radius_of_gyration(coords)
        total_length = sum([len(list(cg.define_residue_num_iterator(d))) for d in sg])

        print(total_length, rmsd)

if __name__ == '__main__':
    main()

