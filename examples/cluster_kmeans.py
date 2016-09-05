#!/usr/bin/python

import copy
import forgi.threedee.model.coarse_grain as ftmc
import operator as oper
import os
import os.path as op
import scipy.cluster.vq as scv
import sys
import numpy as np
from optparse import OptionParser

def main():
    usage = """
    python cluster_kmeans.py struct1.cg struct2.cg ...

    This script assumes that all of the cgs that are passed in
    correspond to the same RNA.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    parser.add_option('-n', '--number-of-clusters', dest='number_of_clusters', default=3, 
                      help="The number of clusters to partition the structures into", type='int')
    parser.add_option('-o', '--output-dir', dest='output_dir', default='/tmp',
                      help='The output directory into which to place the clustered structures', 
                      type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cgs = map(ftmc.CoarseGrainRNA, args)
    points = []
    for cg in cgs:
        points.append(cg.get_coordinates_array().flatten())
    new_cgs = map(copy.deepcopy, [cgs[i] for i in range(options.number_of_clusters)])

    print >>sys.stderr, "len(points)", len(points)

    clusters = scv.kmeans(points, options.number_of_clusters)
    for i,cluster in enumerate(clusters[0]):
        new_cgs[i].load_coordinates_array(cluster.reshape(-1,3))

        if not op.exists(options.output_dir):
            os.makedirs(options.output_dir)
        new_cgs[i].to_file(op.join(options.output_dir, 'cluster_{}.coord'.format(i)))

if __name__ == '__main__':
    main()

