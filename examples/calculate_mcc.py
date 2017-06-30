#!/usr/bin/python

from __future__ import print_function
import sys
from optparse import OptionParser

import forgi.threedee.model.similarity as ftme
import forgi.threedee.model.coarse_grain as ftmc

def main():
    usage = """
    python calculate_mcc.py struct1.cg struct2.cg
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg1 = ftmc.CoarseGrainRNA(args[0])
    cg2 = ftmc.CoarseGrainRNA(args[1])

    confusion_matrix = ftme.confusion_matrix(cg1, cg2)

    
    print("confusion_matrix:", confusion_matrix)
    print("mcc:", ftme.mcc(confusion_matrix))
    print("rmsd:", ftme.cg_rmsd(cg1, cg2))

if __name__ == '__main__':
    main()

