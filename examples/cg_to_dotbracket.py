#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import sys
from optparse import OptionParser

def main():
    usage = """
    python cg_to_ss_fasta.py file.cg

    Convert a coarse-grain RNA file to a fasta file like this:

        >id
        ACCGGGAAA
        (((...)))
    BEWARE: Structures with pseudoknots will give bogus output.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    cg = ftmc.CoarseGrainRNA(args[0])
    print cg.to_fasta_string()

if __name__ == '__main__':
    main()

