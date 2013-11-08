#!/usr/bin/python

import sys
import random
import os.path as op
import forgi.graph.bulge_graph as cgb
import forgi.model.coarse_grain as cmc
from optparse import OptionParser

def main():
    usage = """
    python disturb_multiloop.py

    Add (default) or remove a basepair from a random multiloop in the fast file.
    Print out the resulting structure. An input file must be specified using
    one of the options (fasta or pdb).
    """
    num_args = 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-p', '--pdb', dest='pdb', default='', 
                    help='PDB File to use as input.', type='str')
    parser.add_option('-f', '--fasta', dest='fasta', default='', 
                    help='Fasta file to use as input.', type='str')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    bg = None
    filename = ''

    if len(options.fasta) > 0:
        filename = options.fasta
        with open(options.fastq, 'r'):
            lines = bg.readlines()

            bg = cgb.BulgeGraph()
            bg.from_fasta("".join(lines), False)
    elif len(options.pdb) > 0:
        filename = options.pdb
        cg = cmc.from_pdb(options.pdb)
        bg = cg.bg

    if bg is None:
        parser.print_help()
        sys.exit(1)

    pdb_id = op.basename(op.splitext(filename)[0])
    print ">%s" % (pdb_id)
    print bg.seq
    print bg.to_dotbracket()

    multiloops = bg.find_multiloop_loops()
    for multi in multiloops:
        shortened = set()
        for m in [d for d in multi if d[0] == 'm']:
            connected_stems = list(bg.edges[m])
            for to_shorten in connected_stems:
                if to_shorten in shortened:
                    continue
                shortened.add(to_shorten)

                # find the stems which are connected to this multiloop
                # and pick a random one
                db = list(bg.to_dotbracket())

                # get the side of the stem which is connected to the 
                # multiloop
                (s1b, s1e) = bg.get_sides(to_shorten, m)

                #print to_shorten, s1b, "(", bg.defines[to_shorten], ")"
                # the nucleotides that need to be changed
                to_change = bg.get_side_nucleotides(to_shorten, s1b)
                #print bg.defines[to_shorten], to_change

                db[to_change[0] - 1] = '.'
                db[to_change[1] - 1] = '.'
                print "".join(db)


if __name__ == '__main__':
    main()

