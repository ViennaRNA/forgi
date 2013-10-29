#!/usr/bin/python

import corgy.graph.bulge_graph as cgb
import sys
from optparse import OptionParser

def main():
    usage = """
    ./list_multiloops.py fasta_file

    List all of the multi-loops in each secondary structure listed in the fasta
    file. The fasta file should look like the example below where an id line is
    followed by the sequence and then a number of secondary structures.

        >id
        ACCGCC
        ((()))
        ((..))
        (....)
        
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    with open(args[0], 'r') as f:
        lines = f.readlines()
        # extract the sequence and id, they're not really necessary
        # but what the heck
        id = lines[0].strip('>')
        seq = lines[1].strip()

        # iterate over each secondary structure
        for line in lines[2:]:
            bg = cgb.BulgeGraph()
            bg.from_dotbracket(line.strip())

            for m in bg.find_multiloop_loops():
                loop = [d for d in m if d[0] == 'm']
                for l in loop:
                    print l, bg.defines[l]
                loop.sort(key=lambda x: bg.defines[x][0])
                sizes = [bg.defines[d][1] - bg.defines[d][0] for d in loop]

                print "".join(loop), " ".join(map(str,sizes))

if __name__ == '__main__':
    main()

