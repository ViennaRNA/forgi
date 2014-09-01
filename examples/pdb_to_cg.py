#!/usr/bin/python

import sys
import os.path as op

import forgi.threedee.model.coarse_grain as ftmc

from optparse import OptionParser

def main():
    usage = """
    ./pdb_to_cg_fasta.py pdb_file

    Take a pdb file, extract the secondary structure and print it out
    as a fasta file like this:

        >id
        sequence
        secondary structure

    Where the id will be the part of the filename without the extension.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-p', '--pseudoknots', dest='pseudoknots', default=False,
                      help='Keep pseudoknots in the structure', 
                      action='store_true')
    parser.add_option('-d', '--dump-all', dest='dump_all', 
                      default=None, help='Enter a directory where to dump all of \
                                        temporary and intermediate files.',
                      type = 'str')
    parser.add_option('-c', '--chain', dest='chain', 
                      default=None, help='Specify the chain to coarse-grain',
                      type = 'str')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    pdb_id = op.basename(op.splitext(args[0])[0])
    cg = ftmc.from_pdb(args[0], intermediate_file_dir=options.dump_all, 
                       remove_pseudoknots=not options.pseudoknots,
                      chain_id = options.chain)
    print cg.to_cg_string()

if __name__ == '__main__':
    main()

