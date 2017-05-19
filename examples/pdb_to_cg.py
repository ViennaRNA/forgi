#!python

import sys
import os.path as op

import forgi.threedee.model.coarse_grain as ftmc

from optparse import OptionParser
import logging
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
                      default=None, help='Specify the chain to coarse-grain.',
                      type = 'str')
    parser.add_option('-f', '--to-files', dest='to_files', 
                      default=None, help='Write all connected components to files in the current working directory.',
                      action='store_true' )
    parser.add_option('-v', '--verbose', dest='verbose', default=False,
                      help='Be verbose', 
                      action='store_true')
    parser.add_option( '--debug', dest='debug', default=False,
                      help='Be more verbose', 
                      action='store_true')

    (options, args) = parser.parse_args()

    logging.basicConfig()
    if options.verbose:
        logging.getLogger().setLevel(level=logging.INFO)
    if options.debug:
        logging.getLogger().setLevel(level=logging.DEBUG)

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    if options.to_files:
        if options.chain:
            raise ValueError("Options --to-files and --chain are mutually exclusive!")
        cgs = ftmc.connected_cgs_from_pdb(args[0], remove_pseudoknots=not options.pseudoknots )
        for cg in cgs:
            cg.to_file(cg.name+".cg")
            print("File {}.cg written".format(cg.name))
    else:
        cg = ftmc.from_pdb(args[0], intermediate_file_dir=options.dump_all, 
                       remove_pseudoknots=not options.pseudoknots,
                       chain_id = options.chain)
        print cg.to_cg_string()

if __name__ == '__main__':
    main()

