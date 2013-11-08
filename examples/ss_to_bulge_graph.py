import sys
from optparse import OptionParser

import forgi.graph.bulge_graph as cgb
import forgi.utilities.debug as cud

def main():

    usage = """
        Usage: ./ss_to_bulge_graph.py fastx_file

        Convert a ss file to a structure.

        A ss file contains the following lines:

            >id
            ACCGGGCCC
            ...((..))

        The id can be anything at all, but it must be prefixed with a '>'.
        The sequence must be a sequence of A, C, G, and U, while the 
        dotbracket string must depict the secondary structure of the
        molecule using '.' for unpaired base pairs and '()' for pairsj.
        
        The results are printed to standard out.
        """

    parser = OptionParser(usage = usage)
    parser.add_option('-d', '--dissolve-length-one-stems', dest="dissolve", default=False, action="store_true",  help='Remove any stems that have a length of just one nucleotide')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0])

    lines = f.readlines()

    bg = cgb.BulgeGraph()
    bg.from_dotbracket(lines[-1].strip(), options.dissolve)
    bg.name = lines[0].strip().strip('>')
    bg.seq = lines[-2].strip()

    print bg.to_bg_string()

if __name__ == "__main__":
    main()
