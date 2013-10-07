import sys
from optparse import OptionParser

import corgy.graph.bulge_graph as cgb

def main():

    usage = """
        Usage: ./dotbracket_to_element_string.py file.dotbracket"

        Convert a dotbracket notation to an element-type notation
        where the type of each element is shown.

        For example, the string '..((..))..' will be converted
        to 'ffsshhsstt' to indicate the presence of a fiveprime,
        stem, hairpin, and threeprime section.
        """

    parser = OptionParser(usage = usage)
    parser.add_option('-s', '--secondary-structure', dest='sec_struct', default=False, help="Display the dotbracket represenation of the secondary structure along with the element string.", action='store_true')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0])

    brackets = "".join(f.readlines()).replace('\n', '')
    bg = cgb.BulgeGraph()


    if len(brackets.strip()) == 0:
        print >>sys.stderr, "Invalid input"
        sys.exit(1)

    bg.from_dotbracket(brackets)

    if options.sec_struct:
        print bg.dotbracket_str

    print bg.to_element_string()

if __name__ == "__main__":
    main()
