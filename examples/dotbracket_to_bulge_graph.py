import sys
from optparse import OptionParser

import corgy.graph.bulge_graph as cgb

def main():

    usage = """
        Usage: ./create_bulge_graph.py file.dotbracket [name]"
        
        Creates a graph of the paired and unpaired regions within a
        dotplot. Paired regions are called stems while unpaired regions
        are called bulges.

        The file created contains four sections:
        1. Name (optional):
            The name of this graph. This can be used, for example, to specify which
            PDB file the dotplot was inferred from. The name should not contain any
            spaces.

        2. Length (optional):
            The length of the dotplot that was used to create this graph.

        3. Definitions:
            Each line in the definitions sections begins with the keyword 'define'
            It is followed by the name of the region and the numbers of the 
            nucleotides which define that region.

            The numbers of the nucleotides always come in sets of two for bulges and
            sets of four for stems (which can be interpreted as two sets of two). Each
            set of two numbers indicates the starting and ending nucleotide for that
            region.

            All numbers are 1-based.

        4. Connections:
            This section shows how the different regions are connected to each other.
            The term connected means that they share a common nucleotide.

            Each line begins with the keyword 'connect'. It is followed by the name
            of the region for which the connections are being described. The name
            is followed by a number of different names. Each of these regions is connected
            to the first.

            One region may be defined by more than one connect statement. This section
            of the file is simpy an adjacency list.


        The results are printed to standard out.
        """

    parser = OptionParser(usage = usage)
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
    bg.from_dotbracket(brackets)

    if len(args) == 2:
        bg.name = args[1]

    print bg.to_bg_string()

if __name__ == "__main__":
    main()
