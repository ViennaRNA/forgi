#!/usr/bin/python

import sys, math
import forgi.graph.bulge_graph as cgb

from optparse import OptionParser

def print_neato(bg):
    # The different nodes for different types of bulges
    node_lines = dict()

    node_lines = ''
    connection_lines = ''
    fontsize=20

    print "graph G {"
    print "\tgraph [overlap=false,splines=true];"
    print "\tnode [shape=box];"

    for key2 in bg.defines.keys():
        # Create the nodes with a different color for each type of element
        if key2[0] == 's':
            node_lines += '\t{node [style=filled,fillcolor="#B3E2CD",fontsize=%d,label=\"%s\\n(%d)\"] %s};\n' % (fontsize,key2, bg.stem_length(key2), key2)
            continue
        elif key2[0] == 'i':
            node_lines += '\t{node [style=filled,shape=circle,fillcolor="#FFF2AE",fontsize=%d' % (fontsize)
        elif key2[0] == 'm':
            node_lines += '\t{node [style=filled,shape=circle,fillcolor="#F4CAE4",fontsize=%d' % (fontsize)
        elif key2[0] == 'f':
            node_lines += '\t{node [style=filled,shape=circle,fillcolor="#FDCDAC",fontsize=%d' % (fontsize)
        elif key2[0] == 't':
            node_lines += '\t{node [style=filled,shape=circle,fillcolor="#E6F5C9",fontsize=%d' % (fontsize)
        else:
            node_lines += '\t{node [style=filled,shape=circle,fillcolor="#CBD5E8",fontsize=%d' % (fontsize)

        node_lines += ',label=\"%s \\n' % (key2)

        # figure out the size of the node and use that as a lbel
        node_dims = bg.get_node_dimensions(key2)
        total_bulge = sum(node_dims)
        node_lines += str(node_dims)

        # make bigger interior loops visually bigger
        width = math.sqrt(1.5 * total_bulge / 10.0) 
        height = width

        if key2[0] == 'i':
            node_lines += "\",width=%f,heigh=%f] %s};\n" % (width, height, key2)
        else:
            node_lines += "\"] %s};\n" % (key2)
    
    for key1 in bg.edges:
        if key1[0] == 's':
            for key2 in bg.edges[key1]:
                connection_lines += "\t%s -- %s;\n" % (key1, key2)

    for key1 in bg.longrange.keys():
        for key2 in bg.longrange[key1]:
            connection_lines += "\t%s -- %s [style=dashed]" % (key1, key2)

    print node_lines
    print connection_lines
    print "}"


    #print bg.get_named_define

def main():
    usage = """
        Usage: ./graph_to_neato.py struct.graph

        Convert a BulgeGraph structure to a representation readable by the 'neato' program.
        """
    parser = OptionParser(usage = usage)
    parser.add_option('-d', '--dotbracket', default=False, action='store_true', help='Indicte that the input is in dotbracket format')
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    if args[0] == '-':
        f = sys.stdin
    else:
        f = open(args[0])

    text = "".join(f.readlines())
    brackets = text.replace('\n', '')
    bg = cgb.BulgeGraph()

    if options.dotbracket:
        bg.from_dotbracket(brackets)
    else:
        bg.from_bg_string(text)

    print_neato(bg)

if __name__=="__main__":
    main()
