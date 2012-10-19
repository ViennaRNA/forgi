#!/usr/bin/python

import sys, math
import corgy.graph.bulge_graph as cgb

from optparse import OptionParser

def print_neato(bg):
    stems = bg.get_bulged_stem_names()

    # The different nodes for different types of bulges

    node_defs = dict()
    node_lines = dict()

    # loops
    node_defs[1] = '\t{node [style=filled,shape=circle,fillcolor=red,fontsize=12'
    node_defs[2] = '\t{node [style=filled,shape=circle,fillcolor=yellow,fontsize=12'
    node_defs[3] = '\t{node [style=filled,shape=hexagon,fillcolor=red,fontsize=12'
    node_defs[4] = '\t{node [style=filled,shape=octagon,fillcolor=red,fontsize=12'
    node_defs[5] = '\t{node [style=filled,shape=octagon,fillcolor=red,fontsize=12'
    node_defs[6] = '\t{node [style=filled,shape=octagon,fillcolor=red,fontsize=12'
    node_defs[7] = '\t{node [style=filled,shape=octagon,fillcolor=red,fontsize=12'
    node_defs[8] = '\t{node [style=filled,shape=octagon,fillcolor=red,fontsize=12'
    node_defs[9] = '\t{node [style=filled,shape=octagon,fillcolor=red,fontsize=12'

    node_lines = ''
    connection_lines = ''
    fancy = True

    print "graph G {"
    print "\tgraph [overlap=scale];"
    print "\tnode [shape=box];"

    if fancy:
        for key2 in bg.defines.keys():
            if key2[0] == 's':
                node_lines += '\t{node [style=filled,fillcolor=green,label=\"%s\\n(%d)\"] %s};\n' % (key2, bg.stem_length(key2), key2)
                continue

            if len(bg.edges[key2]) == 2:
                if bg.weights[key2] == 2:
                    node_lines += '\t{node [style=filled,shape=circle,fillcolor=yellow,fontsize=12'

                if bg.weights[key2] == 1:
                    node_lines += '\t{node [style=filled,shape=circle,fillcolor=red,fontsize=12'
            else:
                node_lines += '\t{node [style=filled,shape=circle,fillcolor=blue,fontsize=12'


            #node_lines += node_defs[bg.weights[key2]] 
            node_lines += ',label=\"%s \\n(' % (key2)
            total_bulge = 0

            for j in range(0, len(bg.defines[key2]), 2):
                if j != 0:
                    node_lines += ','

                total_bulge += abs((int(bg.defines[key2][j+1]) - int(bg.defines[key2][j]) + 1))
                node_lines += "%d" % (int(bg.defines[key2][j+1]) - int(bg.defines[key2][j]) + 1)
            j = j / 2

            while j < bg.weights[key2]-1:
                node_lines += ",0"
                j += 1

            width = math.sqrt(1.5 * total_bulge / 10.0) 
            height = width

            if bg.weights[key2] == 2:
                node_lines += ")\",width=%f,heigh=%f] %s};\n" % (width, height, key2)
            else:
                node_lines += ")\"] %s};\n" % (key2)
    else:
        for key2 in bg.defines.keys():
            if key2[0] == 's':
                node_lines += '\t{node [style=filled,fillcolor=green,label=\"%s\"] %s};\n' % (key2, key2)
            else:
                node_lines += node_defs[bg.weights[key2]]
                node_lines += ',label=\"%s\"] %s};\n' % (key2, key2)


    
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
    print_neato(bg)

if __name__=="__main__":
    main()
