#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import sys
from optparse import OptionParser

template = """
<!DOCTYPE html>
<meta charset="utf-8">

This is an RNA container.
<div id='rna_ss'> </div>
This after the RNA container.

    <link rel='stylesheet' type='text/css' href='/css/fornac.css' />
    <script type='text/javascript' src='/js/jquery.js'></script>
    <script type='text/javascript' src='/js/d3.js'></script>
    <script type='text/javascript' src='/js/fornac.js'></script>
    <script type='text/javascript' src='/src/rnagraph.js'></script>
    <script type='text/javascript' src='/src/fornaf.js'></script>

    <script type='text/javascript'>
        var container = new FornaContainer("#rna_ss",
            {{'applyForce': true, 'allowPanningAndZooming': true, "initialSize": [500,800]}});

        var options = {{'structure': '{}',
                        'sequence': '{}'
        }};

        colorStr = "{}" 

        container.addRNA(options.structure, options);
        cs = ColorScheme(colorStr);

        container.addCustomColors(cs.colorsJson);
        container.changeColorScheme('custom'); 
    </script>
"""


def main():
    usage = """
    python ordered_elements_to_forna_colors.py struct.cg element_list

    Output a forna coloring string such that nucleotides in the elements
    are colored according to the order that they are input.
    """
    num_args= 2
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)
    
    cg = ftmc.CoarseGrainRNA(args[0])

    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('Blues')


    out_str = " "
    to_color_nodes = args[1].split(',')
    for i,node in enumerate(to_color_nodes):
        chosen_color = cmap(1 - i / float(len(to_color_nodes)))

        for res in cg.define_residue_num_iterator(node):
            out_str += " {}:rgb({},{},{})".format(res, int(255 * chosen_color[0]),
                                                       int(255 * chosen_color[1]),
                                                       int(255 * chosen_color[2]))

    print template.format(cg.to_dotbracket_string(),
                          cg.seq,
                          out_str)



if __name__ == '__main__':
    main()

