#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import sys
from optparse import OptionParser

output_template = """
<!DOCTYPE html>
<meta charset="utf-8">

This is an RNA container.
<div id='rna_ss'> </div>
This after the RNA container.

    <link rel='stylesheet' type='text/css' href='/css/fornac.css' />
    <script type='text/javascript' src='/js/jquery.js'></script>
    <script type='text/javascript' src='/js/d3.js'></script>
    <script type='text/javascript' src='/js/fornac.js'></script>

    <script type='text/javascript'>
        var container = new FornaContainer("#rna_ss",
            {{'applyForce': true, 'allowPanningAndZooming': true, 'initialSize':[200,500],
            'cssFileLocation':'/css/fornac.css'}});

        var options = {{'structure': '{}',
                        'sequence': '{}',
                        'extraLinks': {}
        }};
        container.addRNA(options.structure, options);
    </script>
"""


def main():
    usage = """
    python cg_to_fornac_html.py file1.cg file2.cg

    Convert coarse grain files to html files using fornac
    to display a 2D version of the structure.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-d', '--distance', dest='distance', default=10000, help="Draw links between elements that are within a certain distance from each other", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    for filename in args:
        cg = ftmc.CoarseGrainRNA(filename)

    sequence_string = cg.seq
    structure_string = cg.to_dotbracket_string()
    extra_links_string = ""
    print output_template.format(sequence_string, structure_string, extra_links_string)


if __name__ == '__main__':
    main()

