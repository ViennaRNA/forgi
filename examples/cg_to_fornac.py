#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud
import forgi.utilities.stuff as fus

import itertools as it
import json

import sys
from optparse import OptionParser

rna_structure_template = """
        var options = {{'sequence': '{}',
                        'structure': '{}',
                        'extraLinks': {}
        }};
"""

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
    <script type='text/javascript' src='/js/d3-grid.js'></script>

    <script type='text/javascript'>
        var svgWidth = 1000;
        var nodeWidth = 300;
        var nodeHeight = 300;
        var padding = [10,10];
        var margin = {{top: 4, left: 4, bottom: 4, right: 4}};

        var seqs = {};
        
        var numCols = Math.floor((svgWidth + padding[0]) / (nodeWidth + padding[0]));
        var svgHeight = Math.ceil(seqs.length / numCols) * (nodeHeight + padding[1]) - padding[1] + margin.bottom;

        var rectGrid = d3.layout.grid()
        .bands()
        .size(svgWidth, svgHeight)
        .cols(numCols)
        .padding(padding)
        .nodeSize([nodeWidth, nodeHeight])


        var rectData = rectGrid( seqs );

        d3.select('#rna_ss')
        .selectAll('.rna-struct')
        .data(rectData)
        .enter()
        .append('div')
        .style('position', 'absolute')
        .attr('id', function(d,i) {{ return "rm" + i; }})
        .style('left', function(d) {{ return d.x + "px"; }})
        .style('top', function(d) {{ return d.y + "px"; }})
        .style('width', function(d) {{ return nodeWidth + "px"; }})
        .style('height', function(d) {{ return nodeHeight + "px"; }})
        .classed('rna-struct', true)
        .each(function(d, i) {{
            console.log('d:', d, i);
            var container = new FornaContainer("#rm" + i,
                {{'applyForce': true, 'allowPanningAndZooming': true, 
                'initialSize':[nodeWidth, nodeHeight],
                'cssFileLocation':'/css/fornac.css'}});
        
                container.addRNA(d.structure, d);

                }});
        
    </script>
"""

def get_residue_num_list(cg, d):
    '''
    Get a list of nucleotides identifying a loop within
    this element. If it's a stem, then pick a cycle
    in the middle.

    @param cg: A CoarseGrainRNA structure.
    @param d: The name of the coarse-grain element
    '''
    if d[0] == 'm':
        return sorted(cg.shortest_bg_loop(d))
    if (d[0] != 's'):
        return list(cg.define_residue_num_iterator(d, adjacent=True))

    stem_length = cg.stem_length(d)
    start_res = cg.defines[d][0] + stem_length / 2


    pt = fus.dotbracket_to_pairtable(cg.to_dotbracket_string())
    other_res = pt[start_res]

    nucleotide_list = [start_res, start_res+1,
                       other_res-1, other_res]

    return nucleotide_list


def extract_extra_links(cg, cutoff_dist=25, bp_distance=sys.maxint):
    '''
    Get a list of extra_links that are within a certain distance of each 
    other.

    @param cg: The coarse grain RNA structure in question.
    @param dist: The distance.
    '''
    links = []

    for e1, e2 in it.combinations(cg.nd_define_iterator(), r=2):
        if (e1[0] == 's' and cg.stem_length(e1) == 1) or (e2[0] == 's' and cg.stem_length(e2) == 1):
            continue

        if cg.connected(e1, e2):
            continue

        bp_dist = cg.min_max_bp_distance(e1, e2)[0]
        #bp_dist = 1000

        if bp_dist < bp_distance:
            continue

        dist = cg.element_physical_distance(e1,e2)

        if dist > cutoff_dist:
            continue

        if dist < cutoff_dist:
            fud.pv('e1,e2,bp_dist')
            links1 = get_residue_num_list(cg, e1)
            links2 = get_residue_num_list(cg, e2)

            links += [[links1, links2]]

    fud.pv('json.dumps(links)')
    return links

def main():
    usage = """
    python cg_to_fornac_html.py file1.cg file2.cg

    Convert coarse grain files to html files using fornac
    to display a 2D version of the structure.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-d', '--distance', dest='distance', default=10000, help="Draw links between elements that are within a certain distance from each other", type='float')
    parser.add_option('-b', '--bp-distance', dest='bp_distance', default=16, help="Draw links only between nucleotides which are so many nucleotides apart", type='int')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    structs = []

    for filename in args:
        cg = ftmc.CoarseGrainRNA(filename)
        seq_struct = {"sequence": cg.seq,
                      "structure": cg.to_dotbracket_string(),
                      "extraLinks": extract_extra_links(cg, options.distance, options.bp_distance)}

        structs += [seq_struct]


    fud.pv('json.dumps(structs)')
    print output_template.format(json.dumps(structs))


if __name__ == '__main__':
    main()

