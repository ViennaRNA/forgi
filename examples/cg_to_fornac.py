#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.comparison as ftme
import forgi.utilities.debug as fud
import forgi.utilities.stuff as fus

import itertools as it
import json

import numpy as np
import os.path as op

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
    <link rel='stylesheet' type='text/css' href='/css/d3-rnaplot.css' />
    <script type='text/javascript' src='/js/jquery.js'></script>
    <script type='text/javascript' src='/js/d3.js'></script>
    <script type='text/javascript' src='/js/d3-grid.js'></script>
    <script type='text/javascript' src='/js/d3-rnaplot.js'></script>
    <script type='text/javascript' src='/js/lib/d3-ForceEdgeBundling.js'></script>

    <script type='text/javascript'>
    // set all of the parameters
    var padding = [10,10];   //the padding between the grid rectangles
    var margin = {{top: 4, left: 4, bottom: 4, right: 4}};
    var svgWidth = 1000 - margin.left - margin.right;  //the width of the svg

    var cellWidth = 200; //the width of each cell
    var cellHeight = 200;  //the height of each cell

    var root = {};

    var minHeight = 500;

    // calculate the number of columns and the height of the SVG,
    // which is dependent on the on the number of data points
    var numCols = Math.floor((svgWidth + padding[0]) / (cellWidth + padding[0]));
    var svgHeight = Math.max(Math.ceil(root.length / numCols) * (cellHeight + padding[1]) - padding[1] + margin.bottom,
                        minHeight);

    var chart = rnaPlot()
    .width(cellWidth)
    .height(cellHeight)
    .bundleExternalLinks(true)

     var rectGrid = d3.layout.grid()
              .bands()
              .size([svgWidth, svgHeight])
              .cols(numCols)
              .padding(padding)
              .nodeSize([cellWidth, cellHeight]);
    var rectData = rectGrid(root)

    // create an svg as a child of the #rna_ss div
    // and then a g for each grid cell
    var svg = d3.select('#rna_ss')
    .append('svg')
    .attr('width', svgWidth + margin.left + margin.right)
    .attr('height', svgHeight + margin.top + margin.bottom)
    .append('g')
    .attr('transform', 'translate(' + margin.left + "," + margin.top + ")")

    var gZoom = svg.append('g')
    .call(d3.behavior.zoom().scaleExtent([0.5, 2 * root.length ]).on("zoom", zoom))

    gZoom.append("rect")
    .attr("class", "overlay")
    .attr("width", svgWidth)
    .attr("height", svgHeight);

    svg = svg.append('g')

    function zoom() {{
        svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
    }}


    svg.selectAll('.rna-cell')
    .data(rectData)
    .enter()
    .append('g')
    .attr('transform', function(d) {{ return 'translate(' + d.x + ',' + d.y + ')'; }})
    .classed('rna-cell', true)
    .call(chart);
    </script>
"""

def reorder_structs(pair_bitmaps):
    '''
    Order the structures according to their first PC as evaluated
    on the bitmap of the pairs.

    @param pair_bitmaps: An array of length n, corresponding to each potentially
        adjacent pair. 1 if that pair is adjacent, 0 if not.
    @return: An array indicating the new ordering of the structures.
    '''
    from sklearn.decomposition import PCA
    pca = PCA(n_components=1)
    one_d_bitmaps = pca.fit_transform(pair_bitmaps)[:,0]

    ix = np.argsort(one_d_bitmaps)
    return ix

def get_residue_num_list(cg, d):
    '''
    Get a list of nucleotides identifying a loop within
    this element. If it's a stem, then pick a cycle
    in the middle.

    @param cg: A CoarseGrainRNA structure.
    @param d: The name of the coarse-grain element
    '''

    return list(cg.define_residue_num_iterator(d, adjacent=True))

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

def extract_extra_links(cg, cutoff_dist=25, bp_distance=sys.maxint,
                       correct_links = None):
    '''
    Get a list of extra_links that are within a certain distance of each 
    other.

    @param cg: The coarse grain RNA structure in question.
    @param dist: The distance.
    @return: A tuple of (links, pairs) where the links contains tuples of 
             the #s of the linked nucleotides and pairs is a bitmap 
             indicating whether a particular pair is linked.
    '''
    links = []
    pair_bitmap = []

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
            pair_bitmap += [0]
            continue
        pair_bitmap += [1]

        if dist < cutoff_dist:
            links1 = get_residue_num_list(cg, e1)
            links2 = get_residue_num_list(cg, e2)
            correct = False

            if correct_links is not None:
                for correct_link in correct_links:
                    if ((links1 == correct_link['from'] and links2 == correct_link['to']) or
                        (links1 == correct_link['from'] and links2 == correct_link['to'])):
                        correct = True
                        break
            else:
                correct = True

            links += [{"from": links1, "to": links2, "linkType": 'correct' if correct else 'incorrect'}]

    return (links, pair_bitmap)

def main():
    usage = """
    python cg_to_fornac_html.py file1.cg file2.cg

    Convert coarse grain files to html files using fornac
    to display a 2D version of the structure.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-d', '--distance', dest='distance', default=25, help="Draw links between elements that are within a certain distance from each other", type='float')
    parser.add_option('-b', '--bp-distance', dest='bp_distance', default=16, help="Draw links only between nucleotides which are so many nucleotides apart", type='int')
    parser.add_option('-s', '--sort-by', dest='sort_by', default='mcc', help="What to sort by (options: mcc, pca)", type='string')
    parser.add_option('-n', '--names', dest='names', default=False, action='store_true', help='Add the name of the structure to the display')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    structs = []
    pair_bitmaps = []
    cgs = []
    all_links = []
    mccs = []

    for filename in args:
        cg = ftmc.CoarseGrainRNA(filename)
        cgs += [cg]
        (links, pair_bitmap) = extract_extra_links(cg, options.distance, options.bp_distance,
                                                  correct_links = None if len(all_links) == 0 else all_links[0])

        all_links += [links]

        pair_bitmaps += [pair_bitmap]
        mcc = ftme.mcc_between_cgs(cgs[0], cg)
        rmsd = ftme.cg_rmsd(cgs[0], cg)

        seq_struct = {"sequence": cg.seq,
                      "structure": cg.to_dotbracket_string(),
                      "extraLinks": links}

        fud.pv('options.names')
        if options.names:
            seq_struct['name'] = op.basename(filename) + " ({:.2f},{:.1f})".format(mcc, rmsd)
        else:
            seq_struct['name'] = ''

        structs += [seq_struct]
        mccs += [mcc]

    if options.sort_by == 'pca':
        print >>sys.stderr, "Sorting by pca"
        ix = reorder_structs(pair_bitmaps) 
    else:
        print >>sys.stderr, "Sorting by mcc"
        ix = np.argsort(-np.array(mccs))

    new_array = [0 for i in range(len(ix))]
    for i,x in enumerate(ix):
        new_array[i] = structs[x]

    print output_template.format(json.dumps(new_array))

if __name__ == '__main__':
    main()

