#!/usr/bin/python

import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud
import forgi.utilities.stuff as fus

import itertools as it
import json

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

        var options = {{'sequence': '{}',
                        'structure': '{}',
                        'extraLinks': {}
        }};
        container.addRNA(options.structure, options);
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
        all_residues = set()
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


def extract_extra_links(cg, cutoff_dist=25):
    '''
    Get a list of extra_links that are within a certain distance of each 
    other.

    @param cg: The coarse grain RNA structure in question.
    @param dist: The distance.
    '''
    links = []

    for e1, e2 in it.combinations(cg.defines, r=2):
        if (e1[0] == 's' and cg.stem_length(e1) == 1) or (e2[0] == 's' and cg.stem_length(e2) == 1):
            continue

        if cg.connected(e1, e2):
            continue

        dist = cg.element_physical_distance(e1,e2)

        if dist > cutoff_dist:
            continue

        fud.pv('dist, cutoff_dist')

        if dist < cutoff_dist:
            links1 = get_residue_num_list(cg, e1)
            links2 = get_residue_num_list(cg, e2)

            links += [[links1, links2]]

    fud.pv('json.dumps(links)')
    return json.dumps(links)

def main():
    usage = """
    python cg_to_fornac_html.py file1.cg file2.cg

    Convert coarse grain files to html files using fornac
    to display a 2D version of the structure.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-d', '--distance', dest='distance', default=10000, help="Draw links between elements that are within a certain distance from each other", type='float')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    for filename in args:
        cg = ftmc.CoarseGrainRNA(filename)

    sequence_string = cg.seq
    structure_string = cg.to_dotbracket_string()
    extra_links_string = extract_extra_links(cg, options.distance)
    print output_template.format(sequence_string, structure_string, extra_links_string)


if __name__ == '__main__':
    main()

