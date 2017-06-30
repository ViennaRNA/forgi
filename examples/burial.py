#!/usr/bin/python

from __future__ import print_function
import Bio.PDB as bpdb
import Bio.KDTree as bkd

import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.pdb as ftup
import forgi.utilities.debug as fud

import numpy as np

import sys
from optparse import OptionParser

def calculate_burial(kd, position, directions, radius=5, step_size=1., max_dist = 500):
    '''
    Calculate how many atoms we have to go by if we travel in all
    directions away from this point)

    :param kd: A Bio.KDTree containing all the other atoms in this
              molecule
    :param position: The position to start at
    :param directions: The directions to travel in
    :param radius: The radius within which to check for neighbors
    :param step_size: How far to step at each point in the simulation
    :param max_dist: The minimum distance to travel from the center
    '''
    distances = []

    for direction in directions:
        # we assume that each direction is a unit vector
        step_vector = step_size * direction
        new_pos = np.array(position)

        while ftuv.magnitude((new_pos - position)) < max_dist:
            kd.search(new_pos, float(radius))
            #kd.all_search(10)
            neighbors = kd.get_indices()

            distance = ftuv.magnitude(new_pos - position)
            if len(neighbors) == 0:
                distances += [ftuv.magnitude(new_pos - position)]
                break

            new_pos += step_vector

    return min(distances)

def main():
    usage = """
    python burial.py pdb_file

    Calculate how buried each nucleotide is in this PDB file.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-r', '--radius', dest='radius', default=5, help="The radius of the search ball", type='float')
    parser.add_option('-s', '--step-size', dest='step_size', default=1, help="The size of each step in the depth search", type='float')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    all_directions = np.array(ftuv.GetPointsEquiAngularlyDistancedOnSphere(numberOfPoints = 45))

    chain = ftup.load_structure(args[0])
    atoms = bpdb.Selection.unfold_entities(chain, 'A')
    atom_coords = np.array([a.get_coord() for a in atoms])

    kd = bkd.KDTree(3)
    kd.set_coords(atom_coords)

    for i,r in enumerate(chain.get_list()):
        distance = calculate_burial(kd, r["C1'"].get_coord(), 
                                    all_directions, radius=options.radius,
                                   step_size = options.step_size)

        print("{}:{}".format(i+1, distance))
        pass


if __name__ == '__main__':
    main()

