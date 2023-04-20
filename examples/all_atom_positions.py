#!/usr/bin/python

from __future__ import print_function
from builtins import map
import collections as c
import itertools as it
import sys
import argparse
import logging

import numpy as np

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.commandline_utils as fuc

log = logging.getLogger(__name__)


def main(parser):
    args = parser.parse_args()

    poss = c.defaultdict(list)
    sources = c.defaultdict(list)

    cgs = fuc.cgs_from_args(
        args, rna_type="pdb", enable_logging=True)
    for i, cg in enumerate(cgs):

        if len(list(cg.stem_iterator())) == 0:
            log.warning("Skipping RNA %s (%s): no stems", i, cg.pdb_name)
            continue

        for d in cg.defines.keys():
            if np.allclose(cg.coords[d][0], cg.coords[d][1]):
                log.warning(
                    "Skipping element %s of RNA %s (%s): degenerate coordinates.", d, i, cg.pdb_name)
                continue

            origin, basis = ftug.element_coord_system(cg, d)

            if d[0] == 'i' or d[0] == 'm':
                conn = cg.connections(d)
                conn_type = cg.connection_type(d, conn)
            else:
                conn_type = 0

            for i, r in zip(it.count(),
                            cg.define_residue_num_iterator(d)):

                # add only the base atoms which are relevant to the calculation
                # of the chi torsion angle
                seq_id = cg.seq.to_resid(r - 1)
                resname = cg.chains[seq_id.chain][seq_id.resid].resname.strip()

                if resname not in ftup.chi_torsion_atoms.keys():
                    print("Unknown nucleotide name:", resname, file=sys.stderr)
                    continue

                atoms = ftup.nonsidechain_atoms + \
                    ftup.chi_torsion_atoms[resname][-2:]
                scatoms = ftup.side_chain_atoms[resname]
                for aname in atoms + scatoms:
                    try:
                        resid = cg.seq.to_resid(r - 1)
                        a = cg.chains[resid.chain][resid.resid][aname]
                    except KeyError as ke:
                        # missing an atom
                        continue

                    # The C1'->B1 and B1->B2 vectors define the plane of the base
                    # The O4'->C1'->B1->B2 sequence defines the torsion
                    # angle chi
                    if aname == ftup.chi_torsion_atoms[resname][-2]:
                        aname = 'B1'
                    elif aname == ftup.chi_torsion_atoms[resname][-1]:
                        aname = 'B2'
                    elif aname in scatoms:
                        aname = resname + "." + aname
                    avec = a.get_vector().get_array()
                    atom_pos = ftuv.change_basis(
                        avec - origin, basis, ftuv.standard_basis)
                    identifier = "%s %s %d %d %s" % (d[0],
                                                     " ".join(
                                                         map(str, cg.get_node_dimensions(d))),
                                                     conn_type,
                                                     i, aname)
                    poss[identifier] += [atom_pos]
                    sources[identifier] += [d]

                    print("{}:{}".format(identifier, ",".join(map(str, atom_pos))))


parser = fuc.get_rna_input_parser("Generate a list of all atom positions found in the pdb files",
                                  nargs="+", rna_type="pdb", enable_logging=True)

if __name__ == '__main__':
    main(parser)
