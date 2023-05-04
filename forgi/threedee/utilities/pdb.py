from __future__ import print_function, unicode_literals

from builtins import zip
from builtins import range
import sys
import re
import warnings
import itertools
import string
import math
import numpy as np
import Bio.PDB as bpdb
from collections import defaultdict

import forgi.utilities.debug as fud
import forgi.threedee.utilities.vector as ftuv
from forgi.threedee.utilities.modified_res import to_4_letter_alphabeth
import forgi.graph.residue as fgr

from logging_exceptions import log_to_exception

import logging
log = logging.getLogger(__name__)


class AtomName(str):
    """
    Like a string, but "C1'" and "C1*" compare equal
    """

    def __eq__(self, other):
        if self.endswith("*"):
            self = self[:-1] + "'"
        if other.endswith("*"):
            other = other[:-1] + "'"
        return str(self) == str(other)

    def __hash__(self):
        if self.endswith("*"):
            self = self[:-1] + "'"
        return hash(str(self))


backbone_atoms = list(map(AtomName, ['P', "O5'", "C5'", "C4'", "C3'", "O3'"]))
ring_atoms = list(map(AtomName, ["C4'", "C3'", "C2'", "C1'", "O4'", "O2'"]))

nonsidechain_atoms = backbone_atoms + ring_atoms

chi_torsion_atoms = dict()
chi_torsion_atoms['A'] = chi_torsion_atoms['G'] = list(
    map(AtomName, ["O4'", "C1'", "N9", "C4"]))
chi_torsion_atoms['C'] = chi_torsion_atoms['U'] = list(
    map(AtomName, ["O4'", "C1'", "N1", "C2"]))

side_chain_atoms = dict()
side_chain_atoms['U'] = list(
    map(AtomName, ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']))
side_chain_atoms['C'] = list(
    map(AtomName, ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']))

side_chain_atoms['A'] = list(
    map(AtomName, ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']))
side_chain_atoms['G'] = list(
    map(AtomName, ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']))

all_side_chains = set(
    side_chain_atoms['U'] + side_chain_atoms['C'] + side_chain_atoms['A'] + side_chain_atoms['G'])

all_rna_atoms = set(nonsidechain_atoms) | all_side_chains

RNA_RESIDUES = ["A", "U", "G", "C", 'rA', 'rC', 'rG', 'rU', 'DU']

interactions = [(AtomName(a), AtomName(b)) for a, b in map(sorted,
                                                           [('P', "O5'"),
                                                            ('P', 'OP1'),
                                                            ('P', 'O1P'),
                                                            ('P', 'OP2'),
                                                            ('P', 'O2P'),
                                                            ("C2'", "O2'"),
                                                            ("O5'", "C5'"),
                                                            ("C5'", "C4'"),
                                                            ("C4'", "O4'"),
                                                            ("C4'", "C3'"),
                                                            ("O4'", "C1'"),
                                                            ("C3'", "C2'"),
                                                            ("C3'", "O3'"),
                                                            ("C2'", "C1'"),
                                                            ("C1'", "N1"),
                                                            ('N1', 'C2'),
                                                            ('N1', 'C6'),
                                                            ('C6', 'C5'),
                                                            ('C5', 'C4'),
                                                            ('C4', 'O4'),
                                                            ('C4', 'N4'),
                                                            ('C4', 'N3'),
                                                            ('N3', 'C2'),
                                                            ('C2', 'O2'),
                                                            ('C2', 'N2'),
                                                            ("C1'", "N9"),
                                                            ('N9', 'C8'),
                                                            ('N9', 'C4'),
                                                            ('C8', 'N7'),
                                                            ('N7', 'C5'),
                                                            ('C6', 'O6'),
                                                            ('C6', 'N6')])]



def rename_chains_for_pdb(chains):
    """
    :param chains: A dict chain_id:chain
    """
    used_chainids = set(chains.keys())
    def get_available_chainid():
        for c in string.ascii_uppercase:
            if c not in used_chainids:
                used_chainids.add(c)
                return c
        raise ValueError("Too many chains. Cannot convert to old PDB format.")
    for chain in chains.values():
        if len(chain.id)>1:
            chain.id = get_available_chainid()
    return {c.id: c for c in chains.values()}

def trim_chain_between(chain, start_res, end_res):
    '''
    Remove all nucleotides between start_res and end_res, inclusive.

    The chain is modified in place so there is no return value.
    '''
    to_detach = []
    for res in chain:
        if start_res <= res.id[1] and res.id[1] <= end_res:
            to_detach += [res]

    for res in to_detach:
        chain.detach_child(res.id)


def extract_subchains_from_seq_ids(all_chains, seq_ids):
    '''
    Extract a portion of one or more pdb chains.
    Creates a list of new chains which contain only
    the specified residues copied from the original chain.

    The chain ids are not modified.

    :param all_chains: A dictionary {chainid:chains}.
    :param seq_ids: An iterable of complete RESIDS.

    :returns: A dictionary chain-id:Bio.PDB.Chain.Chain objects
    '''
    new_chains = {}
    assert isinstance(all_chains, dict)
    for r in seq_ids:
        if r.chain in new_chains:
            chain = new_chains[r.chain]
        else:
            chain = new_chains[r.chain] = bpdb.Chain.Chain(r.chain)
        try:
            chain.add(all_chains[r.chain][r.resid].copy())
        except KeyError:
            log.info(list(sorted(all_chains[r.chain].child_dict.keys())))
            raise
    return new_chains


def is_covalent(contact):
    '''
    Determine if a particular contact is covalent.

    This does not look at the geometric distance but only at the atom names.

    :param contact: A pair of two Atom objects
    :return: `True` if they are covalently bonded
             `False` otherwise
    '''
    r1 = contact[0].parent
    r2 = contact[1].parent

    r1a = (r1, contact[0])
    r2a = (r2, contact[1])

    if contact[0].name.find('H') >= 0 or contact[1].name.find('H') >= 0:
        return True

    ((r1, c1), (r2, c2)) = sorted((r1a, r2a), key=lambda x: x[0].id[1])

    if r1.id == r2.id:
        if tuple(sorted((c1.name, c2.name))) in interactions:
            return True

    if r2.id[1] - r1.id[1] == 1:
        # neighboring residues
        if c1.name == 'O3*' and c2.name == 'P':
            return True

    return False


def num_noncovalent_clashes(chain):
    '''
    Check if a chain has non-covalent clashes. Non-covalent clashes are found
    when two atoms that aren't covalently linked are within 1.8 A of each other.

    :param chain: The chain to evaluate
    :param return: The number of non-covalent clashes.
    '''
    all_atoms = bpdb.Selection.unfold_entities(chain, 'A')
    ns = bpdb.NeighborSearch(all_atoms)

    contacts = ns.search_all(1.9)

    return len([c for c in contacts if not is_covalent(c)])


def noncovalent_distances(chain, cutoff=0.3):
    '''
    Print out the distances between all non-covalently bonded atoms
    which are closer than cutoff to each other.

    :param chain: The Bio.PDB chain.
    :param cutoff: The maximum distance
    '''
    all_atoms = bpdb.Selection.unfold_entities(chain, 'A')
    ns = bpdb.NeighborSearch(all_atoms)

    contacts = ns.search_all(cutoff)

    return [ftuv.magnitude(c[1] - c[0]) for c in contacts if not is_covalent(c)]



def pdb_rmsd(c1, c2, sidechains=False, superimpose=True ):
    '''
    Calculate the all-atom rmsd between two RNA chains.

    :param c1: A Bio.PDB.Chain
    :param c2: Another Bio.PDB.Chain
    :return: The rmsd between the locations of all the atoms in the chains.
    '''
    c1_list = [cr for cr in c1.get_list() if cr.resname.strip()
               in RNA_RESIDUES]
    c2_list = [cr for cr in c2.get_list() if cr.resname.strip()
               in RNA_RESIDUES]
    return residuelist_rmsd(c1_list, c2_list, sidechains, superimpose )

def residuelist_rmsd(c1_list, c2_list, sidechains=False, superimpose=True):
    import forgi.threedee.model.similarity as ftms

    if len(c1_list) != len(c2_list):
        log.error("c1_list (len %s): %s", c1_list)
        log.error("c2_list (len %s): %s", c2_list)
        raise Exception(
            "Chains of different length. (Maybe an RNA-DNA hybrid?)")

    #c1_list.sort(key=lambda x: x.id[1])
    #c2_list.sort(key=lambda x: x.id[1])
    to_residues = []
    crds1 = []
    crds2 = []
    all_atoms1 = []
    all_atoms2 = []
    for r1, r2 in zip(c1_list, c2_list):
        if sidechains:
            anames = nonsidechain_atoms + \
                side_chain_atoms[r1.resname.strip()]
        else:
            anames = nonsidechain_atoms
        #anames = a_5_names + a_3_names

        for a in anames:
            try:
                at1 = r1[a]
                at2 = r2[a]
            except:
                continue
            else:
                all_atoms1.append(at1)
                all_atoms2.append(at2)
                crds1.append(at1.coord)
                crds2.append(at2.coord)
                to_residues.append(r1)

    diff_vecs = ftms._pointwise_deviation(crds1, crds2)
    dev_per_res = defaultdict(list)
    for i, res in enumerate(to_residues):
        dev_per_res[res].append(diff_vecs[i])

    if superimpose:
        sup = bpdb.Superimposer()
        sup.set_atoms(all_atoms1, all_atoms2)

        return (len(all_atoms1), sup.rms, sup.rotran, dev_per_res)
    else:
        return (len(all_atoms1), ftuv._vector_set_rmsd(crds1, crds2), None, dev_per_res)


def get_first_chain(filename):
    '''
    Load a PDB file using the Bio.PDB module and return the first chain.

    :param filename: The path to the pdb file
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = bpdb.PDBParser(PERMISSIVE=False).get_structure('t', filename)
        return list(s.get_chains())[0]


def pdb_file_rmsd(fn1, fn2):
    '''
    Calculate the RMSD of all the atoms in two pdb structures.

    :param fn1: The first filename.
    :param fn2: The second filename.
    :return: The rmsd between the two structures.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        s1 = bpdb.PDBParser().get_structure('t', fn1)
        s2 = bpdb.PDBParser().get_structure('t', fn2)

    c1, _, _ = get_biggest_chain(fn1)
    c2, _, _ = get_biggest_chain(fn2)

    rmsd = pdb_rmsd(c1, c2)

    return rmsd


def renumber_chain(chain, resids=None):
    '''
    Renumber all the residues in this chain so that they start at 1 and end at
    len(chain)

    :param chain: A Bio.PDB.Chain object
    :return: The same chain, but with renamed nucleotides
    '''

    counter = 1

    if resids is None:
        resids = [(' ', i + 1, ' ') for i in range(len(chain))]

    new_child_dict = dict()
    new_child_list = []

    for res, r_new in zip(chain, resids):
        res.id = r_new
        new_child_dict[res.id] = res
        new_child_list.append(res)

    chain.child_dict = new_child_dict
    chain.child_list = new_child_list

    return chain


def output_chain(chain, filename, fr=None, to=None):
    '''
    Dump a chain to an output file. Remove the hydrogen atoms.

    :param chain: The Bio.PDB.Chain to dump.
    :param filename: The place to dump it.
    '''
    class HSelect(bpdb.Select):
        def accept_atom(self, atom):
            if atom.name.find('H') >= 0:
                return False
            else:
                return True
    m = bpdb.Model.Model(' ')
    s = bpdb.Structure.Structure(' ')

    m.add(chain)
    s.add(m)

    io = bpdb.PDBIO()
    io.set_structure(s)
    io.save(filename, HSelect())


def output_multiple_chains(chains, filename, file_type="pdb"):
    '''
    Dump multiple chains to an output file. Remove the hydrogen atoms.

    :param chains: An iterable of Bio.PDB.Chain to dump.
    :param filename: The place to dump it.
    '''
    class HSelect(bpdb.Select):
        def accept_atom(self, atom):
            if atom.name.find('H') >= 0:
                return False
            else:
                return True
    m = bpdb.Model.Model(0)
    s = bpdb.Structure.Structure('stru')
    for chain in chains:
        log.debug("Adding chain %s with %s residues", chain.id, len(chain))
        m.add(chain)
        if file_type=="pdb" and len(chain.id)!=1:
            raise ValueError("Cannot save chain with name %s (not a single character) "
                             "in PDB format. Use cif format instead!")
    s.add(m)
    if file_type == "pdb":
        io = bpdb.PDBIO()
    else:
        io = bpdb.MMCIFIO()
    io.set_structure(s)
    try:
        io.save(filename, HSelect())
    except Exception as e:
        with log_to_exception(log, e):
            log.error("Could not output PDB with chains and residues:")
            for chain in s[0]:
                log.error("%s: %s", chain.id, [r.id for r in chain])
        raise


def get_particular_chain(in_filename, chain_id, parser=None):
    '''
    Load a PDB file and return a particular chain.

    :param in_filename: The name of the pdb file.
    :param chain_id: The id of the chain.
    :return: A Bio.PDB.Chain object containing that particular chain.
    '''
    chains, mr, ir = get_all_chains(in_filename, parser)
    chain, = [c for c in chains if c.id == chain_id]
    return chain, mr, ir


def get_biggest_chain(in_filename, parser=None):
    '''
    Load the PDB file located at filename, select the longest
    chain and return it.

    :param in_filename: The location of the original file.
    :return: A Bio.PDB chain structure corresponding to the longest
             chain in the structure stored in in_filename
    '''
    chains, mr, ir = get_all_chains(in_filename, parser)
    biggest = 0
    biggest_len = 0

    for i in range(len(chains)):
        c = chains[i]

        # Only count RNA residues
        num_residues = 0
        for res in c:
            if (res.resname.strip() == 'A' or
                res.resname.strip() == 'C' or
                res.resname.strip() == 'G' or
                    res.resname.strip() == 'U'):
                num_residues += 1

        if num_residues > biggest_len:
            biggest = i
            biggest_len = num_residues

    # sys.exit(1)

    orig_chain = chains[biggest]
    return orig_chain, mr, ir

def _extract_symmetrymatrices_from_cif_dict(cif_dict):
    """
    Code originally by G. Entzian, modified

    Extract matrices for symmetry operations from the cif-dict to
    convert the assymetric unit to the biological unit

    :returns: A tuple matrices, vectors.
              matrices is a dictionary {operation_id: rot_matrix}
              vectors a dictionary {operation_id: vector}
    """
    symmetry_ids = cif_dict["_pdbx_struct_oper_list.id"]
    matrices = defaultdict(lambda: [[0.0 for x in range(3)] for y in range(3)])
    vectors = defaultdict(lambda: [0.0 for x in range(3)])
    for i in range(1,4):
        v_key = '_pdbx_struct_oper_list.vector['+str(i)+']'
        value_vectors = cif_dict[v_key]
        if not isinstance(value_vectors, list):
            value_vectors = [value_vectors]
        for k, v in enumerate(value_vectors):
            vectors[symmetry_ids[k]][i-1] = float(v)

        for j in range(1,4):
            key = '_pdbx_struct_oper_list.matrix['+str(i)+']['+str(j)+']'
            value_per_matrix = cif_dict[key]
            if not isinstance(value_per_matrix, list):
                value_per_matrix = [value_per_matrix]
            for k, v in enumerate(value_per_matrix):
                matrices[symmetry_ids[k]][i-1][j-1] = float(v)
    return matrices, vectors

def _convert_cif_operation_id_expression(expression):
    """
    By Gregor Entzian
    """
    tmp_expr = expression
    op_id_list = []
    match = re.findall("\([\w\-\,\d]+\)", tmp_expr)#("(\d+)\-(\d+)", id)
    if match:
        for m in match:
            match_s = re.search("\((\d+)\-(\d+)\)", m)
            if match_s:
                f = match_s.group(1)
                t = match_s.group(2)
                f = int(f)
                t = int(t)
                for n in range(f,t+1):
                    op_id_list.append(str(n))
                tmp_expr = tmp_expr.replace(m,",")
                continue
            else:
                #split ',' in case of more entries. In case of one entry use the same code.
                l = m[1:-1] # remove parentheses
                parts = l.split(',')
                for p in parts:
                    op_id_list.append(str(p))
                tmp_expr = tmp_expr.replace(m,",")
                continue
    if tmp_expr != None and tmp_expr != "":
        tmp_expr = tmp_expr.replace("(","").replace(")","")
        parts = tmp_expr.split(',')
        for p in parts:
            if ',' not in p and p != '':
                match_s = re.search("(\d+)\-(\d+)", p)
                if match_s:
                    f = match_s.group(1)
                    t = match_s.group(2)
                    f = int(f)
                    t = int(t)
                    for n in range(f,t+1):
                        op_id_list.append(str(n))
                    continue
                op_id_list.append(str(p))
    return tuple(op_id_list)

def _extract_assembly_gen(cif_dict):
    """
    Extracts the information, how assemblies are generated out
    of the chains and symmetry operations, from the cifdict.

    :returns: A dictionary {assembly_id: [(chainid, operationids)...]}
              Where operation_ids is a list
    """
    assembly_ids = cif_dict['_pdbx_struct_assembly_gen.assembly_id']       #1
    operation_ids = cif_dict['_pdbx_struct_assembly_gen.oper_expression']   #1,2
    chain_id_lists = cif_dict['_pdbx_struct_assembly_gen.asym_id_list']      #A,B

    if not isinstance(assembly_ids, list):
        assembly_ids=[assembly_ids]
        assert not isinstance(operation_ids, list)
        operation_ids=[operation_ids]
        assert not isinstance(chain_id_lists, list)
        chain_id_lists=[chain_id_lists]

    assembly_components = defaultdict(lambda:defaultdict(list)) # assembly: operation-lists: chains
    for i, aid in enumerate(assembly_ids):
        op_ids=_convert_cif_operation_id_expression(operation_ids[i])
        chain_ids=chain_id_lists[i].split()
        assembly_components[0][operation_ids].extend(chain_id_lists)
    return assembly_components

def _get_new_chainid(chain_id, op_id, taken_ids):
    new_id= chain_id+"op"+str(op_id)
    while new_id in taken_ids:
        new_id+="a"
    return new_id

def _extract_asym_auth_id_map(cif_dict):
    label_ids = cif_dict['_atom_site.label_asym_id']
    auth_ids = cif_dict['_atom_site.auth_asym_id']
    l2a=defaultdict(lambda: None)
    a2l=defaultdict(list)
    for i, lid in enumerate(label_ids):
        if lid in l2a:
            if l2a[lid]!=auth_ids[i]:
                log.error("lid %s, authid %s but before %s", lid, auth_ids[i], l2a[lid])
                raise ValueError("Inconsistent cif.")
        else:
            l2a[lid]=auth_ids[i]
            a2l[auth_ids[i]].append(lid)
    log.debug(l2a)
    return l2a, a2l

def _get_assemblies(chains, cifdict):
    assemblies = _extract_assembly_gen(cifdict)
    id2chainid, chainid2ids = _extract_asym_auth_id_map(cifdict)
    if assembly_nr is None:
        assembly_nr = list(sorted(assemblies.keys()))[0]
    assembly_gen = assemblies[assembly_nr]
    old_chains = { c.id:c for c in chains}
    new_chains = {}
    log.debug("Original chains are %s. Now performing symmetry "
              "for assembly %s", old_chains, assembly_nr)
    log.debug("AG: %s",assembly_gen)
    for operations, labelids in assembly_gen[0].items():
        chainids = set(id2chainid[lid] for lid in labelids)
        lids_back = set( x for x in a2l for cid in chainids for a2l in chainid2ids[cid])
        if lids_back!=set(labelids):
            log.error("%s, %s, extra: %s", lids_back, labelids, lids_back-set(labelids))
            raise ValueError("Operation on part of an author designated assymetric unit "
                             "(auth_asym_id) nt supported.")
        for chain_id in chainids:
            if chain_id not in old_chains:
                log.debug ("Skipping chain %s: not RNA? chains: %s", chain_id, old_chains.keys())
                continue
            for op_id in op_ids:
                log.debug("Applying op %s to %s", op_id, chain_id)
                if chain_id in new_chains:
                    chain = old_chains[chain_id].copy()
                    newid = _get_new_chainid(chain_id, op_id,
                                             set(old_chains.keys())| set(new_chains.keys()))
                    log.debug("Setting id %s", newid)
                    chain.id = newid
                    chain.transform(operation_mat[op_id], operation_vec[op_id])
                new_chains[chain.id]=chain
    log.debug("new_chains: %s", new_chains)
    chains = list(new_chains.values())
    return chains

def get_all_chains(in_filename, parser=None, no_annotation=False, assembly_nr=None):
    '''
    Load the PDB file located at filename, read all chains and return them.

    :param in_filename: The location of the original file.
    :param assembly_nr: Which assembly to return. Default: The first.
    :return: a tuple chains, missing_residues

             * chains: A list of Bio.PDB chain structures corresponding to all
                       RNA structures stored in in_filename
             * missing_residues: A list of dictionaries, describing the missing residues.
             * interacting residues: A list of residues
    '''
    if parser is None:
        if in_filename.endswith(".pdb"):
            parser = bpdb.PDBParser()
        elif in_filename.endswith(".cif"):
            parser = bpdb.MMCIFParser()
        else:  # Cannot determine filetype by extention. Try to read first line.
            with open(in_filename) as pdbfile:
                line = pdbfile.readline(20)
                # According to
                # page 10 of ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf
                # a HEADER entry is mandatory. Biopython sometime starts directly with ATOM
                if line.startswith("HEADER") or line.startswith("ATOM"):
                    parser = bpdb.PDBParser()
                else:
                    parser = bpdb.MMCIFParser()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = parser.get_structure('temp', in_filename)
    if len(s) > 1:
        warnings.warn("Multiple models in file. Using only the first model")

    # Let's detach all H2O, to speed up processing.
    for chain in s[0]:
        log.debug("Before detaching water from %s: chain has %s residues", chain.id, len(chain))
        for r in chain.child_list[:]: # We need a copy here, because we are modifying it during iteration
            if r.resname.strip() == "HOH":
                chain.detach_child(r.id)
        log.debug("After detaching water from %s: chain has %s residues", chain.id, len(s[0][chain.id]))

    # Rename residues from other programs
    for chain in s[0]:
        for r in chain:
            # rename rosetta-generated structures
            # Old: with whitespace
            if r.resname == ' rA':
                r.resname = '  A'
            elif r.resname == ' rC':
                r.resname = '  C'
            elif r.resname == ' rG':
                r.resname = '  G'
            elif r.resname == ' rU':
                r.resname = '  U'
            # newer biopython: No whitespace
            if r.resname == 'rA':
                r.resname = 'A'
            elif r.resname == 'rC':
                r.resname = 'C'
            elif r.resname == 'rG':
                r.resname = 'G'
            elif r.resname == 'rU':
                r.resname = 'U'
            # rename iFoldRNA-generated structures
            if r.resname == 'ADE':
                r.resname = '  A'
            elif r.resname == 'CYT':
                r.resname = '  C'
            elif r.resname == 'GUA':
                r.resname = '  G'
            elif r.resname == 'URI':
                r.resname = '  U'
    # Now search for protein interactions.
    if not no_annotation:
        interacting_residues = enumerate_interactions_kdtree(s[0])
    else:
        interacting_residues = set()

    # The chains containing RNA
    chains = list(chain for chain in s[0] if contains_rna(chain))
    
    if isinstance(parser, bpdb.MMCIFParser):
        cifdict = parser._mmcif_dict # We read a private attribute here, because parsing the mmcif dictionary a second time would cause a performance penalty.
        # Generate an assembly
        try:
            operation_mat, operation_vec = _extract_symmetrymatrices_from_cif_dict(cifdict)
        except KeyError:
            pass
        else:
            if False: # Still experimental and not working correctly.
                chains = _get_assemblies(chains, cifdict)
        mr = []
        try:
            mask = np.array(
                cifdict["_pdbx_poly_seq_scheme.pdb_mon_id"], dtype=str) == "?"
            int_seq_ids = np.array(
                cifdict["_pdbx_poly_seq_scheme.pdb_seq_num"], dtype=int)[mask]
            cs = np.array(
                cifdict["_pdbx_poly_seq_scheme.pdb_strand_id"], dtype=str)[mask]
            insertions = np.array(
                cifdict["_pdbx_poly_seq_scheme.pdb_ins_code"], dtype=str)[mask]
            insertions[insertions == "."] = " "
            symbol = np.array(
                cifdict["_pdbx_poly_seq_scheme.mon_id"], dtype=str)[mask]
        except KeyError:
            pass
        else:
            if not no_annotation:
                for i, sseq in enumerate(int_seq_ids):
                    mr.append({
                        "model": None,
                        "res_name": symbol[i],
                        "chain": cs[i],
                        "ssseq": sseq,
                        "insertion": insertions[i]
                    })

    else:
        log.debug("PDB header %s", parser.header)
        try:
            mr = parser.header["missing_residues"]
            if assembly_nr is not None:
                warnings.warn("Getting an assembly is not supported for the old PDB format.")
        except KeyError:
            mr = []
            with open(in_filename) as f:
                for wholeline in f:
                    if wholeline.startswith("REMARK 465"):
                        line = wholeline[10:].strip()
                        mr_info = _parse_remark_465(line)
                        if mr_info is not None:
                            mr.append(mr_info)
                    else:
                        continue

    if mr:
        log.info("This PDB has missing residues")
    elif not no_annotation:
        log.info("This PDB has no missing residues")
    log.debug("LOADING DONE: chains %s, mr %s, ir: %s",
             chains, mr, interacting_residues)
    return chains, mr, interacting_residues



def _parse_remark_465(line):
    """Parse missing residue remarks.
    Returns a dictionary describing the missing residue.
    The specification for REMARK 465 at
    http://www.wwpdb.org/documentation/file-format-content/format33/remarks2.html#REMARK%20465
    only gives templates, but does not say they have to be followed.
    So we assume that not all pdb-files with a REMARK 465 can be understood.
    Returns a dictionary with the following keys:
    "model", "res_name", "chain", "ssseq", "insertion"
    """
    if line:
        # Note that line has been stripped.
        assert line[0] != " " and line[-1] not in "\n ", "line has to be stripped"
    pattern = (r"""
                (\d+\s[\sA-Z][\sA-Z][A-Z] |   # Either model number + residue name
                 [A-Z]?[A-Z]?[A-Z])           # Or only residue name with
                                              # 1 (RNA) to 3 letters
                \s ([A-Za-z0-9])              # A single character chain
                \s+(\d+[A-Za-z]?)$            # Residue number: A digit followed
                                              # by an optional insertion code
                                              # (Hetero-flags make no sense in
                                              # context with missing res)
                """)
    match = re.match(pattern, line, re.VERBOSE)
    if match is not None:
        residue = {}
        if " " in match.group(1):
            model, residue["res_name"] = match.group(1).split(" ")
            residue["model"] = int(model)
        else:
            residue["model"] = None
            residue["res_name"] = match.group(1)
        residue["chain"] = match.group(2)
        try:
            residue["ssseq"] = int(match.group(3))
        except ValueError:
            residue["insertion"] = match.group(3)[-1]
            residue["ssseq"] = int(match.group(3)[:-1])
        else:
            residue["insertion"] = None
        return residue
    return None

def enumerate_interactions_kdtree(model):
    relevant_atoms = [a for a in model.get_atoms() if a.name[0] in [
        "C", "N", "O"]]
    if not relevant_atoms:
        return set()
    kdtree = bpdb.NeighborSearch(relevant_atoms)
    pairs = kdtree.search_all(6, "A")
    res_pair_list = set()
    for a1, a2 in pairs:
        if a1.name not in all_side_chains and a2.name not in all_side_chains:
            continue
        p1 = a1.get_parent()
        p2 = a2.get_parent()
        if p1.id == p2.id:
            continue
        elif p1 < p2:
            res_pair_list.add((p1, p2))
        else:
            res_pair_list.add((p2, p1))
    interacting_residues = set()
    for res1, res2 in res_pair_list:
        rna_res = None
        other_res = None
        if res1.resname.strip() in RNA_RESIDUES and not res1.id[0].startswith("H_"):
            rna_res = res1
        else:
            other_res = res1
        if res2.resname.strip() in RNA_RESIDUES and not res2.id[0].startswith("H_"):
            rna_res = res2
        else:
            other_res = res2
        if rna_res is None or other_res is None:
            continue
        log.debug("%s(chain %s) and %s(chain %s, resname %s) are close", rna_res,rna_res.parent.id,other_res, other_res.parent.id, other_res.resname.strip())
        # Only consider C and N. So no ions etc
        if any(a.name[0] in ["C", "N"] for a in other_res.get_atoms()):
            interacting_residues.add(rna_res)
        else:
            log.debug("but %s has wrong atoms %s", other_res,
                      list(a.name for a in other_res.get_atoms()))
    log.debug("Interacting: %s", interacting_residues)
    return interacting_residues


"""def residues_interact(rna_res, other_res):
    for rna_atom in rna_res:
        if rna_atom.get_name() in all_side_chains:
            for other_atom in other_res:
                atom_symbol="".join(s for s in other_atom.get_name() if not s.isdigit())
                if atom_symbol in ["C", "N"]:
                    d=ftuv.vec_distance(rna_atom.coord, other_atom.coord)
                    if d<6:
                        return True
    return False"""

HBOND_CUTOFF = 4.5  # 4.5 and 0.9 are values optimized against DSSR for 5T5H_A-B-C
OOP_CUTOFF = 0.9


def _get_points(res1, res2):
    labels = {res1.resname.strip(), res2.resname.strip()}
    if labels == {"A", "U"}:
        return _points_AU(res1, res2)
    elif labels == {"G", "C"}:
        return _points_GC(res1, res2)
    elif labels == {"G", "U"}:
        return _points_GC(res1, res2)
    else:
        return None


def _points_AU(res1, res2):
    if res1.resname.strip() == "A":
        resA = res1
        resU = res2
    else:
        resA = res2
        resU = res1
    a = resA["N6"].coord
    b = resU["O4"].coord
    c = resA["N1"].coord
    d = resU["N3"].coord
    return (resA["C8"].coord, resU["C6"].coord), (a, b), (c, d)


def _points_GU(res1, res2):
    if res1.resname.strip() == "G":
        resG = res1
        resU = res2
    else:
        resG = res2
        resU = res1
    a = resG["O6"].coord
    b = resU["N3"].coord
    c = resG["N1"].coord
    d = resU["O2"].coord
    return (resG["C8"].coord, resU["C6"].coord), (a, b), (c, d)


def _points_GC(res1, res2):
    if res1.resname.strip() == "G":
        resG = res1
        resC = res2
    else:
        resG = res2
        resC = res1
    a = resG["O6"].coord
    c = resG["N1"].coord
    e = resG["N2"].coord
    b = resC["N4"].coord
    d = resC["N3"].coord
    f = resC["O2"].coord
    return (resC["C6"].coord, resG["C8"].coord), (a, b), (c, d), (e, f)


def is_basepair_pair(res1, res2):

    pairs = _get_points(res1, res2)
    if not pairs:
        return False
    for pair in pairs[1:]:  # pairs[0] is only for coplanarity]]
        d = ftuv.vec_distance(pair[0], pair[1])
        if d >= HBOND_CUTOFF:
            return False
    if is_almost_coplanar(*[point for pair in pairs for point in pair]):
        return True
    return False


def _coplanar_point_indices(*points):
    """ Thanks to https://stackoverflow.com/a/18968498"""
    from numpy.linalg import svd
    points = np.array(points).T
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(
        points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)  # Could also use np.cov(x) here.
    normal = svd(M)[0][:, -1]

    out = []
    for i, p in enumerate(points.T):
        w = p - ctr
        oop_distance = ftuv.magnitude(
            np.dot(w, normal)) / ftuv.magnitude(normal)
        if oop_distance <= OOP_CUTOFF:
            out.append(i)
    return out, ctr, normal


def is_almost_coplanar(*points):
    indices, c, n = _coplanar_point_indices(*points)
    return len(indices) == len(points)


def annotate_fallback(chain_list):
    """
    If neither DSSR nor MC-Annotate are available, we use an ad-hoc implementation of canonical
    basepair detection as fallback.
    This does not work well for missing atoms or modified residues.
    """
    kdtree = bpdb.NeighborSearch(
        [atom for chain in chain_list for atom in chain.get_atoms()])
    pairs = kdtree.search_all(10, "R")
    basepairs = {}
    # Sorted, so conflicting basepairs are deterministically solved
    for res1, res2 in sorted(pairs):
        if res1.resname.strip() not in RNA_RESIDUES or res1.id[0].startswith("H_"):
            continue
        if res2.resname.strip() not in RNA_RESIDUES or res2.id[0].startswith("H_"):
            continue
        labels = {res1.resname.strip(), res2.resname.strip()}
        try:
            is_bp = is_basepair_pair(res1, res2)
            if is_bp:
                res1_id = fgr.resid_from_biopython(res1)
                res2_id = fgr.resid_from_biopython(res2)
                if res1_id in basepairs:
                    warnings.warn("More than one basepair detected for {}."
                                  " Ignoring {}-{} because {}-{} is already"
                                  " part of the structure".format(res1_id, res1_id, res2_id, res1_id, basepairs[res1_id]))
                    continue
                if res2_id in basepairs:
                    warnings.warn("More than one basepair detected for {}."
                                  " Ignoring {}-{} because {}-{} is already"
                                  " part of the structure".format(res2_id, res2_id, res1_id, res2_id, basepairs[res2_id]))
                    continue
                basepairs[res1_id] = res2_id
                basepairs[res2_id] = res1_id
        except KeyError as e:
            log.debug("Missing atom %s. %s has atoms %s, %s has atoms %s",
                      e, res1, res1.child_dict, res2, res2.child_dict)
            pass

    seq_ids = []
    for chain in sorted(chain_list, key=lambda x: x.id):
        for residue in chain:
            seq_ids.append(fgr.resid_from_biopython(residue))
    bpseq = ""
    chain_dict = {c.id: c for c in chain_list}
    for i, seqid in enumerate(seq_ids):
        if seqid in basepairs:
            bp = seq_ids.index(basepairs[seqid]) + 1
        else:
            bp = 0

        bpseq += "{} {} {}\n".format(i + 1,
                                     chain_dict[seqid.chain][seqid.resid].resname.strip(
                                     ),
                                     bp)
    return bpseq, seq_ids


def rename_rosetta_atoms(chain):
    '''
    Rosetta names all the backbone atoms with an asterisk rather than an
    apostrophe. All that needs to be reversed.

    :param chain. A Bio.PDB.Chain structure generated by Rosetta
    :return: The same chain with renamed atoms
    '''
    for a in bpdb.Selection.unfold_entities(chain, 'A'):
        oldid = a.id
        a.name = a.name.replace('*', "'")
        a.fullname = a.name.replace('*', "'")
        a.id = a.id.replace('*', "'")
        #: Not needed with newer biopython versions
        #: Seems to be needed again?
        del a.parent.child_dict[oldid]
        a.parent.child_dict[a.id] = a
    # log.debug("Replaced rosetta atoms. \n%s\n%s",
    #                    chain.child_list[0].child_list,
    #                    chain.child_list[0].child_dict
    #                    )

    return chain


def remove_disordered(chain):
    for i, residue in enumerate(chain):
        if hasattr(residue, "selected_child"):
            new_res = residue.selected_child
            chain.detach_child(residue.id)
            chain.insert(i, new_res)
            residue = new_res
        for j, atom in enumerate(residue):
            if hasattr(atom, "selected_child"):
                new_atom = atom.selected_child
                new_atom.altloc = " "
                new_atom.occupancy = 1.0
                new_atom.disordered_flag = 0
                residue.detach_child(atom.id)
                residue.insert(j, new_atom)
    return chain


def load_structure(pdb_filename):
    '''
    Load a Bio.PDB.Structure object and return the largest chain.
    This chain will be modified so that all hetatms are removed, modified
    residues will be renamed to regular residues, etc...

    '''
    chain, mr, ir = get_biggest_chain(pdb_filename)
    return clean_chain(chain, True)[0]


def clean_chain(chain, query_PDBeChem=False):
    """
    Clean a pdb chain for further use with forgi.

    It will be modified so that all hetatms are removed, modified
    residues will be renamed to regular residues, residue ids will be positive integers, ...

    :param chain: A Bio.PDB.Chain object
    :param query_PDBeChem: If true, query the PDBeChem database whenever a
                          modified residue with unknown 3-letter code
                          is encountered.

    :returns: A modified version of this chain
    """
    chain, modifications = to_4_letter_alphabeth(chain, query_PDBeChem)
    chain = rename_rosetta_atoms(chain)
    chain = remove_disordered(chain)
    return chain, modifications


def interchain_contacts(struct):
    all_atoms = bpdb.Selection.unfold_entities(struct, 'A')

    ns = bpdb.NeighborSearch(all_atoms)
    pairs = ns.search_all(2.8)

    ic_pairs = []

    for (a1, a2) in pairs:
        if a1.parent.parent != a2.parent.parent:
            ic_pairs += [(a1, a2)]
    return ic_pairs


def contains_rna(chain):
    '''
    Determine if a Bio.PDB.Chain structure corresponds to an RNA
    molecule.

    :param chain: A Bio.PDB.Chain molecule
    :return: True if it is an RNA molecule, False if at least one residue is not an RNA.
    '''
    for res in chain:
        if res.resname.strip() in RNA_RESIDUES:
            return True
    return False


def is_protein(chain):
    '''
    Determine if a Bio.PDB.Chain structure corresponds to an protein
    molecule.

    :param chain: A Bio.PDB.Chain molecule
    :return: True if it is a protein molecule, False otherwise
    '''
    for res in chain:
        if res.resname in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
            return True
    return False
