from __future__ import print_function

import sys, warnings
import numpy as np
import Bio.PDB as bpdb

import forgi.utilities.debug as fud
import forgi.threedee.utilities.vector as ftuv

import logging
log=logging.getLogger(__name__)

backbone_atoms_real = ['P', "O5'", "C5'", "C4'", "C3'", "O3'"] 
backbone_atoms = ['P', 'O5*', 'C5*', 'C4*', 'C3*', 'O3*']
backbone_atoms += ['P', "O5'", "C5'", "C4'", "C3'", "O3'"]
ring_atoms = ['C4*', 'C3*', 'C2*', 'C1*', 'O4*']
ring_atoms_real = ["C4'", "C3'", "C2'", "C1'", "O4'"]

nonsidechain_atoms = backbone_atoms_real + ring_atoms_real

chi_torsion_atoms = dict()
chi_torsion_atoms['A'] = ["O4'", "C1'", "N9", "C4"]
chi_torsion_atoms['G'] = chi_torsion_atoms['A']
chi_torsion_atoms['C'] = ["O4'", "C1'", "N1", "C2"]
chi_torsion_atoms['U'] = chi_torsion_atoms['C']

side_chain_atoms = dict()
side_chain_atoms['U'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
side_chain_atoms['C'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']

side_chain_atoms['A'] = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']
side_chain_atoms['G'] = ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']

all_side_chains = set(side_chain_atoms['U'] + side_chain_atoms['C'] + side_chain_atoms['A'] + side_chain_atoms['G'])

all_rna_atoms = backbone_atoms_real + ring_atoms_real
for v in side_chain_atoms.values():
    all_rna_atoms += v
all_rna_atoms = set(all_rna_atoms)

interactions = [('P', 'O5*'),
                ('P', "O5'"),
                ('P', 'OP1'),
                ('P', 'O1P'),
                ('P', 'OP2'),
                ('P', 'O2P'),
                ('C2*', 'O2*'),
                ("C2'", "O2'"),
                       ('O5*', 'C5*'),
                       ("O5'", "C5'"),
                       ('C5*', 'C4*'),
                       ("C5'", "C4'"),
                       ('C4*', 'O4*'),
                       ("C4'", "O4'"),
                       ('C4*', 'C3*'),
                       ("C4'", "C3'"),
                       ('O4*', 'C1*'),
                       ("O4'", "C1'"),
                       ('C3*', 'C2*'),
                       ("C3'", "C2'"),
                       ('C3*', 'O3*'),
                       ("C3'", "O3'"),
                       ('C2*', 'C1*'),
                       ("C2'", "C1'"),
                       ('C1*', 'N1'),
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
                       ('C1*', 'N9'),
                       ("C1'", "N9"),
                       ('N9', 'C8'),
                       ('N9', 'C4'),
                       ('C8', 'N7'),
                       ('N7', 'C5'),
                       ('C6', 'O6'),
                       ('C6', 'N6')]

interactions_set = [tuple(sorted(i)) for i in interactions]

def trim_chain(chain, start_res, end_res):
    '''
    Remove all residues that are not between start_res and end_res, inclusive.
    '''
    to_detach = []
    for res in chain:
        if res.id[1] <= start_res or end_res <= res.id[1]:
            to_detach += [res]

    for res in to_detach:
        chain.detach_child(res.id)

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

def extract_subchain(chain, start_res, end_res):
    '''
    Extract a portion of a particular chain. The new chain
    will contain residues copied from the original chain.

    :param chain: The source chain.
    :param start_res: The number of the first nucleotide to extract
    :param last_res: The number of the last nucleotide to extract
    '''
    new_chain = bpdb.Chain.Chain(' ')
    for r in chain:
        if start_res <= r.id and r.id <= end_res:
            new_chain.add(r.copy())

    return new_chain

def extract_subchain_from_res_list(chain, res_list):
    '''
    Extract a portion of a particular chain. The new chain
    will contain residues copied from the original chain.

    :param chain: The source chain.
    :param res_list: The list of residue identifiers of the nucleotides
                     to extract
    '''
    new_chain = bpdb.Chain.Chain(' ')
    for r in res_list:
        new_chain.add(chain[r].copy())

    return new_chain

def is_covalent(contact):
    '''
    Determine if a particular contact is covalent.

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
        if tuple(sorted((c1.name, c2.name))) in interactions_set:
            return True

    if r2.id[1] - r1.id[1] == 1:
        #neighboring residues
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

def pdb_rmsd(c1, c2, sidechains=False, superimpose=True, apply_sup=False):
    '''
    Calculate the all-atom rmsd between two RNA chains.

    :param c1: A Bio.PDB.Chain
    :param c2: Another Bio.PDB.Chain
    :return: The rmsd between the locations of all the atoms in the chains.
    '''

    a_5_names = ['P', 'O5*', 'C5*', 'C4*', 'O4*', 'O2*']
    a_5_names += ['P', "O5'", "C5'", "C4'", "O4'", "O2'"]
    a_3_names = ["C1*", "C2*", "C3*", "O3*"]
    a_3_names += ["C1'", "C2'", "C3'", "O3'"]

    a_names = dict()
    a_names['U'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'] + a_3_names
    a_names['C'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'] + a_3_names

    a_names['A'] = a_5_names + ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9'] + a_3_names
    a_names['G'] = a_5_names + ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9'] + a_3_names

    a_names['U'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'] + a_3_names
    a_names['C'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'] + a_3_names

    a_names['A'] = a_5_names + ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9'] + a_3_names
    a_names['G'] = a_5_names + ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9'] + a_3_names

    all_atoms1 = []
    all_atoms2 = []

    acceptable_residues = ['A','C','G','U','rA','rC','rG','rU','DG']
    c1_list = [cr for cr in c1.get_list() if cr.resname.strip() in acceptable_residues]
    c2_list = [cr for cr in c2.get_list() if cr.resname.strip() in acceptable_residues]

    if len(c1_list) != len(c2_list):
        #print >>sys.stderr, "Chains of different length", len(c1.get_list()), len(c2.get_list())
        raise Exception("Chains of different length. (Maybe an RNA-DNA hybrid?)")

    #c1_list.sort(key=lambda x: x.id[1])
    #c2_list.sort(key=lambda x: x.id[1])
    
    for r1,r2 in zip(c1_list, c2_list):
        if sidechains:
            anames = backbone_atoms + a_names[c1[i].resname.strip()]
        else:
            anames = backbone_atoms
        #anames = a_5_names + a_3_names

        for a in anames:
	    try:
		at1 = r1[a]
		at2 = r2[a]

		all_atoms1 += [at1]
		all_atoms2 += [at2]
	    except:
	        continue

    #print "rmsd len:", len(all_atoms1), len(all_atoms2)
    if superimpose:
        sup = bpdb.Superimposer()
        sup.set_atoms(all_atoms1, all_atoms2)

        if apply_sup:
            sup.apply(c2.get_atoms())

        return (len(all_atoms1), sup.rms, sup.rotran)
    else:
        crvs1 = np.array([a.get_vector().get_array() for a in all_atoms1])
        crvs2 = np.array([a.get_vector().get_array() for a in all_atoms2])

        return (len(all_atoms1), ftuv.vector_set_rmsd(crvs1, crvs2), None)

def get_first_chain(filename):
    '''
    Load a PDB file using the Bio.PDB module and return the first chain.

    :param filename: The path to the pdb file
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = bpdb.PDBParser().get_structure('t', filename)
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

        s1= bpdb.PDBParser().get_structure('t', fn1)
        s2= bpdb.PDBParser().get_structure('t', fn2)

    c1 = get_biggest_chain(fn1)
    c2 = get_biggest_chain(fn2)

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
        resids = [(' ', i+1, ' ') for i in range(len(chain))]

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

def output_multiple_chains(chains, filename):
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
    for chain in chains:
        m.add(chain)

    s.add(m)    

    io = bpdb.PDBIO()
    io.set_structure(s)
    io.save(filename, HSelect()) 
    
def get_particular_chain(in_filename, chain_id, parser=None):
    '''
    Load a PDB file and return a particular chain.

    :param in_filename: The name of the pdb file.
    :param chain_id: The id of the chain.
    :return: A Bio.PDB.Chain object containing that particular chain.
    '''
    if parser is None:
        parser = bpdb.PDBParser()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = parser.get_structure('temp', in_filename)

    # always take the first model
    m = s.get_list()[0]

    return m[chain_id]

def get_biggest_chain(in_filename, parser=None):
    '''
    Load the PDB file located at filename, select the longest
    chain and return it.

    :param in_filename: The location of the original file.
    :return: A Bio.PDB chain structure corresponding to the longest
             chain in the structure stored in in_filename
    '''
    if parser is None:
        parser = bpdb.PDBParser()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = parser.get_structure('temp', in_filename)

    chains = list(s.get_chains())
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

        #print c, num_residues
    #sys.exit(1)

    orig_chain = chains[biggest]
    return orig_chain
    
def get_all_chains(in_filename, parser=None):
    '''
    Load the PDB file located at filename, select the longest
    chain and return it.
    
    :param in_filename: The location of the original file.
    :return: A Bio.PDB chain structure corresponding to the longest
             chain in the structure stored in in_filename
    '''
    if parser is None:
        #print("in_filename is {}".format(in_filename), file=sys.stderr)
        if in_filename.endswith(".pdb"):
            parser = bpdb.PDBParser()
        elif in_filename.endswith(".cif"):
            parser = bpdb.MMCIFParser()
        else: #Cannot determine filetype by extention. Try to read first line.
            with open(in_filename) as pdbfile:
                line = pdbfile.readline(20)
                # According to 
                #page 10 of ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf
                # a HEADER entry is mandatory.
                if line.startswith("HEADER"):
                    #print("HEADER found", file=sys.stderr)
                    parser = bpdb.PDBParser()
                else:
                    parser = bpdb.MMCIFParser()
                    #print("HEADER NOT found", file=sys.stderr)


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = parser.get_structure('temp', in_filename)

    chains = list(chain for chain in s.get_chains() if is_rna(chain))
    return chains
    
def change_residue_id(residue, new_id):
    chain = residue.parent
    if new_id in chain:
        raise ValueError("Cannot change id")
    old_id = residue.id
    chain.child_dict[new_id] = residue
    del chain.child_dict[old_id]
    residue.id = new_id
    
def rename_modified_ress(chain):
    '''
    Rename the modified residues so that they have the same
    names as unmodified residues

    :param chain: A Bio.PDB.Chain structure
    :return: The same chain, but with modified residue names.
    '''
    for r in chain:
        # rename rosetta-generated structures
        if r.resname == ' rA':
            r.resname = '  A'
        elif r.resname == ' rC':
            r.resname = '  C'
        elif r.resname == ' rG':
            r.resname = '  G'
        elif r.resname == ' rU':
            r.resname = '  U'

        # rename modified residues
        if r.id[0] == 'H_PSU':
            r.resname = '  U'
            change_residue_id(r, (' ', r.id[1], r.id[2]))
        elif r.id[0] == 'H_5MU':
            r.resname = '  U'
            change_residue_id(r, (' ', r.id[1], r.id[2]))
        elif r.id[0] == 'H_5MC':
            r.resname = '  C'
            change_residue_id(r, (' ', r.id[1], r.id[2]))
        elif r.id[0] == 'H_1MG':
            r.resname = '  G'
            change_residue_id(r, (' ', r.id[1], r.id[2]))
        elif r.id[0] == 'H_H2U':
            r.resname = '  U'
            change_residue_id(r, (' ', r.id[1], r.id[2]))

        # treat deoxyuridine as uridine
        # TODO: is this acceptable?
        if r.resname == ' DU':
            r.resname = '  U'

    return chain

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

        del a.parent.child_dict[oldid]
        a.parent.child_dict[a.id] = a

    return chain

def remove_hetatm(chain):
    '''
    Remove all the hetatms in the chain.

    :param chain: A Bio.PDB.Chain
    :return: The same chain, but missing all hetatms
    '''
    to_detach = []
    for r in chain:
        if len(r.id[0].strip()) > 0:
            to_detach += [r.id]

    for r in to_detach:
        chain.detach_child(r)

    return chain

def load_structure(pdb_filename):
    '''
    Load a Bio.PDB.Structure object and return the largest chain.
    This chain will be modified so that all hetatms are removed, modified
    residues will be renamed to regular residues, etc...
    '''
    chain = get_biggest_chain(pdb_filename) 
    return clean_chain(chain)

def clean_chain(chain):
    """
    Clean a pdb chain for further use with forgi.
    
    It will be modified so that all hetatms are removed, modified
    residues will be renamed to regular residues, residue ids will be positive integers, ...
    
    :param chaion: A Bio.PDB.Chain object
    :returns: A modified version of this chain
    """
    chain = rename_modified_ress(chain)
    chain = rename_rosetta_atoms(chain)
    chain = remove_hetatm(chain)
    # chain = renumber_chain(chain) #We cannot do this because we need to parse FR3D output!
    return chain
    
def interchain_contacts(struct):
    all_atoms = bpdb.Selection.unfold_entities(struct, 'A')

    ns = bpdb.NeighborSearch(all_atoms)
    pairs = ns.search_all(2.8)

    ic_pairs = []

    for (a1, a2) in pairs:
        if a1.parent.parent != a2.parent.parent:
            ic_pairs += [(a1,a2)]

    return ic_pairs

def is_rna(chain):
    '''
    Determine if a Bio.PDB.Chain structure corresponds to an RNA
    molecule.

    :param chain: A Bio.PDB.Chain molecule
    :return: True if it is an RNA molecule, False otherwise
    '''
    for res in chain:
        if res.resname.strip() in ['A', 'C', 'G', 'U']:
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
