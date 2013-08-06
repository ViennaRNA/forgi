import sys, warnings
import numpy as np
import Bio.PDB as bpdb
import corgy.utilities.rmsd as brmsd
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

backbone_atoms = ['P', 'O5*', 'C5*', 'C4*', 'C3*', 'O3*']
ring_atoms = ['C4*', 'C3*', 'C2*', 'C1*', 'O4*']

side_chain_atoms = dict()
side_chain_atoms['U'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
side_chain_atoms['C'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']

side_chain_atoms['A'] = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']
side_chain_atoms['G'] = ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']

all_rna_atoms = backbone_atoms + ring_atoms
for v in side_chain_atoms.values():
    all_rna_atoms += v
all_rna_atoms = set(all_rna_atoms)

interactions = [('P', 'O5*'),
                ('P', 'OP1'),
                ('P', 'O1P'),
                ('P', 'OP2'),
                ('P', 'O2P'),
                ('C2*', 'O2*'),
                       ('O5*', 'C5*'),
                       ('C5*', 'C4*'),
                       ('C4*', 'O4*'),
                       ('C4*', 'C3*'),
                       ('O4*', 'C1*'),
                       ('C3*', 'C2*'),
                       ('C3*', 'O3*'),
                       ('C2*', 'C1*'),
                       ('C1*', 'N1'),
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

def is_covalent(contact):
    '''
    Determine if a particular contact is covalent.

    @param contact: A pair of two Atom objects
    @return True if they are covalently bonded
            False otherwise
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

    #cud.pv('((r1.id[1], c1.name), (r2.id[1], c2.name))')
    #cud.pv('r2.id[1] - r1.id[1]')

    return False

def num_noncovalent_clashes(chain):
    '''
    Check if a chain has non-covalent clashes. Non-covalent clashes are found
    when two atoms that aren't covalently linked are within 1.8 A of each other.

    @param chain: The chain to evaluate
    @param return: The number of non-covalent clashes.
    '''
    all_atoms = bpdb.Selection.unfold_entities(chain, 'A')
    ns = bpdb.NeighborSearch(all_atoms)

    contacts = ns.search_all(1.9)

    return len([c for c in contacts if not is_covalent(c)])

def noncovalent_distances(chain, cutoff=0.3):
    '''
    Print out the distances between all non-covalently bonded atoms
    which are closer than cutoff to each other.

    @param chain: The Bio.PDB chain.
    @param cutoff: The maximum distance
    '''
    all_atoms = bpdb.Selection.unfold_entities(chain, 'A')
    ns = bpdb.NeighborSearch(all_atoms)

    contacts = ns.search_all(cutoff)

    return [cuv.magnitude(c[1] - c[0]) for c in contacts if not is_covalent(c)]

def pdb_rmsd(c1, c2, sidechains=False, superimpose=True, apply_sup=False):
    '''
    Calculate the all-atom rmsd between two RNA chains.

    @param c1: A Bio.PDB.Chain
    @param c2: Another Bio.PDB.Chain
    @return: The rmsd between the locations of all the atoms in the chains.
    '''

    a_5_names = ['P', 'O5*', 'C5*', 'C4*', 'O4*', 'O2*']
    a_3_names = ['C1*', 'C2*', 'C3*', 'O3*']

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

    if len(c1.get_list()) != len(c2.get_list()):
        cud.pv('len(c1.get_list())')
        cud.pv('len(c2.get_list())')
        print >>sys.stderr, "Chains of different length"
        raise Exception("Chains of different length.")

    c1_list = c1.get_list()
    c2_list = c2.get_list()

    c1_list.sort(key=lambda x: x.id[1])
    c2_list.sort(key=lambda x: x.id[1])

    for r1,r2 in zip(c1_list, c2_list):
        if sidechains:
            anames = backbone_atoms + a_names[c1[i].resname.strip()]
        else:
            anames = backbone_atoms
        #anames = a_5_names + a_3_names

        for a in anames:
            if a in r1 and a in r2:
                all_atoms1 += [r1[a]]
                all_atoms2 += [r2[a]]
                '''
                except KeyError:
                    print >>sys.stderr, "Residue number %d is missing an atom, continuing with the rest." % (i)
                    continue
                '''

        '''
        if len(atoms1) != len(atoms2):
            print >>sys.stderr, "Number of atoms differs in the two chains."
            raise Exception("Missing atoms.")
        '''
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

        #cud.pv('cuv.vector_set_rmsd(crvs1, crvs2)')
        #return (len(all_atoms1), brmsd.rmsd(crvs1, crvs2), None)
        return (len(all_atoms1), cuv.vector_set_rmsd(crvs1, crvs2), None)

def get_first_chain(filename):
    '''
    Load a PDB file using the Bio.PDB module and return the first chain.

    @param filename: The path to the pdb file
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = bpdb.PDBParser().get_structure('t', filename)
        return list(s.get_chains())[0]

def pdb_file_rmsd(fn1, fn2):
    '''
    Calculate the RMSD of all the atoms in two pdb structures.

    @param fn1: The first filename.
    @param fn2: The second filename.
    @return: The rmsd between the two structures.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        s1= bpdb.PDBParser().get_structure('t', fn1)
        s2= bpdb.PDBParser().get_structure('t', fn2)

    c1 = list(s1.get_chains())[0]
    c2 = list(s2.get_chains())[0]

    rmsd = pdb_rmsd(c1, c2)

    return rmsd

def get_biggest_chain(in_filename, out_filename):
    '''
    Load the PDB file located at filename, select the longest
    chain, and save the file to out_filename.

    @param in_filename: The location of the original file.
    @param out_filename: The location of the file containing just
                         the largest chain.
    '''

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = bpdb.PDBParser().get_structure('temp', in_filename)

    chains = list(s.get_chains())
    biggest = 0
    biggest_len = 0

    for i in range(len(chains)):
        c = chains[i]
        res_list = list(c.get_list())
        
        #print >> sys.stderr, res_list[0].resname
        rna = False

        # Only count RNA residues
        num_residues = 0
        for res in res_list:
            if (res.resname.strip() == 'A' or
                res.resname.strip() == 'C' or
                res.resname.strip() == 'G' or
                res.resname.strip() == 'U'):
                num_residues += 1
                break

        if num_residues > biggest_len:
            biggest = i
            biggest_len = num_residues

    orig_chain = chains[biggest]

    # only output the chain we deemed biggest
    class FirstChainSelect(bpdb.Select):
        def accept_chain(self, chain):
            if chain.get_id() == orig_chain.get_id():
                return 1
            else:
                return 0

    io=bpdb.PDBIO()
    io.set_structure(s)

    io.save(out_filename, FirstChainSelect())
