"""
A module for dealing with modified residues.

For now the only functionality is to find out the corresponding
unmodified residue.
"""
from __future__ import print_function

from builtins import str
from builtins import object
try:
    from urllib.request import urlopen  # python3
except ImportError:
    from urllib import urlopen  # python2
import re
import logging
import warnings
import sys
from pprint import pprint


from ...graph import residue as fgr

#from .modified_res_lookup import RESIDUE_DICT

log = logging.getLogger(__name__)


def _parse_table_row(tr):
    tds = [td for td in tr.contents if td.name == "td"]
    if len(tds) == 2:
        key = tds[0].h3.text.strip()
        value = tds[1].text.strip()
        return key, value
    else:
        return None


def _html_to_info_dict(html):
    try:
        from bs4 import BeautifulSoup
    except ImportError:
        warnings.warn("Could not query PDBeChem, because BeautifulSoup is not installed"
                      "(install with `pip install forgi[pdbechem]` or `pip install forgi[all]`)")
        return {}
    with warnings.catch_warnings():
        # We don't care what xml parser beautiful soup uses
        warnings.simplefilter("ignore")
        parsed_html = BeautifulSoup(html)
    for h1 in parsed_html.body.find_all("h1"):
        if h1.text == "wwPDB Information":
            table = h1.parent.parent.parent  # (h1_>td->tr->tbody)
            rows = table.find_all("tr")
            info = {}
            for row in rows:
                data = _parse_table_row(row)
                if data:
                    info[data[0]] = data[1]
            return info
    return {}


def query_PDBeChem(three_letter_code):
    """
    Query the PDBeChem and return a dictionary with key-value pairs.

    Unfortunately, this information does not seem to be provided in any
    machine-readable format. We have to parse the html.
    """
    if len(three_letter_code) > 3:
        raise ValueError("Illegal 3-letter code")
    if not re.match("^[A-Z0-9][A-Z0-9][A-Z0-9]$", three_letter_code):
        raise ValueError("Illegal 3-letter code")
    html = urlopen(
        "http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{}".format(three_letter_code))
    residue_info = _html_to_info_dict(html)
    return residue_info


def change_residue_id(residue, new_id):
    chain = residue.parent
    if chain is not None and new_id in chain:
        raise ValueError("Cannot change id {old} to {new}. {new} already exists".format(
            old=residue.id, new=new_id))
    old_id = residue.id
    residue.id = new_id
    if chain is not None:
        try:
            del chain.child_dict[old_id]
        except KeyError:
            pass  # New version of biopython. Id is a property
        else:
            chain.child_dict[new_id] = residue


def to_4_letter_alphabeth(chain, query_db=False):
    '''
    Rename the modified residues so that they have the same
    names as unmodified residues.

    If a residue's name is unknown, query PDBeChem.

    If the residue is not a ribonucleotide, remove it from the chain.

    If it is a modified ribonucleotde, replace it by the unmodified "standard parent"


    :param chain: A Bio.PDB.Chain structure
    :param query_db: If true, query PDBeChem whenever a 3-letter code is unknown.
    :return: The same chain, but with only AUGC residues.
    '''
    import forgi.threedee.utilities.pdb as ftup
    modifications = {}
    i = 0
    while i < len(chain):
        r = chain.child_list[i]
        # Non-standard residues
        if r.resname.strip() == "I":
            # "I" has no standart parent (AUGC) to replace it with.
            warnings.warn("Inosine will be replaced by G")
            r.resname = "G"
            modifications[fgr.RESID(chain=chain.id, resid=(
                " ", r.id[1], r.id[2]))] = r.resname.strip()
        else:
            if r.resname.strip() not in ftup.RNA_RESIDUES:
                res_info = ModifiedResidueLookup(query_db)[r.resname]
                if not res_info:
                    # Unknown code. Remove residue
                    log.info(
                        "Detaching %s, because we do not understand the code %s.", r, r.resname)
                    chain.detach_child(r.id)
                    continue  # Continue with same i (now different residue)
                if res_info["Polymer type"] != "Ribonucleotide":
                    # DNA/ Amino Acid. Remove residue
                    chain.detach_child(r.id)
                    continue  # Continue with same i (now different residue)
                if res_info["Type description"] == "NON-POLYMER":
                    # e.g. GTP in 3DIR
                    log.warning(
                        "Detaching %s, because %s is classified as non-polymeric.", r, r.resname)
                    chain.detach_child(r.id)
                    continue  # Continue with same i (now different residue)
                if res_info["Standard parent"] not in ["A", "U", "G", "C"]:
                    log.warning("Detaching %s, because %s has standard parent '%s'.",
                                r, r.resname, res_info["Standard parent"])
                    chain.detach_child(r.id)
                    continue  # Continue with same i (now different residue)

                modifications[fgr.RESID(chain=chain.id, resid=(
                    " ", r.id[1], r.id[2]))] = r.resname.strip()
                r.resname = res_info["Standard parent"]

        # rename modified residues
        if r.id[0].strip():
            log.debug("Modified residue: %s", r.id)
            # The following was added because of 1NYI.pdb. It includes a single G-residue denoted as HETATOM in chain A.
            # It forms a base-pair, but is probably not connected to the backbone.
            # Instead of removing it/ in addition to removing it,
            # introducing cutpoint whereever the PDB Chain
            # contains gaps would be another option (TODO).
            # A plain AUGC as a HETATOM means it is a ligand.
            # 1Y26 has H_ADE instead of H_A
            if r.id[0].strip() in ["H_A", "H_U", "H_G", "H_C", "H_  G", "H_  C", "H_  A", "H_  U", "H_ADE", "H_CYT", "H_GUA", "H_URI"]:
                log.debug("Detaching %s, because it is a ligand", r.id)
                chain.detach_child(r.id)
                continue  # Continue with same i (now different residue)

            change_residue_id(r, (' ', r.id[1], r.id[2]))

        i += 1  # Go to next residue.
    return chain, modifications


class ModifiedResidueLookup(object):
    """
    Convenience wrapper to access modified_res_lookup.RESIDUE_DICT.
    If a key does not exist, query the database to add it."""

    def __init__(self, query_db=False):
        """
        :param query_db: If true, query PDBeChem whenever a 3-letter code is not understood.
        """
        from . import modified_res_lookup
        self._dict = modified_res_lookup.RESIDUE_DICT
        self._query_db= query_db
    def __getitem__(self, key):
        key = key.strip()
        if key not in self._dict:
            if not self._query_db:
                return None
            try:
                self._dict[key] = query_PDBeChem(key)
                log.debug("Successfully looked up %s", key)
            except:
                log.warning("Could not look-up modified residue key %s", key)
                self._dict[key] = None
        return self._dict[key]

    def clean_failed(self):
        for k, v in self._dict.items():
            if v is None:
                del self._dict[k]

    def __str__(self):
        return str(self._dict)


def dict_from_pdbs(filenames):
    # Just look at the SEQRES headers.
    ignore = set("AUGC")
    out_dict = {}
    for filename in filenames:
        with open(filename) as f:
            log.info("Opened file %s", filename)
            for line in f:
                if line.startswith("SEQRES"):
                    codes = line[20:].split()
                    for code in codes:
                        code = code.strip()
                        if code not in out_dict and code not in ignore:
                            try:
                                res_info = query_PDBeChem(code)
                                out_dict[code] = res_info
                            except (ValueError, LookupError) as e:
                                log.warning(
                                    "3-letter code '%s' not found: %s", code, e)
                                out_dict[code] = None
    return {k: v for k, v in out_dict.items()}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("# Created by `python {}`".format(" ".join(sys.argv)))
    print("RESIDUE_DICT = ", end="")
    d = dict_from_pdbs(sys.argv[1:])
    for ion in ["MG", "NA", "K", "CL"]:
        if ion not in d:
            d[ion] = None
    pprint({k: v for k, v in d.items() if v is not None})
