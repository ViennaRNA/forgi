"""
A module for dealing with modified residues.

For now the only functionality is to find out the corresponding 
unmodified residue.
"""

try:
    from urllib.request import urlopen #python3
except ImportError:
    from urllib import urlopen #python2
import urllib
import re
try: 
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup

import logging
import sys
from pprint import pprint
log = logging.getLogger(__name__)
def query_PDBeChem(three_letter_code):
    """
    Query the PDBeChem for the corresponding 1-letter code
    
    Unfortunately, this information does not seem to be provided in any
    machine-readable format. We have to parse the html.
    """
    if len(three_letter_code)!=3:
        raise ValueError("Illegal 3-letter code")
    if not re.match("^[A-Z0-9][A-Z0-9][A-Z0-9]$", three_letter_code):
        raise ValueError("Illegal 3-letter code")
    html = urlopen("http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{}".format(three_letter_code))
    parsed_html = BeautifulSoup(html, "lxml")
    content_area = parsed_html.body.find("td", attrs = {"class":"contentsarea"})
    #log.info(content_area)
    olc = content_area.find_all("h3")#, string=re.compile("One-letter code"))
    for h3 in olc:
        for stri in h3.stripped_strings:
            if stri.startswith("One-letter code"):
                for sibling in h3.parent.next_siblings:
                    if sibling.name == "td":
                        return sibling.string.strip()

    return None
    
def dict_from_pdbs(filenames):
    #Just look at the SEQRES headers.
    out_dict = { "A":"A", "U":"U", "G":"G", "C":"C"}
    for filename in filenames:
        with open(filename) as f:
            log.info("Opened file %s", filename)
            for line in f:
                if line.startswith("SEQRES"):
                    codes = line[20:].split()
                    for code in codes:
                        code = code.strip()
                        if code not in out_dict:
                            try:
                                one_letter_c = query_PDBeChem(code)
                                log.info("Found code for %r: %r", code, one_letter_c)
                                out_dict[code] = str(one_letter_c)
                            except (ValueError, LookupError) as e:
                                log.warning("3-letter code '%s' not found: %s", code, e)
    return out_dict
                                
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("# Created by `{}`".format(" ".join(sys.argv)))
    pprint(dict_from_pdbs(sys.argv[1:]))
