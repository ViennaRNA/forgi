from builtins import str
import collections
import logging

log = logging.getLogger(__name__)
Resid_base = collections.namedtuple("complete_resid", ["chain", "resid"])


class RESID(Resid_base):
    # Thanks to https://stackoverflow.com/a/42146452 for hints on extending named_tuples
    def __repr__(self):
        return "Res(\"" + resid_to_str(self) + "\")"

    def __new__(cls, chain, resid=None):
        if resid is None:
            r = resid_from_str(chain)
            return r
        else:
            resid = (resid[0], int(resid[1]), str(resid[2]))
            if chain is not None:
                chain = str(chain)
            return Resid_base.__new__(cls, chain, resid)


Res = RESID  # So the repr can be used to create a RESID


def resid_to_str(resid):
    if resid.chain is not None:
        out = "{}:{}".format(resid.chain, resid.resid[1])
    else:
        out = str(resid.resid[1])
    if resid.resid[2] != " ":
        out += ".{}".format(resid.resid[2])
    return out


def resid_from_str(resstr):
    resstr = str(resstr)  # Make sure we use future's newstring on python2
    if ":" in resstr:
        chain, resid = resstr.split(":")
    else:
        resid = resstr
        log.debug("No chain given in string {!r}".format(resstr))
        chain = None
    if chain is not None:
        chain = str(chain)
    idparts = resid.split(".")
    if len(idparts) == 1:
        idparts.append(" ")
    return RESID(chain, (' ', int(idparts[0]), str(idparts[1])))


def resid_from_biopython(residue):
    if residue.parent is not None:
        return RESID(residue.parent.id, residue.id)
    else:
        return RESID(None, residue.id)
