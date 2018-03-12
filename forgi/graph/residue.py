from builtins import str
import collections
import logging

log = logging.getLogger(__name__)

RESID = collections.namedtuple("complete_resid", ["chain", "resid"])


def resid_to_str(resid):
    if resid.chain is not None:
        out="{}:{}".format(resid.chain, resid.resid[1])
    else:
        out=str(resid.resid[1])
    if resid.resid[2]!=" ":
        out+=".{}".format(resid.resid[2])
    return out

def resid_from_str(resstr):
    resstr = str(resstr) # Make sure we use future's newstring on python2
    if ":" in resstr:
        chain, resid = resstr.split(":")
    else:
        resid=resstr
        log.debug("No chain given in string {!r}".format(resstr))
        chain=str('A')
    idparts=resid.split(".")
    if len(idparts)==1:
        idparts.append(" ")
    return RESID(chain, (' ', int(idparts[0]), idparts[1]))

RESID.__repr__ = lambda x: "res"+resid_to_str(x)
