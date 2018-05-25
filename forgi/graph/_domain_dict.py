"""
Dictionaries that use element names as regular keys (i.e. "s0"),
but also support lists of elements (corresponding to domains) and
capital letter keys as aliases of such lists.


While the BulgeGraph holds these dictionaries, they contain a reference to it,
which is ONLY used for complex lookup of domains.
"""

import logging
try:
    from collections.abc import Mapping
except ImportError:
    from colelctions import Mapping

log=logging.getLogger(__name__)

class _DomainSupportingDictionary(Mapping):
    """
    A abstract base class for a dictionary-like class that supports
    domain names ("H0",...) and element lists (["s0", "i0", "s1"]) for lookup.
    """
    def __init__(self, d, cg):
        """
        :param d: The underlying dictionary.
        """
        if not isinstance(d, dict):
            raise TypeError("Expecting first argument's type to be a dictionary (sub)class")
        self._dict = d
        self._helices = None
        self._cg = cg

    def __getattr__(self, attr):
        """
        If no implementation is provided, delegate to the underlying dictionary.
        """
        return getattr(self._dict, attr)

    def __getitem__(self, key):
        log.debug("Getitem called for %s", key)
        try:
            return self._dict[key]
        except KeyError:
            log.debug("%s not in %s", key, self._dict)
            if isinstance(key, list):
                return self.lookup(key)
            elif key[0]=="H":
                # Helix lookup
                if self._helices is None:
                    self._set_helices()
                    return self.lookup(self._helices[key], _type="helix")
            else:
                log.debug("Raising KeyError")
                raise
    def __len__(self):
        return len(self._dict)
    def __iter__(self):
        return iter(self._dict)

    def helices(self):
        if self._helices is None:
            self._set_helices()
        return sorted(self._helices.keys())

    def _set_helices(self):
        domains = self._cg.get_domains()
        for i,helix in enumerate(self.domains["rods"]):
            self._helices["H{}".format(i)] = helix

    def lookup(self, elements, _type=None):
        raise NotImplementedError("Must be implemented by subclass")


class _DefineDictionary(_DomainSupportingDictionary):
    def lookup(self, elements, _type=None):
        """
        :param elements: A list of coarse grained element names
        :param _type: Internal, for speedup
        """
        if _type=="helix" or self._is_helix(elements):
            return self._helix_lookup(elements)
        else:
            raise KeyError("For defines, only helices can be used as keys.")
    def _is_helix(self, elements):
        """
        Return whether or not the elements are a helix or a continuous sub-helix.
        """
        if any(elem[0] not in "si" for elem in elemets):
            return False
        # Add one element to connected component and try to extend it.
        connected_component = [elements[0]]
        i=0
        while True:
            try:
                elem = connected_component[i]
            except IndexError:
                break
            neighbors = self._cg.edges[elem]
            for n in neighbors:
                if n in elements and n not in connected_component:
                    connected_component.append(n)
        return set(connected_component)==set(elements)

    def _helix_lookup(self, elements):
        """
        Lookup where the elements are a helix.

        Has undefined bahaviour if elements are not a helix.
        """
        stems = [e for e in elements if e[0]=="s"]
        defines = np.nans((4, len(stems)))
        for i in range(3):
            defines[i,:] = [self[s][0] for s in stems]
        return min(defines[0,:]),max(defines[1,:]), min(defines[2,:]),max(defines[3,:])

class _EdgeDictionary(_DomainSupportingDictionary):
    def lookup(self, elements, _type=None):
        # ignore type
        all_edges = { n for n in self[e] for e in elements}
        internal_edges = set(elements)
        return all_edges - internal_edges
