from collections.abs import Mapping

class Ensemble(object):
    def __init__(self, cgs, sort_key = None):
        """
        An Ensemble is a sequence of Coarse grained RNAs, all of which must correspond 
                    to the same RNA 2D structure.

        :param cgs: An iterable of coarse grain RNAs or a mapping key->cg.
                    The Ensemble instance may keep references to the cg objects in cgs.
                    It modifies them by centering them on their coords's centroid.
        :param sort_key: Optional. A function that takes a cg (if cgs is a sequence)
                         or a key (if cgs is a mapping)
        """
        self._cgs=[]
        self._cg_lookup={}
        self._cg_rev_lookup={}
        if isinstance(cgs, Mapping):
            for key in sorted(cgs.keys(), key=sort_key):
                self._add_to_cg_list(cgs[key], key)
        else:
            if sort_key is not None:
                ids = list(sorted(range(len(cgs))), key=lambda x: sort_key(cgs[x]))
            else:
                ids = list(range(len(cgs)))
            for key in ids:
                self._add_to_cg_list(cgs[key], key)
        self._rmsd = np.ones((len(self._cgs), len(self._cgs)))*np.nan
    def _add_to_cg_list(self, cg, key):
        #In the future, we might have to check the rmsd in addition, 
        # if the MCMC will allow rotation/translation
        if not self._cgs or cg.coords != self._cgs[-1].coords or cg.twists != self._cgs[-1].twists:
            self._cgs.append(cg)
            self._cgs[-1].coords.center()
        #Else: only store lookup entry pointing to previous cg
        self._cg_lookup[key]=len(self._cgs)-1
        self._cg_rev_lookup[len(self._cgs)-1]=key
            
    def _calculate_rmsd_matrix(self):
        for i in range(len(self._cgs)):
            self._rmsd[i,i]=0.
        for i,j in it.combinations(range(len(self._cgs,2))):
            self._rmsd[i,j]=self._rmsd[j,i] = cg1.coords.rmsd_to(cg2.coords)
