from forgi.graph.transform_graphs import _GCDummy, BGTransformer


class CGTransformer(BGTransformer):
    def condensed(self):
        """
        Return a condensed copy of the CoarseGrainRNA.

        See superclass docstring for definition of condensed.

        In the case of CoarseGrainRNA objects, the 3D coordinates are not
        affected by condensing the RNA. This means that stems which has multiple
        base-pairs have the same length as 1-bp stems in the condensed version.
        """
        new_cg = super(CGTransformer, self).condended()
        new_cg.coords = self.bg.coords
        new_cg.twists = self.bg.twists
        new_cg.sampled = self.bg.sampled
        new_cg.longrange = self.bg.longrange
        new_cg.chains = self.bg.chains
        return new_cg
