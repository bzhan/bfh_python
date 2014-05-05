"""Unit test for extendbyid.py"""

from extendbyid import *
from dastructure import SimpleDAGenerator, SimpleDAStructure
from localpmc import LocalIdempotent
from localpmc import restrictPMC
from pmc import Idempotent
from pmc import splitPMC
import unittest

class ExtendedDAStructureTest(unittest.TestCase):
    def setUp(self):
        self.pmc = splitPMC(2)
        local_pmc, mapping = restrictPMC(self.pmc, [(0, 2)])
        outer_pmc, outer_mapping = restrictPMC(self.pmc, [(3, 7)])

        # Construct the local generators
        local_da = SimpleDAStructure(
            F2, local_pmc.getAlgebra(), local_pmc.getAlgebra())
        idems = {"x1" : ([], []),
                 "x2" : ([0], [0]),
                 "x3" : ([0], [1]),
                 "x4" : ([1], [1]),
                 "x5" : ([0, 1], [0, 1])}
        gens = {}
        for name, (l_idem, r_idem) in idems.items():
            gens[name] = SimpleDAGenerator(
                local_da, LocalIdempotent(local_pmc, l_idem),
                LocalIdempotent(local_pmc, r_idem), name)
            local_da.addGenerator(gens[name])

        for name_from, name_to, algs_a, alg_d in [
                # Example from _short_underslide_down_middle in arcslideda.py
                # (subtract one from there since the bottom is closed).
                ("x3", "x4", [], [(0, 1)]),
                ("x2", "x2", [[(0, 2)]], [(0, 2)]),
                ("x5", "x5", [[1, (0, 2)]], [1, (0, 2)]),
                ("x2", "x1", [[(2, 3)]], [(2, 3)]),
                ("x5", "x4", [[1, (2, 3)]], [1, (2, 3)]),
                ("x4", "x3", [[(1, 2)], [(0, 1)]], [(1, 2)])]:
            local_da.addDelta(gens[name_from], gens[name_to],
                              local_pmc.sd(alg_d),
                              [local_pmc.sd(alg_a) for alg_a in algs_a], 1)

        self.extended_da = ExtendedDAStructure(
            local_da, outer_pmc, self.pmc, self.pmc, mapping, outer_mapping,
            mapping, outer_mapping)
        mod_gens = self.extended_da.getGenerators()
        # Set up short names for the extended generators
        extended_idems = {"y1" : ([0, 1], [0, 1]),
                          "y2" : ([0, 2], [0, 2]),
                          "y3" : ([0, 3], [0, 3]),
                          "y4" : ([1, 2], [1, 2]),
                          "y5" : ([1, 3], [1, 3]),
                          "y6" : ([2, 3], [2, 3]),
                          "y7" : ([0, 2], [1, 2]),
                          "y8" : ([0, 3], [1, 3])}
        self.extended_gens = {}
        for name, (l_idem, r_idem) in extended_idems.items():
            for gen in mod_gens:
                if gen.idem1 == Idempotent(self.pmc, l_idem) and \
                   gen.idem2 == Idempotent(self.pmc, r_idem):
                    self.extended_gens[name] = gen

    def testGetGenerators(self):
        self.assertEquals(len(self.extended_da.getGenerators()), 8)
        self.assertEquals(len(self.extended_gens), 8)
        
    def testDelta(self):
        for x, algs_a, alg_d, y in [
                ("y7", [], [4, (0, 1)], "y4"),
                ("y8", [], [5, (0, 1)], "y5"),
                ("y2", [[4, (0, 2)]], [4, (0, 2)], "y2"),
                ("y2", [[4, (2, 3)]], [4, (2, 3)], "y4"),
                ("y1", [[1, (0, 2)]], [1, (0, 2)], "y1"),
                ("y1", [[1, (2, 4)]], [1, (2, 4)], "y4"),
                ("y4", [[(1, 2), (4, 5)], [(0, 1), (5, 6)]],
                 [(1, 2), (4, 6)], "y7")]:
            self.assertEquals(self.extended_da.delta(
                self.extended_gens[x], [self.pmc.sd(a) for a in algs_a]),
                              self.pmc.sd(alg_d) * self.extended_gens[y])

    def testDeltaPrefix(self):
        # Restriction of one of the testDelta cases
        self.assertTrue(self.extended_da.deltaPrefix(
            self.extended_gens["y4"], [self.pmc.sd([(1, 2), (4, 5)])]))
        # Haven't stored such a local arrow
        self.assertFalse(self.extended_da.deltaPrefix(
            self.extended_gens["y4"], [self.pmc.sd([(1, 3), (4, 5)])]))
        
if __name__ == "__main__":
    unittest.main()
