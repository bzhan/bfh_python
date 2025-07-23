"""Unit test for extendbyid.py"""

from ..extendbyid import *
from ..dastructure import SimpleDAGenerator
from ..ddstructure import identityDD
from ..localpmc import LocalIdempotent
from ..pmc import Idempotent
from ..pmc import linearPMC, splitPMC
import unittest

class ExtendedDAStructureTest(unittest.TestCase):
    def setUp(self):
        self.pmc = splitPMC(2)
        self.splitting = PMCSplitting(self.pmc, [(0, 2)])
        local_pmc = self.splitting.local_pmc
        outer_pmc = self.splitting.outer_pmc

        # Construct the local generators
        local_da = LocalDAStructure(
            F2, local_pmc.getAlgebra(), local_pmc.getAlgebra())
        idems = {"x1" : ([], []),
                 "x2" : ([0], [0]),
                 "x3" : ([0], [1]),
                 "x4" : ([1], [1]),
                 "x5" : ([0, 1], [0, 1])}
        gens = {}
        for name, (l_idem, r_idem) in list(idems.items()):
            gens[name] = SimpleDAGenerator(
                local_da, LocalIdempotent(local_pmc, l_idem),
                LocalIdempotent(local_pmc, r_idem), name)
            local_da.addGenerator(gens[name])

        local_da.auto_u_map()

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
            local_da, self.splitting, self.splitting)
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
        for name, (l_idem, r_idem) in list(extended_idems.items()):
            for gen in mod_gens:
                if gen.idem1 == Idempotent(self.pmc, l_idem) and \
                   gen.idem2 == Idempotent(self.pmc, r_idem):
                    self.extended_gens[name] = gen

    def testGetGenerators(self):
        self.assertEqual(len(self.extended_da.getGenerators()), 8)
        self.assertEqual(len(self.extended_gens), 8)
        
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
            self.assertEqual(self.extended_da.delta(
                self.extended_gens[x], [self.pmc.sd(a) for a in algs_a]),
                              self.pmc.sd(alg_d) * self.extended_gens[y])

    def testDeltaPrefix(self):
        # Restriction of one of the testDelta cases
        self.assertTrue(self.extended_da.deltaPrefix(
            self.extended_gens["y4"], [self.pmc.sd([(1, 2), (4, 5)])]))
        # Haven't stored such a local arrow
        self.assertFalse(self.extended_da.deltaPrefix(
            self.extended_gens["y4"], [self.pmc.sd([(1, 3), (4, 5)])]))

class AntiBraidTest(unittest.TestCase):
    # A test of the local situation of anti-braid resolution, involving two
    # unpaired points. Note this is a simple case before making an isotopy. It
    # cannot be used in actual calculations because the DA action would be
    # infinitely large.
    def setUp(self):
        self.pmc = linearPMC(2)
        self.splitting = PMCSplitting(self.pmc, [(1, 4)])

        # Correspondence between points in pmc and local / outer pmc:
        #       pmc - 0      1 2 3 4      5 6 7
        # local_pmc -     0* 1 2 3 4 5*
        # outer_pmc - 0 1*             2* 3 4 5
        local_pmc = self.splitting.local_pmc
        outer_pmc = self.splitting.outer_pmc

        # Construct the local generators. In local_pmc, idempotent 0 is (1, 4),
        # idempotent 1 is (2,), idempotent 2 is (3,).
        local_da = LocalDAStructure(
            F2, local_pmc.getAlgebra(), local_pmc.getAlgebra(),
            single_idems1 = [1, 2], single_idems2 = [1, 2])
        idems = {"x1" : ([0], [1]),
                 "x2" : ([0], [2]),
                 "x3" : ([0, 2], [1, 2]),
                 "x4" : ([0, 1], [1, 2])}
        gens = {}
        for name, (l_idem, r_idem) in list(idems.items()):
            gens[name] = SimpleDAGenerator(
                local_da, LocalIdempotent(local_pmc, l_idem),
                LocalIdempotent(local_pmc, r_idem), name)
            local_da.addGenerator(gens[name])

        local_da.auto_u_map()

        for name_from, name_to, algs_a, alg_d in [
            # Some examples of local DA actions
            ("x4", "x3", [], [1, (2, 3)]),
            ("x1", "x2", [[(2, 3)]], [1]),
            ]:
            local_da.addDelta(gens[name_from], gens[name_to],
                              local_pmc.sd(alg_d),
                              [local_pmc.sd(alg_a) for alg_a in algs_a], 1)

        self.extended_da = ExtendedDAStructure(
            local_da, self.splitting, self.splitting)
        mod_gens = self.extended_da.getGenerators()
        # Set up short names for the extended generators. Here y_i is the unique
        # extension of x_i.
        extended_idems = {"y1" : ([1, 3], [0, 3]),
                          "y2" : ([1, 3], [2, 3]),
                          "y3" : ([1, 2], [0, 2]),
                          "y4" : ([0, 1], [0, 2])}
        self.extended_gens = {}
        for name, (l_idem, r_idem) in list(extended_idems.items()):
            for gen in mod_gens:
                if gen.idem1 == Idempotent(self.pmc, l_idem) and \
                   gen.idem2 == Idempotent(self.pmc, r_idem):
                    self.extended_gens[name] = gen

    def testGetGenerators(self):
        self.assertEqual(len(self.extended_da.getGenerators()), 4)
        self.assertEqual(len(self.extended_gens), 4)

    def testDelta(self):
        for x, algs_a, alg_d, y in [
            ("y4", [], [1, (2, 3)], "y3"),
            ("y1", [[5, (2, 3)]], [1, 5], "y2")]:
            self.assertEqual(self.extended_da.delta(
                    self.extended_gens[x], [self.pmc.sd(a) for a in algs_a]),
                              self.pmc.sd(alg_d) * self.extended_gens[y])

    def testDeltaPrefix(self):
        # TODO: add test cases
        pass

class IdentityDDLocalTest(unittest.TestCase):
    def testIdentityDDLocal(self):
        splitting = PMCSplitting(linearPMC(2), [(1, 4)])
        local_pmc = splitting.local_pmc
        id_local = identityDALocal(local_pmc)
        extended_id = ExtendedDAStructure(id_local, splitting, splitting)
        id_dd = extended_id.tensorDD(identityDD(linearPMC(2)))
        self.assertTrue(id_dd.compareDDStructures(identityDD(linearPMC(2))))

if __name__ == "__main__":
    unittest.main()
