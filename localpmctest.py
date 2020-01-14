"""Unit test for localpmc.py"""

from localpmc import *
from pmc import splitPMC
import unittest

class LocalPMCTest(unittest.TestCase):
    def testLocalPMC(self):
        # One piece: (0-1-2-3*), with 0 paired with 2. Appears in short
        # underslide at bottom of full PMC.
        pmc1 = LocalPMC(4, [(0, 2),(1,)], [3])
        self.assertEqual(pmc1.num_pair, 2)
        self.assertEqual(pmc1.otherp, [2, 1, 0, -1])
        self.assertEqual(pmc1.pairid, [0, 1, 0, -1])
        self.assertEqual(pmc1.pairs, [(0, 2),(1,)])

        # One piece: (0*-1-2-3), with 1 paired with 3. Appears in short
        # underslide at top of full PMC.
        pmc2 = LocalPMC(4, [(1, 3),(2,)], [0])
        self.assertEqual(pmc2.num_pair, 2)
        self.assertEqual(pmc2.otherp, [-1, 3, 2, 1])
        self.assertEqual(pmc2.pairid, [-1, 0, 1, 0])
        self.assertEqual(pmc2.pairs, [(1, 3),(2,)])

        # Two pieces: (0*-1-2-3*) (4*-5-6*), with 2 paired with 5. Appears in
        # general arcslides.
        pmc3 = LocalPMC(7, [(1,),(2, 5)], [0, 3, 4, 6])
        self.assertEqual(pmc3.num_pair, 2)
        self.assertEqual(pmc3.otherp, [-1, 1, 5, -1, -1, 2, -1])
        self.assertEqual(pmc3.pairid, [-1, 0, 1, -1, -1, 1, -1])
        self.assertEqual(pmc3.pairs, [(1,),(2, 5)])

    def testGetLocalStrandDiagrams(self):
        # Tests number of local strand diagrams.
        pmc1 = LocalPMC(4, [(0, 2),(1,)], [3])
        self.assertEqual(len(pmc1.getStrandDiagrams()), 17)
        pmc2 = LocalPMC(4, [(1, 3),(2,)], [0])
        self.assertEqual(len(pmc2.getStrandDiagrams()), 17)
        pmc3 = LocalPMC(5, [(1, 3),(2,)], [0, 4])
        self.assertEqual(len(pmc3.getStrandDiagrams()), 41)
        pmc4 = LocalPMC(7, [(1,),(2, 5)], [0, 3, 4, 6])
        self.assertEqual(len(pmc4.getStrandDiagrams()), 78)

class LocalStrandsTest(unittest.TestCase):
    def setUp(self):
        # One piece: (0*-1-2-3-4*), with 1 paired with 3. Appears in short
        # underslide in the middle of PMC.
        self.pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])

    def testConstructStrands(self):
        strands1 = LocalStrands(self.pmc1, [(0, 1)])
        self.assertEqual(strands1.multiplicity, [1, 0, 0, 0])
        strands2 = LocalStrands(self.pmc1, [(1, 4)])
        self.assertEqual(strands2.multiplicity, [0, 1, 1, 1])

    def testPropagateRight(self):
        strands1 = LocalStrands(self.pmc1, [(0, 1)])
        strands2 = LocalStrands(self.pmc1, [(2, 4)])
        test_cases = [
            # Introduce pair (1, 3).
            (strands1, [], [0]),
            # Does not interfere with the single point (2,) that already exists.
            (strands1, [1], [0, 1]),
            # However, returns None (not compatible) if pair (1, 3) already
            # exists.
            (strands1, [0], None),

            # Removes single point (2,)
            (strands2, [1], []),
            # Does not interfere with pair (1, 3)
            (strands2, [0, 1], [0]),
            # Returns none if single point (2,) does not appear in left_idem
            (strands2, [0], None),
        ]
        for strand, left_idem, right_idem in test_cases:
            if right_idem is None:
                self.assertEqual(strand.propagateRight(left_idem), None)
            else:
                self.assertEqual(strand.propagateRight(left_idem),
                                 LocalIdempotent(self.pmc1, right_idem))

class LocalStrandDiagramTest(unittest.TestCase):
    def testMultiply(self):
        # 0*-1-2-3-4*, with 1 and 3 paired
        pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])
        sd1 = pmc1.sd([(0, 1)])
        sd2 = pmc1.sd([(1, 2)])
        sd3 = pmc1.sd([2, (0, 1)])
        sd4 = pmc1.sd([(0, 2)])
        sd5 = pmc1.sd([(0, 1),(1, 2)])
        sd6 = pmc1.sd([(3, 4)])
        self.assertEqual(sd1 * sd2, 1*sd4)
        self.assertEqual(sd2 * sd1, 0)
        self.assertEqual(sd2 * sd3, 1*sd5)
        self.assertEqual(sd3 * sd4, 0)
        self.assertEqual(sd1 * sd6, 0)

        # 0-1-2-3*-4*-5-6-7, with 0 and 2, 5 and 7 paired.
        pmc2 = LocalPMC(8, [(0, 2),(1,),(5, 7),(6,)], [3, 4])
        sd1 = pmc2.sd([(0, 3)])
        sd2 = pmc2.sd([(4, 6)])
        sd3 = pmc2.sd([(0, 3),(4, 6)])
        self.assertEqual(sd1 * sd2, 1*sd3)

    def testDiff(self):
        # 0*-1-2-3-4*, with 1 and 3 paired
        pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])
        sd1 = pmc1.sd([(0, 4)])
        sd2 = pmc1.sd([1, (0, 4)])
        sd3 = pmc1.sd([2, (0, 4)])
        sd4 = pmc1.sd([(0, 1),(1, 4)])
        sd5 = pmc1.sd([(0, 2),(2, 4)])
        sd6 = pmc1.sd([(0, 3),(3, 4)])
        self.assertEqual(sd1.diff(), 0)
        self.assertEqual(sd2.diff(), 1*sd4 + 1*sd6)
        self.assertEqual(sd3.diff(), 1*sd5)

    def testAntiDiff(self):
        # 0*-1-2-3-4*, with 1 and 3 paired
        pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])

        sd1 = pmc1.sd([(0, 1),(1, 2)])
        self.assertEqual(len(sd1.antiDiff()), 1)
        sd2 = pmc1.sd([(0, 1),(1, 2),(2, 4)])
        self.assertEqual(len(sd2.antiDiff()), 2)

    def testFactor(self):
        # 0*-1-2-3-4*, with 1 and 3 paired
        pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])

        # (0->2)*idem, (0->1)*(1->2), idem*(0->2)
        sd1 = pmc1.sd([(0, 2)])
        self.assertEqual(len(sd1.factor()), 3)

        # (1,0->2)*idem, idem*(1,0->2)
        sd2 = pmc1.sd([1, (0, 2)])    
        self.assertEqual(len(sd2.factor()), 2)

class PMCSplittingTest(unittest.TestCase):
    def testRestrictPMC(self):
        pmc1 = splitPMC(1)
        pmc2 = splitPMC(2)
        self.assertEqual(PMCSplitting.restrictPMC(pmc1, [(0, 2)]),
                         (LocalPMC(4, [(0, 2), (1,)], [3]), {0:0, 1:1, 2:2}))
        self.assertEqual(PMCSplitting.restrictPMC(pmc1, [(1, 3)]),
                         (LocalPMC(4, [(1, 3), (2,)], [0]), {1:1, 2:2, 3:3}))
        self.assertEqual(PMCSplitting.restrictPMC(pmc2, [(5, 7)]),
                         (LocalPMC(4, [(1, 3), (2,)], [0]), {5:1, 6:2, 7:3}))
        # Restriction is (0*-1-2*) (3*-4-5*), with 2 and 4 paired
        self.assertEqual(PMCSplitting.restrictPMC(pmc2, [(4, 4), (6, 6)]),
                         (LocalPMC(6, [(1, 4)], [0, 2, 3, 5]), {4:1, 6:4}))

    def testComplementIntervals(self):
        self.assertEqual(
            PMCSplitting.complementIntervals(splitPMC(2), [(0, 2), (6, 6)]),
            [(3, 5), (7, 7)])

    def testJoin(self):
        pmc = splitPMC(2)

        # Local restriction is (0-1-2-3*), (4*-5-6*), where 0 and 2 are paired
        # (pair-id = 0), and 1, 5 have pair-id 1 and 2, respectively.
        # Outer restriction is (0*-1-2-3-4*), (5*-6), where 3 and 6 are paired
        # (pair-id = 2), and 1, 2 have pair-id 0 and 1, respectively.

        # Correspondence between points in pmc and local_pmc1/2:
        #       pmc - 0 1 2      3 4 5      6      7
        # local_pmc - 0 1 2 3*           4* 5 6*
        # outer_pmc -         0* 1 2 3 4*       5* 6
        splitting = PMCSplitting(pmc, [(0, 2), (6, 6)])

        local_pmc = splitting.local_pmc
        outer_pmc = splitting.outer_pmc

        # Format is as follows:
        # - input to local_pmc.sd
        # - input to outer_pmc.sd
        # - input to pmc.sd. None indicates join should fail.
        test_data = [
            # One strand
            ([(0, 1)], [], [(0, 1)]),
            ([(0, 3)], [(0, 2)], [(0, 4)]),
            ([(0, 3),(4, 5)], [(0, 4)], [(0, 6)]),
            ([(0, 3),(4, 6)], [(0, 4),(5, 6)], [(0, 7)]),
            ([(1, 3)], [(0, 3)], [(1, 5)]),
            ([(1, 3)], [(0, 3)], [(1, 5)]),
            # Mismatch cases
            ([(4, 5)], [], None),
            ([(0, 2)], [(0, 3)], None),
            ([(0, 3),(4, 5)], [], None),
            ([(0, 3)], [(5, 6)], None),
            ([1], [(1, 4)], None),
            ([(0, 3), 5], [(2, 3)], None),
            # Several strands
            ([1,(2, 3)], [(0, 1),(1, 2)], [(2, 3),(3, 4)]),
            ([(4, 5),(5, 6)], [2,(3, 4),(5, 6)], [(5, 6),(6, 7)]),
            # Double horizontal
            ([0, 1, 5], [1, 2], [0, 1, 4]),
            ([0, 1], [1, 2], [0, 1, 4]),
            # Strands starting and ending at paired points
            ([(0, 1),(2, 3)], [(0, 2)], None),
            ([(0, 1),(1, 3)], [(0, 1)], None)
        ]
        for sd1, sd2, sd_total in test_data:
            joined_sd = splitting.joinStrandDiagram(local_pmc.sd(sd1),
                                                    outer_pmc.sd(sd2))
            if sd_total == None:
                self.assertEqual(joined_sd, None)
            else:
                self.assertEqual(joined_sd, pmc.sd(sd_total))

    def testRestrictStrandDiagram(self):
        # First test case: one local PMC at the boundary
        # local_pmc1 is 0-1-2-3*, where 0 and 2 are paired.
        pmc = splitPMC(2)
        splitting1 = PMCSplitting(pmc, [(0, 2)])
        local_pmc1 = splitting1.local_pmc
        for full_sd, local_sd in [
            ([(0, 1)], [(0, 1)]),
            ([(0, 5)], [(0, 3)]),
            ([3, (0, 5)], [1, (0, 3)]),
            ([3, (4, 5)], [1]),
            ([(2, 5)], [(2, 3)]),
            ([0, (1, 3)], [0, (1, 3)]),
            ]:
            self.assertEqual(
                splitting1.restrictStrandDiagramLocal(pmc.sd(full_sd)),
                local_pmc1.sd(local_sd))

        # Second test case: one local PMC's in the middle
        # local_pmc2 is 0*-1-2-3-4*, where no points are paired.
        # so point 1 has pair-id 0, point 2 has pair-id 1, and point3 has
        # pair-id 2.
        splitting2 = PMCSplitting(pmc, [(3, 5)])
        local_pmc2 = splitting2.local_pmc
        for full_sd, local_sd in [
            ([(0, 1)], []),
            ([(0, 6)], [(0, 4)]),
            ([(1, 6)], [(0, 4)]),
            ([1], [1]),
            ([(0, 4)], [(0, 2)]),
            ([(4, 6)], [(2, 4)]),
            ]:
            self.assertEqual(
                splitting2.restrictStrandDiagramLocal(pmc.sd(full_sd)),
                local_pmc2.sd(local_sd))

        # Third test case: local PMC has two separated intervals.
        # local_pmc3 is (0-1-2-3*), (4*-5-6*), where 0 and 2 are paired
        # (pair-id = 0), 1 and 5 are single (pair-id = 1 and 2).
        splitting3 = PMCSplitting(pmc, [(0, 2), (6, 6)])
        local_pmc3 = splitting3.local_pmc
        for full_sd, local_sd in [
            ([(0, 7)], [(0, 3), (4, 6)]),
            ([(0, 6)], [(0, 3), (4, 5)]),
            ([(0, 5)], [(0, 3)]),
            ([(2, 6)], [(2, 3), (4, 5)]),
            ([(3, 6)], [(4, 5)]),
            ([0, 4], [0, 5]),
            ]:
            self.assertEqual(
                splitting3.restrictStrandDiagramLocal(pmc.sd(full_sd)),
                local_pmc3.sd(local_sd))

if __name__ == "__main__":
    unittest.main()
