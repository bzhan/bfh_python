"""Unit test for localpmc.py"""

from localpmc import *
from pmc import splitPMC
import unittest

class LocalPMCTest(unittest.TestCase):
    def testLocalPMC(self):
        pmc1 = LocalPMC(4, [(0, 2),(1,)], [3])
        self.assertEqual(pmc1.num_pair, 2)
        self.assertEqual(pmc1.otherp, [2, 1, 0, -1])
        self.assertEqual(pmc1.pairid, [0, 1, 0, -1])
        self.assertEqual(pmc1.pairs, [(0, 2),(1,)])

        pmc2 = LocalPMC(4, [(1, 3),(2,)], [0])
        self.assertEqual(pmc2.num_pair, 2)
        self.assertEqual(pmc2.otherp, [-1, 3, 2, 1])
        self.assertEqual(pmc2.pairid, [-1, 0, 1, 0])
        self.assertEqual(pmc2.pairs, [(1, 3),(2,)])

        # Two pieces: (0*-1-2-3*) (4*-5-6*), with 2 paired with 5
        pmc3 = LocalPMC(7, [(1,),(2, 5)], [0, 3, 4, 6])
        self.assertEqual(pmc3.num_pair, 2)
        self.assertEqual(pmc3.otherp, [-1, 1, 5, -1, -1, 2, -1])
        self.assertEqual(pmc3.pairid, [-1, 0, 1, -1, -1, 1, -1])
        self.assertEqual(pmc3.pairs, [(1,),(2, 5)])

class LocalStrandsTest(unittest.TestCase):
    def setUp(self):
        # One piece: (0*-1-2-3-4*), with 1 paired with 3
        self.pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])

    def testConstructStrands(self):
        strands1 = LocalStrands(self.pmc1, [(0, 1)])
        self.assertEqual(strands1.multiplicity, [1, 0, 0, 0])
        strands2 = LocalStrands(self.pmc1, [(1, 4)])
        self.assertEqual(strands2.multiplicity, [0, 1, 1, 1])

    def testPropagateRight(self):
        strands1 = LocalStrands(self.pmc1, [(0, 1)])
        # Introduce pair (1, 3)
        self.assertEqual(strands1.propagateRight([]), [0])
        # Does not interfere with the single point (2,) that already exists
        self.assertEqual(strands1.propagateRight([1]), [0, 1])
        # However, returns None (not compatible) if pair (1, 3) already exists
        self.assertEqual(strands1.propagateRight([0]), None)

        strands2 = LocalStrands(self.pmc1, [(2, 4)])
        # Removes single point (2,)
        self.assertEqual(strands2.propagateRight([1]), [])
        # Does not interfere with pair (1, 3)
        self.assertEqual(strands2.propagateRight([0, 1]), [0])
        # Returns none if single point (2,) does not appear in left_idem
        self.assertEqual(strands2.propagateRight([0]), None)

class LocalStrandDiagramTest(unittest.TestCase):
    def testMultiply(self):
        # 0*-1-2-3-4*, with 1 and 3 paired
        pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])
        sd1 = LocalStrandDiagram(pmc1, [], [(0, 1)])
        sd2 = LocalStrandDiagram(pmc1, [0], [(1, 2)])
        sd3 = LocalStrandDiagram(pmc1, [1], [(0, 1)])
        sd4 = LocalStrandDiagram(pmc1, [], [(0, 2)])
        sd5 = LocalStrandDiagram(pmc1, [0], [(0, 1), (1, 2)])
        sd6 = LocalStrandDiagram(pmc1, [0], [(3, 4)])
        self.assertEqual(sd1.multiply(sd2), sd4)
        self.assertEqual(sd2.multiply(sd1), sd5)
        self.assertEqual(sd2.multiply(sd3), sd5)
        self.assertEqual(sd3.multiply(sd4), None)
        self.assertEqual(sd1.multiply(sd6), None)

        # 0-1-2-3*-4*-5-6-7, with 0 and 2, 5 and 7 paired.
        pmc2 = LocalPMC(8, [(0, 2),(1,),(5, 7),(6,)], [3, 4])
        sd1 = LocalStrandDiagram(pmc2, [0], [(0, 3)])
        sd2 = LocalStrandDiagram(pmc2, [], [(4, 6)])
        sd3 = LocalStrandDiagram(pmc2, [0], [(0, 3), (4, 6)])
        self.assertEqual(sd1.multiply(sd2), sd3)

    def testJoin(self):
        pmc = splitPMC(2)
        # Restriction is (0-1-2-3*), (4*-5-6*), where 0 and 2 are paired
        # (pair-id = 0), and 1, 5 have pair-id 1 and 2, respectively.
        local_pmc1, mapping1 = restrictPMC(pmc, [(0, 2), (6, 6)])
        # Restriction is (0*-1-2-3-4*), (5*-6), where 3 and 6 are paired
        # (pair-id = 2), and 1, 2 have pair-id 0 and 1, respectively.
        local_pmc2, mapping2 = restrictPMC(pmc, [(3, 5), (7, 7)])

        # Format is as follows:
        # - idempotent for local_pmc1
        # - strands for local_pmc1
        # - idempotent for local_pmc2
        # - strands for local_pmc2
        # - input to pmc.sd. None indicates join should fail.
        test_data = [
            # One strand
            ([0], [(0, 1)], [], [], [(0, 1)]),
            ([0], [(0, 3)], [], [(0, 2)], [(0, 4)]),
            ([0], [(0, 3),(4, 5)], [], [(0, 4)], [(0, 6)]),
            ([0], [(0, 3),(4, 6)], [], [(0, 4),(5, 6)], [(0, 7)]),
            ([1], [(1, 3)], [], [(0, 3)], [(1, 5)]),
            ([1], [(1, 3)], [0], [(0, 3)], [(1, 5)]),
            # Mismatch cases
            ([], [(4, 5)], [], [], None),
            ([0], [(0, 2)], [], [(0, 3)], None),
            ([0], [(0, 3),(4, 5)], [], [], None),
            ([0], [(0, 3)], [], [(5, 6)], None),
            ([1], [], [0], [(1, 4)], None),
            ([0, 2], [(0, 3)], [1], [(2, 3)], None),
            # Several strands
            ([0, 1], [(2, 3)], [0], [(0, 1),(1, 2)], [(2, 3),(3, 4)]),
            ([2], [(4, 5),(5, 6)], [1, 2], [(3, 4),(5, 6)], [(5, 6),(6, 7)]),
            # Double horizontal
            ([0, 1, 2], [], [0, 1], [], [0, 1, 4]),
            ([0, 1], [], [0, 1], [], [0, 1, 4]),
            # Strands starting and ending at paired points
            ([0], [(0, 1),(2, 3)], [], [(0, 2)], None),
            ([0, 1], [(0, 1),(1, 3)], [0], [(0, 1)], None)
        ]
        for idem1, sd1, idem2, sd2, sd_total in test_data:
            # print idem1, sd1, idem2, sd2, sd_total
            joined_sd = LocalStrandDiagram(local_pmc1, idem1, sd1).join(
                LocalStrandDiagram(local_pmc2, idem2, sd2),
                pmc, mapping1, mapping2)
            if sd_total == None:
                self.assertEqual(joined_sd, None)
            else:
                self.assertEqual(joined_sd, pmc.sd(sd_total))

class RestrictionsTest(unittest.TestCase):
    def testRestrictPMC(self):
        pmc1 = splitPMC(1)
        pmc2 = splitPMC(2)
        self.assertEqual(restrictPMC(pmc1, [(0, 2)]),
                         (LocalPMC(4, [(0, 2), (1,)], [3]), {0:0, 1:1, 2:2}))
        self.assertEqual(restrictPMC(pmc1, [(1, 3)]),
                         (LocalPMC(4, [(1, 3), (2,)], [0]), {1:1, 2:2, 3:3}))
        self.assertEqual(restrictPMC(pmc2, [(5, 7)]),
                         (LocalPMC(4, [(1, 3), (2,)], [0]), {5:1, 6:2, 7:3}))
        # Restriction is (0*-1-2*) (3*-4-5*), with 2 and 4 paired
        self.assertEqual(restrictPMC(pmc2, [(4, 4), (6, 6)]),
                         (LocalPMC(6, [(1, 4)], [0, 2, 3, 5]), {4:1, 6:4}))

    def testRestrictStrandDiagram(self):
        # First test case: one local PMC at the boundary
        # local_pmc1 is 0-1-2-3*, where 0 and 2 are paired.
        pmc1 = splitPMC(2)
        local_pmc1, mapping1 = restrictPMC(pmc1, [(0, 2)])
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 1)]),
                                               local_pmc1, mapping1),
                         LocalStrandDiagram(local_pmc1, [0], [(0, 1)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 5)]),
                                               local_pmc1, mapping1),
                         LocalStrandDiagram(local_pmc1, [0], [(0, 3)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([3, (0, 5)]),
                                               local_pmc1, mapping1),
                         LocalStrandDiagram(local_pmc1, [0, 1], [(0, 3)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([3, (4, 5)]),
                                               local_pmc1, mapping1),
                         LocalStrandDiagram(local_pmc1, [1], []))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(2, 5)]),
                                               local_pmc1, mapping1),
                         LocalStrandDiagram(local_pmc1, [0], [(2, 3)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([0, (1, 3)]),
                                               local_pmc1, mapping1),
                         LocalStrandDiagram(local_pmc1, [0, 1], [(1, 3)]))

        # Second test case: one local PMC's in the middle
        # local_pmc2 is 0*-1-2-3-4*, where no points are paired.
        # so point 1 has pair-id 0, point 2 has pair-id 1, and point3 has
        # pair-id 2.
        local_pmc2, mapping2 = restrictPMC(pmc1, [(3, 5)])
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 1)]),
                                               local_pmc2, mapping2),
                         LocalStrandDiagram(local_pmc2, [], []))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 6)]),
                                               local_pmc2, mapping2),
                         LocalStrandDiagram(local_pmc2, [], [(0, 4)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(1, 6)]),
                                               local_pmc2, mapping2),
                         LocalStrandDiagram(local_pmc2, [], [(0, 4)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([1]),
                                               local_pmc2, mapping2),
                         LocalStrandDiagram(local_pmc2, [0], []))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 4)]),
                                               local_pmc2, mapping2),
                         LocalStrandDiagram(local_pmc2, [], [(0, 2)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(4, 6)]),
                                               local_pmc2, mapping2),
                         LocalStrandDiagram(local_pmc2, [1], [(2, 4)]))

        # Third test case: local PMC has two separated intervals.
        # local_pmc3 is (0-1-2-3*), (4*-5-6*), where 0 and 2 are paired
        # (pair-id = 0), 1 and 5 are single (pair-id = 1 and 2).
        local_pmc3, mapping3 = restrictPMC(pmc1, [(0, 2), (6, 6)])
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 7)]),
                                               local_pmc3, mapping3),
                         LocalStrandDiagram(local_pmc3, [0], [(0, 3), (4, 6)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 6)]),
                                               local_pmc3, mapping3),
                         LocalStrandDiagram(local_pmc3, [0], [(0, 3), (4, 5)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(0, 5)]),
                                               local_pmc3, mapping3),
                         LocalStrandDiagram(local_pmc3, [0], [(0, 3)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(2, 6)]),
                                               local_pmc3, mapping3),
                         LocalStrandDiagram(local_pmc3, [0], [(2, 3), (4, 5)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([(3, 6)]),
                                               local_pmc3, mapping3),
                         LocalStrandDiagram(local_pmc3, [], [(4, 5)]))
        self.assertEqual(restrictStrandDiagram(pmc1, pmc1.sd([0, 4]),
                                               local_pmc3, mapping3),
                         LocalStrandDiagram(local_pmc3, [0, 2], []))

if __name__ == "__main__":
    unittest.main()
