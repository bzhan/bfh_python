"""Unit test for pmc.py"""

from grading import averageRefinement
from pmc import *
import unittest

class PMCTest(unittest.TestCase):
    def testPMC(self):
        pmc1 = PMC([(0,2),(1,3)])
        for p, q in [(0,2),(2,0),(1,3),(3,1)]:
            self.assertEqual(pmc1.otherp[p], q)
        for p, i in [(0,0),(2,0),(1,1),(3,1)]:
            self.assertEqual(pmc1.pairid[p], i)
        for i, (p,q) in [(0, (0,2)), (1, (1,3))]:
            self.assertEqual(pmc1.pairs[i], (p,q))

    def testPMCEqual(self):
        pmc1 = PMC([(0,2),(1,3)])
        pmc2 = PMC([(1,3),(0,2)])
        pmc3 = PMC([(2,0),(3,1)])
        self.assertEqual(pmc1, pmc2)
        self.assertEqual(pmc1, pmc3)
        self.assertEqual(hash(pmc1), hash(pmc2))
        self.assertEqual(hash(pmc1), hash(pmc3))

    def testPMCOpp(self):
        pmc1 = PMC([(0,2),(1,3)])
        self.assertEqual(pmc1.opp(), pmc1)
        pmc2 = PMC([(0,2),(1,6),(3,5),(4,7)])
        pmc3 = PMC([(0,3),(1,6),(2,4),(5,7)])
        self.assertEqual(pmc2.opp(), pmc3)
        self.assertEqual(pmc3.opp(), pmc2)
        self.assertTrue(pmc2 != pmc3)

    def testSplitPMC(self):
        self.assertEqual(splitPMC(1), PMC([(0,2),(1,3)]))
        self.assertEqual(splitPMC(2), PMC([(0,2),(1,3),(4,6),(5,7)]))

    def testLinearPMC(self):
        self.assertEqual(linearPMC(1), PMC([(0,2),(1,3)]))
        self.assertEqual(linearPMC(2), PMC([(0,2),(1,4),(3,6),(5,7)]))

    def testAntipodalPMC(self):
        self.assertEqual(antipodalPMC(1), PMC([(0,2),(1,3)]))
        self.assertEqual(antipodalPMC(2), PMC([(0,4),(1,5),(2,6),(3,7)]))

    def testConnectSumPMC(self):
        self.assertEqual(connectSumPMC(linearPMC(2), splitPMC(1)),
                         PMC([(0,2),(1,4),(3,6),(5,7),(8,10),(9,11)]))

    def testUnconnectSumPMC(self):
        self.assertEqual(unconnectSumPMC(splitPMC(2), 1),
                         (splitPMC(1), splitPMC(1)))
        self.assertEqual(
            unconnectSumPMC(connectSumPMC(linearPMC(2), splitPMC(1)), 2),
            (linearPMC(2), splitPMC(1)))

    def testGetIdempotents(self):
        pmc = splitPMC(2)
        idems = pmc.getIdempotents()
        self.assertEqual(len(idems), 6)
        for idem in idems:
            self.assertTrue(isinstance(idem, Idempotent))
            self.assertEqual(idem.pmc, pmc)

    def testGetStrandDiagrams(self):
        sds = splitPMC(1).getStrandDiagrams()
        self.assertEqual(len(sds), 8)
        pmc = splitPMC(2)
        sds2 = pmc.getStrandDiagrams()
        for sd in sds2:
            self.assertTrue(isinstance(sd, StrandDiagram))
            self.assertEqual(sd.pmc, pmc)
            self.assertEqual(sd.mult_one, False)

    def testGetMultOneStrandDiagrams(self):
        pmc = splitPMC(1)
        m1_sds = pmc.getMultOneStrandDiagrams(idem_size = 2)
        for sd in m1_sds:
            self.assertTrue(isinstance(sd, StrandDiagram))
            self.assertEqual(sd.pmc, pmc)
            self.assertEqual(sd.mult_one, True)
            self.assertTrue(all([x<=1 for x in sd.multiplicity]))
        self.assertEqual(len(m1_sds), 5)

class IdempotentTest(unittest.TestCase):
    def testIdempotent(self):
        pmc = antipodalPMC(2)
        idem1 = Idempotent(pmc, [0,1])
        self.assertEqual(idem1.opp(), Idempotent(pmc, [2,3]))
        self.assertEqual(idem1.comp(), Idempotent(pmc, [2,3]))
        idem2 = Idempotent(pmc, [0,3])
        self.assertEqual(idem2.opp(), Idempotent(pmc, [0,3]))
        self.assertEqual(idem2.comp(), Idempotent(pmc, [1,2]))
        pmc2 = PMC([(0,2),(1,6),(3,5),(4,7)])
        # Note pmc2.opp() = PMC([(0,3),(1,6),(2,4),(5,7)])
        idem3 = Idempotent(pmc2, [0,1])
        self.assertEqual(idem3.opp(), Idempotent(pmc2.opp(), [1,3]))
        self.assertEqual(idem3.comp(), Idempotent(pmc2, [2,3]))

    def testUnconnectSumIdem(self):
        pmc1, pmc2, pmc3 = splitPMC(1), splitPMC(2), splitPMC(3)
        idem1 = Idempotent(pmc3, [0, 2, 4])
        self.assertEqual(unconnectSumIdem(idem1, 1),
                         (Idempotent(pmc1, [0]), Idempotent(pmc2, [0, 2])))
        self.assertEqual(unconnectSumIdem(idem1, 2),
                         (Idempotent(pmc2, [0, 2]), Idempotent(pmc1, [0])))
        idem2 = Idempotent(pmc3, [0, 1, 2])
        self.assertEqual(unconnectSumIdem(idem2, 1),
                         (Idempotent(pmc1, [0, 1]), Idempotent(pmc2, [0])))
        self.assertEqual(unconnectSumIdem(idem2, 2),
                         (Idempotent(pmc2, [0, 1, 2]), Idempotent(pmc1, [])))

class StrandsTest(unittest.TestCase):
    def testStrands(self):
        pmc = PMC([(0,2),(1,6),(3,5),(4,7)])
        sd1 = Strands(pmc, [(0,1)])
        self.assertEqual(sd1.opp(), Strands(pmc.opp(), [(6,7)]))

    def testIdemCompatible(self):
        pmc = splitPMC(2)
        sd1 = Strands(pmc, [(0,1)])
        self.assertTrue(sd1.leftCompatible(pmc.idem([0,4])))
        self.assertFalse(sd1.leftCompatible(pmc.idem([0,1])))
        self.assertFalse(sd1.leftCompatible(pmc.idem([1,4])))
        self.assertFalse(sd1.leftCompatible(pmc.idem([])))
        self.assertTrue(sd1.rightCompatible(pmc.idem([3])))
        self.assertFalse(sd1.rightCompatible(pmc.idem([2])))

    def testUnconnectSumStrands(self):
        pmc1, pmc2 = splitPMC(1), splitPMC(2)
        sd1 = Strands(pmc2, [(0,3),(4,5)])
        self.assertEqual(unconnectSumStrands(sd1, 1),
                         (Strands(pmc1, [(0,3)]), Strands(pmc1, [(0,1)])))
        sd2 = Strands(pmc2, [(5,7)])
        self.assertEqual(unconnectSumStrands(sd2, 1),
                         (Strands(pmc1, []), Strands(pmc1, [(1,3)])))

class StrandDiagramTest(unittest.TestCase):
    def setUp(self):
        self.pmc = splitPMC(2)
        self.sd1 = self.pmc.sd([(0,3),(1,4)], mult_one = False)
        self.sd2 = self.pmc.sd([1,(0,4)], mult_one = True)
        self.sd3 = self.pmc.sd([(0,1),(1,2)], mult_one = False)
        self.sd4 = self.pmc.sd([(0,4),(1,3)], mult_one = False)

    def testStrandDiagramInit(self):
        self.assertEqual(self.sd1.multiplicity, [1,2,2,1,0,0,0])
        self.assertEqual(self.sd1.parent,
                         self.pmc.getAlgebra(mult_one = False))
        self.assertEqual(self.sd1.double_hor, ())
        self.assertEqual(self.sd2.multiplicity, [1,1,1,1,0,0,0])
        self.assertEqual(self.sd2.parent, self.pmc.getAlgebra(mult_one = True))
        self.assertEqual(self.sd2.double_hor, (1,))

    def testStrandDiagramOpp(self):
        self.assertEqual(self.sd1.opp(),
                         self.pmc.opp().sd([(4,7),(3,6)], mult_one = False))
        self.assertEqual(self.sd2.opp(),
                         self.pmc.opp().sd([6,(3,7)], mult_one = True))

    def testStrandDiagramEqual(self):
        pmc2 = splitPMC(2)
        pmc2alg = pmc2.getAlgebra(mult_one = False)
        sd10 = StrandDiagram(pmc2alg, left_idem = [1,0],
                             strands = [(1,4),(0,3)])
        self.assertEqual(self.sd1, sd10)
        self.assertTrue(self.sd1 != self.sd4)
        self.assertEqual(hash(self.sd1), hash(sd10))

    def testPropagateRight(self):
        self.assertEqual(self.sd1.right_idem, self.pmc.idem([1,4]))
        self.assertEqual(self.sd2.right_idem, self.pmc.idem([1,4]))
        self.assertEqual(self.sd3.right_idem, self.pmc.idem([0,1]))

    def testNumCrossing(self):
        self.assertEqual(self.sd1.numCrossing(), 0)
        self.assertEqual(self.sd3.numCrossing(), 0)
        self.assertEqual(self.sd4.numCrossing(), 1)

    def testMaslov(self):
        self.assertEqual(self.sd1.maslov(), Fraction(-2))
        self.assertEqual(self.sd2.maslov(), Fraction(-1,2))
        self.assertEqual(self.sd3.maslov(), Fraction(-3,2))
        self.assertEqual(self.sd4.maslov(), Fraction(-1))

    def testGrading(self):
        self.assertEqual(self.sd1.getBigGrading(),
                         self.pmc.big_gr(-2, [1,2,2,1,0,0,0]))

    def testDiff(self):
        self.assertEqual(self.sd1.diff(), 0)
        self.assertEqual(self.sd4.diff(), 1*self.sd1)
        sd2d1 = self.pmc.sd([(0,1),(1,4)], True)
        sd2d2 = self.pmc.sd([(0,3),(3,4)], True)
        self.assertEqual(self.sd2.diff(), 1*sd2d1 + 1*sd2d2)

        # Test for double crossing
        pmc2 = antipodalPMC(2)
        sd10 = pmc2.sd([(0,7),(1,6),(2,5)], False)
        sd10d1 = pmc2.sd([(0,6),(1,7),(2,5)], False)
        sd10d2 = pmc2.sd([(0,7),(1,5),(2,6)], False)
        self.assertEqual(sd10.diff(), 1*sd10d1 + 1*sd10d2)

    def testMultiply(self):
        # TODO: improve structure of this test
        sd1 = self.pmc.sd([4,(0,1)], False)
        sd2 = self.pmc.sd([4,(1,2)], False)
        sd12 = self.pmc.sd([4,(0,2)], False)
        self.assertEqual(sd1 * sd2, 1*sd12)
        self.assertEqual(sd2 * sd1, 0)

        sd3 = self.pmc.sd([1,(4,5)], False)
        sd4 = self.pmc.sd([5,(3,4)], False)
        sd34 = self.pmc.sd([(3,4),(4,5)], False)
        self.assertEqual(sd3 * sd4, 1*sd34)
        self.assertEqual(sd4 * sd3, 0)

        # Test for double crossing
        sd5 = self.pmc.sd([0,(1,3)], False)
        sd6 = self.pmc.sd([1,(2,4)], False)
        self.assertEqual(sd5 * sd6, 0)

        sd7 = self.pmc.sd([0,(1,3)], False)
        sd8 = self.pmc.sd([3,(0,2)], False)
        sd9 = self.pmc.sd([(0,2),(1,3)], False)
        self.assertEqual(sd7 * sd8, 1*sd9)
        self.assertEqual(sd8 * sd7, 0)

        sd10 = self.pmc.sd([0,(1,3)], True)
        sd11 = self.pmc.sd([3,(0,2)], True)
        self.assertEqual(sd10 * sd11, 0)

    def testAntiDiff(self):
        self.assertEqual(self.sd1.antiDiff(), 1*self.sd4)
        self.assertEqual(self.sd2.antiDiff(), 0)
        sd3_ans = self.pmc.sd([1,(0,2)], False)
        self.assertEqual(self.sd3.antiDiff(), 1*sd3_ans)
        self.assertEqual(self.sd4.antiDiff(), 0)

    def testFactor(self):
        # Number after each strand diagram is the number of factorizations
        tests = [(self.sd1, 8), (self.sd2, 3), (self.sd3, 2), (self.sd4, 9)]
        for sd, n in tests:
            self.assertEqual(len(sd.factor()), n)

    def testStrandDiagramOppGrading(self):
        # Verify the relations gr'(a).opp() = gr'(a.opp()) and
        # gr(a).opp() = gr(a.opp()).
        sd_to_test = [self.sd1, self.sd2, self.sd3, self.sd4]
        for sd in sd_to_test:
            self.assertEqual(sd.getBigGrading().opp(),
                             sd.opp().getBigGrading())
            # This only works if DEFAULT_REFINEMENT is set to averageRefinement.
            self.assertEqual(
                sd.getSmallGrading(refinement = averageRefinement).opp(),
                sd.opp().getSmallGrading(refinement = averageRefinement))

    def testUnconnectSumStrandDiagram(self):
        sd1 = self.pmc.sd([4,(0,3)], False)
        pmc1 = splitPMC(1)
        self.assertEqual(unconnectSumStrandDiagram(sd1, 1),
                         (pmc1.sd([(0,3)], False), pmc1.sd([0], False)))
        sd2 = self.pmc.sd([(1,2),(6,7)], False)
        self.assertEqual(unconnectSumStrandDiagram(sd2, 1),
                         (pmc1.sd([(1,2)], False), pmc1.sd([(2,3)], False)))

if __name__ == "__main__":
    unittest.main()
