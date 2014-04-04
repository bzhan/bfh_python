"""Unit test for grading.py"""

from grading import *
from pmc import *
import unittest

class BigGradingTest(unittest.TestCase):
    def testMultiply(self):
        pmc = splitPMC(1)
        bgrp = BigGradingGroup(pmc)
        elt1 = BigGradingElement(bgrp, 0, [1,0,0])
        elt2 = BigGradingElement(bgrp, 0, [0,1,0])
        elt12 = BigGradingElement(bgrp, Fraction(1,2), [1,1,0])
        elt21 = BigGradingElement(bgrp, Fraction(-1,2), [1,1,0])
        self.assertEqual(elt1 * elt2, elt12)
        self.assertEqual(elt2 * elt1, elt21)
        self.assertEqual(elt12, elt12.toSmallGrading().toBigGrading())

class SmallGradingTest(unittest.TestCase):
    def testMultiply(self):
        pmc = splitPMC(1)
        sgrp = SmallGradingGroup(pmc)
        elt1 = SmallGradingElement(sgrp, 0, [1,0])
        elt2 = SmallGradingElement(sgrp, 0, [0,1])
        elt12 = SmallGradingElement(sgrp, 1, [1,1])
        elt21 = SmallGradingElement(sgrp, -1, [1,1])
        self.assertEqual(elt1 * elt2, elt12)
        self.assertEqual(elt2 * elt1, elt21)
        self.assertEqual(elt1, elt1.toBigGrading().toSmallGrading())
        self.assertEqual(elt2, elt2.toBigGrading().toSmallGrading())

class TestCommonRefinement(unittest.TestCase):
    def testStandardRefinement(self):
        pmc_to_test = [splitPMC(1), splitPMC(2), antipodalPMC(2)]
        for pmc in pmc_to_test:
            refinement = standardRefinement(pmc)

    def testLowerRefinement(self):
        pmc_to_test = [splitPMC(1), splitPMC(2), antipodalPMC(2)]
        for pmc in pmc_to_test:
            refinement = lowerRefinement(pmc)

    def testAverageRefinement(self):
        pmc_to_test = [splitPMC(1), splitPMC(2), antipodalPMC(2)]
        for pmc in pmc_to_test:
            refinement = averageRefinement(pmc)

class TestSimpleGradingSet(unittest.TestCase):
    def testSimpleGradingSet(self):
        pmc = splitPMC(2)
        periodic_domains = [pmc.small_gr(0, [1,0,0,0]),
                            pmc.small_gr(0, [0,0,1,0])]
        gr_set = SimpleGradingSet(SmallGradingGroup(pmc), ACTION_LEFT,
                                  periodic_domains)
        elt1 = gr_set.zero()
        elt2 = elt1 * pmc.small_gr(0, [0,1,0,0])
        elt3 = elt1 * pmc.small_gr(0, [1,0,0,0])
        elt4 = SimpleGradingSetElement(gr_set, pmc.small_gr(0, [0,1,0,0]))
        elt5 = elt4 * pmc.small_gr(-2, [1,0,0,0])
        elt6 = elt4 * pmc.small_gr(0, [1,0,0,0])
        self.assertNotEqual(elt1, elt2)
        self.assertEqual(elt1, elt3)
        self.assertEqual(elt4, elt5)
        self.assertNotEqual(elt4, elt6)

class TestSimpleDbGradingSet(unittest.TestCase):
    def testSimpleDbGradingSet(self):
        pmc = splitPMC(1)
        periodic_domains = [(pmc.big_gr(0, [1,1,0]), pmc.big_gr(0, [-1,-1,0])),
                            (pmc.big_gr(0, [0,1,1]), pmc.big_gr(0, [0,-1,-1]))]
        gr_set = SimpleDbGradingSet(BigGradingGroup(pmc), ACTION_LEFT,
                                    BigGradingGroup(pmc), ACTION_LEFT,
                                    periodic_domains)
        elt1 = gr_set.zero()
        elt2 = elt1 * [pmc.big_gr(0, [-1,-1,0]), pmc.big_gr(0, [0,0,0])]
        elt3 = elt1 * [pmc.big_gr(0, [-1,-1,0]), pmc.big_gr(0, [1,1,0])]
        self.assertNotEqual(elt1, elt2)
        self.assertEqual(elt1, elt3)

class TestGeneralGradingSet(unittest.TestCase):
    def testGeneralGradingSet(self):
        pmc = splitPMC(1)
        periodic_domains1 = [
            (pmc.big_gr(0, [1,1,0]), pmc.big_gr(0, [-1,-1,0])),
            (pmc.big_gr(0, [0,1,1]), pmc.big_gr(0, [0,-1,-1]))]
        gr_set1 = SimpleDbGradingSet(BigGradingGroup(pmc), ACTION_LEFT,
                                     BigGradingGroup(pmc), ACTION_RIGHT,
                                     periodic_domains1)
        periodic_domains2 = [pmc.big_gr(0, [1,1,0])]
        gr_set2 = SimpleGradingSet(BigGradingGroup(pmc), ACTION_LEFT,
                                   periodic_domains2)
        gr_set = GeneralGradingSet([gr_set1, gr_set2])
        elt1 = gr_set.zero()
        l2 = gr_set1.zero() * [pmc.big_gr(0, [1,1,0]), pmc.big_gr(0, [0,0,0])]
        l3 = gr_set1.zero() * [pmc.big_gr(0, [0,1,1]), pmc.big_gr(0, [0,0,0])]
        elt2 = GeneralGradingSetElement(gr_set, [l2, gr_set2.zero()])
        elt3 = GeneralGradingSetElement(gr_set, [l3, gr_set2.zero()])
        self.assertEqual(elt1, elt2)
        self.assertNotEqual(elt1, elt3)

    def testSimplifiedGradingSet(self):
        pmc = splitPMC(1)
        periodic_domains1 = [(pmc.small_gr(0, [1,1]), pmc.small_gr(0, [1,0])),
                             (pmc.small_gr(0, [2,3]), pmc.small_gr(0, [0,1]))]
        gr_set1 = SimpleDbGradingSet(SmallGradingGroup(pmc), ACTION_LEFT,
                                     SmallGradingGroup(pmc), ACTION_RIGHT,
                                     periodic_domains1)
        periodic_domains2 = [pmc.small_gr(0, [0,1])]
        gr_set2 = SimpleGradingSet(SmallGradingGroup(pmc), ACTION_LEFT,
                                   periodic_domains2)
        gr_set = GeneralGradingSet([gr_set1, gr_set2])
        gr_set_short = gr_set.simplifiedSet()

        r1 = gr_set2.zero() * [pmc.small_gr(0, [1,-1])]
        elt1 = GeneralGradingSetElement(gr_set, [gr_set1.zero(), r1])
        elt1_short = gr_set.simplifiedElt(elt1)

        # Actually simplified
        self.assertTrue(isinstance(gr_set_short, SimpleGradingSet))
        self.assertTrue(isinstance(elt1_short, SimpleGradingSetElement))

        # elt1_short should be the 'same' element as before
        l2 = gr_set1.zero() * [elt1_short.data, pmc.small_gr(0, [0,0])]
        elt2 = GeneralGradingSetElement(gr_set, [l2, gr_set2.zero()])
        self.assertEqual(elt1, elt2)

    def testNoSimplification(self):
        pmc = splitPMC(1)
        # Does not form automorphism
        periodic_domains1 = [(pmc.small_gr(0, [1,0]), pmc.small_gr(0, [0,0])),
                             (pmc.small_gr(0, [0,0]), pmc.small_gr(0, [0,1]))]
        gr_set1 = SimpleDbGradingSet(SmallGradingGroup(pmc), ACTION_LEFT,
                                     SmallGradingGroup(pmc), ACTION_RIGHT,
                                     periodic_domains1)
        periodic_domains2 = [pmc.small_gr(0, [0,1])]
        gr_set2 = SimpleGradingSet(SmallGradingGroup(pmc), ACTION_LEFT,
                                   periodic_domains2)
        gr_set = GeneralGradingSet([gr_set1, gr_set2])
        gr_set_short = gr_set.simplifiedSet()
        self.assertEqual(gr_set, gr_set_short)

        elt1 = GeneralGradingSetElement(
            gr_set, [gr_set1.zero(), gr_set2.zero()])
        elt1_short = gr_set.simplifiedElt(elt1)
        self.assertEqual(elt1, elt1_short)

if __name__ == "__main__":
    unittest.main()
