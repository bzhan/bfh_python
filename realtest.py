"""Unit test for real.py"""

from real import realAZ, realAZda, selfGluingHB, connectSumPMC
from arcslide import Arcslide
from arcslideda import ArcslideDA
from pmc import linearPMC, PMC
from braid import Braid
import unittest

class RealTest(unittest.TestCase):
    def testRealAZ(self):
        "Test realAZ satisfies type D structure relation"
        P = realAZ(linearPMC(2))
        self.assertTrue(P.testDelta())

    def testSelfGluingHB(self):
        "Test self gluing handlebody satisfies type D structure relation"
        P = selfGluingHB(linearPMC(2))
        self.assertTrue(P.testDelta())

    def testRealAZda(self):
        "Test realAZda satisfies type DA structure relation"
        M = realAZda(linearPMC(1)) #Testing it for 2 takes a long time.
        self.assertTrue(M.testDelta())

    def testTrefoil(self):
        "Test that HFR of branched double cover of trefoil has rank 3, as does HF^. Uses only type D modules."
        hb = selfGluingHB(linearPMC(1))
        adehntwist = Arcslide(hb.algebra.pmc,2,1)
        bdehntwist = Arcslide(hb.algebra.pmc,1,0)
        ada = ArcslideDA(adehntwist)
        bda = ArcslideDA(bdehntwist)
        twisted = ada.tensorD(bda.tensorD(hb))
        az = realAZ(linearPMC(1))
        cx = az.morToD(twisted)
        cx.simplify()
        self.assertEqual(len(cx),3) #HFR(Sigma(trefoil))
        moretwisted = ada.tensorD(bda.tensorD(twisted))
        cx2 = hb.morToD(moretwisted)
        cx2.simplify()
        self.assertEqual(len(cx2),3) #HF^(Sigma(trefoil))

    def testFigure8(self):
        "Test that HFR of branched double cover of figure 8 knot has rank 3, as does HF^. Uses only type D modules."
        hb = selfGluingHB(linearPMC(1))
        adehntwist = Arcslide(hb.algebra.pmc,1,0)
        binvdehntwist = Arcslide(hb.algebra.pmc,2,3)
        ada = ArcslideDA(adehntwist)
        binvda = ArcslideDA(binvdehntwist)
        twisted = ada.tensorD(binvda.tensorD(hb))
        az = realAZ(linearPMC(1))
        cx = az.morToD(twisted)
        cx.simplify()
        self.assertEqual(len(cx),5) #HFR(Sigma(4_1))
        moretwisted = ada.tensorD(binvda.tensorD(twisted))
        cx2 = hb.morToD(moretwisted)
        cx2.simplify()
        self.assertEqual(len(cx2),5) #HF^(Sigma(4_1))

    def testTrefoilDA(self):
        "Test that HFR of branched double cover of trefoil has rank 3, as does HF^. Uses da module for real AZ piece."
        hb = selfGluingHB(linearPMC(1))
        adehntwist = Arcslide(hb.algebra.pmc,2,1)
        bdehntwist = Arcslide(hb.algebra.pmc,1,0)
        ada = ArcslideDA(adehntwist)
        bda = ArcslideDA(bdehntwist)
        twisted = ada.tensorD(bda.tensorD(hb))
        az = realAZda(linearPMC(1))
        cx = az.tensorD(twisted)
        cx.simplify()
        self.assertEqual(len(cx),3) #HFR(Sigma(trefoil))

    def testFigure8DA(self):
        "Test that HFR of branched double cover of figure 8 has rank 5. Uses da module for real AZ piece."
        hb = selfGluingHB(linearPMC(1))
        adehntwist = Arcslide(hb.algebra.pmc,1,0)
        binvdehntwist = Arcslide(hb.algebra.pmc,2,3)
        ada = ArcslideDA(adehntwist)
        binvda = ArcslideDA(binvdehntwist)
        twisted = ada.tensorD(binvda.tensorD(hb))
        az = realAZda(linearPMC(1))
        cx = az.tensorD(twisted)
        cx.simplify()
        self.assertEqual(len(cx),5) #HFR(Sigma(4_1))

    def testMoreKnots(self):
        "Compute HFR of branched double cover of some more knots, using D structures. Data is list of pairs (name, braid word whose branched double cover is monodromy of the knot, expected rank of HFR). This takes a long time."
        data = [("10_145",[2,2,3,-2,4,1],5)]#,("12n121",[2,2,3,2,-4,-5,1],3),("12n642",[2,2,3,3,4,-3,-2,1],29)]
        braid_group = Braid(6)
        pmc: PMC = linearPMC(2)
        pmc2: PMC = connectSumPMC(pmc,pmc.opp())
        for (name, word, rank) in data:
            slides = braid_group.getArcslides(word)
            actual_slides: list[Arcslide] = list()
            cur_pmc = pmc2
            for s in slides:
                actual_slides.append(Arcslide(cur_pmc,s.b1,s.c1))
                cur_pmc = actual_slides[-1].end_pmc
            Q = selfGluingHB(pmc)
            actual_slides = [slide.inverse() for slide in actual_slides]
            for s in actual_slides:
                Q = ArcslideDA(s).tensorD(Q)
                Q.reindex()
                Q.simplify()
            cx = Q.morToD(realAZ(pmc))
            cx.simplify()
            self.assertEqual(len(cx),rank)

    def testMoreKnotsDA(self):
        "Compute HFR of branched double cover of some more knots, using A AZ at end. Data is list of pairs (name, braid word whose branched double cover is monodromy of the knot, expected rank of HFR). This takes a long time."
        data = [("10_145",[2,2,3,-2,4,1],5),("12n121",[2,2,3,2,-4,-5,1],3),("12n642",[2,2,3,3,4,-3,-2,1],29)]
        braid_group = Braid(6)
        pmc: PMC = linearPMC(2)
        pmc2: PMC = connectSumPMC(pmc,pmc.opp())
        for (name, word, rank) in data:
            slides = braid_group.getArcslides(word)
            actual_slides: list[Arcslide] = list()
            cur_pmc = pmc2
            for s in slides:
                actual_slides.append(Arcslide(cur_pmc,s.b1,s.c1))
                cur_pmc = actual_slides[-1].end_pmc
            Q = selfGluingHB(pmc)
            actual_slides = [slide.inverse() for slide in actual_slides]
            for s in actual_slides:
                Q = ArcslideDA(s).tensorD(Q)
                Q.reindex()
                Q.simplify()
            cx = realAZda(pmc).tensorD(Q)
            cx.simplify()
            self.assertEqual(len(cx),rank)

if __name__ == "__main__":
    unittest.main()
