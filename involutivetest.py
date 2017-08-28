"""Unit test for involutive.py"""

from involutive import *
from dstructure import zeroTypeD, infTypeD
from arcslide import Arcslide
from arcslideda import ArcslideDA
from pmc import splitPMC
import unittest

class InvolutiveTest(unittest.TestCase):
    def testInvolutiveRankS3(self):
        "Check that the involutive Floer homology of S^3 is F_2^2."
        P = zeroTypeD(1)
        Q = infTypeD(1)
        cx = involutiveCx(P,Q)
        self.assertTrue(len(cx)==2)

    def testInvolutiveRankS2S1(self):
        "Check that the involutive Floer homology of S^2 x S^1 is F_2^4."
        P = zeroTypeD(1)
        cx = involutiveCx(P,P)
        self.assertTrue(len(cx)==4)

    def testInvolutiveRankL31(self):
        "Check that the involutive Floer homology of L(3,1) is F_2^4."
        P = infTypeD(1)
        Q = infTypeD(1)
        R = ArcslideDA(Arcslide(splitPMC(1),1,0))
        Q2 = R.tensorD(Q)
        Q3 = R.tensorD(Q2)
        Q4 = R.tensorD(Q3)
        othcx = P.morToD(Q4)
        othcx.simplify()
        cx = involutiveCx(P,Q3)
        self.assertTrue(len(othcx)==3)
        self.assertTrue(len(cx)==4)

    def testInvolutiveRankS2S1S2S1(self):
        "Check that the involutive Floer homology of a connect sum of S^2 x S^1 with itself is SOMETHING."
        P = zeroTypeD(2)
        cx = involutiveCx(P,P)
        self.assertTrue(len(cx)==8)

    def testIFHDCov(self):
        "Check the involutive Floer homology of some branched double covers."
        self.assertTrue(len(invOfDCov("4_1 3 1 0 5 4 3 2 4 1 2 3 2 1 2 1 0 5 4 3 2 1 0"))==6)


if __name__ == "__main__":
    unittest.main()
