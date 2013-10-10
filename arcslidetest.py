"""Unit test for arcslide.py"""

from arcslide import *
from utility import DEFAULT_GRADING, SMALL_GRADING
import unittest

class ArcslideTest(unittest.TestCase):
    def testArcslide(self):
        slide1 = Arcslide(splitPMC(1), 0, 1)
        self.assertEqual(slide1.c2, 3)
        self.assertEqual(slide1.slide_type, OVER_SLIDE)
        self.assertEqual(slide1.end_pmc, splitPMC(1))
        self.assertEqual(slide1.inverse(), Arcslide(splitPMC(1), 3, 2))

        slide2 = Arcslide(splitPMC(2), 3, 4)
        self.assertEqual(slide2.end_pmc, PMC([(0,2),(1,6),(3,5),(4,7)]))
        self.assertEqual(slide2.inverse().end_pmc, splitPMC(2))

        slide3 = Arcslide(PMC([(0,2),(1,6),(3,5),(4,7)]), 5, 6)
        self.assertEqual(slide3.slide_type, UNDER_SLIDE)
        self.assertEqual(slide3.end_pmc, PMC([(0,3),(1,6),(2,4),(5,7)]))

    def testUnderslideDDStructure(self):
        slide1 = Arcslide(splitPMC(1), 1, 0)
        ddstr = slide1.getDDStructure()
        self.assertTrue(ddstr.testDelta())
        slide2 = Arcslide(splitPMC(2), 1, 0)
        ddstr2 = slide2.getDDStructure()
        self.assertTrue(ddstr2.testDelta())
        slide3 = Arcslide(PMC([(0,2),(1,6),(3,5),(4,7)]), 5, 6)
        ddstr3 = slide3.getDDStructure()
        self.assertTrue(ddstr3.testDelta())

    def testOverslideDDStructure(self):
        slide1 = Arcslide(splitPMC(1), 0, 1)
        ddstr = slide1.getDDStructure()
        self.assertTrue(ddstr.testDelta())
        slide2 = Arcslide(splitPMC(2), 1, 0)
        ddstr2 = slide2.getDDStructure()
        self.assertTrue(ddstr2.testDelta())
        slide3 = Arcslide(splitPMC(2), 3, 4)
        ddstr3 = slide3.getDDStructure()
        self.assertTrue(ddstr3.testDelta())
        slide4 = Arcslide(splitPMC(2), 4, 3)
        ddstr4 = slide4.getDDStructure()
        self.assertTrue(ddstr4.testDelta())

    def testGenus3DDStructure(self):
        slide1 = Arcslide(splitPMC(3), 0, 1)
        ddstr1 = slide1.getDDStructure()

    def testArcslideMatchDiagram(self):
        slide_to_test = [Arcslide(PMC([(0,2),(1,6),(3,5),(4,7)]), 5, 6),
                         Arcslide(splitPMC(2), 4, 3)]
        for slide in slide_to_test:
            ddstr = slide.getDDStructure()
            if DEFAULT_GRADING == SMALL_GRADING:
                self.assertTrue(ddstr.gr_set.isAutomorphism())

if __name__ == "__main__":
    unittest.main()
