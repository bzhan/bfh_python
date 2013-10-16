"""Unit test for arcslideda.py"""

from arcslideda import *
from arcslide import Arcslide
from pmc import PMC
from pmc import linearPMC, splitPMC
import unittest

class ArcslideDATest(unittest.TestCase):
    def testShortArcslide(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 2),
            Arcslide(splitPMC(1), 2, 3),
            Arcslide(splitPMC(2), 1, 2),
            Arcslide(splitPMC(2), 6, 7),
            Arcslide(linearPMC(2), 1, 2),
            Arcslide(linearPMC(2), 6, 7),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 1, 2),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 6, 7),
            Arcslide(splitPMC(2), 2, 3),
            Arcslide(splitPMC(2), 5, 6),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 4, 5),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 3, 4),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).getDAStructure()
            self.assertTrue(dastr.testDelta())

    def testShortArcslideLocal(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 2),
            Arcslide(splitPMC(1), 2, 3),
            Arcslide(splitPMC(2), 2, 3),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testLength3Arcslide(self):
        slides_to_test = [
            # Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 3),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

if __name__ == "__main__":
    unittest.main()
