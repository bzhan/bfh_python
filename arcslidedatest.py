"""Unit test for arcslideda.py"""

from arcslideda import *
import unittest

class ArcslideDATest(unittest.TestCase):
    def testShortArcslide(self):
        dastr1 = ArcslideDA.getDAStructure(Arcslide(splitPMC(1), 1, 2))
        self.assertTrue(dastr1.testDelta())
        dastr2 = ArcslideDA.getDAStructure(Arcslide(splitPMC(2), 1, 2))
        self.assertTrue(dastr2.testDelta())

if __name__ == "__main__":
    unittest.main()
