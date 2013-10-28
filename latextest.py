"""Unit test for latex.py"""

from latex import *
from arcslide import Arcslide
from arcslideda import ArcslideDA
from localpmc import LocalPMC, LocalStrandDiagram
from pmc import PMC
from pmc import splitPMC
import unittest

class LatexTest(unittest.TestCase):
    def testPrintLocalStrandDiagram(self):
        pmc1 = LocalPMC(5, [(1, 3),(2,)], [0, 4])
        sd1 = pmc1.sd([2, (0, 1)])
        sd2 = pmc1.sd([1, (0, 4)])
        print beginDoc()
        print showLocalStrandDiagram(sd1)
        print showLocalStrandDiagram(sd2)
        print endDoc()

    def testPrintDAStructure(self):
        # slide = Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 1, 0)
        # slide = Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 2, 1)
        slide = Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1)
        local_dastr = ArcslideDA(slide).getLocalDAStructure()
        f = open("latex_output.txt", "w")
        f.write(beginDoc())
        # For middle length 3
        # f.write(showDAStructure(local_dastr, [0,1,2,3,5,6], [0,1,3,4,5,6]))
        # For general underslide
        f.write(showDAStructure(local_dastr, [0,1,2,3,4,6,7], [0,1,3,4,5,6,7]))
        f.write(endDoc())
        f.close()

    def testPrintArrows(self):
        slide = Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 1, 0)
        patterns = ArcslideDA(slide).arrow_patterns
        f = open("latex_output.txt", "w")
        f.write(beginDoc())
        for coeffs_a, lst_d in patterns.items():
            for coeff_d in lst_d:
                f.write(showArrow(coeff_d, coeffs_a, [0,1,2,4,5], [0,2,3,4,5]))
        f.write(endDoc())
        f.close()

if __name__ == "__main__":
    unittest.main()
