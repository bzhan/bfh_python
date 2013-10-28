"""Unit test for arcslideda.py"""

from arcslideda import *
from arcslide import Arcslide
from latex import beginDoc, endDoc, showArrow
from pmc import PMC
from pmc import linearPMC, splitPMC
import unittest

class ArcslideDATest(unittest.TestCase):
    def testShortArcslide(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 0),
            Arcslide(splitPMC(1), 2, 1),
            Arcslide(splitPMC(2), 1, 0),
            Arcslide(splitPMC(2), 6, 5),
            Arcslide(linearPMC(2), 1, 0),
            Arcslide(linearPMC(2), 6, 5),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 1, 0),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 6, 5),
            Arcslide(splitPMC(2), 2, 1),
            Arcslide(splitPMC(2), 5, 4),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 4, 3),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 3, 2),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).getDAStructure()
            self.assertTrue(dastr.testDelta())

    def testShortArcslideLocal(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 0),
            Arcslide(splitPMC(1), 2, 1),
            Arcslide(splitPMC(2), 2, 1),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testLength3Arcslide(self):
        slides_to_test = [
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 1, 0),
            Arcslide(PMC([(0, 3), (1, 4), (2, 6), (5, 7)]), 1, 0),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 1, 0),
            Arcslide(PMC([(0, 3), (1, 7), (2, 5), (4, 6)]), 1, 0),
            Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 2, 1),
            Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 4, 3),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).getDAStructure()
            self.assertTrue(dastr.testDelta())

    def testLength3ArcslideLocal(self):
        slides_to_test = [
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 1, 0),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testUnderslide(self):
        slides_to_test = [
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 2, 1),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 2, 1),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).getDAStructure()
            self.assertTrue(dastr.testDelta())

    def testUnderslideLocal(self):
        slides_to_test = [
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())        

    def testAutoCompleteArcslide(self):
        slide = Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1)
        local_dastr = ArcslideDA(slide).getLocalDAStructure()
        self.assertTrue(local_dastr.testDelta())

        arrows = []
        patterns = ArcslideDA(slide).arrow_patterns
        for coeffs_a, lst_d in patterns.items():
            for coeff_d in lst_d:
                arrows.append((coeff_d, coeffs_a))
        # Order is: 5 for top, 0 for bottom, 4 for upper middle, 2 for lower middle
        arrows_base = [(coeff_d, coeffs_a) for coeff_d, coeffs_a in arrows
                       if coeff_d.multiplicity[2] == 0]
        arrows_new = [(coeff_d, coeffs_a) for coeff_d, coeffs_a in arrows
                      if coeff_d.multiplicity[2] != 0]
        add_strs = []
        single_idems = [(1, 1)]
        arrows_to_extend = ArcslideDA.autoCompleteByLinAlg(
            arrows_base, arrows_new, single_idems,
            local_dastr.getGenerators())
        for coeff_d, coeffs_a in arrows_to_extend:
            add_strs.append("(%s)," % ", ".join(coeff.inputForm() for coeff in
                                                list(coeffs_a) + [coeff_d]))
            arrows.append((coeff_d, coeffs_a))

        while True:
            result = ArcslideDA.autoCompleteArrows(
                arrows, single_idems, local_dastr.getGenerators())
            if result[0] is True:
                print "Finished"
                for add_str in add_strs:
                    print add_str
                break
            elif result[0] is False:
                print "Cannot continue"
                for add_str in add_strs:
                    print add_str
                f = open("tocancel_arrows.txt", "w")
                f.write(beginDoc())
                for coeff_d, coeffs_a in result[1]:
                    print "To cancel: ", coeff_d, coeffs_a
                    f.write(showArrow(coeff_d, coeffs_a, [0,1,2,4,5], [0,2,3,4,5]))
                f.write(endDoc())
                f.close()
                suggestions = []
                for (coeff_d, coeffs_a), times in result[2].items():
                    suggestions.append((times, coeff_d, coeffs_a))
                for times, coeff_d, coeffs_a in reversed(sorted(suggestions)):
                    print times, coeffs_a, coeff_d
                break
            else:
                coeff_d, coeffs_a = result
                add_str = "(%s)," % ", ".join(coeff.inputForm() for coeff in
                                              list(coeffs_a) + [coeff_d])
                print "To add:", add_str
                add_strs.append(add_str)
                arrows.append((coeff_d, coeffs_a))

if __name__ == "__main__":
    unittest.main()
