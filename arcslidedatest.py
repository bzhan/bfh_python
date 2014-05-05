"""Unit test for arcslideda.py"""

from arcslideda import *
from arcslide import Arcslide
from dstructure import zeroTypeD
from latex import beginDoc, endDoc, showArrow
from pmc import PMC
from pmc import antipodalPMC, linearPMC, splitPMC
import unittest

class ArcslideDATest(unittest.TestCase):
    def testShortUnderslideDown(self):
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
            dastr = ArcslideDA(slide).toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())

    def testShortUnderslideDownLocal(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 0),
            Arcslide(splitPMC(1), 2, 1),
            Arcslide(splitPMC(2), 2, 1),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testGeneralUnderslideDown(self):
        slides_to_test = [
            Arcslide(antipodalPMC(2), 1, 0),
            Arcslide(antipodalPMC(2), 2, 1),
            Arcslide(antipodalPMC(2), 3, 2),
            Arcslide(antipodalPMC(2), 4, 3),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 1, 0),
            Arcslide(PMC([(0, 5), (1, 3), (2, 6), (4, 7)]), 1, 0),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 1, 0),
            Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 2, 1),
            Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 4, 3),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 2, 1),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 2, 1),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 5, 4),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 3, 2),
            Arcslide(PMC([(0, 5), (1, 3), (2, 6), (4, 7)]), 5, 4),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())

    def testGeneralUnderslideDownLocal(self):
        slides_to_test = [
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 1, 0),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 5, 4),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testGeneralOverslideDown(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 3, 2),
            Arcslide(splitPMC(2), 3, 2),
            Arcslide(splitPMC(2), 4, 3),
            Arcslide(splitPMC(2), 7, 6),
            Arcslide(linearPMC(2), 3, 2),
            Arcslide(linearPMC(2), 5, 4),
            Arcslide(linearPMC(2), 7, 6),
            Arcslide(antipodalPMC(2), 5, 4),
            Arcslide(antipodalPMC(2), 6, 5),
            Arcslide(antipodalPMC(2), 7, 6),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())

    def testGeneralOverslideDownLocal(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 3, 2),
            Arcslide(splitPMC(2), 3, 2),
            Arcslide(splitPMC(2), 4, 3),
            Arcslide(splitPMC(2), 7, 6),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testShortUnderslideUp(self):
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
            dastr = ArcslideDA(slide).toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())

    def testShortUnderslideUpLocal(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 2),
            Arcslide(splitPMC(1), 2, 3),
            Arcslide(splitPMC(2), 2, 3),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testGeneralUnderslideUp(self):
        slides_to_test = [
            Arcslide(antipodalPMC(2), 3, 4),
            Arcslide(antipodalPMC(2), 4, 5),
            Arcslide(antipodalPMC(2), 5, 6),
            Arcslide(antipodalPMC(2), 6, 7),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 3),
            Arcslide(PMC([(0, 5), (1, 3), (2, 6), (4, 7)]), 4, 5),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 2, 3),
            Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 3, 4),
            Arcslide(PMC([(0, 2), (1, 4), (3, 6), (5, 7)]), 5, 6),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 5, 6),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 5, 6),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 4, 5),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 6, 7),
            Arcslide(PMC([(0, 3), (1, 5), (2, 7), (4, 6)]), 6, 7),
            Arcslide(PMC([(0, 5), (1, 3), (2, 6), (4, 7)]), 6, 7),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())

    def testGeneralUnderslideUpLocal(self):
        slides_to_test = [
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 3),
            Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 5, 6),
            Arcslide(PMC([(0, 2), (1, 6), (3, 5), (4, 7)]), 6, 7),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testGeneralOverslideUp(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 0, 1),
            Arcslide(splitPMC(2), 0, 1),
            Arcslide(splitPMC(2), 3, 4),
            Arcslide(splitPMC(2), 4, 5),
            Arcslide(linearPMC(2), 0, 1),
            Arcslide(linearPMC(2), 2, 3),
            Arcslide(linearPMC(2), 4, 5),
            Arcslide(antipodalPMC(2), 0, 1),
            Arcslide(antipodalPMC(2), 1, 2),
            Arcslide(antipodalPMC(2), 2, 3),
        ]
        for slide in slides_to_test:
            print slide
            dastr = ArcslideDA(slide).toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())

    def testGeneralOverslideUpLocal(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 0, 1),
            Arcslide(splitPMC(2), 0, 1),
            Arcslide(splitPMC(2), 3, 4),
            Arcslide(splitPMC(2), 4, 5),
        ]
        for slide in slides_to_test:
            local_dastr = ArcslideDA(slide).getLocalDAStructure()
            self.assertTrue(local_dastr.testDelta())

    def testAutoCompleteArcslide(self):
        def run_test(slide, d_side_order, intervals, single_idems):
            local_dastr = ArcslideDA(slide).getLocalDAStructure()

            # Recover list of arrows.
            arrows = []
            patterns = ArcslideDA(slide).arrow_patterns
            for coeffs_a, lst_d in patterns.items():
                for coeff_d in lst_d:
                    arrows.append((coeff_d, coeffs_a))

            # Auto-complete by linear algebra proceeds stages as specified by
            # d_side_order.
            for i in range(len(d_side_order)):
                arrows_base, arrows_new, arrows_expected = [], [], []
                int_d, int_a = intervals[i]
                for coeff_d, coeffs_a in arrows:
                    if all([coeff_d.multiplicity[p] == 0
                            for p in d_side_order[i:]]):
                        # already exists
                        arrows_base.append((coeff_d, coeffs_a))
                    elif coeff_d.strands == (int_d,) and \
                         len(coeffs_a) == 1 and coeffs_a[0].strands == (int_a,):
                        # seed for this round
                        arrows_new.append((coeff_d, coeffs_a))
                # Triggers assert failure if it does not succeed
                arrows_to_extend = ArcslideDA.autoCompleteByLinAlg(
                    arrows_base, arrows_new, single_idems,
                    local_dastr.getGenerators())

                # Uncomment the following to see the added arrows.
                # print "Added at this stage:", len(arrows_to_extend)
                # for coeff_d, coeffs_a in arrows_to_extend:
                #     print "(%s)," % ", ".join(
                #         coeff.inputForm()
                #         for coeff in list(coeffs_a) + [coeff_d])

        # Short underslide down
        run_test(slide = Arcslide(splitPMC(2), 2, 1),
                 d_side_order = (3, 0),
                 intervals = [((3, 4), (3, 4)), ((0, 1), (0, 1))],
                 single_idems = [(1, 1)])

        # General underslide down
        run_test(slide = Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 2, 1),
                 d_side_order = (5, 0, 4, 2),
                 intervals = [((5, 6), (5, 6)), ((0, 1), (0, 1)),
                              ((4, 5), (3, 4)), ((2, 3), (1, 2))],
                 single_idems = [(1, 1)])

        # Short underslide up
        run_test(slide = Arcslide(splitPMC(2), 2, 3),
                 d_side_order = (3, 0),
                 intervals = [((3, 4), (3, 4)), ((0, 1), (0, 1))],
                 single_idems = [(1, 1)])

        # General underslide up
        run_test(slide = Arcslide(PMC([(0, 3), (1, 6), (2, 4), (5, 7)]), 5, 6),
                 d_side_order = (5, 0, 3, 1),
                 intervals = [((5, 6), (5, 6)), ((0, 1), (0, 1)),
                              ((3, 4), (4, 5)), ((1, 2), (2, 3))],
                 single_idems = [(1, 1)])

    def testGrading(self):
        slides_to_test = [
            Arcslide(splitPMC(1), 1, 0),
            Arcslide(splitPMC(2), 4, 3),
            Arcslide(linearPMC(2), 1, 0),
        ]
        for slide in slides_to_test:
            dastr = ArcslideDA(slide)
            dastr.toSimpleDAStructure().checkGrading()

class TensorTest(unittest.TestCase):
    def testDATensorD(self):
        # So far mostly checking that it will run in a reasonable amount of
        # time.
        slide = Arcslide(splitPMC(5), 2, 1)  # will change zeroTypeD
        dastr = ArcslideDA(slide)
        dstr = zeroTypeD(5)
        dstr_result = dastr.tensorD(dstr)
        dstr_result.reindex()
        self.assertEqual(len(dstr_result), 2)

if __name__ == "__main__":
    unittest.main()
