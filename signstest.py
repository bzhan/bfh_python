"""Unit test for signs.py."""

from signs import *
from grading import DEFAULT_REFINEMENT, lowerRefinement
from pmc import PMC
from pmc import antipodalPMC, linearPMC, splitPMC
from utility import ZZ
import unittest

class AbsZ2GradingTest(unittest.TestCase):
    def testAbsGrading(self):
        def testOneAlgebra(alg, test_op = True):
            abs_gr = AbsZ2Grading(alg)
            for gen in alg.getGenerators():
                # Test asserts in getAbsGrading
                abs_gr.getAbsGrading(gen)
            if not test_op:
                return
            # Test differential and multiplication
            for a in alg.getGenerators():
                for term in a.diff():
                    a_gr, da_gr = [abs_gr.getAbsGrading(gen) for gen in (a, term)]
                    assert (a_gr - 1) % 2 == da_gr
            for a in alg.getGenerators():
                for b in alg.getGenerators():
                    if a * b != 0:
                        a_gr, b_gr, ab_gr = [abs_gr.getAbsGrading(gen)
                                             for gen in (a, b, (a*b).getElt())]
                        assert (a_gr + b_gr) % 2 == ab_gr

        for pmc in [splitPMC(1), splitPMC(2), linearPMC(2)]:
            testOneAlgebra(pmc.getAlgebra())
        for pmc in [antipodalPMC(2), splitPMC(3)]:
            testOneAlgebra(pmc.getAlgebra(), test_op = False)
        for (pmc, idem_size) in [(splitPMC(1), 1), (splitPMC(1), 2),
                                 (splitPMC(2), 1), (splitPMC(2), 2)]:
            testOneAlgebra(PreStrandAlgebra(F2, pmc, idem_size))
        for (pmc, idem_size) in [(splitPMC(2), 3), (splitPMC(2), 4),
                                 (splitPMC(3), 2)]:
            testOneAlgebra(PreStrandAlgebra(F2, pmc, idem_size),
                           test_op = False)

class PreStrandAlgebraTest(unittest.TestCase):
    def testGetGenerators(self):
        for pmc, idem_size, n in [(splitPMC(1), 1, 10),
                                  (splitPMC(1), 2, 25),
                                  (splitPMC(2), 1, 36),
                                  (splitPMC(2), 2, 462),
                                  (splitPMC(2), 3, 2646),
                                  (splitPMC(2), 4, 6951),
                                  (splitPMC(3), 2, 2431),
                                  (splitPMC(3), 3, 39325)]:
            # Further numbers:
            # splitPMC(3), 4 --> 359502
            algebra = PreStrandAlgebra(F2, pmc, idem_size)
            self.assertEqual(n, len(algebra.getGenerators()))

    def testDiff(self):
        algebra = PreStrandAlgebra(F2, splitPMC(1), 2)
        for sd, sd_diff in [([(1, 3), (2, 2)], [(1, 2), (2, 3)]),
                            ([(0, 3), (1, 2)], [(0, 2), (1, 3)])]:
            self.assertEqual(PreStrandDiagram(algebra, sd).diff(),
                              PreStrandDiagram(algebra, sd_diff).elt())
        for sd in [[(0, 2), (1, 3)], [(1, 2), (2, 3)]]:
            self.assertEqual(PreStrandDiagram(algebra, sd).diff(), E0)

    def testSignedDiff(self):
        for pmc, idem_size in [(splitPMC(2), 3),
                               (splitPMC(2), 4),
                               (splitPMC(3), 3)]:
            algebra = PreStrandAlgebra(ZZ, pmc, idem_size)
            # Test d^2 = 0.
            for gen in algebra.getGenerators():
                assert gen.diff().diff() == 0

    def testMultiply(self):
        algebra = PreStrandAlgebra(F2, splitPMC(1), 2)
        for sd1, sd2, prod in [([(0, 0), (1, 3)], [(0, 1), (3, 3)],
                                [(0, 1), (1, 3)]),
                               ([(0, 3), (1, 1)], [(1, 2), (3, 3)],
                                [(0, 3), (1, 2)]),
                               ([(1, 3), (0, 0)], [(0, 2), (3, 3)],
                                [(0, 2), (1, 3)])]:
            self.assertEqual(PreStrandDiagram(algebra, sd1) *
                              PreStrandDiagram(algebra, sd2),
                              PreStrandDiagram(algebra, prod).elt())
        for sd1, sd2 in [([(0, 2), (1, 1)], [(1, 3), (2, 2)])]:
            self.assertEqual(PreStrandDiagram(algebra, sd1) *
                              PreStrandDiagram(algebra, sd2), E0)

    # def testSignedMultiply(self):
    #     for pmc, idem_size in [(splitPMC(1), 1),
    #                            (splitPMC(1), 2),
    #                            (splitPMC(2), 1),
    #                            (splitPMC(2), 2),
    #                            (splitPMC(2), 3),
    #                            (splitPMC(2), 4),
    #                            (splitPMC(3), 2)]:
    #         algebra = PreStrandAlgebra(ZZ, pmc, idem_size)
    #         print(algebra)
    #         for gen1 in algebra.getGenerators():
    #             for gen2 in algebra.getGeneratorsForPtIdem(
    #                 l_pt_idem = gen1.right_pt_idem):
    #                 if gen1 * gen2 != E0:
    #                     # Test d(ab) = (da)*b + (-1)^gr(a)*a*(db)
    #                     self.assertEqual(
    #                         (gen1 * gen2).diff(), gen1.diff() * gen2 + \
    #                             algebra.grSign(gen1) * gen1 * gen2.diff())
    #                 for gen3 in algebra.getGeneratorsForPtIdem(
    #                     l_pt_idem = gen2.right_pt_idem):
    #                     if gen2 * gen3 != E0:
    #                         # Tests associativity of multiplication
    #                         self.assertEqual(
    #                             (gen1 * gen2) * gen3, gen1 * (gen2 * gen3))

class SignLinAlgTest(unittest.TestCase):
    def testCreateRowSystem(self):
        sign = SignLinAlg(StrandAlgebra(F2, antipodalPMC(2), idem_size = 2,
                                        mult_one = True))
        sign.createRowSystem()

if __name__ == "__main__":
    unittest.main()
