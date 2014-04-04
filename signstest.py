"""Unit test for signs.py."""

from signs import *
from grading import DEFAULT_REFINEMENT, lowerRefinement
from pmc import antipodalPMC, linearPMC, splitPMC
import unittest

class AbsZ2GradingTest(unittest.TestCase):
    def testAbsGrading(self):
        def testOnePMC(pmc, test_alg = True):
            alg = pmc.getAlgebra()
            abs_gr = AbsZ2Grading(pmc)
            for gen in alg.getGenerators():
                # Test asserts in getAbsGrading
                abs_gr.getAbsGrading(gen)
            if not test_alg:
                return
            # Test differential and multiplication
            for a in alg.getGenerators():
                for term in a.diff():
                    a_gr, da_gr = [abs_gr.getAbsGrading(gen) for gen in a, term]
                    assert (a_gr - 1) % 2 == da_gr
            for a in alg.getGenerators():
                for b in alg.getGenerators():
                    if a * b != 0:
                        a_gr, b_gr, ab_gr = [abs_gr.getAbsGrading(gen)
                                             for gen in a, b, (a*b).getElt()]
                        assert (a_gr + b_gr) % 2 == ab_gr

        for pmc in [splitPMC(1), splitPMC(2)]:
            testOnePMC(pmc)
        for pmc in [linearPMC(2), antipodalPMC(2), splitPMC(3)]:
            testOnePMC(pmc, test_alg = False)

if __name__ == "__main__":
    unittest.main()
