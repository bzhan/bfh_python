"""Unit test for ddstructure.py"""

from ddstructure import *
from dstructure import infTypeD, zeroTypeD
from pmc import PMC
from pmc import linearPMC, splitPMC
from utility import DEFAULT_GRADING, SMALL_GRADING
import unittest

class DDStructureTest(unittest.TestCase):
    def testCommonIdentityDD(self):
        id2 = identityDD(splitPMC(2))
        self.assertEqual(len(id2), 6)
        self.assertTrue(id2.testDelta())

    def testMorToD(self):
        id2 = identityDD(splitPMC(2))
        d = infTypeD(2)
        d2 = id2.morToD(d)
        self.assertTrue(d2.testDelta())
        d2.simplify()
        self.assertEqual(len(d2), 1)
        d3 = zeroTypeD(2)
        d4 = id2.morToD(d3)
        self.assertTrue(d4.testDelta())
        d4.simplify()
        self.assertEqual(len(d4), 1)

    def testToDStructure(self):
        ddstr = identityDD(linearPMC(1))
        dstr = ddstr.toDStructure()
        hochchild = dstr.morToD(dstr)
        self.assertEqual(len(hochchild), 18)
        hochchild.simplify(find_homology_basis = True)
        self.assertEqual(len(hochchild), 4)
        meaning_len = [len(gen.prev_meaning)
                       for gen in hochchild.getGenerators()]
        self.assertEqual(sorted(meaning_len), [1,1,1,2])

    def testIdentityMatchDiagram(self):
        pmc_to_test = [splitPMC(2), PMC([(0,2),(1,6),(3,5),(4,7)])]
        for pmc in pmc_to_test:
            for idem_size in [0,1,2,3,4]:
                ddstr = identityDD(pmc, idem_size)
                # Special check for the identity diagram: all gradings should be
                # zero
                if DEFAULT_GRADING == SMALL_GRADING:
                    for gen in ddstr.generators:
                        self.assertEqual(ddstr.grading[gen],
                                         ddstr.gr_set.zero())

    def testDual(self):
        pmc = PMC([(0,2),(1,6),(3,5),(4,7)])
        ddstr = identityDD(pmc)
        ddstr_dual = ddstr.dual()
        self.assertEqual(ddstr_dual.algebra1, pmc.getAlgebra())
        self.assertEqual(ddstr_dual.algebra2, pmc.opp().getAlgebra())

    def testHochschildCohomology(self):
        # Check HH^*(A)
        # Rank 4 is right for HH^* of A in genus 1 case
        pmc = splitPMC(1)
        ddstr = identityDD(pmc)
        cx = ddstr.hochschildCochains()
        cx.simplify()
        self.assertEqual(len(cx), 4)

        # Rank 1 in the extremal strands grading
        pmc = splitPMC(2)
        ddstr = identityDD(pmc, 0)
        cx = ddstr.hochschildCochains()
        cx.simplify()
        self.assertEqual(len(cx), 1)

        # Rank should be indepdendent of PMC (equals 16)
        for pmc in [splitPMC(2), linearPMC(2), PMC([(0,2),(1,6),(3,5),(4,7)])]:
            ddstr = identityDD(pmc)
            cx = ddstr.hochschildCochains()
            cx.simplify()
            self.assertEqual(len(cx), 16)

class DDMorphismTest(unittest.TestCase):
    def testMappingCone(self):
        pmc = splitPMC(1)
        alg = pmc.getAlgebra()
        tensor_alg = TensorDGAlgebra((alg, alg))
        ddstr1 = SimpleDDStructure(F2, alg, alg)
        ddstr2 = SimpleDDStructure(F2, alg, alg)

        gens = dict()
        gens["x"] = SimpleDDGenerator(ddstr1, pmc.idem([0]), pmc.idem([0]), "x")
        gens["y"] = SimpleDDGenerator(ddstr2, pmc.idem([1]), pmc.idem([1]), "y")
        ddstr1.addGenerator(gens["x"])
        ddstr2.addGenerator(gens["y"])
        ddstr1.addDelta(
            gens["x"], gens["x"], pmc.sd([(0, 2)]), pmc.sd([(0, 2)]), 1)
        ddstr2.addDelta(
            gens["y"], gens["y"], pmc.sd([(1, 3)]), pmc.sd([(1, 3)]), 1)

        morphism_cx = MorDDtoDDComplex(F2, ddstr1, ddstr2)
        morphism = 1*MorDDtoDDGenerator(
            morphism_cx, gens["x"],
            TensorGenerator((pmc.sd([(0, 1)]), pmc.sd([(0, 1)])), tensor_alg),
            gens["y"])
        morphism += 1*MorDDtoDDGenerator(
            morphism_cx, gens["x"],
            TensorGenerator((pmc.sd([(2, 3)]), pmc.sd([(2, 3)])), tensor_alg),
            gens["y"])
        # Need diff of morphism to be zero for the mapping cone to satisfy
        # d^2 = 0.
        self.assertEqual(morphism.diff(), 0)
        mapping_cone = morphism_cx.getMappingCone(morphism)

        self.assertEqual(len(mapping_cone), 2)
        mc_gens = dict()
        for gen in mapping_cone.getGenerators():
            mc_gens[gen.name] = gen
        self.assertTrue("S_x" in mc_gens and "T_y" in mc_gens)
        self.assertEqual(len(mapping_cone.delta(mc_gens["S_x"])), 3)
        self.assertEqual(len(mapping_cone.delta(mc_gens["T_y"])), 1)
        self.assertTrue(mapping_cone.testDelta())

    def testMultiply(self):
        pmc = splitPMC(2)
        alg = pmc.getAlgebra()
        tensor_alg = TensorDGAlgebra((alg, alg))
        ddstr = SimpleDDStructure(F2, alg, alg)
        gens = dict()
        gens["x"] = SimpleDDGenerator(ddstr, pmc.idem([0]), pmc.idem([0]), "x")
        gens["y"] = SimpleDDGenerator(ddstr, pmc.idem([1]), pmc.idem([1]), "y")
        ddstr.addGenerator(gens["x"])
        ddstr.addGenerator(gens["y"])
        
        morphism_cx = MorDDtoDDComplex(F2, ddstr, ddstr)
        morphism1 = 1*MorDDtoDDGenerator(
            morphism_cx, gens["x"],
            TensorGenerator((pmc.sd([(0, 1)]), pmc.sd([(0, 1)])), tensor_alg),
            gens["y"])
        morphism1 += 1*MorDDtoDDGenerator(
            morphism_cx, gens["x"],
            TensorGenerator((pmc.sd([(2, 3)]), pmc.sd([(2, 3)])), tensor_alg),
            gens["y"])
        morphism2 = 1*MorDDtoDDGenerator(
            morphism_cx, gens["y"],
            TensorGenerator((pmc.sd([(1, 2)]), pmc.sd([(1, 2)])), tensor_alg),
            gens["x"])
        self.assertEqual(len(morphism1 * morphism2), 1)
        self.assertEqual(len(morphism2 * morphism1), 1)
        self.assertNotEqual(morphism1 * morphism2, morphism2 * morphism1)

if __name__ == "__main__":
    unittest.main()
