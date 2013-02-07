"""Unit test for dstructure.py"""

from dstructure import *
from hdiagram import *
import unittest

class DStructureTest(unittest.TestCase):
    def testDStructure(self):
        pmc = splitPMC(1)
        dstr = SimpleDStructure(F2, pmc.getAlgebra())
        genx = SimpleDGenerator(dstr, pmc.idem([0]), "x")
        dstr.addGenerator(genx)
        dstr.addDelta(genx, genx, pmc.sd([(0,2)]), 1)
        self.assertEqual(genx.delta(), 1*TensorGenerator(
                (pmc.sd([(0,2)]), genx), dstr.AtensorM))
        dstr.reindex()
        self.assertEqual(len(dstr), 1)

    def testCommonDStructures(self):
        dsZero = [zeroTypeD(i) for i in range(1,4)]
        for dstr in dsZero:
            dstr.testDelta()
        # Uncomment to see printout of structures
        # print dsZero[0], dsZero[1], dsZero[2]
        dsInf = [infTypeD(i) for i in range(1,4)]
        for dstr in dsInf:
            dstr.testDelta()
        # Uncomment to see printout of structures
        # print dsInf[0], dsInf[1], dsInf[2]
        dsPlat = [platTypeD(i) for i in range(1,4)]
        for dstr in dsPlat:
            dstr.testDelta()
        # Uncomment to see printout of structures
        # print dsPlat[0], dsPlat[1], dsPlat[2]

    def testMorToD(self):
        cx = zeroTypeD(2).morToD(zeroTypeD(2))
        self.assertEqual(len(cx), 4)
        cx2 = infTypeD(2).morToD(infTypeD(2))
        self.assertEqual(len(cx2), 4)
        cx3 = zeroTypeD(2).morToD(infTypeD(2))
        self.assertEqual(len(cx3), 9)
        cx3.simplify()
        self.assertEqual(len(cx3), 1)

        cx = zeroTypeD(1).morToD(zeroTypeD(1))
        cx2 = infTypeD(1).morToD(infTypeD(1))
        cx3 = zeroTypeD(1).morToD(infTypeD(1))
        cx4 = infTypeD(1).morToD(zeroTypeD(1))
        self.assertEqual(cx3.getGradingInfo(), {(1, 2) : 1})
        # Uncomment to see morphism complexes for the torus
        # print cx, cx2, cx3, cx4

    def testCommonDStrMatchDiagram(self):
        dstr1 = infTypeD(2)
        dstr2 = zeroTypeD(2)
        dstr3 = platTypeD(1)
        dstr4 = platTypeD(2)

    def testAdmDStructures(self):
        dstr = zeroTypeDAdm(1)
        dstr2 = zeroTypeDAdm(2)
        dstr3 = zeroTypeDAdm(3)
        dstr3.simplify()
        # Uncomment to see printout of structures
        # print dstr, dstr2, dstr3

    def testDual(self):
        dstr1 = zeroTypeD(1).dual()
        dstr1.testDelta()
        dstr2 = zeroTypeDAdm(1).dual()
        dstr2.testDelta()
        dstr1.checkGrading()
        # Uncomment to see printout of structures
        # print dstr1, dstr2

if __name__ == "__main__":
    unittest.main()
