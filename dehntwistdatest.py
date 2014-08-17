"""Unit test for dehntwistda.py"""

from ddstructure import identityDD
from dehntwist import AntiBraid
from dehntwistda import *
import unittest

class AntiBraidDATest(unittest.TestCase):
    def testLocalDA(self):
        da = AntiBraidDA(2, 1)
        self.assertTrue(da.local_da.testDelta())

    def testAntiBraidDA(self):
        for genus, c_pair in [(2, 1), (2, 2)]:
            da = AntiBraidDA(genus, c_pair)
            self.assertTrue(da.toSimpleDAStructure().testDelta())

    def testAntiBraidAgreesWithDD(self):
        for genus, c_pair in [(2, 1), (2, 2), (3, 2)]:
            dastr = AntiBraidDA(genus, c_pair)
            ddstr = dastr.tensorDD(identityDD(dastr.pmc))
            ddstr.simplify()
            ori_ddstr = AntiBraid(genus, c_pair).getDDStructure()
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

    def testLocalDAShort(self):
        for genus, c_pair in [(2, 0), (2, 3)]:
            da = AntiBraidDA(genus, c_pair)
            self.assertTrue(da.local_da.testDelta())

    def testAntiBraidDAShort(self):
        for genus, c_pair in [(2, 0), (2, 3)]:
            da = AntiBraidDA(genus, c_pair)
            self.assertTrue(da.toSimpleDAStructure().testDelta())

    def testAntiBraidAgreesWithDDShort(self):
        for genus, c_pair in [(2, 0), (3, 0), (2, 3)]:
            dastr = AntiBraidDA(genus, c_pair)
            ddstr = dastr.tensorDD(identityDD(dastr.pmc))
            ddstr.simplify()
            ori_ddstr = AntiBraid(genus, c_pair).getDDStructure()
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

if __name__ == "__main__":
    unittest.main()
