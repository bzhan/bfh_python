"""Unit test for dehntwistda.py"""

from ddstructure import identityDD
from dehntwist import AntiBraid, DehnSurgery
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

class DehnSurgeryDATest(unittest.TestCase):
    def testLocalDA(self):
        for genus, c_pair, orientation in [(2, 1, NEG), (2, 1, POS)]:
            ds = DehnSurgeryDA(genus, c_pair, orientation)
            morphism = ds.getLocalMorphism()
            self.assertEqual(ds.getLocalMorphism().diff(), 0)
            self.assertTrue(ds.getLocalMappingCone().testDelta())

    def testDehnSurgeryDA(self):
        for genus, c_pair, orientation in [(2, 1, NEG), (2, 1, POS)]:
            ds = DehnSurgeryDA(genus, c_pair, orientation)
            self.assertTrue(
                ds.getMappingCone().toSimpleDAStructure().testDelta())

    def testDehnSurgeryAgreesWithDD(self):
        for genus, c_pair, orientation in [(2, 1, NEG), (2, 1, POS)]:
            ds = DehnSurgeryDA(genus, c_pair, orientation)
            dastr = ds.getMappingCone()
            ddstr = dastr.tensorDD(identityDD(dastr.algebra1.pmc))
            ddstr.simplify(
                cancellation_constraint = lambda x, y: (
                    x.filtration == y.filtration))

            ori_mor = DehnSurgery(genus, c_pair, orientation).getMorphism()
            ori_mor_cx = ori_mor.getElt().parent
            ori_ddstr = ori_mor_cx.getMappingCone(ori_mor)
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

    def testLocalDAShort(self):
        for genus, c_pair, orientation in [
                (2, 0, NEG), (2, 0, POS), (2, 3, POS)]:
            ds = DehnSurgeryDA(genus, c_pair, orientation)
            morphism = ds.getLocalMorphism()
            self.assertEqual(ds.getLocalMorphism().diff(), 0)
            self.assertTrue(ds.getLocalMappingCone().testDelta())

    def testDehnSurgeryDAShort(self):
        for genus, c_pair, orientation in [
                (2, 0, NEG), (2, 0, POS), (2, 3, POS)]:
            ds = DehnSurgeryDA(genus, c_pair, orientation)
            self.assertTrue(
                ds.getMappingCone().toSimpleDAStructure().testDelta())

    def testDehnSurgeryAgreesWithDDShort(self):
        for genus, c_pair, orientation in [
                (2, 0, NEG), (2, 0, POS), (2, 3, POS)]:
            ds = DehnSurgeryDA(genus, c_pair, orientation)
            dastr = ds.getMappingCone()
            ddstr = dastr.tensorDD(identityDD(dastr.algebra1.pmc))
            ddstr.simplify(
                cancellation_constraint = lambda x, y: (
                    x.filtration == y.filtration))

            ori_mor = DehnSurgery(genus, c_pair, orientation).getMorphism()
            ori_mor_cx = ori_mor.getElt().parent
            ori_ddstr = ori_mor_cx.getMappingCone(ori_mor)
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

if __name__ == "__main__":
    unittest.main()
