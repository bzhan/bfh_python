"""Unit test for cobordismtest.py"""

from cobordism import *
from dstructure import platTypeD
import unittest

class CobordismTest(unittest.TestCase):
    def testCobordism(self):
        for genus, c_pair, side in [
                (2, 1, RIGHT), (2, 2, RIGHT), (3, 1, RIGHT),
                (3, 2, RIGHT), (3, 3, RIGHT), (3, 4, RIGHT)]:
            c = Cobordism(genus, c_pair, side)
            c.getDDStructure()  # verifies testDelta()

    def testCobordismShort(self):
        for genus, c_pair, side in [
                (2, 0, RIGHT), (2, 3, RIGHT), (3, 0, RIGHT),
                (3, 5, RIGHT)]:
            c = Cobordism(genus, c_pair, side)
            c.getDDStructure()  # verifies testDelta()

    def testMorToPlatRight(self):
        for genus, c_pair, side in [
                (2, 1, RIGHT), (2, 3, RIGHT), (3, 1, RIGHT),
                (3, 3, RIGHT), (3, 5, RIGHT)]:
            dd = Cobordism(genus, c_pair, side).getDDStructure()
            plat_d = platTypeD(genus-1)
            plat_d2 = dd.morToD(plat_d)
            plat_d2.simplify()
            self.assertTrue(plat_d2.compareDStructures(platTypeD(genus)))

        for genus, c_pair, side in [
                (2, 0, RIGHT), (2, 2, RIGHT), (3, 0, RIGHT),
                (3, 2, RIGHT), (3, 4, RIGHT)]:
            dd = Cobordism(genus, c_pair, side).getDDStructure()
            plat_d = platTypeD(genus-1).dual()
            plat_d2 = dd.morToD(plat_d)
            plat_d2.simplify()
            self.assertTrue(plat_d2.compareDStructures(platTypeD(genus).dual()))

    def testMorToPlatLeft(self):
        for genus, c_pair, side in [
                (2, 0, LEFT), (2, 2, LEFT), (3, 0, LEFT),
                (3, 2, LEFT), (3, 4, LEFT)]:
            dd = Cobordism(genus, c_pair, side).getDDStructure()
            plat_d = platTypeD(genus)
            plat_d2 = dd.morToD(plat_d)
            plat_d2.simplify()
            self.assertTrue(plat_d2.compareDStructures(platTypeD(genus-1)))

        for genus, c_pair, side in [
                (2, 1, LEFT), (2, 3, LEFT), (3, 1, LEFT),
                (3, 3, LEFT), (3, 5, LEFT)]:
            dd = Cobordism(genus, c_pair, side).getDDStructure()
            plat_d = platTypeD(genus).dual()
            plat_d2 = dd.morToD(plat_d)
            plat_d2.simplify()
            self.assertTrue(
                plat_d2.compareDStructures(platTypeD(genus-1).dual()))

if __name__ == "__main__":
    unittest.main()
