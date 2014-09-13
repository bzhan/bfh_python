"""Unit test for cobordismda.py"""

from cobordismda import *
from ddstructure import identityDD
import unittest

class CobordismDATest(unittest.TestCase):
    def testLeftCobordismDA(self):
        for genus, c_pair in [(2, 0), (2, 1), (2, 2), (2, 3),
                              (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5)]:
            c = Cobordism(genus, c_pair, LEFT)
            c_da = CobordismDA(c)
            dastr = c_da.toSimpleDAStructure()
            self.assertTrue(dastr.testDelta())
            ddstr = dastr.tensorDD(identityDD(c.end_pmc))
            ori_ddstr = c.getDDStructure()
            self.assertTrue(ddstr.compareDDStructures(ori_ddstr))

if __name__ == "__main__":
    unittest.main()
