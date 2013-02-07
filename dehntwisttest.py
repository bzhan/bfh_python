"""Unit test for dehntwist.py"""

from dehntwist import *
import unittest

class DehnTwistTest(unittest.TestCase):
    def testDehnTwist(self):
        twist = DehnTwist(3, 1, POS)
        twist_dd = twist.getDDStructure() 

if __name__ == "__main__":
    unittest.main()
