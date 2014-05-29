"""Unit test for dehntwist.py"""

from dehntwist import *
import unittest

class DehnTwistTest(unittest.TestCase):
    def testDehnTwist(self):
        twist = DehnTwist(3, 1, POS)
        twist_dd = twist.getDDStructure()

class AntiBraidTest(unittest.TestCase):
    def testAntiBraid(self):
        for genus, c_pair in [(3, 1), (3, 0), (3, 5)]:
            antibraid = AntiBraid(genus, c_pair)
            antibraid_dd = antibraid.getDDStructure()

class DehnSurgeryTest(unittest.TestCase):
    def testDehnSurgery(self):
        for genus, c_pair, orientation in [(3, 1, NEG), (3, 2, POS),
                                           (3, 0, NEG), (3, 5, POS)]:
            surgery = DehnSurgery(genus, c_pair, orientation)
            morphism = surgery.getMorphism()
            self.assertEqual(morphism.diff(), 0)

    def testDehnSurgeryMappingCone(self):
        surgery = DehnSurgery(2, 1, NEG)
        morphism = surgery.getMorphism()
        morphism_cx = morphism.getElt().parent
        mapping_cone = morphism_cx.getMappingCone(morphism)
        # Six generators in identity, four in anti-braid
        self.assertEqual(len(mapping_cone), 10)
        self.assertTrue(mapping_cone.testDelta())

if __name__ == "__main__":
    unittest.main()
