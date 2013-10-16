"""Unit test for utility.py"""

from utility import *
import unittest

class ToListTest(unittest.TestCase):
    def testToList(self):
        """Testing various cases for the tolist function."""
        self.assertEqual(tolist(3), [3])
        self.assertEqual(tolist("3"), ["3"])
        self.assertEqual(tolist(()), [()])
        self.assertEqual(tolist([]), [])
        self.assertEqual(tolist([3, 4]), [3, 4])

class SubsetTest(unittest.TestCase):
    def testSubset(self):
        self.assertEqual(sorted(subset([])), [()])
        self.assertEqual(sorted(subset([0])), [(), (0,)])
        self.assertEqual(sorted(subset(['a', 'b'])),
                         [(), ('a',), ('a', 'b'), ('b',)])

class SummableDictTest(unittest.TestCase):
    def setUp(self):
        self.a = SummableDict({"a" : 5, "b" : 3, "d" : 3})
        self.b = SummableDict({"a" : 5, "c" : 4, "d" : -3})
        self.c = SummableDict({})
        self.sumab = SummableDict({"a" : 10, "b" : 3, "c" : 4})
        self.a2 = SummableDict({"a" : 10, "b" : 6, "d" : 6})
        self.refa = self.a.copy()
        self.refb = self.b.copy()

    def testDictAdd(self):
        """Testing __add__ operation in SummableDict."""
        self.assertEqual(self.a + self.b, self.sumab)
        self.assertEqual(self.a + self.a, self.a2)
        self.assertEqual(self.a, self.refa)
        self.assertEqual(self.b, self.refb)

    def testDictIAdd(self):
        """Testing __iadd__ and accumulate operation in SummableDict."""
        self.a += self.c
        self.assertEqual(self.a, self.refa)
        self.assertEqual(self.a.accumulate([self.b, self.c]), self.sumab)
        self.a = self.refa
        self.a += self.a
        self.assertEqual(self.a, self.a2)
        self.assertEqual(self.b, self.refb)
        self.assertEqual(self.c, {})

    def testDictMul(self):
        """Testing __mul__ and __rmul__ operation in SummableDict."""
        self.assertEqual(self.a * 2, self.a2)
        self.assertEqual(2 * self.a, self.a2)
        self.assertEqual(0 * self.a, {})

    def testDictEqual(self):
        self.assertTrue(self.a == self.refa)
        self.assertTrue(self.a != self.b)
        self.assertTrue(self.a != 0)
        self.assertTrue(self.c == 0)

    def testDictTranslate(self):
        d1 = SummableDict({"a" : 5, "b" : 4, "c" : 3})
        d2 = d1.translateKey(dict({"a" : "aa", "b" : "bb", "c" : "cc"}))
        self.assertTrue(d2 != d1)
        self.assertEqual(d2, SummableDict({"aa" : 5, "bb" : 4, "cc" : 3}))

class FiniteRingTest(unittest.TestCase):
    def testFiniteRing(self):
        self.assertTrue(F2.one != 0)
        self.assertEqual(F2.one + F2.one, 0)

if __name__ == "__main__":
    unittest.main()
