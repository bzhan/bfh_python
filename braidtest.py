"""Unit test for braid.py"""

from braid import *
import unittest

class BraidTest(unittest.TestCase):
    def testGetArcslide(self):
        br2 = Braid(6)
        pos_size = [1,2,2,1,4]
        for i in range(1, 6):
            self.assertTrue(len(br2.getArcslides(i)), pos_size[i-1])
            self.assertTrue(len(br2.getArcslides(-i)), pos_size[i-1])
        self.assertTrue(len(br2.getArcslides(range(1,6))), sum(pos_size))

class BraidCapTest(unittest.TestCase):
    def testBraidCap(self):
        dcap_len = [len(BraidCap(end3).getHandlebody())
                    for end3 in BraidCap.ends3]
        self.assertEqual(dcap_len, [2,1,2,1,1])

    def testPlatTypeD2(self):
        self.assertEqual(len(platTypeD2(3, True)), 1)

    def testGenus2Algebra(self):
        dstrs = dict()
        for end3 in BraidCap.ends3:
            dstrs[end3] = BraidCap(end3).getHandlebody()
        algs = dict()
        total_gen = 0
        for end3, dstr in dstrs.items():
            algs[end3] = dict()
            for end3to, dstrto in dstrs.items():
                algs[end3][end3to] = dstr.morToD(dstrto)
                algs[end3][end3to].simplify()
                total_gen += len(algs[end3][end3to])
        self.assertEqual(total_gen, 52)

class HFTest(unittest.TestCase):
    def testHFPretzel(self):
        std_cap = [6,3,2,5,4,1]
        br = BridgePresentation("pretzel_-2_3_5",
                                 std_cap, 5*[1]+3*[3]+2*[-5], std_cap)
        self.assertEqual(len(br.getHF()), 1)

    def testHFFromFile(self):
        to_test = ["3_1", "4_1",
                   # "12n_0210", # 1*[3 2 2 2 2 2]
                   # "12n_0292", # 1*[1 0 2 2 0 0 2 2]
                   # "11n_6", "11n_9", "11n_24", # 3-bridge
                   # "11a_14", "12n_0055", "12n_0056", # 4-bridge
                   ]
        with open('data/input_12_FL.txt', 'r') as input_file:
            with open('data/output_12.txt', 'r') as check_file:
                while True:
                    line = input_file.readline()
                    if len(line) == 0:
                        break
                    expected_hf = int(check_file.readline().split()[1])
                    cur_br = readBridgePresentation(line)
                    if cur_br.name in to_test:
                        cx = cur_br.getHF()
                        self.assertEqual(len(cx), expected_hf)
                        print cur_br.name, cx.getGradingInfo()

if __name__ == "__main__":
    unittest.main()
