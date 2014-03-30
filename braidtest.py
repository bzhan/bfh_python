"""Unit test for braid.py"""

from fractions import gcd
from braid import *
import unittest
import time  # benchmarking

class BraidTest(unittest.TestCase):
    def testGetArcslide(self):
        br2 = Braid(6)
        pos_size = [1,2,2,1,4]
        for i in range(1, 6):
            self.assertTrue(len(br2.getArcslides(i)), pos_size[i-1])
            self.assertTrue(len(br2.getArcslides(-i)), pos_size[i-1])
        self.assertTrue(len(br2.getArcslides(range(1, 6))), sum(pos_size))

class BraidCapTest(unittest.TestCase):
    def testVerifyEnds(self):
        BraidCap.verifyEnds()

    def testBraidCap3(self):
        ends3 = reversed(sorted(
                [end for end in BraidCap.ends if len(end) == 6]))
        dcap_len = [len(BraidCap(end).getHandlebody()) for end in ends3]
        self.assertEqual(dcap_len, [2, 1, 2, 1, 1])

    def testBraidCap3LocalDA(self):
        ends3 = reversed(sorted(
                [end for end in BraidCap.ends if len(end) == 6]))
        dcap_len = [len(BraidCap(end).getHandlebodyByLocalDA())
                    for end in ends3]
        self.assertEqual(dcap_len, [2, 1, 2, 1, 1])

    def testBraidCap4LocalDA(self):
        ends4 = reversed(sorted(
                [end for end in BraidCap.ends if len(end) == 8]))
        dcap_len = [len(BraidCap(end).getHandlebodyByLocalDA())
                    for end in ends4]
        self.assertEqual(dcap_len, [4, 3, 2, 2, 1, 4, 3, 2, 2, 2, 1, 2, 1, 1])

    def testPlatTypeD2(self):
        self.assertEqual(len(platTypeD2(3, True)), 1)

    def testGenus2Algebra(self):
        dstrs = dict()
        for end3 in [end for end in BraidCap.ends if len(end) == 6]:
            dstrs[end3] = BraidCap(end3).getHandlebodyByLocalDA()
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
        # Test both methods of finding HF
        self.assertEqual(len(br.getHF()), 1)
        self.assertEqual(len(br.getHFByLocalDA()), 1)

    def testHFFromFile(self):
        to_test = ["3_1", "4_1",
                   "12n_0210", # 1*[3 2 2 2 2 2]
                   "12n_0292", # 1*[1 0 2 2 0 0 2 2]
                   "11n_6", "11n_9", "11n_24", # 3-bridge
                   "11a_14", "12n_0055", "12n_0056", # 4-bridge
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
                        print "Testing: %s (genus %d)" % \
                            (cur_br.name, cur_br.num_strands / 2 - 1)
                        start_time = time.time()
                        cx = cur_br.getHFByLocalDA()
                        # cx = cur_br.getHF()
                        print cx.getGradingInfo()
                        self.assertEqual(len(cx), expected_hf)
                        print "Time elapsed (s): ", time.time() - start_time

    def testTorus(self):
        def singleTest(m, n):
            start_time = time.time()
            cap = list(range(2*m, 0, -1))
            half_twist = list(range(m-1, 0, -1))
            br = BridgePresentation("T%d_%d" % (m, n), cap, n*half_twist, cap)
            cx = br.getHFByLocalDA()
            print br, len(cx)
            print "Time elapsed (s): ", time.time() - start_time

        for n in [1,2,4,5,7,8,10,20,50,100]:
            singleTest(3, n)
        for n in [1,3,5,7,9,11,13,15,17,19]:
            singleTest(4, n)
        for n in [1,2,3,4,6]:
            singleTest(5, n)
        for n in [1,5]:
            singleTest(6, n)

    def testGenus5FromFile(self):
        # Empty means test all
        to_test = ["14n_6302"]
        with open('data/input_14_FL.txt', 'r') as input_file:
            while True:
                line = input_file.readline()
                if len(line) == 0:
                    break
                cur_br = readBridgePresentation(line)
                if len(to_test) == 0 or cur_br.name in to_test:
                    print "Testing:", cur_br.name
                    start_time = time.time()
                    cx = cur_br.getHFByLocalDA()
                    print len(cx)
                    print "Time elapsed (s): ", time.time() - start_time

if __name__ == "__main__":
    # Can use this when running profiler.
    # suite = unittest.TestLoader().loadTestsFromTestCase(HFTest)
    # unittest.TextTestRunner(verbosity=2).run(suite)    
    unittest.main()
