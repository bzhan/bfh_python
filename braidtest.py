"""Unit test for braid.py"""

from fractions import gcd
from braid import *
import unittest
import time  # benchmarking
import cProfile
import pstats

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

    def testGetCobordismSequence(self):
        for matching, result in [
                ((2, 1), []), ((4, 3, 2, 1), [1]), ((2, 1, 4, 3), [0]),
                ((6, 5, 4, 3, 2, 1), [2, 1]),
                ((6, 3, 2, 5, 4, 1), [1, 1]),
                ((2, 1, 6, 5, 4, 3), [0, 1]),
                ((10, 3, 2, 7, 6, 5, 4, 9, 8, 1), [1, 2, 1, 1])]:
            self.assertEquals(
                BraidCap(matching).getCobordismSequence(), result)

    def testCobordisms(self):
        # 6 strands (genus 2)
        ends3 = reversed(sorted(
                [end for end in BraidCap.ends if len(end) == 6]))
        dcap_len = [len(BraidCap(end).openCap()) for end in ends3]
        self.assertEqual(dcap_len, [2, 1, 2, 1, 1])

        # 8 strands (genus 3)
        ends4 = reversed(sorted(
                [end for end in BraidCap.ends if len(end) == 8]))
        dcap_len = [len(BraidCap(end).getHandlebodyByLocalDA())
                    for end in ends4]
        self.assertEqual(dcap_len, [4, 3, 2, 2, 1, 4, 3, 2, 2, 2, 1, 2, 1, 1])

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
                        if hasattr(cx, "grading"):
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
        to_test = []
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
                    print cx.getGradingInfo()
                    print "Time elapsed (s): ", time.time() - start_time

    def testGetSpecSeq(self):
        to_test = ["3_1", "4_1",
                   "11n_9", "11n_12", "11n_19", # 3-bridge
                   "12n_0475", # four pages
                   "12n_0553", # 4-bridge
        ]
        with open('data/input_12_FL.txt', 'r') as input_file:
            with open('data/output_12_sseq.txt', 'r') as check_file:
                while True:
                    # Read input
                    line = input_file.readline()
                    if len(line) == 0:
                        break
                    cur_br = readBridgePresentation(line)
                    # Read expected output
                    output_header = check_file.readline().split()
                    output_header[0] == cur_br.name
                    num_pages = int(output_header[1])
                    output_lines = [check_file.readline()
                                    for i in range(num_pages)]
                    if cur_br.name in to_test:
                        print "Testing: %s (genus %d)" % \
                            (cur_br.name, cur_br.num_strands / 2 - 1)
                        start_time = time.time()
                        filt_grs = cur_br.getSpecSeq()
                        self.assertEqual(len(filt_grs), len(output_lines))
                        for i in range(len(output_lines)):
                            self.assertEqual(filt_grs[i], [
                                int(s) for s in output_lines[i].split()])
                        print "Time elapsed (s): ", time.time() - start_time

    def testTorusSpecSeq(self):
        def singleTest(m, n):
            start_time = time.time()
            cap = list(range(2*m, 0, -1))
            half_twist = list(range(m-1, 0, -1))
            br = BridgePresentation("T%d_%d" % (m, n), cap, n*half_twist, cap)
            filt_grs = br.getSpecSeq()
            print
            for filt_gr in filt_grs:
                print filt_gr
            print "Time elapsed (s): ", time.time() - start_time

        for n in [1,2,4,5,7,8,10,20,50,100]:
            singleTest(3, n)
        for n in [1,3,5,7,9,11,13,15,17,19]:
            singleTest(4, n)
        for n in [4,6]:
            singleTest(5, n)
        for n in [1,5]:
            singleTest(6, n)

    def testGetSpecSeqProfile(self):
        cProfile.runctx('self.testGetSpecSeq()', globals(), locals(), 'restats')
        p = pstats.Stats('restats')
        p.sort_stats('cumulative').print_stats(50)

if __name__ == "__main__":
    unittest.main()
