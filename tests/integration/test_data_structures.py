import random
import timeit
import unittest
from random import shuffle

from pde_utils.classes.data_structures import (
        BloomFilter, BloomFilterStack, CountMinSketch)

NUCLEOTIDE_MAP = {1: "A", 2: "T", 3: "C", 4: "G"}


def generate_random_nucleotide_seq():
    nucleotides = []
    for i in range(20):
        nucleotides.append(NUCLEOTIDE_MAP[random.randint(1, 4)])

    return "".join(nucleotides)


class DataStructuresTest(unittest.TestCase):
    def setUp(self):
        self.sequences_present = []
        self.sequences_absent = []

        for i in range(280):
            self.sequences_present.append(generate_random_nucleotide_seq())

        absent_seqs = 0
        while absent_seqs < 5:
            random_seq = generate_random_nucleotide_seq()

            if random_seq not in self.sequences_present:
                absent_seqs += 1
                self.sequences_absent.append(random_seq)

    def test_bloom_filter_1(self):
        start = timeit.default_timer()

        for i in range(1000):
            with self.subTest(trial=i+1):
                bloomf = BloomFilter(len(self.sequences_present))

                for seq in self.sequences_present:
                    bloomf.add(seq)

                shuffle(self.sequences_present)
                shuffle(self.sequences_absent)

                test_seqs = self.sequences_present[:5] + self.sequences_absent
                shuffle(test_seqs)

                false_positives = 0
                for seq in test_seqs:
                    if not bloomf.check(seq):
                        self.assertTrue(seq not in self.sequences_present)
                    else:
                        if seq in self.sequences_absent:
                            false_positives += 1

                self.assertTrue(false_positives <= 1)

        stop = timeit.default_timer()
        print(f"Bloom filter test 1: {stop-start}")

    def test_bloom_filter_2(self):
        start = timeit.default_timer()

        for i in range(1000):
            with self.subTest(trial=i+1):
                bloomf = BloomFilter(10, fpp=1)

                for seq in self.sequences_present:
                    bloomf.add(seq)

                shuffle(self.sequences_present)
                shuffle(self.sequences_absent)

                test_seqs = self.sequences_present[:5] + self.sequences_absent
                shuffle(test_seqs)

                false_positives = 0
                for seq in test_seqs:
                    if not bloomf.check(seq):
                        self.assertTrue(seq not in self.sequences_present)
                    else:
                        if seq in self.sequences_absent:
                            false_positives += 1

                self.assertTrue(false_positives == 5)

        stop = timeit.default_timer()
        print(f"Bloom filter test 2: {stop-start}")

    def test_bloom_filter_stack_1(self):
        start = timeit.default_timer()

        target_sequence = generate_random_nucleotide_seq()

        for i in range(1000):
            with self.subTest(trial=i+1):
                bfstack = BloomFilterStack(len(self.sequences_present) + 10,
                                           10)

                for seq in self.sequences_present:
                    bfstack.add(seq)

                for i in range(10):
                    bfstack.add(target_sequence)

                self.assertTrue(bfstack.check(target_sequence))

        stop = timeit.default_timer()
        print(f"Bloom filter stack test 1: {stop-start}")

    def test_count_min_sketch_1(self):
        start = timeit.default_timer()

        target_sequence = generate_random_nucleotide_seq()

        count_overshoots = 0
        for i in range(1000):
            with self.subTest(trial=i+1):
                cmsketch = CountMinSketch(len(self.sequences_present) + 10)

                for seq in self.sequences_present:
                    cmsketch.add(seq)

                for i in range(10):
                    cmsketch.add(target_sequence)

                target_count = cmsketch.check(target_sequence)

                self.assertTrue(target_count >= 10)

                if target_count != 10:
                    count_overshoots += 1

        if count_overshoots > 0:
            print("TEST WARNING: CountMinSketch inaccurately counted"
                  f"{round(count_overshoots/1000, 3)*100}% "
                  "of trial iterations")

        stop = timeit.default_timer()
        print(f"Count min-sketch test 1: {stop-start}")


if __name__ == "__main__":
    unittest.main()
