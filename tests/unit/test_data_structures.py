import unittest
import random
from random import shuffle

from pde_utils.classes.data_structures import BloomFilter

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

        for i in range(100):
            self.sequences_present.append(generate_random_nucleotide_seq())

        absent_seqs = 0
        while absent_seqs < 5:
            random_seq = generate_random_nucleotide_seq()

            if random_seq not in self.sequences_present:
                absent_seqs += 1
                self.sequences_absent.append(random_seq)

    def test_bloom_filter_1(self):
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

    def test_bloom_filter_2(self):
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


if __name__ == "__main__":
    unittest.main()
