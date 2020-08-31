import math
import random

import mmh3
from bitarray import bitarray

RANDOM_MAX = 4294967290


class BloomFilter:
    def __init__(self, count, hashes=None, fpp=0.0001):
        self.fpp = fpp
        self.count = count

        self.set_size()
        if hashes is None:
            self.set_hashes()
        else:
            self.hashes = hashes

        self.seed = random.randint(1, RANDOM_MAX)

        self.bitarray = bitarray(self.size)
        self.bitarray.setall(0)

    def check(self, item):
        for i in range(self.hashes):
            hash_seed = self.seed + i
            digest = mmh3.hash(item, hash_seed) % self.size
            if not self.bitarray[digest]:
                return False

        return True

    def add(self, item):
        for i in range(self.hashes):
            hash_seed = self.seed + i
            digest = mmh3.hash(item, hash_seed) % self.size
            self.bitarray[digest] = True

    def set_size(self):
        self.size = self.calc_size(self.count, self.fpp)

    def set_hashes(self):
        self.hashes = self.calc_hashes(self.size, self.count)

    @classmethod
    def calc_size(self, count, fpp):
        size = round(-(count * math.log(fpp))/(math.log(2) ** 2), 0)
        return int(size)

    @classmethod
    def calc_hashes(self, size, count):
        hashes = round((size/count) * (math.log(2)), 0)
        return int(hashes)
