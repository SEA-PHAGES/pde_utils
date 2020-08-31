import math
import random

import mmh3
from bitarray import bitarray

RANDOM_MAX = 4294967290

class BloomFilter:
    def  __init__(self, count, hashes=None, fpp=0.0001):
        self.fpp = fpp

        self.size = self.calc_size(count)
        if hashes is None:
            self.calc_hashes(size, count)

        self.seed = random.randint(1, RANDOM_MAX)

        self.bitarray = bitarray(self.size)  
        self.bitarray.setall(0)


    def check(self, item):
        for i in range(self.hash_count):
            hash_seed = self.seed + i
            digest = mmh3.hash(item, hash_seed)
            if not self.bitarray[digest]:
                return False

        return True


    def add(self, item):
        for i in range(self.hashes):
            digest = mmh3.hash(item, self.seed) % self.size
            self.bit_array[digest] = True
 

    @classmethod
    def calc_size(self, count, fpp):
        size = round(-(count * math.log(fpp))/(math.log(2) ** 2), 0)
        return int(size)


    @classmethod
    def calc_hashes(self, size, count):
        hashes = round((size/count) * (math.log(2)), 0)
        return int(hashes)
