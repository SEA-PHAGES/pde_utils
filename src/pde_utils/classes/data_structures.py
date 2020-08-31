import math
import random

import mmh3
from bitarray import bitarray

RANDOM_MAX = 4294967290


class BloomFilter:
    def __init__(self, item_count, hash_count=None, fpp=0.0001):
        self.fpp = fpp
        self.count = item_count

        self.set_size()
        if hash_count is None:
            self.set_hash_count()
        else:
            self.hash_count = hash_count
        self.set_salts()

        self.bitarray = bitarray(self.size)
        self.bitarray.setall(0)

    def check(self, item):
        for salt in self.salts:
            digest = mmh3.hash(item, salt) % self.size
            if not self.bitarray[digest]:
                return False

        return True

    def add(self, item):
        for salt in self.salts:
            digest = mmh3.hash(item, salt) % self.size
            self.bitarray[digest] = True

    def set_size(self):
        self.size = self.calc_size(self.count, self.fpp)

    def set_hash_count(self):
        self.hash_count = self.calc_hash_count(self.size, self.count)

    def set_salts(self):
        salts = []
        for i in range(self.hash_count):
            salts.append(random.randint(1, RANDOM_MAX))

        self.salts = salts

    @classmethod
    def calc_size(self, count, fpp):
        size = round(-(count * math.log(fpp))/(math.log(2) ** 2), 0)
        return int(size)

    @classmethod
    def calc_hash_count(self, size, count):
        hash_count = round((size/count) * (math.log(2)), 0)
        return int(hash_count)


class BloomFilterStack:
    def __init__(self, item_count, stack_size, hash_count=None, fpp=0.0001):
        self.blooms = []
        for i in range(stack_size):
            self.blooms.append(BloomFilter(item_count, hash_count=hash_count,
                                           fpp=fpp))

    def check(self, item):
        return self.blooms[-1].check(item)

    def add(self, item):
        for bfilter in self.blooms:
            if not bfilter.check(item):
                bfilter.add(item)
                break


class CountMinSketch(BloomFilter):
    def __init__(self, item_count, hash_count=None, fpp=0.0001):
        self.fpp = fpp
        self.count = item_count

        self.set_size()
        if hash_count is None:
            self.set_hash_count()
        else:
            self.hash_count = hash_count
        self.set_salts()

        self.numarray = [0] * self.size

    def check(self, item):
        count = None
        for salt in self.salts:
            digest = mmh3.hash(item, salt) % self.size
            curr_bucket = self.numarray[digest]

            if count is None:
                count = curr_bucket
            elif count > curr_bucket:
                count = curr_bucket

        return count

    def add(self, item):
        for salt in self.salts:
            digest = mmh3.hash(item, salt) % self.size
            self.numarray[digest] += 1
