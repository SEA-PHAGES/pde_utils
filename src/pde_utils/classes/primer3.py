import re

from Bio import SeqUtils
from Bio.Seq import Seq
from primer3 import (
        calcTm, calcHairpin, calcHomodimer, calcHeterodimer)
from primer3.bindings import calcEndStability


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
class Clamp:
    def __init__(self, oligomer, target):
        self.oligomer = oligomer
        self.target = target

        thermo = calcEndStability(self.oligomer, self.target)

        self.tm = thermo.tm
        self.dg = thermo.dg
        self.dh = thermo.dh
        self.ds = thermo.ds


class Heterodimer:
    def __init__(self, oligomer, target):
        self.oligomer = oligomer
        self.target = target

        thermo = calcHeterodimer(self.oligomer, self.target,
                                 output_structure=True)

        self.tm = thermo.tm
        self.dg = thermo.dg
        self.dh = thermo.dh
        self.ds = thermo.ds

        self.structure = thermo.ascii_structure
        self.structure_lines = thermo.ascii_structure_lines


class Homodimer:
    def __init__(self, oligomer):
        self.oligomer = oligomer

        thermo = calcHomodimer(self.oligomer, output_structure=True)

        self.tm = thermo.tm
        self.dg = thermo.dg
        self.dh = thermo.dh
        self.ds = thermo.ds

        self.structure = thermo.ascii_structure
        self.structure_lines = thermo.ascii_structure_lines


class Hairpin:
    def __init__(self, oligomer):
        self.oligomer = oligomer

        thermo = calcHairpin(self.oligomer, output_structure=True)

        self.tm = thermo.tm
        self.dg = thermo.dg
        self.dh = thermo.dh
        self.ds = thermo.ds

        self.structure = thermo.ascii_structure
        self.structure_lines = thermo.ascii_structure_lines


class Oligomer:
    def __init__(self, oligomer, start=None, max_runs=4):
        base_runs_format = re.compile("".join(["A{", str(max_runs), "}|"]) +
                                      "".join(["T{", str(max_runs), "}|"]) +
                                      "".join(["G{", str(max_runs), "}|"]) +
                                      "".join(["C{", str(max_runs), "}"]))

        self.seq = oligomer

        self.Tm = calcTm(oligomer)
        self.hairpin = Hairpin(oligomer)
        self.homodimer = Homodimer(oligomer)
        self.GC = SeqUtils.GC(oligomer)

        self.base_run = (re.match(base_runs_format, oligomer) is not None)

        self.start = start

        self._rating = None

    @property
    def rating(self):
        return self.get_rating()

    def set_rating(self):
        self._rating = self.calc_rating(self)

    def get_rating(self):
        if self._rating is None:
            self.set_rating()

        return self._rating

    @classmethod
    def calc_rating(self, oligomer):
        rating = 0

        rating += -4 * (10**(oligomer.hairpin.dg/-1400) / 10**(-2000/-1400))
        rating += -4 * (10**(oligomer.homodimer.dg/-1400) /
                        10**(-5000/-1400))

        return rating


class PrimerPair:
    def __init__(self, fwd, rvs, start=None, end=None, genome=None):
        if isinstance(fwd, str):
            self.fwd = Oligomer(fwd)
        elif isinstance(fwd, Oligomer):
            self.fwd = fwd
        else:
            raise TypeError

        if isinstance(rvs, str):
            self.rvs = Oligomer(rvs)
        elif isinstance(fwd, Oligomer):
            self.rvs = rvs
        else:
            raise TypeError

        self.start = start
        self.end = end

        if self.start is None:
            self.start = self.fwd.start
        if self.end is None and self.start is not None:
            self.end = (self.rvs.start - self.start) + 1

        self.Tm_gap = abs(self.fwd.Tm - self.rvs.Tm)

        self.unst_primer = self.fwd
        if self.rvs.Tm < self.fwd.Tm:
            self.unst_primer = self.rvs
        self.heterodimer = Heterodimer(self.fwd.seq, self.rvs.seq)

        self._genome = genome
        self._product = None
        self._rvs_product = None
        self._annealing_Tm = None
        self._annealing_Tm_gap = None
        self._fwd_clamp = None
        self._rvs_clamp = None
        self._fwd_antisense_heterodimer = None
        self._rvs_antisense_heterodimer = None
        self._rvs_sense_heterodimer = None
        self._fwd_sense_heterodimer = None
        self._rating = None

    @property
    def genome(self):
        return self._genome

    @genome.setter
    def genome(self, genome):
        self._genome = genome

        self._product = None
        self._anti_product = None
        self._annealing_Tm = None
        self._annealing_Tm_gap = None
        self._fwd_clamp = None
        self._rvs_clamp = None
        self._fwd_antisense_heterodimer = None
        self._rvs_antisense_heterodimer = None
        self._fwd_sense_heterodimer = None
        self._rvs_sense_heteroimer = None

    @property
    def product(self):
        return self.get_product()

    @property
    def anti_product(self):
        return self.get_anti_product()

    @property
    def annealing_Tm(self):
        return self.get_annealing_Tm()

    @property
    def annealing_Tm_gap(self):
        return self.get_annealing_Tm_gap()

    @property
    def fwd_clamp(self):
        return self.get_fwd_clamp()

    @property
    def rvs_clamp(self):
        return self.get_rvs_clamp()

    @property
    def fwd_antisense_heterodimer(self):
        return self.get_fwd_antisense_heterodimer()

    @property
    def rvs_antisense_heterodimer(self):
        return self.get_rvs_antisense_heterodimer()

    @property
    def fwd_sense_heterodimer(self):
        return self.get_fwd_sense_heterodimer()

    @property
    def rvs_sense_heterodimer(self):
        return self.get_rvs_sense_heterodimer()

    @property
    def rating(self):
        return self.get_rating()

    def set_product(self):
        if self._genome is None:
            raise AttributeError("Genome sequence is required to extract "
                                 "a primer product.")

        rvs_compl = str(Seq(self.rvs.seq).reverse_complement())

        product_format = re.compile("".join([self.fwd.seq, "[ATGC]*",
                                             rvs_compl]))

        matches = re.findall(product_format, self._genome)

        if len(matches) == 0:
            raise ValueError("PCR product cannot be found in genome.")
        if len(matches) > 1:
            raise ValueError("Given primers have multiple products.")

        product = matches[0]
        self._product = product

    def get_product(self):
        if self._product is None:
            self.set_product()

        return self._product

    def set_anti_product(self):
        product = self.get_product()

        self._anti_product = str(Seq(product).reverse_complement())

    def get_anti_product(self):
        if self._anti_product is None:
            self.set_anti_product()

        return self._anti_product

    def set_annealing_Tm(self):
        product = self.get_product()
        anneal_tm = 0.3*self.unst_primer.Tm + 0.7*calcTm(product) - 14.9

        self._annealing_Tm = anneal_tm

    def get_annealing_Tm(self):
        if self._annealing_Tm is None:
            self.set_annealing_Tm()

        return self._annealing_Tm

    def set_annealing_Tm_gap(self):
        annealing_Tm = self.get_annealing_Tm()

        fwd_gap = abs(annealing_Tm - self.fwd.Tm)
        rvs_gap = abs(annealing_Tm - self.fwd.Tm)

        annealing_Tm_gap = fwd_gap
        if rvs_gap > fwd_gap:
            annealing_Tm_gap = rvs_gap

        self._annealing_Tm_gap = annealing_Tm_gap

    def get_annealing_Tm_gap(self):
        if self._annealing_Tm_gap is None:
            self.set_annealing_Tm_gap()

        return self._annealing_Tm_gap

    def set_fwd_clamp(self):
        anti_product = self.get_anti_product()

        fwd_clamp = Clamp(self.fwd.seq, anti_product)
        self._fwd_clamp = fwd_clamp

    def get_fwd_clamp(self):
        if self._fwd_clamp is None:
            self.set_fwd_clamp()

        return self._fwd_clamp

    def set_rvs_clamp(self):
        product = self.get_product()

        rvs_clamp = Clamp(self.rvs.seq, product)
        self._rvs_clamp = rvs_clamp

    def get_rvs_clamp(self):
        if self._rvs_clamp is None:
            self.set_rvs_clamp()

        return self._rvs_clamp

    def set_fwd_antisense_heterodimer(self):
        internal_anti_product = self.anti_product[:-len(self.fwd.seq)]

        self._fwd_antisense_heterodimer = Heterodimer(self.fwd.seq,
                                                      internal_anti_product)

    def get_fwd_antisense_heterodimer(self):
        if self._fwd_antisense_heterodimer is None:
            self.set_fwd_antisense_heterodimer()

        return self._fwd_antisense_heterodimer

    def set_rvs_antisense_heterodimer(self):
        self._rvs_antisense_heterodimer = Heterodimer(self.rvs.seq,
                                                      self.anti_product)

    def get_rvs_antisense_heterodimer(self):
        if self._rvs_antisense_heterodimer is None:
            self.set_rvs_antisense_heterodimer()

        return self._rvs_antisense_heterodimer

    def set_fwd_sense_heterodimer(self):
        self._fwd_sense_heterodimer = Heterodimer(self.fwd.seq,
                                                  self.product)

    def get_fwd_sense_heterodimer(self):
        if self._fwd_sense_heterodimer is None:
            self.set_fwd_sense_heterodimer()

        return self._fwd_sense_heterodimer

    def set_rvs_sense_heterodimer(self):
        internal_product = self.product[:-len(self.rvs.seq)]

        self._rvs_sense_heterodimer = Heterodimer(self.rvs.seq,
                                                  internal_product)

    def get_rvs_sense_heterodimer(self):
        if self._rvs_sense_heterodimer is None:
            self.set_rvs_sense_heterodimer()

        return self._rvs_sense_heterodimer

    def set_rating(self):
        self._rating = self.calc_rating(self)

    def get_rating(self):
        if self._rating is None:
            self.set_rating()

        return self._rating

    @classmethod
    def calc_rating(self, primer_pair):
        rating = 0

        rating += primer_pair.fwd.rating
        rating += primer_pair.rvs.rating

        rating += -4 * (10**(primer_pair.heterodimer.dg/-1400) /
                        10**(-5000/-1400))

        rating += -4 * (1 / ((10**(primer_pair.fwd_clamp.dg / -1400)) /
                             (10**(-3000 / -1400))))
        rating += -4 * (1 / ((10**(primer_pair.rvs_clamp.dg / -1400)) /
                             (10**(-3000 / -1400))))

        rating += -4 * (10**(primer_pair.fwd_antisense_heterodimer.dg/-1400)
                        / 10**(-5000/-1400))
        rating += -4 * (10**(primer_pair.rvs_antisense_heterodimer.dg/-1400)
                        / 10**(-5000/-1400))
        rating += -4 * (10**(primer_pair.fwd_sense_heterodimer.dg/-1400)
                        / 10**(-5000/-1400))
        rating += -4 * (10**(primer_pair.rvs_sense_heterodimer.dg/-1400)
                        / 10**(-5000/-1400))

        return rating
