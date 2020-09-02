import re

from Bio import SeqUtils
from Bio.Seq import Seq
from primer3 import (
        calcTm, calcHairpin, calcHomodimer, calcHeterodimer)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
class OligomerTDynamics:
    def __init__(self, oligomer, max_runs=4):
        base_runs_format = re.compile("".join(["A{", str(max_runs), "}|"]) +
                                      "".join(["T{", str(max_runs), "}|"]) +
                                      "".join(["G{", str(max_runs), "}|"]) +
                                      "".join(["C{", str(max_runs), "}"]))

        self.seq = oligomer

        self.Tm = calcTm(oligomer)
        self.hairpin = calcHairpin(oligomer)
        self.homodimer = calcHomodimer(oligomer)
        self.GC = SeqUtils.GC(oligomer)

        self.base_run = (re.match(base_runs_format, oligomer) is not None)


class PrimerPairTDynamics:
    def __init__(self, fwd, rvs, start=None, genome=None):
        self.fwd = OligomerTDynamics(fwd)
        self.rvs = OligomerTDynamics(rvs)
        self.Tm_gap = abs(self.fwd.Tm - self.rvs.Tm)

        self.unst_primer = self.fwd
        if self.rvs.Tm < self.fwd.Tm:
            self.unst_primer = self.rvs
        self.heterodimer = calcHeterodimer(self.fwd.seq, self.rvs.seq)

        self.start = start

        self._genome = genome
        self._product = None
        self._annealing_Tm = None
        self._annealing_Tm_gap = None

    @property
    def genome(self):
        return self._genome

    @property
    def product(self):
        return self.get_product()

    @property
    def annealing_Tm(self):
        return self.get_annealing_Tm()

    @property
    def annealing_Tm_gap(self):
        return self.get_annealing_Tm_gap()

    def set_product(self):
        if self._genome is None:
            raise AttributeError("Genome sequence is required to extract "
                                 "a primer product.")

        rvs_compl = str(Seq(self.rvs.seq).reverse_complement())

        product_format = re.compile("".join([self.fwd.seq, "[ATGC]*",
                                             rvs_compl]))

        matches = re.findall(product_format, self._genome)

        if len(matches) == 0:
            raise Exception("PCR product cannot be found in genome.")
        if len(matches) > 1:
            raise Exception("Given primers have multiple products.")

        product = matches[0]
        self._product = product

    def get_product(self):
        if self._product is None:
            self.set_product()

        return self._product

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
