from Bio import SeqUtils
from primer3 import (
        calcTm, calcHairpin, calcHomodimer, calcHeterodimer)


def filter_oligomers(oligomers, tmMin=52, tmMax=58, GC_max=60,
                     hpn_dG_min=-2000, homo_dG_min=-5000):
    stable_oligomers = []
    for oligomer in oligomers:
        if check_oligomer(oligomer, tmMin=tmMin, tmMax=tmMax, GC_max=GC_max,
                          hpn_dG_min=hpn_dG_min, homo_dG_min=homo_dG_min):
            stable_oligomers.append(oligomer)

    return stable_oligomers


def check_oligomer(oligomer, tmMin=52, tmMax=58, hpn_dG_min=-2000,
                   homo_dG_min=-5000, GC_max=60):
    stable = (check_stable_tm(oligomer, tmMin=tmMin, tmMax=tmMax) and
              (not check_hairpin(oligomer, delG_min=hpn_dG_min)) and
              check_optimal_GC(oligomer, GC_max=GC_max) and
              (not check_homodimer(oligomer, delG_min=homo_dG_min)))

    return stable


def check_stable_tm(oligomer, tmMin=52, tmMax=58):
    tm = calcTm(oligomer)

    return (tm >= tmMin and tm <= tmMax)


def check_hairpin(oligomer, delG_min=-2000):
    hairpin_result = calcHairpin(oligomer)

    if hairpin_result.dg < delG_min:
        return True

    return False


def check_optimal_GC(oligomer, GC_max=60):
    gc = SeqUtils.GC(oligomer)

    return (gc <= GC_max)


def check_homodimer(oligomer, delG_min=-5000):
    homodimer_result = calcHomodimer(oligomer)

    if homodimer_result.dg < delG_min:
        return True

    return False


def check_heterodimer(oligomer_1, oligomer_2, delG_min=-5000):
    heterodimer_result = calcHeterodimer(oligomer_1, oligomer_2)

    if heterodimer_result.dg < delG_min:
        return True

    return False
