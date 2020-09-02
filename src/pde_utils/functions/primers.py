from pde_utils.classes import primer_TD


def filter_oligomers(oligomers, tmMin=52, tmMax=58, GC_max=60,
                     hpn_dG_min=-2000, homo_dG_min=-5000, max_runs=4):
    stable_oligomers = []
    for oligomer in oligomers:
        if check_oligomer(oligomer, tmMin=tmMin, tmMax=tmMax, GC_max=GC_max,
                          hpn_dG_min=hpn_dG_min, homo_dG_min=homo_dG_min):
            stable_oligomers.append(oligomer)

    return stable_oligomers


def check_oligomer(oligomer, tmMin=52, tmMax=58, hpn_dG_min=-2000,
                   homo_dG_min=-5000, GC_max=60, max_runs=4):
    oligomer_TD = primer_TD.OligomerTDynamics(oligomer)

    return ((oligomer_TD.Tm >= tmMin and oligomer_TD.Tm <= tmMax) and
            (oligomer_TD.hairpin.dg > hpn_dG_min) and
            (oligomer_TD.GC < GC_max) and
            (oligomer_TD.homodimer.dg > homo_dG_min) and
            (not oligomer_TD.base_run))
