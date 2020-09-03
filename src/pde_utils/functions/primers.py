from pde_utils.classes import primer_TD


def get_stable_oligomer(oligomer, tmMin=52, tmMax=58, hpn_dG_min=-2000,
                        homo_dG_min=-5000, GC_max=60, max_runs=4):
    oligomer_TD = primer_TD.OligomerTDynamics(oligomer)

    stable = ((oligomer_TD.Tm >= tmMin and oligomer_TD.Tm <= tmMax) and
              (oligomer_TD.hairpin.dg > hpn_dG_min) and
              (oligomer_TD.GC < GC_max) and
              (oligomer_TD.homodimer.dg > homo_dG_min) and
              (not oligomer_TD.base_run))

    if stable:
        return oligomer_TD
    else:
        return None


def write_primer_txt_file(primer_pair, file_path):
    with file_path.open(mode="w") as filehandle:
        filehandle.write(f"Annealing Temp: {primer_pair.annealing_Tm}\u00BAC")
        filehandle.write("\n\n")

        filehandle.write(f"~{primer_pair.start}......\n"
                         f"{primer_pair.fwd.seq}......\n"
                         f"Forward Melting Temp: {primer_pair.fwd.Tm}\u00BAC\n"
                         "Forward Homodimerization:\n"
                         f"\t\u0394G: {primer_pair.fwd.homodimer.dg}\n"
                         f"\t\u0394H: {primer_pair.fwd.homodimer.dh}\n"
                         f"\t\u0394S: {primer_pair.fwd.homodimer.ds}\n"
                         "Forward 3' end stability:\n"
                         f"\t\u0394G: {primer_pair.fwd_clamp.dg}\n"
                         f"\t\u0394H: {primer_pair.fwd_clamp.dh}\n"
                         f"\t\u0394S: {primer_pair.fwd_clamp.ds}\n")

        if primer_pair.fwd.hairpin.dg != 0:
            filehandle.write(
                     "Forward Hairpin Formation:\n"
                     f"\t\u0394G: {primer_pair.fwd.hairpin.dg}\n"
                     f"\t\u0394H: {primer_pair.fwd.hairpin.dh}\n"
                     f"\t\u0394S: {primer_pair.fwd.hairpin.ds}\n")

        filehandle.write("\n")

        filehandle.write(
                     "......~"
                     f"{primer_pair.start + len(primer_pair.product)}\n"
                     f"......{primer_pair.rvs.seq}\n"
                     f"Reverse Melting Temp: {primer_pair.rvs.Tm}\u00BAC\n"
                     "Reverse Homodimerization:\n"
                     f"\t\u0394G: {primer_pair.rvs.homodimer.dg}\n"
                     f"\t\u0394H: {primer_pair.rvs.homodimer.dh}\n"
                     f"\t\u0394S: {primer_pair.rvs.homodimer.ds}\n"
                     "Reverse 3' end stability:\n"
                     f"\t\u0394G: {primer_pair.rvs_clamp.dg}\n"
                     f"\t\u0394H: {primer_pair.rvs_clamp.dh}\n"
                     f"\t\u0394S: {primer_pair.rvs_clamp.ds}\n")

        if primer_pair.rvs.hairpin.dg != 0:
            filehandle.write(
                     "Reverse Hairpin Formation:\n"
                     f"\t\u0394G: {primer_pair.rvs.hairpin.dg}\n"
                     f"\t\u0394H: {primer_pair.rvs.hairpin.dh}\n"
                     f"\t\u0394S: {primer_pair.rvs.hairpin.ds}\n")

        filehandle.write("\n")

        filehandle.write(
                     "Forward/Reverse Dimerization:\n"
                     f"\t\u0394G: {primer_pair.heterodimer.dg}\n"
                     f"\t\u0394H: {primer_pair.heterodimer.dh}\n"
                     f"\t\u0394S: {primer_pair.heterodimer.ds}")
