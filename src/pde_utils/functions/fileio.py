import textwrap

from Bio.Seq import Seq


def write_primer_txt_file(primer_pair, file_path):
    with file_path.open(mode="w") as filehandle:
        filehandle.write(
                "<pde_utils find_primers analysis>\n"
                f"Primer penalty rating: {round(primer_pair.rating, 2)}\n\n\n")

        filehandle.write(
                     "PRIMER PRODUCT\n"
                     "-----------------------------------------"
                     "----------------------------------------\n"
                     f"Annealing Temp: {primer_pair.annealing_ta}\u00BAC\n")

        split_product = textwrap.wrap(primer_pair.product, 60)
        split_anti_product = textwrap.wrap(primer_pair.anti_product, 60)

        for i in range(len(split_product)):
            product_line = split_product[i]
            anti_product_line = split_anti_product[i]

            filehandle.write(" ".join(["5'", product_line, "3'\n"]))
            filehandle.write(" ".join(["3'", anti_product_line, "5'\n"]))
            filehandle.write("\n")

        filehandle.write("\n\n")

        filehandle.write(
                     "FORWARD PRIMER\n"
                     "----------------------------------------"
                     "---------------------------------------\n"
                     f"Melting Temp: {primer_pair.fwd.Tm}\u00BAC\n\n"
                     f"~{primer_pair.start}...\n"
                     f"5'[{primer_pair.fwd.seq}]      3'\n"
                     f"3' {Seq(primer_pair.fwd.seq).reverse_complement()}"
                     "...... 5'\n\n"
                     "3' GC clamp:\n"
                     f"\t\u0394G: {primer_pair.fwd_clamp.dg}\n"
                     f"\t\u0394H: {primer_pair.fwd_clamp.dh}\n"
                     f"\t\u0394S: {primer_pair.fwd_clamp.ds}\n\n"
                     "Homodimerization:\n"
                     f"\t\u0394G: {primer_pair.fwd.homodimer.dg}\n"
                     f"\t\u0394H: {primer_pair.fwd.homodimer.dh}\n"
                     f"\t\u0394S: {primer_pair.fwd.homodimer.ds}\n\n"
                     f"{primer_pair.fwd.homodimer.structure}\n")

        if primer_pair.fwd.hairpin.dg != 0:
            filehandle.write(
                     "\nHairpin Formation:\n"
                     f"\t\u0394G: {primer_pair.fwd.hairpin.dg}\n"
                     f"\t\u0394H: {primer_pair.fwd.hairpin.dh}\n"
                     f"\t\u0394S: {primer_pair.fwd.hairpin.ds}\n\n"
                     f"{primer_pair.fwd.hairpin.structure}\n")

        filehandle.write("\n\n")

        filehandle.write(
                     "REVERSE PRIMER\n"
                     "----------------------------------------"
                     "---------------------------------------\n"
                     f"Melting Temp: {primer_pair.rvs.Tm}\u00BAC\n\n"
                     f"...~{primer_pair.start + len(primer_pair.product)}\n"
                     f"5'     [{primer_pair.rvs.seq}] 3'\n"
                     f"3'......{Seq(primer_pair.rvs.seq).reverse_complement()}"
                     "  5\n\n"
                     "3' GC clamp:\n"
                     f"\t\u0394G: {primer_pair.rvs_clamp.dg}\n"
                     f"\t\u0394H: {primer_pair.rvs_clamp.dh}\n"
                     f"\t\u0394S: {primer_pair.rvs_clamp.ds}\n\n"
                     "Homodimerization:\n"
                     f"\t\u0394G: {primer_pair.rvs.homodimer.dg}\n"
                     f"\t\u0394H: {primer_pair.rvs.homodimer.dh}\n"
                     f"\t\u0394S: {primer_pair.rvs.homodimer.ds}\n\n"
                     f"{primer_pair.rvs.homodimer.structure}\n")

        if primer_pair.rvs.hairpin.dg != 0:
            filehandle.write(
                     "\n"
                     "Hairpin Formation:\n"
                     f"\t\u0394G: {primer_pair.rvs.hairpin.dg}\n"
                     f"\t\u0394H: {primer_pair.rvs.hairpin.dh}\n"
                     f"\t\u0394S: {primer_pair.rvs.hairpin.ds}\n\n"
                     f"{primer_pair.rvs.hairpin.structure}\n")

        filehandle.write("\n\n")

        filehandle.write(
                    "HETERODIMER FORMATION\n"
                    "----------------------------------------"
                    "---------------------------------------\n"
                    "Forward/Reverse Dimerization:\n"
                    f"\t\u0394G: {primer_pair.heterodimer.dg}\n"
                    f"\t\u0394H: {primer_pair.heterodimer.dh}\n"
                    f"\t\u0394S: {primer_pair.heterodimer.ds}\n\n"
                    f"{primer_pair.heterodimer.structure}\n\n")

        if primer_pair.fwd_antisense_heterodimer.structure_lines:
            heterodimer = primer_pair.fwd_antisense_heterodimer
            filehandle.write(
                    "Forward/Antisense Internal Dimerization:\n"
                    f"\t\u0394G: {heterodimer.dg}\n"
                    f"\t\u0394H: {heterodimer.dh}\n"
                    f"\t\u0394S: {heterodimer.ds}\n"
                    "\n"
                    f"{heterodimer.reduced_structure}")

            for i in range(len(fwd_antisense_lines[0])):
                for split_lines in fwd_antisense_lines:
                    filehandle.write("".join([split_lines[i], "\n"]))

            filehandle.write("\n")

        if primer_pair.rvs_antisense_heterodimer.structure_lines:
            filehandle.write(
                    "Reverse/Antisense Dimerization:\n"
                    f"\t\u0394G: {primer_pair.rvs_antisense_heterodimer.dg}\n"
                    f"\t\u0394H: {primer_pair.rvs_antisense_heterodimer.dh}\n"
                    f"\t\u0394S: {primer_pair.rvs_antisense_heterodimer.ds}\n"
                    "\n")
            rvs_antisense_lines = []
            for line in primer_pair.rvs_antisense_heterodimer.structure_lines:
                rvs_antisense_lines.append(textwrap.wrap(
                                                line, 60,
                                                replace_whitespace=False,
                                                drop_whitespace=False))

            for i in range(len(rvs_antisense_lines[0])):
                for split_lines in rvs_antisense_lines:
                    filehandle.write("".join([split_lines[i], "\n"]))

            filehandle.write("\n")

        if primer_pair.fwd_sense_heterodimer.structure_lines:
            filehandle.write(
                    "Forward/Sense Dimerization:\n"
                    f"\t\u0394G: {primer_pair.fwd_sense_heterodimer.dg}\n"
                    f"\t\u0394H: {primer_pair.fwd_sense_heterodimer.dh}\n"
                    f"\t\u0394S: {primer_pair.fwd_sense_heterodimer.ds}\n"
                    "\n")
            fwd_sense_lines = []
            for line in primer_pair.fwd_sense_heterodimer.structure_lines:
                fwd_sense_lines.append(textwrap.wrap(
                                                line, 60,
                                                replace_whitespace=False,
                                                drop_whitespace=False))

            for i in range(len(fwd_sense_lines[0])):
                for split_lines in fwd_sense_lines:
                    filehandle.write("".join([split_lines[i], "\n"]))

            filehandle.write("\n")

        if primer_pair.rvs_sense_heterodimer.structure_lines:
            filehandle.write(
                    "Reverse/Sense Internal Dimerization:\n"
                    f"\t\u0394G: {primer_pair.rvs_sense_heterodimer.dg}\n"
                    f"\t\u0394H: {primer_pair.rvs_sense_heterodimer.dh}\n"
                    f"\t\u0394S: {primer_pair.rvs_sense_heterodimer.ds}\n")
            rvs_sense_lines = []
            for line in primer_pair.rvs_sense_heterodimer.structure_lines:
                rvs_sense_lines.append(textwrap.wrap(
                                                line, 60,
                                                replace_whitespace=False,
                                                drop_whitespace=False))

            for i in range(len(rvs_sense_lines[0])):
                for split_lines in rvs_sense_lines:
                    filehandle.write("".join([split_lines[i], "\n"]))
