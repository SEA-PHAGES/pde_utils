import json
import textwrap

from Bio.Seq import Seq
from networkx import readwrite


def write_graph(graph, file_format, export_path, file_name):
    file_path = export_path.joinpath(f"{file_name}.{file_format}")

    if file_format == "csv":
        readwrite.edgelist.write_edgelist(graph, file_path, delimiter=",",
                                          data=EDGE_WEIGHTS)
    elif file_format == "gexf":
        readwrite.gexf.write_gexf(graph, file_path)
    elif file_format == "gml":
        readwrite.gml.write_gml(graph, file_path)
    elif file_format == "gpickle":
        readwrite.gpickle.write_gpickle(graph, file_path)
    elif file_format == "graphml":
        readwrite.graphml.write_graphml(graph, file_path)
    elif file_format == "json" or file_format == "cyjs":
        if file_format == "json":
            json_data = readwrite.json_graph.node_link_data(graph)
        else:
            json_data = readwrite.json_graph.cytoscape_data(graph)

        file_path.touch()
        file_handle = file_path.open(mode="w")
        json.dump(json_data, file_handle)
        file_handle.close()
    elif file_format == "yaml":
        readwrite.nx_yaml.write_yaml(graph, file_path)
    elif file_format == "pajek":
        readwrite.pajek.write_pajek(graph, file_path)
    elif file_format == "shp":
        readwrite.nx_shp.write_shp(graph, export_path)
    elif file_format == "al":
        readwrite.adjlist.write_adjlist(graph, file_path)
    elif file_format == "nl":
        file_path.touch()
        file_handle = file_path.open(mode="w")
        for node in list(graph.nodes()):
            file_handle.write(f"{node}\n")
        file_handle.close()
    else:
        raise ValueError("Graph output format is not recognized.")


def write_primer_txt_file(primer_pair, file_path):
    with file_path.open(mode="w") as filehandle:
        filehandle.write(
                "<pde_utils find_primers analysis>\n"
                f"Primer penalty rating: {round(primer_pair.rating, 2)}\n"
                f"Annealing Temp: {primer_pair.annealing_ta}\u00BAC\n\n\n")

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
                     f"{primer_pair.fwd.homodimer.formatted_structure}\n")

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
                     f"{primer_pair.rvs.homodimer.formatted_structure}\n")

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
                    f"{primer_pair.heterodimer.formatted_structure}\n\n")

        if primer_pair.fwd_antisense_heterodimer.structure_lines:
            heterodimer = primer_pair.fwd_antisense_heterodimer
            filehandle.write(
                    "Forward/Antisense Internal Mispriming:\n"
                    f"\t\u0394G: {heterodimer.dg}\n"
                    f"\t\u0394H: {heterodimer.dh}\n"
                    f"\t\u0394S: {heterodimer.ds}\n"
                    "\n"
                    f"{heterodimer.reduced_structure}\n")

        if primer_pair.rvs_antisense_heterodimer.structure_lines:
            heterodimer = primer_pair.rvs_antisense_heterodimer
            filehandle.write(
                    "Reverse/Antisense Mispriming:\n"
                    f"\t\u0394G: {heterodimer.dg}\n"
                    f"\t\u0394H: {heterodimer.dh}\n"
                    f"\t\u0394S: {heterodimer.ds}\n"
                    "\n"
                    f"{heterodimer.reduced_structure}\n")

        if primer_pair.fwd_sense_heterodimer.structure_lines:
            heterodimer = primer_pair.fwd_sense_heterodimer
            filehandle.write(
                    "Forward/Sense Mispriming:\n"
                    f"\t\u0394G: {heterodimer.dg}\n"
                    f"\t\u0394H: {heterodimer.dh}\n"
                    f"\t\u0394S: {heterodimer.ds}\n"
                    "\n"
                    f"{heterodimer.reduced_structure}\n")

        if primer_pair.rvs_sense_heterodimer.structure_lines:
            heterodimer = primer_pair.fwd_sense_heterodimer
            filehandle.write(
                    "Reverse/Sense Internal Mispriming:\n"
                    f"\t\u0394G: {heterodimer.dg}\n"
                    f"\t\u0394H: {heterodimer.dh}\n"
                    f"\t\u0394S: {heterodimer.ds}\n"
                    "\n"
                    f"{heterodimer.reduced_structure}\n")

        filehandle.write(
                    "SAMPLE PRIMER PRODUCT\n"
                    "----------------------------------------"
                    "---------------------------------------\n\n")

        split_product = textwrap.wrap(primer_pair.product, 60)
        split_anti_product = textwrap.wrap(primer_pair.anti_product, 60)

        for i in range(len(split_product)):
            product_line = split_product[i]
            anti_product_line = split_anti_product[i]

            filehandle.write(" ".join(["5'", product_line, "3'\n"]))
            filehandle.write(" ".join(["3'", anti_product_line, "5'\n"]))
            filehandle.write("\n")

        filehandle.write("\n\n")
