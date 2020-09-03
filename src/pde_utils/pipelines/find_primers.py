"""Pipeline to find potential primers from given sequences quickly
   using phamilies as an potential indicator of conserved nucleotide
   regions.
   """

import argparse
import sys
import time
from pathlib import Path

from Bio.Seq import Seq
from numpy import std
from pdm_utils.functions import annotation
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import pipelines_basic
from pdm_utils.pipelines import export_db

from pde_utils.classes import primer_TD
from pde_utils.functions import primers
from pde_utils.functions import seq

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_find_primers")

PHAGE_QUERY = "SELECT * FROM phage"
GENE_QUERY = "SELECT * FROM gene"
TRNA_QUERY = "SELECT * FROM trna"
TMRNA_QUERY = "SELECT * FROM tmrna"


def main(unparsed_args):
    args = parse_find_primers(unparsed_args)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)

    execute_find_primers(alchemist, folder_path=args.folder_path,
                         folder_name=args.folder_name, values=values,
                         filters=args.filters, groups=args.groups,
                         verbose=args.verbose)


def parse_find_primers(unparsed_args):
    DATABASE_HELP = """Name of the MySQL database to pull sequences from."""

    CONFIG_FILE_HELP = """
        Find_primers option that enables use of a config file for sourcing
        credentials.
            Follow selection argument with the path to the config file.
        """
    VERBOSE_HELP = """
        Find_primers option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Find_primers option to change the path of the directory where the
        exported files are stored.
            Follow selection argument with the path to the desired export
            directory.
        """
    FOLDER_NAME_HELP = """
        Export option to change the name of the directory where the exported
        files are stored.
            Follow selection argument with the desired name.
        """

    IMPORT_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selection argument with path to the
            csv file containing the names of each genome in the first column.
        """
    SINGLE_GENOMES_HELP = """
        Selection input option that imports values from cmd line input.
            Follow selection argument with space separated names of genomes
            in the database.
        """
    WHERE_HELP = """
        Data filtering option that filters data by the inputted expressions.
            Follow selection argument with formatted filter expression:
                {Table}.{Column}={Value}
        """
    GROUP_BY_HELP = """
        Data selection option that groups data by the inputted columns.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """

    parser = argparse.ArgumentParser()

    parser.add_argument("database", type=str, help=DATABASE_HELP)

    parser.add_argument("-c", "--config_file",
                        type=pipelines_basic.convert_file_path,
                        help=CONFIG_FILE_HELP)
    parser.add_argument("-m", "--folder_name",
                        type=str,  help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path",
                        type=Path, help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)
    parser.add_argument("-th", "--threads")

    parser.add_argument("-if", "--import_file", dest="input",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP)
    parser.add_argument("-in", "--import_names", nargs="*", dest="input",
                        help=SINGLE_GENOMES_HELP)
    parser.add_argument("-w", "--where", nargs="?", dest="filters",
                        help=WHERE_HELP)
    parser.add_argument("-g", "--group_by", nargs="*", dest="groups",
                        help=GROUP_BY_HELP)

    parser.add_argument("-prc")
    parser.add_argument("-minD")
    parser.add_argument("-maxD")
    parser.add_argument("-hpn_min")
    parser.add_argument("-ho_min")
    parser.add_argument("-het_min")
    parser.add_argument("-GC")
    parser.add_argument("-max_std")
    parser.add_argument("-len")
    parser.add_argument("-ex")

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, verbose=False, input=[],
                        filters="", groups=[])

    parsed_args = parser.parse_args(unparsed_args[2:])
    return parsed_args


def execute_find_primers(alchemist, folder_path=None,
                         folder_name=DEFAULT_FOLDER_NAME, values=None,
                         filters="", groups=[], verbose=False,
                         threads=4, prc=0.8, max_std=3000, len_oligomer=20,
                         minD=900, maxD=1100, tmMin=52, tmMax=58,
                         hpn_dG_min=-2000, homo_dG_min=-5000, GC_max=60,
                         hetero_dG_min=-5000, exclude=True):
    db_filter = pipelines_basic.build_filter(alchemist, "phage", filters)

    working_path = pipelines_basic.create_working_path(
                                                    folder_path, folder_name)

    conditionals_map = pipelines_basic.build_groups_map(
                                                    db_filter, working_path,
                                                    groups=groups,
                                                    verbose=verbose)
    results_map = {}

    if verbose:
        print("Prepared query and path structure, beginning primer search...")

    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        db_filter.reset()
        db_filter.key = "phage"
        db_filter.values = values

        conditionals = conditionals_map[mapped_path]

        db_filter.values = db_filter.build_values(where=conditionals)

        if db_filter.hits() == 0:
            print("No database entries received from phage "
                  f" for '{mapped_path.name}'.")

        genome_map = {}
        for genome_id in db_filter.values:
            export_db.get_single_genome(alchemist, genome_id,
                                        data_cache=genome_map)
        pham_gene_map = build_pham_gene_map(db_filter, conditionals,
                                            verbose=verbose)

        if verbose:
            print(f"Identifying primer pairs for '{mapped_path}'...")
        primer_pairs = find_primer_pairs(
                                    alchemist, pham_gene_map, genome_map,
                                    verbose=verbose, threads=threads, prc=prc,
                                    max_std=max_std, minD=minD, maxD=maxD,
                                    len_oligomer=len_oligomer, tmMin=tmMin,
                                    tmMax=tmMax, hpn_dG_min=hpn_dG_min,
                                    homo_dG_min=homo_dG_min, GC_max=GC_max)

        primer_pairs = test_primer_pairs(primer_pairs, genome_map,
                                         verbose=verbose,
                                         hetero_dG_min=hetero_dG_min)

        if verbose:
            print(f"Selecting primer_pair from {len(primer_pairs)} "
                  f"primer pairs for '{mapped_path}'...")

        if primer_pairs:
            results_map[mapped_path] = (primer_pairs, genome_map)

    if exclude:
        results_map = select_primer_pairs(results_map)

    for mapped_path, primer_pairs in results_map.items():
        pipelines_basic.create_working_dir(mapped_path)
        file_path = mapped_path.joinpath("primer.txt")
        primers.write_primer_txt_file(primer_pairs[0][0], file_path)


def find_primer_pairs(alchemist, pham_gene_map, genome_map, verbose=False,
                      threads=4, prc=0.8, max_std=3000, len_oligomer=20,
                      minD=900, maxD=1100, tmMin=52, tmMax=58,
                      hpn_dG_min=-2000, homo_dG_min=-5000, GC_max=60):
    pham_histogram = {}
    for pham, genes in pham_gene_map.items():
        pham_histogram[pham] = len(genes)
    pham_histogram = basic.sort_histogram(pham_histogram)

    F_pos_oligomer_map = {}
    R_pos_oligomer_map = {}
    for pham, count in pham_histogram.items():
        pham_per = count/len(genome_map)

        if pham_per < prc:
            break

        if verbose:
            print(f"...Analyzing pham {pham}...")
            print(f"......Pham is represented in {round(pham_per, 3)*100}%"
                  " of currently viewed genomes...")

        cds_list = export_db.parse_feature_data(alchemist,
                                                values=pham_gene_map[pham])

        starts = []
        avg_start = 0
        for cds in cds_list:
            starts.append(cds.start)
            avg_start += cds.start
        avg_start = int(round(avg_start/len(cds_list), 0))

        start_std = std(starts)
        if start_std > max_std:
            if verbose:
                print(f"......Pham {pham} has unstable synteny...")
            continue

        avg_orientation = 0
        for cds in cds_list:
            if cds.orientation == "F":
                avg_orientation += 1
        avg_orientation = avg_orientation/len(cds_list)
        if avg_orientation != 0 and avg_orientation != 1:
            if verbose:
                print(f"......Pham {pham} has unstable direction...")
            continue

        cds_to_seq = seq.map_cds_to_seq(alchemist, cds_list,
                                        data_cache=genome_map)
        conserved_kmer_data = seq.find_conserved_kmers(cds_to_seq,
                                                       len_oligomer)

        num_conserved_kmers = len(conserved_kmer_data)
        if num_conserved_kmers == 0:
            if verbose:
                print(f"......Pham {pham} has no conserved kmers...")
            continue

        if verbose:
            print(f"......Pham {pham} has {num_conserved_kmers} "
                  "conserved kmers...")

        F_oligomers = get_stable_oligomers(
                                        conserved_kmer_data, avg_start, "F",
                                        tmMin=tmMin, tmMax=tmMax,
                                        hpn_dG_min=hpn_dG_min,
                                        homo_dG_min=homo_dG_min, GC_max=GC_max)
        R_oligomers = get_stable_oligomers(
                                        conserved_kmer_data, avg_start, "R",
                                        tmMin=tmMin, tmMax=tmMax,
                                        hpn_dG_min=hpn_dG_min,
                                        homo_dG_min=homo_dG_min, GC_max=GC_max)

        num_stable_oligomers = len(F_oligomers) + len(R_oligomers)

        if num_stable_oligomers == 0:
            if verbose:
                print(f"......Pham {pham} has no stable oligomers...")
            continue

        if verbose:
            print(f"......Pham {pham} has {len(F_oligomers)} "
                  "stable forward putative primer oligomers.")
            print(f"......Pham {pham} has {len(R_oligomers)} "
                  "stable reverse putative primer oligomers.")

        pham_F_pos_oligomer_map = map_pos_to_oligomer(F_oligomers,
                                                      verbose=verbose)

        pham_R_pos_oligomer_map = map_pos_to_oligomer(R_oligomers,
                                                      verbose=verbose)

        F_pos_oligomer_map.update(pham_F_pos_oligomer_map)
        R_pos_oligomer_map.update(pham_R_pos_oligomer_map)

    if verbose:
        print("Matching oligomers to create primer pairs...")
    primer_pairs = link_primer_oligomers(
                                        F_pos_oligomer_map, R_pos_oligomer_map,
                                        minD=minD, maxD=maxD)
    return primer_pairs


def test_primer_pairs(primer_pairs, genome_map, verbose=False,
                      hetero_dG_min=-5000, Tm_gap_max=5, anneal_gap_max=5,
                      minD=900, maxD=1100):
    tested_primer_pairs = []

    num_pairs = len(primer_pairs)
    for i in range(num_pairs):
        primer_pair = primer_pairs[i]

        if verbose:
            print(f"...Testing primer pair {i} out of "
                  f"{num_pairs} putative primers...")

        valid_primers = True
        for genome_id, genome in genome_map.items():
            primer_pair.genome = str(genome.seq)

            try:
                valid_primers = (
                        (primer_pair.Tm_gap < Tm_gap_max) and
                        (len(primer_pair.product) > minD) and
                        (len(primer_pair.product) < maxD) and
                        (primer_pair.heterodimer.dg >= hetero_dG_min) and
                        (primer_pair.annealing_Tm_gap < anneal_gap_max))
            except AttributeError:
                valid_primers = False

            if not valid_primers:
                break

        if valid_primers:
            tested_primer_pairs.append(primer_pair)

    if verbose:
        print("...Rating tested primer pairs...")
    tested_primer_pairs = sorted(tested_primer_pairs,
                                 key=lambda pair: pair.rating,
                                 reverse=True)

    return tested_primer_pairs


def select_primer_pairs(results_map):
    selected_results_map = {}

    for mapped_path, data in results_map.items():

        valid_primer_pair = None
        for primer_pair in data[0]:

            valid_pair = True
            for inner_mapped_path, inner_data in results_map.items():
                if mapped_path == mapped_path:
                    continue

                for genome_id, genome in inner_data[1]:
                    try:
                        primer_pair.genome = str(genome.seq)
                        primer_pair.product
                        valid_pair = False
                    except AttributeError:
                        pass

                if not valid_pair:
                    break

            if valid_pair:
                valid_primer_pair = primer_pair
                break

        if valid_primer_pair is not None:
            selected_results_map[mapped_path] = ([valid_primer_pair],)
        else:
            print(f"All primer pairs for '{mapped_path}' generate products "
                  "against non-grouped genomes")
            sys.exit(1)

    return selected_results_map


# TO FIX IN BASIC
def invert_dictionary(dictionary):
    """Inverts a dictionary, where the values and keys are swapped.

    :param dictionary: A dictionary to be inverted.
    :type dictionary: dict
    :returns: Returns an inverted dictionary of the given dictionary.
    :rtype: dict
    """
    new_dict = {}
    for key, value in dictionary.items():
        if isinstance(value, list):
            for subvalue in value:
                if new_dict.get(subvalue) is None:
                    new_dict[subvalue] = [key]
                else:
                    new_dict[subvalue].append(key)
        else:
            if new_dict.get(value) is None:
                new_dict[value] = [key]
            else:
                new_dict[value].append(key)

    return new_dict


def get_stable_oligomers(conserved_kmer_data, avg_start, orientation,
                         tmMin=52, tmMax=58, hpn_dG_min=-2000,
                         homo_dG_min=-5000, GC_max=60):
    oligomers = []

    for kmer, kmer_data in conserved_kmer_data.items():
        if kmer_data[0][0].orientation != orientation:
            kmer = str(Seq(kmer).reverse_complement())

        oligomer = primers.get_stable_oligomer(
                                       kmer, tmMin=tmMin, tmMax=tmMax,
                                       GC_max=GC_max, hpn_dG_min=hpn_dG_min,
                                       homo_dG_min=homo_dG_min)
        if oligomer is None:
            continue

        avg_kmer_start = 0
        for data in kmer_data:
            avg_kmer_start += data[1]
        avg_kmer_start = int(round(avg_kmer_start/len(kmer_data), 0))
        oligomer.start = avg_kmer_start

        oligomers.append(oligomer)

    return oligomers


def map_pos_to_oligomer(oligomers, verbose=False):
    pos_oligomer_map = {}
    for oligomer in oligomers:
        mapped_oligomers = pos_oligomer_map.get(oligomer.start, [])
        mapped_oligomers.append(oligomer)
        pos_oligomer_map[oligomer.start] = mapped_oligomers

    return pos_oligomer_map


def link_primer_oligomers(forward_position_map, reverse_position_map,
                          minD=900, maxD=1100):
    primer_pairs = []

    for pos, F_oligomers in forward_position_map.items():
        for i in range(minD, maxD+1):
            R_oligomers = reverse_position_map.get(pos+i)
            if R_oligomers:
                for F_oligomer in F_oligomers:
                    for R_oligomer in R_oligomers:
                        primer_pair = primer_TD.PrimerPairTDynamics(
                                                F_oligomer, R_oligomer)
                        primer_pairs.append(primer_pair)

    primer_pairs = sorted(primer_pairs, key=lambda pair: pair.start)
    return primer_pairs


# TO RELOCATE
def build_pham_gene_map(db_filter, conditionals, verbose=False):
    gene_col = db_filter.get_column("gene.GeneID")
    pham_col = db_filter.get_column("gene.PhamID")

    phams = db_filter.transpose("gene.PhamID")

    pham_gene_map = {}
    for pham in phams:
        pham_conditionals = conditionals + [(pham_col == pham)]
        pham_gene_map[pham] = db_filter.build_values(column=gene_col,
                                                     where=pham_conditionals)

    return pham_gene_map


# TO FIX IN ANNOTATION
def get_count_phams_in_genes(alchemist, geneids, incounts=None):
    phams = annotation.get_phams_from_genes(alchemist, geneids)

    pham_histogram = {}
    if incounts is not None:
        pham_histogram = incounts

    basic.increment_histogram(phams, pham_histogram)
    return pham_histogram


if __name__ == "__main__":
    main(sys.argv)
