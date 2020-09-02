"""Pipeline to find potential primers from given sequences quickly
   using phamilies as an potential indicator of conserved nucleotide
   regions.
   """

import argparse
import sys
import time
from collections import OrderedDict
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

    parser.add_argument("-if", "--import_file", dest="input",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP)
    parser.add_argument("-in", "--import_names", nargs="*", dest="input",
                        help=SINGLE_GENOMES_HELP)
    parser.add_argument("-w", "--where", nargs="?", dest="filters",
                        help=WHERE_HELP)
    parser.add_argument("-g", "--group_by", nargs="*", dest="groups",
                        help=GROUP_BY_HELP)

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
                         hpn_dG_min=-500, homo_dG_min=-1000, GC_max=60,
                         hetero_dG_min=-1000):
    db_filter = pipelines_basic.build_filter(alchemist, "phage", filters)

    working_path = pipelines_basic.create_working_path(
                                                    folder_path, folder_name)

    conditionals_map = pipelines_basic.build_groups_map(
                                                    db_filter, working_path,
                                                    groups=groups,
                                                    verbose=verbose)
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
            print("Identifying primer pairs for '{mapped_path}'...")
        primer_pairs = find_primer_pairs(
                                    alchemist, pham_gene_map, genome_map,
                                    verbose=verbose, threads=threads, prc=prc,
                                    max_std=max_std, minD=minD, maxD=maxD,
                                    len_oligomer=len_oligomer, tmMin=tmMin,
                                    tmMax=tmMax, hpn_dG_min=hpn_dG_min,
                                    homo_dG_min=homo_dG_min, GC_max=GC_max)

        primer_pair_TDs = test_primer_pairs(primer_pairs, genome_map,
                                            verbose=verbose,
                                            hetero_dG_min=hetero_dG_min)

        if verbose:
            print(f"Selecting primer_pair from {len(primer_pair_TDs)} "
                  f"primer pairs for '{mapped_path}'...")

        opt_primer_pair = select_primer_pair(primer_pair_TDs)
        print(opt_primer_pair)
        sys.exit(1)
        pipelines_basic.create_working_dir(working_path)


def find_primer_pairs(alchemist, pham_gene_map, genome_map, verbose=False,
                      threads=4, prc=0.8, max_std=3000, len_oligomer=20,
                      minD=900, maxD=1100, tmMin=52, tmMax=58,
                      hpn_dG_min=-500, homo_dG_min=-1000, GC_max=60):
    pham_histogram = {}
    for pham, genes in pham_gene_map.items():
        pham_histogram[pham] = len(genes)
    pham_histogram = basic.sort_histogram(pham_histogram)

    F_oligomer_pos_map = OrderedDict()
    R_oligomer_pos_map = OrderedDict()
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

        F_oligomers_cds_map = get_oriented_oligomers(
                                        conserved_kmer_data, "F")
        R_oligomers_cds_map = get_oriented_oligomers(
                                        conserved_kmer_data, "R")

        F_oligomers = primers.filter_oligomers(
                                    F_oligomers_cds_map.keys(), tmMin=tmMin,
                                    tmMax=tmMax, hpn_dG_min=hpn_dG_min,
                                    homo_dG_min=homo_dG_min, GC_max=GC_max)
        R_oligomers = primers.filter_oligomers(
                                    R_oligomers_cds_map.keys(), tmMin=tmMin,
                                    tmMax=tmMax, hpn_dG_min=hpn_dG_min,
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

        pham_F_oligomer_pos_map = map_oligomer_to_pos(F_oligomers_cds_map,
                                                      F_oligomers,
                                                      avg_start,
                                                      verbose=verbose)

        pham_R_oligomer_pos_map = map_oligomer_to_pos(R_oligomers_cds_map,
                                                      R_oligomers,
                                                      avg_start,
                                                      verbose=verbose)

        F_oligomer_pos_map.update(pham_F_oligomer_pos_map)
        R_oligomer_pos_map.update(pham_R_oligomer_pos_map)

    primer_pairs = []

    F_pos_oligomer_map = invert_dictionary(F_oligomer_pos_map)
    R_pos_oligomer_map = invert_dictionary(R_oligomer_pos_map)

    for pos, F_oligomers in F_pos_oligomer_map.items():
        for i in range(minD, maxD+1):
            R_oligomers = R_pos_oligomer_map.get(pos+i)
            if R_oligomers:
                for F_oligomer in F_oligomers:
                    for R_oligomer in R_oligomers:
                        pair_dict = {}
                        pair_dict["forward"] = F_oligomer
                        pair_dict["reverse"] = R_oligomer
                        pair_dict["length"] = i
                        pair_dict["start"] = pos
                        primer_pairs.append(pair_dict)

    primer_pairs = sorted(primer_pairs, key=lambda pair: pair['start'])
    return primer_pairs


def test_primer_pairs(primer_pairs, genome_map, verbose=False,
                      hetero_dG_min=-1000, Tm_gap_max=5, anneal_gap_max=5):
    tested_primer_pair_TDs = []

    num_pairs = len(primer_pairs)
    for i in range(num_pairs):
        primer_pair = primer_pairs[i]

        for genome_id, genome in genome_map.items():
            primer_pair_TD = primer_TD.PrimerPairTDynamics(
                                                    primer_pair["forward"],
                                                    primer_pair["reverse"],
                                                    genome=str(genome.seq))
        if verbose:
            print(f"...Testing primer pair {i} out of "
                  f"{num_pairs} putative primers...")

        if primer_pair_TD.Tm_gap > Tm_gap_max:
            continue

        if primer_pair_TD.heterodimer.dg < hetero_dG_min:
            continue

        if primer_pair_TD.annealing_Tm_gap > anneal_gap_max:
            continue

        tested_primer_pair_TDs.append(primer_pair_TD)

    return tested_primer_pair_TDs


def select_primer_pair(primer_pair_TDs):
    return primer_pair_TDs[0]


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


def get_oriented_oligomers(conserved_kmer_data, orientation):
    oligomers_cds_map = {}

    for kmer, kmer_data in conserved_kmer_data.items():
        if kmer_data[0][0].orientation != orientation:
            kmer = str(Seq(kmer).reverse_complement())

        oligomers_cds_map[kmer] = kmer_data

    return oligomers_cds_map


def map_oligomer_to_pos(oligomers_cds_map, oligomers, avg_start,
                        verbose=False):
    oligomer_pos_map = OrderedDict()
    for oligomer in oligomers:
        data = oligomers_cds_map[oligomer]
        avg_oligomer_start = 0
        for kmer_data in data:
            avg_oligomer_start += kmer_data[1]
        avg_oligomer_start = int(round(avg_oligomer_start/len(data),
                                       0))

        avg_oligomer_start = avg_start + avg_oligomer_start

        oligomer_pos_map[oligomer] = avg_oligomer_start

    oligomer_pos_map = basic.sort_histogram(oligomer_pos_map)

    return oligomer_pos_map


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
