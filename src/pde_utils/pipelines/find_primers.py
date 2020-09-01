"""Pipeline to find potential primers from given sequences quickly
   using phamilies as an potential indicator of conserved nucleotide
   regions.
   """

import argparse
import sys
import time
from pathlib import Path

from numpy import std
from pdm_utils.functions import annotation
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import mysqldb
from pdm_utils.functions import pipelines_basic
from pdm_utils.pipelines import export_db
from primer3 import calcTm

from pde_utils.classes import kmers

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
                         threads=4, prc=0.8, max_std=500):
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

        parent_genomes = {}
        for genome_id in db_filter.values:
            export_db.get_single_genome(alchemist, genome_id,
                                        data_cache=parent_genomes)

        total_genomes = db_filter.hits()
        if total_genomes == 0:
            print("No database entries received from phage "
                  f" for '{mapped_path.name}'.")

        pham_gene_map = build_pham_gene_map(db_filter, conditionals,
                                            verbose=verbose)

        pham_histogram = {}
        for pham, genes in pham_gene_map.items():
            pham_histogram[pham] = len(genes)
        pham_histogram = basic.sort_histogram(pham_histogram)

        if verbose:
            print(f"Identifying prevelant phams for '{mapped_path.name}'...")

        for pham, count in pham_histogram.items():
            pham_per = count/total_genomes

            if pham_per < prc:
                break

            if verbose:
                print(f"...Analyzing pham {pham}...")
                print(f"......Pham is represented in {round(pham_per, 3)*100}%"
                      " of currently viewed genomes...")

            pham_cds = export_db.parse_feature_data(alchemist, 
                                                    values=pham_gene_map[pham])

            starts = []
            avg_start = 0
            for cds in pham_cds:
                starts.append(cds.start)
                avg_start += cds.start
            avg_start = int(round(avg_start/len(pham_cds), 0))
            start_std = std(starts)

            if start_std > max_std:
                print(f"......Pham {pham} has unstable synteny...")
                continue

            gs_to_seq = get_cds_nucleotide_seqs(alchemist, pham_cds, 
                                                data_cache=parent_genomes)
            conserved_kmer_data = find_conserved_kmers(gs_to_seq)

            num_conserved_kmers = len(conserved_kmer_data)
            if num_conserved_kmers == 0:
                print(f"......Pham {pham} has no conserved kmers...")
                continue

            if verbose:
                print(f"......Pham {pham} has {num_conserved_kmers} "
                      "conserved kmers...")
                
            thermostable_oligomers = find_thermostabilize_oligomers(
                                                        conserved_kmer_data)
            num_stable_oligomers = len(thermostable_oligomers)

            if num_stable_oligomers == 0:
                print(f"......Pham {pham} has no stable oligomers...")
                continue
            
            if verbose:
                print(f"......Pham {pham} has {num_stable_oligomers} "
                      "stable putative primer oligomers.")  

            for oligomer, data in thermostable_oligomers.items():
                avg_oligomer_start = 0
                for kmer_data in data:
                    avg_oligomer_start += kmer_data[1]
                avg_oligomer_start = int(round(avg_oligomer_start/len(data), 
                                               0))

                avg_oligomer_start = avg_start + avg_oligomer_start
            
                print(f".........Oligomer '{oligomer}' is located at "
                      f"{avg_oligomer_start}")
        
        sys.exit(1)
        pipelines_basic.create_working_dir(working_path)


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


def find_conserved_kmers(seq_id_to_sequence_map, kmer_length=20,
                         hash_count=None, fpp=0.0001):
    sequences = list(seq_id_to_sequence_map.values())
    num_sequences = len(sequences)
    approx_count = (len(sequences[0]) - kmer_length) * num_sequences

    cmsketch = kmers.CountMinSketch(approx_count,
                                    hash_count=hash_count, fpp=fpp)

    conserved_kmer_data = {}
    for seq_id, seq in seq_id_to_sequence_map.items():
        for i in range(len(seq) - kmer_length):
            subseq = seq[i:i+kmer_length]
            cmsketch.add(subseq)
            if cmsketch.check(subseq) >= num_sequences:
                kmer_data = conserved_kmer_data.get(subseq, [])
                kmer_data.append((seq_id, i))
                conserved_kmer_data[subseq] = kmer_data

    return conserved_kmer_data


def find_thermostabilize_oligomers(conserved_kmer_data, tmMin=52, tmMax=58):
    thermostable_oligomers = {}
    for oligomer, data in conserved_kmer_data.items(): 
        if not (calcTm(oligomer) < tmMin or calcTm(oligomer) > tmMax):
            thermostable_oligomers[oligomer] = data

    return thermostable_oligomers


def get_count_phams_in_genes(alchemist, geneids, incounts=None):
    phams = annotation.get_phams_from_genes(alchemist, geneids)

    pham_histogram = {}
    if incounts is not None:
        pham_histogram = incounts

    basic.increment_histogram(phams, pham_histogram)
    return pham_histogram


def get_cds_nucleotide_seqs(alchemist, cds_list, data_cache=None):
    if data_cache is None:
        data_cache = {}

    pham_nuc_genes = {}
    for cds in cds_list:
        parent_genome = data_cache.get(cds.genome_id)

        if parent_genome is None:
            parent_genome = export_db.get_single_genome(
                                        alchemist, cds.genome_id,
                                        data_cache=data_cache)

        cds.set_seqfeature()
        cds.set_nucleotide_sequence(
                            parent_genome_seq=parent_genome.seq)

        pham_nuc_genes[cds.id] = str(cds.seq)

    return pham_nuc_genes


if __name__ == "__main__":
    main(sys.argv)
