"""Pipeline to find potential primers from given sequences quickly
   using phamilies as an potential indicator of conserved nucleotide
   regions.
   """

import argparse
import sys
import time
from pathlib import Path

from pdm_utils.functions import annotation
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import fileio
from pdm_utils.functions import pipelines_basic
from pdm_utils.pipelines import export_db

from pde_utils.functions import alignment

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_find_primers")


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
                         threads=4, prc=0.8):
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

        genomes = db_filter.values
        total_genomes = db_filter.hits()
        if total_genomes == 0:
            print("No database entries received from phage "
                  f" for '{mapped_path.name}'.")

        pipelines_basic.create_working_dir(working_path)

        fasta_dir = working_path.joinpath("fasta")
        fasta_dir.mkdir()
        aln_dir = working_path.joinpath("aln")
        aln_dir.mkdir()
        hmmp_dir = working_path.joinpath("hmmp")
        hmmp_dir.mkdir()

        genes = db_filter.transpose("gene.GeneID")
        pham_histogram = get_count_phams_in_genes(alchemist, genes)
        pham_histogram = basic.sort_histogram(pham_histogram)
        pham_col = db_filter.get_column("gene.PhamID")
        phage_col = db_filter.get_column("gene.PhageID")

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

            db_filter.reset()
            db_filter.key = "gene"
            db_filter.values = db_filter.build_values(
                                            where=[(pham_col == pham),
                                                   phage_col.in_(genomes)])

            gs_to_seq = get_pham_nucleotide_genes(alchemist, db_filter.values)

            fasta_path = fasta_dir.joinpath(f"{pham}.fasta")
            fileio.write_fasta(gs_to_seq, fasta_path)

            aln_path = aln_dir.joinpath(f"{pham}.aln")
            alignment.clustalo(fasta_path, aln_path, outfmt="fasta",
                               threads=threads)

            hmmp_path = hmmp_dir.joinpath(f"{pham}.hmm")
            alignment.hmmbuild(aln_path, hmmp_path)


def get_count_phams_in_genes(alchemist, geneids, incounts=None):
    phams = annotation.get_phams_from_genes(alchemist, geneids)

    pham_histogram = {}
    if incounts is not None:
        pham_histogram = incounts

    basic.increment_histogram(phams, pham_histogram)
    return pham_histogram


def get_pham_nucleotide_genes(alchemist, genes):
    cds_list = export_db.parse_feature_data(alchemist, values=genes)

    parent_genomes = {}
    pham_nuc_genes = {}
    for cds in cds_list:
        parent_genome = parent_genomes.get(cds.genome_id)

        if parent_genome is None:
            parent_genome = export_db.get_single_genome(
                                        alchemist, cds.genome_id,
                                        data_cache=parent_genomes)

        cds.set_seqfeature()
        cds.set_nucleotide_sequence(
                            parent_genome_seq=parent_genome.seq)

        pham_nuc_genes[cds.id] = str(cds.seq)

    return pham_nuc_genes


if __name__ == "__main__":
    main(sys.argv)
