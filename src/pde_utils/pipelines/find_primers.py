"""Pipeline to find potential primers from given sequences quickly
                ionserved_kmer_data.append(subseq)
   using phamilies as an potential indicator of conserved nucleotide
   regions.
   """

import argparse
import heapq
import multiprocessing
import pickle
import sys
import time
import math
from pathlib import Path

from Bio.Seq import Seq
from pdm_utils.functions import annotation
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import pipelines_basic
from pdm_utils.pipelines import export_db

from pde_utils.classes import primer3
from pde_utils.functions import fileio
from pde_utils.functions import seq

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
TEMP_DIR = Path("/tmp/pde_utils_find_primers_cache")
PICKLED_FILE_NAME = "PICKLED_PRIMERS"

DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_find_primers")

PHAGE_QUERY = "SELECT * FROM phage"
GENE_QUERY = "SELECT * FROM gene"
TRNA_QUERY = "SELECT * FROM trna"
TMRNA_QUERY = "SELECT * FROM tmrna"


def main(unparsed_args):
    """Uses parsed args to run the entirety of the find primers pipeline.

    :param unparsed_args: Input a list of command line args.
    :type unparsed_args: list[str]
    """
    args = parse_find_primers(unparsed_args)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)

    execute_find_primers(alchemist, folder_path=args.folder_path,
                         folder_name=args.folder_name, values=values,
                         filters=args.filters, groups=args.groups,
                         verbose=args.verbose, threads=args.threads,
                         prc=args.prc, minD=args.minD, maxD=args.maxD,
                         hpn_min=args.hpn_min, ho_min=args.ho_min,
                         het_min=args.het_min, GC_max=args.GC,
                         len_oligomer=args.oligomer_length, tm_min=args.tm_min,
                         tm_max=args.tm_max, tm_gap=args.tm_gap,
                         ta_min=args.ta_min, ta_max=args.ta_max,
                         mode=args.mode, soft_cap=args.soft_cap,
                         phams_in=args.phams_in,
                         fwd_in=args.fwd_in, rvs_in=args.rvs_in)


def parse_find_primers(unparsed_args):
    """Parses find primers arguments and stores them with an argparse object.
    """
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
    NUMBER_THREADS_HELP = """
        Pipeline option that allows for multithreading of workflows.
            Follow selection argument with number of threads to be used
        """

    PHAMS_IN_HELP = """Find primer selection option that narrows the scope
        of the primer pair search to the nucleotide regions of inputted phams
            Follow selection argument with the desired phams"""

    MODE_HELP = """
        Find primer parameter option that changes the characteristics of
        the primer pairs outputted.
        """
    OLIGOMER_LENGTH_HELP = """
        Find primer parameter option that sets the desired length of the
        oligomers used to generate given primers
            Follow parameter argument with the desired length in single bases.
        """
    PHAM_REPRESENTATION_CUTOFF_HELP = """
        Find primer parameter option that changes the pham representation
        prefilter threshold.
            Follow parameter argument with cutoff value [0-1]
        """
    MINIMUM_PRODUCT_LENGTH_HELP = """
        Find primer parameter option that sets the minimum length of the
        PCR product created from a pair of given primers
            Follow parameter argument with the desired length in single bases.
        """
    MAXIMUM_PRODUCT_LENGTH_HELP = """
        Find primer parameter option that sets the minimum length of the
        PCR product created from a pair of given primers
            Follow parameter argument with the desired length in single bases.
        """
    HAIRPIN_GIBBS_FREE_ENERGY_MINIMUM_HELP = """
        Find primer parameter option that sets the minimum threshold for the
        Gibbs free energy of hairpin formation for any oligomer in a pair
        of given primers.
            Follow parameter argument with the desired minimum kcal/mol.
        """
    HOMODIMER_GIBBS_FREE_ENERGY_MINIMUM_HELP = """
        Find primer parameter option that sets the minimum threshold for the
        Gibbs free energy of homodimer formation for any oligomer in a pair
        of given primers.
            Follow parameter argument with the desired minimum kcal/mol.
        """
    HETERODIMER_GIBBS_FREE_ENERGY_MINIMUM_HELP = """
        Find primer parameter option that sets the minimum threshold for the
        Gibbs free energy of heterodimer formation for the oligomers in a pair
        of given primers.
            Follow parameter argument with the desired minimum kcal/mol.
        """
    MAXIMUM_GC_CONTENT_HELP = """
        Find primer parameter option that sets the maximum threshold for the
        GC content percentage of any oligomer in a pair of given primers.
            Follow parameter argument with the desired percentage [0-100].
        """
    MINIMUM_MELTING_TEMPERATURE_HELP = """
        Find primer parameter option that sets the minimum melting temperature
        threshold of any oligomer in a pair of given primers.
            Follow parameter argument with the desired temperature in Celsius.
        """
    MAXIMUM_MELTING_TEMPERATURE_HELP = """
        Find primer parameter option that sets the maximum melting temperature
        threshold of any oligomer in a pair of given primers.
            Follow parameter argument with the desired temperature in Celsius.
        """
    MAXIMUM_MELTING_TEMPERATURE_GAP_HELP = """
        Find primer parameter option that sets the maximum melting temperature
        gap threshold between the oligomers in a pair of given primers.
            Follow parameter argument with the desired temperature gap in
            Celsius.
        """
    MAXIMUM_OPTIMAL_ANNEALING_TEMPERATURE_HELP = """
        Find primer parameter option that sets the maximum optimal annealing
        temperature threshold of a any pair of given primers.
            Follow parameter argument with the desired temperature in Celsius.
        """
    MINIMUM_OPTIMAL_ANNEALING_TEMPERATURE_HELP = """
        Find primer parameter option that sets the maximum optimal annealing
        temperature threshold of a any pair of given primers.
            Follow parameter argument with the desired temperature in Celsius.
        """
    START_DEVIATION_NET_HELP = """
        Find primer parameter option that expands the prefiltered search radius
        during primer matching steps
            Follow parameter argument with the desired net length in base pairs
        """
    SOFT_CAP_HELP = """
        Find primer parameter option that restricts the amount
        of primer pairs to be evaluated after testing, to limit memory usage
            Follow parameter argument with the desired number of pairs
        """
    FORWARD_IN_HELP = """
        Find primer parameter option that allows for the manual selection
        of the forward primer oligomer sequence
            Follow parameter argument with the desired sequence
        """
    REVERSE_IN_HELP = """
        Find primer parameter option that allows for the manual selection
        of the reverse primer oligomer sequence
            Follow parameter argument with the desired sequence
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
    parser.add_argument("-th", "--threads", type=int, help=NUMBER_THREADS_HELP)

    parser.add_argument("-if", "--import_file", dest="input",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP)
    parser.add_argument("-in", "--import_names", nargs="*", dest="input",
                        help=SINGLE_GENOMES_HELP)
    parser.add_argument("-w", "--where", nargs="?", dest="filters",
                        help=WHERE_HELP)
    parser.add_argument("-g", "--group_by", nargs="*", dest="groups",
                        help=GROUP_BY_HELP)

    parser.add_argument("-phin", "--phams_in", nargs="*", type=int,
                        help=PHAMS_IN_HELP)

    parser.add_argument("-mo", "--mode", type=int, help=MODE_HELP)
    parser.add_argument("-prc", type=float,
                        help=PHAM_REPRESENTATION_CUTOFF_HELP)
    parser.add_argument("-minD", type=int, help=MINIMUM_PRODUCT_LENGTH_HELP)
    parser.add_argument("-maxD", type=int, help=MAXIMUM_PRODUCT_LENGTH_HELP)
    parser.add_argument("-hpn_min", type=int,
                        help=HAIRPIN_GIBBS_FREE_ENERGY_MINIMUM_HELP)
    parser.add_argument("-ho_min", type=int,
                        help=HOMODIMER_GIBBS_FREE_ENERGY_MINIMUM_HELP)
    parser.add_argument("-het_min", type=int,
                        help=HETERODIMER_GIBBS_FREE_ENERGY_MINIMUM_HELP)
    parser.add_argument("-GC", type=float, help=MAXIMUM_GC_CONTENT_HELP)
    parser.add_argument("-dev_net", type=int, help=START_DEVIATION_NET_HELP)
    parser.add_argument("-len", "--oligomer_length", type=int,
                        help=OLIGOMER_LENGTH_HELP)
    parser.add_argument("-tm_min", type=float,
                        help=MINIMUM_MELTING_TEMPERATURE_HELP)
    parser.add_argument("-tm_max", type=float,
                        help=MAXIMUM_MELTING_TEMPERATURE_HELP)
    parser.add_argument("-tm_gap", type=float,
                        help=MAXIMUM_MELTING_TEMPERATURE_GAP_HELP)
    parser.add_argument("-ta_min", type=float,
                        help=MAXIMUM_OPTIMAL_ANNEALING_TEMPERATURE_HELP)
    parser.add_argument("-ta_max", type=float,
                        help=MINIMUM_OPTIMAL_ANNEALING_TEMPERATURE_HELP)
    parser.add_argument("-sc", "--soft_cap", type=int, help=SOFT_CAP_HELP)
    parser.add_argument("--fwd_in", type=str, help=FORWARD_IN_HELP)
    parser.add_argument("--rvs_in", type=str, help=REVERSE_IN_HELP)

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, verbose=False, input=[], threads=1,
                        filters="", groups=[], prc=0.70, minD=900, maxD=1100,
                        hpn_min=-2000, ho_min=-5000, het_min=-5000, GC=60.0,
                        oligomer_length=20, tm_min=52.0, tm_max=58, dev_net=0,
                        tm_gap=5.0, ta_min=48.0, ta_max=68.0, mode=0,
                        soft_cap=None, phams_in=[], fwd_in=None, rvs_in=None)

    parsed_args = parser.parse_args(unparsed_args[2:])
    return parsed_args


def execute_find_primers(alchemist, folder_path=None,
                         folder_name=DEFAULT_FOLDER_NAME, values=None,
                         filters="", groups=[], verbose=False,
                         threads=4, prc=0.7, dev_net=0, len_oligomer=20,
                         minD=900, maxD=1100, tm_min=52.0, tm_max=58.0,
                         hpn_min=-2000, ho_min=-5000, GC_max=60.0,
                         het_min=-5000, tm_gap=5.0, ta_min=48.0,
                         fwd_in=None, rvs_in=None,
                         ta_max=68.0, mode=0, full_genome=False,
                         soft_cap=None, phams_in=[]):
    """Executes the entirety of the file export pipeline.

    :param alchemist: A connected and fully build AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param folder_path: Path
    :type folder_path: Path
    :param folder_name: A name for the working directory folder
    :type folder_name: str
    :param values: List of values to filter database results
    :type values: list[str]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param filters: A pseudo-SQL WHERE clause string to filter values.
    :type filters: str
    :param groups: A list of SQL column names to filter values.
    :type groups: list[str]
    :param threads: Number of child process workers to utilize
    :type threads: int
    :param prc: Percentage of genomes a pham must exist in to pass prefiltering
    :type prc: float
    :param dev_net: Allowance for the primer positions to pass prefiltering
    :type dev_net: int
    :param len_oligomer: Length of the oligomers used to create the primers
    :type len_oligomer: int
    :param minD: Minimum primer product length to pass primer testing
    :type minD: int
    :param maxD: Maximum primer product length to pass primer testing
    :type maxD: int
    :param tm_min: Minimum primer melting temperature to pass primer testing
    :type tm_min: float
    :param tm_max: Maximum primer melting temperature to pass primer testing
    :type tm_max: float
    :param hpn_min: Minimum hairpin Gibbs free energy to pass primer testing
    :type hpn_min: int
    :param ho_min: Minimum homodimer Gibbs free energy to pass primer testing
    :type ho_min: int
    :param GC_max: Maximum GC content percentage allowed for an oligomer
    :type GC_max: float
    :param het_min: Minimum heterodimer Gibbs free energy to pass testing
    :type het_min: int
    :param tm_gap: Maximum allowed melting temperature gap between oligomers
    :type tm_gap: float
    :param ta_min: Minimum allowed optimal annealing temperature
    :type ta_min: float
    :param ta_max: Maximum allowed optimal annealing temperature
    :type ta_max: float
    :param fwd_in: Fixed forward sequence to find primer pairs for
    :type fwd_in: str
    :param rvs_in: Fixed reverse sequence to find primer pairs for
    :type rvs_in: str
    :param mode: Run mode for find primers analysis
    :type mode: int
    :param soft_cap: Cap limit on number of pairs evaluated after testing
    :type soft_cap: int
    :param phams_in: Phams to evaluate during count min sketch eval of kmers
    :type phams_in: list[str]
    """
    db_filter = pipelines_basic.build_filter(alchemist, "phage", filters)

    working_path = pipelines_basic.create_working_path(
                                                    folder_path, folder_name)

    conditionals_map = pipelines_basic.build_groups_map(
                                                    db_filter, working_path,
                                                    groups=groups,
                                                    verbose=verbose)

    if verbose:
        print("Prepared query and path structure, beginning primer search...")

    if not TEMP_DIR.is_dir():
        TEMP_DIR.mkdir()
    pickled_results_file = TEMP_DIR.joinpath(PICKLED_FILE_NAME)
    if pickled_results_file.is_file():
        pickled_results_file.unlink()

    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        results_map = {}
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

        if verbose:
            print(f"...Identifying primer pairs for '{mapped_path}'...")

        if full_genome:
            F_results, R_results = find_full_genome_oligomers(
                                    genome_map,
                                    verbose=verbose, threads=threads, prc=prc,
                                    minD=minD, maxD=maxD,
                                    len_oligomer=len_oligomer, tm_min=tm_min,
                                    tm_max=tm_max, hpn_min=hpn_min,
                                    ho_min=ho_min, GC_max=GC_max)
        else:
            pham_gene_map = build_pham_gene_map(db_filter, conditionals,
                                                phams_in=phams_in,
                                                verbose=verbose)
            if not pham_gene_map:
                print(f"No valid phams found for '{mapped_path}' with current "
                      "settings")

            F_results, R_results = find_oligomers(
                                    alchemist, pham_gene_map, genome_map,
                                    verbose=verbose, threads=threads, prc=prc,
                                    minD=minD, maxD=maxD,
                                    len_oligomer=len_oligomer, tm_min=tm_min,
                                    tm_max=tm_max, hpn_min=hpn_min,
                                    ho_min=ho_min, GC_max=GC_max,
                                    fwd_in=fwd_in, rvs_in=rvs_in)

        if (not F_results) or (not R_results):
            if verbose:
                print(f"No valid oligomers found for '{mapped_path.name}'")
            continue

        if verbose:
            print("...Matching oligomers to create primer pairs...")

        primer_pairs = match_oligomers(F_results, R_results,
                                       minD=minD, maxD=maxD, dev_net=dev_net,
                                       threads=threads)

        if not primer_pairs:
            print(f"No valid primer pairs found for '{mapped_path}' with "
                  "current parameters...")
            continue

        if verbose:
            print(f"...Identified {len(primer_pairs)} valid primer pairs.")

        if verbose:
            print(f"...Testing primer pairs for '{mapped_path}'...")
        primer_pairs = test_primer_pairs(primer_pairs, genome_map,
                                         threads=threads, verbose=verbose,
                                         het_min=het_min,
                                         ta_min=ta_min, ta_max=ta_max,
                                         tm_gap_max=tm_gap)

        if verbose:
            print(f"...{len(primer_pairs)} passed primer testing.")

        if soft_cap is not None:
            if len(primer_pairs) > soft_cap:
                primer_pairs = primer_pairs[:soft_cap]

        if pickled_results_file.is_file():
            with pickled_results_file.open(mode="rb") as filehandle:
                results_map = pickle.load(filehandle)

        if primer_pairs:
            results_map[mapped_path] = (primer_pairs, genome_map)

        with pickled_results_file.open(mode="wb") as filehandle:
            pickle.dump(results_map, filehandle)

    if pickled_results_file.is_file():
        pickled_results_file.unlink()

    if not results_map:
        print("No primer pairs found with current parameters...")

    results_map = select_primer_pairs(results_map, verbose=verbose, mode=mode,
                                      het_min=het_min)

    for mapped_path, primer_pairs in results_map.items():
        pipelines_basic.create_working_dir(mapped_path)
        file_path = mapped_path.joinpath("primer.txt")
        fileio.write_primer_txt_file(primer_pairs[0][0], file_path)


def find_oligomers(alchemist, pham_gene_map, genome_map,
                   verbose=False, threads=4, prc=0.8,
                   len_oligomer=20, minD=900, maxD=1100, tm_min=52, tm_max=58,
                   hpn_min=-2000, ho_min=-5000, GC_max=60, fwd_in=None,
                   rvs_in=None):
    """

    :param alchemist: A connected and fully build AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param threads: Number of child process workers to utilize
    :type threads: int
    :param prc: Percentage of genomes a pham must exist in to pass prefiltering
    :type prc: float
    :param len_oligomer: Length of the oligomers used to create the primers
    :type len_oligomer: int
    :param minD: Minimum primer product length to pass primer testing
    :type minD: int
    :param maxD: Maximum primer product length to pass primer testing
    :type maxD: int
    :param tm_min: Minimum primer melting temperature to pass primer testing
    :type tm_min: float
    :param tm_max: Maximum primer melting temperature to pass primer testing
    :type tm_max: float
    :param hpn_min: Minimum hairpin Gibbs free energy to pass primer testing
    :type hpn_min: int
    :param ho_min: Minimum homodimer Gibbs free energy to pass primer testing
    :type ho_min: int
    :param GC_max: Maximum GC content percentage allowed for an oligomer
    :type GC_max: float
    :param fwd_in: Fixed forward sequence to find primer pairs for
    :type fwd_in: str
    :param rvs_in: Fixed reverse sequence to find primer pairs for
    :type rvs_in: str
    :returns: Returns forward and reverse oligomers mapped to a position
    :rtype: list[(int, list[Oligomer])]
    """
    pham_histogram = {}
    for pham, genes in pham_gene_map.items():
        pham_histogram[pham] = len(genes)
    pham_histogram = basic.sort_histogram(pham_histogram)

    thread_pool = multiprocessing.Pool(processes=threads)
    thread_manager = multiprocessing.Manager()

    managed_genome_map = thread_manager.dict()
    managed_genome_map.update(genome_map)

    work_items = []
    F_pos_oligomer_map = {}
    R_pos_oligomer_map = {}

    for pham, count in pham_histogram.items():
        pham_per = count/len(genome_map)

        if pham_per < prc:
            break

        if verbose:
            print(f".........Pham {pham} is represented in "
                  f"{round(pham_per, 3)*100}% of currently viewed genomes...")

        cds_list = export_db.parse_feature_data(alchemist,
                                                values=pham_gene_map[pham])
        work_items.append((pham, cds_list))

    results = []
    for work_bundle in work_items:
        results.append(thread_pool.apply_async(
                process_find_oligomers, args=(work_bundle, managed_genome_map,
                                              len_oligomer, minD,
                                              maxD, tm_min, tm_max, hpn_min,
                                              ho_min, GC_max, fwd_in,
                                              rvs_in, verbose)))

    for result in results:
        F_pos_map, R_pos_map = result.get()
        F_pos_oligomer_map.update(F_pos_map)
        R_pos_oligomer_map.update(R_pos_map)

    thread_pool.close()
    thread_pool.join()

    if (fwd_in is not None) or (rvs_in is not None):
        genome_sequences = []
        for genome_id, genome_obj in genome_map.items():
            genome_sequences.append(str(genome_obj.seq))

        if fwd_in is not None:
            try:
                fwd_start = find_oligomer_pos(fwd_in, genome_sequences, "F")
            except ValueError:
                F_pos_oligomer_map = {}
                print("Forward in sequence has an undesired number of "
                      "positions within the given genome sequences")
            else:
                F_pos_oligomer_map = {fwd_start: [primer3.Oligomer(
                                                    fwd_in, start=fwd_start)]}

        if rvs_in is not None:
            try:
                rvs_start = find_oligomer_pos(rvs_in, genome_sequences, "R")
            except ValueError:
                R_pos_oligomer_map = {}
                print("Reverse in sequence has an undesired number of "
                      "positions within the given genome sequences")
            else:
                R_pos_oligomer_map = {rvs_start: [primer3.Oligomer(
                                                    rvs_in, start=rvs_start)]}

    F_results = []
    for pos, oligomers in F_pos_oligomer_map.items():
        F_results.append((pos, list(oligomers)))

    R_results = []
    for pos, oligomers in R_pos_oligomer_map.items():
        R_results.append((pos, list(oligomers)))

    return F_results, R_results


def find_full_genome_oligomers(
                   genome_map, verbose=False, threads=4, prc=0.8,
                   len_oligomer=20, minD=900, maxD=1100, tm_min=52, tm_max=58,
                   hpn_min=-2000, ho_min=-5000, GC_max=60):
    raise NotImplementedError(
                    "Full genome oligomer search is not yet supported")
    thread_pool = multiprocessing.Pool(processes=threads)
    thread_manager = multiprocessing.Manager()

    managed_genome_map = thread_manager.dict()
    managed_genome_map.update(genome_map)

    work_items = []
    F_pos_oligomer_map = {}
    R_pos_oligomer_map = {}

    work_items.append([])

    results = []
    for work_bundle in work_items:
        results.append(thread_pool.apply_async(
                process_find_oligomers, args=(work_bundle, managed_genome_map,
                                              len_oligomer, minD,
                                              maxD, tm_min, tm_max, hpn_min,
                                              ho_min, GC_max, verbose)))

    for result in results:
        F_pos_map, R_pos_map = result.get()
        F_pos_oligomer_map.update(F_pos_map)
        R_pos_oligomer_map.update(R_pos_map)

    thread_pool.close()
    thread_pool.join()

    F_results = []
    for pos, oligomers in F_pos_oligomer_map.items():
        F_results.append((pos, list(oligomers)))

    R_results = []
    for pos, oligomers in R_pos_oligomer_map.items():
        R_results.append((pos, list(oligomers)))

    return F_results, R_results
    pass


def match_oligomers(F_oligomer_results, R_oligomer_results,
                    minD=900, maxD=1100, dev_net=0, threads=4):
    thread_pool = multiprocessing.Pool(processes=threads)
    thread_manager = multiprocessing.Manager()

    managed_R_pos_oligomer_map = thread_manager.dict()
    managed_R_pos_oligomer_map.update(dict(R_oligomer_results))

    chunk_size = math.ceil(math.sqrt(len(F_oligomer_results)))
    work_chunks = basic.partition_list(F_oligomer_results, chunk_size)

    results = []
    for work_items in work_chunks:
        results.append(thread_pool.apply_async(
            process_match_oligomers, args=(
                                    work_items, managed_R_pos_oligomer_map,
                                    minD, maxD, dev_net)))

    primer_pairs = []
    for result in results:
        primer_pairs = primer_pairs + result.get()

    thread_pool.close()
    thread_pool.join()

    return primer_pairs


def test_primer_pairs(primer_pairs, genome_map, verbose=False,
                      threads=4, het_min=-5000, tm_gap_max=5.0,
                      ta_min=48.0, ta_max=68.0, minD=900, maxD=1100):
    thread_pool = multiprocessing.Pool(processes=threads)
    thread_manager = multiprocessing.Manager()

    tested_primer_pairs = thread_manager.list()

    managed_genome_map = thread_manager.dict()
    managed_genome_map.update(genome_map)

    chunk_size = math.ceil(math.sqrt(len(primer_pairs)))
    work_chunks = basic.partition_list(primer_pairs, chunk_size)

    results = []
    for work_items in work_chunks:
        results.append(thread_pool.apply_async(
            process_test_primer_pairs, args=(
                                 work_items, managed_genome_map, minD,
                                 maxD, tm_gap_max, het_min, ta_min, ta_max)))
    total_run_info = [0] * 4
    tested_primer_pairs = []
    for result in results:
        pair_results, run_info = result.get()

        tested_primer_pairs = tested_primer_pairs + pair_results
        for i in range(len(run_info)):
            total_run_info[i] += run_info[i]

    thread_pool.close()
    thread_pool.join()

    if verbose:
        print(f"......{total_run_info[0]} primer "
              "pairs formed no products on at least one genome")
        print(f"......{total_run_info[1]} primer "
              "pairs formed multiple products on at least one genome")
        print(f"......{total_run_info[2]} primer pairs "
              "formed product with incorrect lengths on at least one genome")
        print(f"......{total_run_info[3]} primer pairs "
              "failed thermodynamic checks")

    tested_primer_pairs = heapq.merge(tested_primer_pairs,
                                      key=lambda pair: pair.rating,
                                      reverse=True)

    return list(tested_primer_pairs)


def select_primer_pairs(results_map, verbose=False, mode=0, het_min=-5000):
    if mode == 0:
        return results_map

    if mode >= 0:
        working_map = results_map
        primer_solution = False
        while not primer_solution:
            if mode >= 1:
                temp_map = {}
                for mapped_path, data in working_map.items():
                    if verbose:
                        print(f"Selecting primers for '{mapped_path}'...")

                    valid_primer_num = None
                    for i in range(len(data[0])):
                        primer_pair = data[0][i]

                        valid_pair = True
                        for i_mapped_path, i_data in working_map.items():
                            if mapped_path == mapped_path:
                                continue

                            for genome_id, genome in i_data[1]:
                                try:
                                    primer_pair.genome = str(genome.seq)
                                    primer_pair.product
                                    valid_pair = False
                                except ValueError:
                                    pass

                            if not valid_pair:
                                break

                        if valid_pair:
                            valid_primer_num = i
                            break

                    if valid_primer_num is not None:
                        temp_map[mapped_path] = (data[0][valid_primer_num:],
                                                 data[1])
                    else:
                        print(f"All primer pairs for '{mapped_path}' generate "
                              "products with other nucleotide reagents")
                        sys.exit(1)

            working_map = temp_map
            primer_solution = True
            if mode >= 2:
                for mapped_path, data in working_map.items():
                    curr_primer = data[0][0]

                    for i_mapped_path, i_data in working_map.items():
                        if mapped_path == i_mapped_path:
                            continue

                        i_primer = data[0][0]

                        for curr_oligo in [curr_primer.fwd, curr_primer.rvs]:
                            for i_oligo in [i_primer.fwd, i_primer.rvs]:
                                primer_dimer = primer3.Heterodimer(
                                            curr_oligo.seq, i_oligo.seq)

                                if primer_dimer.dg < het_min:
                                    primer_solution = False
                                    break

                            if not primer_solution:
                                break

                    if not primer_solution:
                        break

        return working_map


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


def find_oligomer_pos(subsequence, sequences, orientation):
    if orientation == "R":
        match_seq = str(Seq(subsequence).reverse_complement())
    else:
        match_seq = subsequence

    avg_start = 0
    for sequence in sequences:
        starts = seq.find_subsequence_starts(match_seq, sequence,
                                             zero_index=False)

        if len(starts) != 1:
            raise ValueError

        avg_start += starts[0]

    avg_start = int(round((avg_start/len(sequences)), 0))

    return avg_start


def process_find_oligomers(work_bundle, work_data_cache, len_oligomer,
                           minD, maxD, tm_min, tm_max, hpn_min, ho_min,
                           GC_max, fwd_in, rvs_in, verbose):
    pham = work_bundle[0]
    cds_list = work_bundle[1]

    avg_orientation = 0
    for cds in cds_list:
        if cds.orientation == "F":
            avg_orientation += 1
    avg_orientation = avg_orientation/len(cds_list)
    if avg_orientation != 0 and avg_orientation != 1:
        if verbose:
            print(f"...Pham {pham} has inconsistent orientation")
        return {}, {}

    cds_to_seq = seq.map_cds_to_seq(None, cds_list,
                                    data_cache=work_data_cache)
    conserved_kmer_data = seq.find_conserved_kmers(cds_to_seq,
                                                   len_oligomer)

    num_conserved_kmers = len(conserved_kmer_data)
    if num_conserved_kmers == 0:
        if verbose:
            print(f"...Pham {pham} has no conserved kmers")
        return {}, {}

    if verbose:
        print(f"...Pham {pham} has {num_conserved_kmers} "
              "conserved kmers...")

    if fwd_in is not None:
        F_oligomers = []
    else:
        F_oligomers = get_stable_oligomers(
                                    conserved_kmer_data, work_data_cache, "F",
                                    tm_min=tm_min, tm_max=tm_max,
                                    hpn_min=hpn_min,
                                    ho_min=ho_min,
                                    GC_max=GC_max)

    if rvs_in is not None:
        R_oligomers = []
    else:
        R_oligomers = get_stable_oligomers(
                                    conserved_kmer_data, work_data_cache, "R",
                                    tm_min=tm_min, tm_max=tm_max,
                                    hpn_min=hpn_min,
                                    ho_min=ho_min,
                                    GC_max=GC_max)

    num_stable_oligomers = len(F_oligomers) + len(R_oligomers)

    if num_stable_oligomers == 0:
        if verbose:
            print(f"...Pham {pham} has no stable oligomers")
        return {}, {}

    if verbose:
        if fwd_in is None:
            print(f"......Pham {pham} has {len(F_oligomers)} "
                  "stable forward putative primer oligomers.")
        if rvs_in is None:
            print(f"......Pham {pham} has {len(R_oligomers)} "
                  "stable reverse putative primer oligomers.")

    pham_F_pos_oligomer_map = map_pos_to_oligomer(F_oligomers,
                                                  verbose=verbose)
    pham_R_pos_oligomer_map = map_pos_to_oligomer(R_oligomers,
                                                  verbose=verbose)

    return pham_F_pos_oligomer_map, pham_R_pos_oligomer_map


def process_match_oligomers(work_items, reverse_position_map, minD, maxD,
                            dev_net):
    primer_pairs = []
    for work_bundle in work_items:
        pos = work_bundle[0]
        F_oligomers = work_bundle[1]

        for i in range((minD-dev_net), (maxD+1+dev_net)):
            R_oligomers = reverse_position_map.get(pos+i)
            if R_oligomers:
                for F_oligomer in F_oligomers:
                    for R_oligomer in R_oligomers:
                        primer_pair = primer3.PrimerPair(
                                                F_oligomer, R_oligomer)
                        primer_pairs.append(primer_pair)

    return primer_pairs


def process_test_primer_pairs(work_items, genome_map, minD, maxD,
                              tm_gap_max, het_min, ta_min, ta_max):
    failed_no_product = 0
    failed_many_product = 0
    failed_length = 0
    failed_thermo = 0

    pair_results = []
    for primer_pair in work_items:
        valid_primers = True
        for genome_id, genome in genome_map.items():
            primer_pair.genome = str(genome.seq)

            try:
                primer_pair.set_product()
            except primer3.NoProductError:
                valid_primers = False
                failed_no_product += 1
                break
            except primer3.MultipleProductError:
                valid_primers = False
                failed_many_product += 1

            product_len = len(primer_pair.product)
            if product_len < minD or product_len > maxD:
                valid_primers = False
                failed_length += 1
                break

            valid_primers = (
                    (primer_pair.Tm_gap < tm_gap_max) and
                    (primer_pair.heterodimer.dg >= het_min) and
                    (primer_pair.annealing_ta > ta_min) and
                    (primer_pair.annealing_ta < ta_max))

            if not valid_primers:
                failed_thermo += 1
                break

        if valid_primers:
            pair_results.append(primer_pair)

    pair_results = sorted(pair_results, key=lambda pair: pair.rating,
                          reverse=True)

    return (pair_results, (failed_no_product, failed_many_product,
                           failed_length, failed_thermo))


def get_stable_oligomers(conserved_kmer_data, genome_map, orientation,
                         tm_min=52, tm_max=58, hpn_min=-2000,
                         ho_min=-5000, GC_max=60):
    oligomers = []

    genome_sequences = []
    for genome_id, genome_obj in genome_map.items():
        genome_sequences.append(str(genome_obj.seq))

    for kmer, kmer_data in conserved_kmer_data.items():
        if kmer_data[0][0].orientation != orientation:
            kmer = str(Seq(kmer).reverse_complement())

        oligomer = primer3.Oligomer(kmer)

        stable = ((oligomer.Tm >= tm_min and oligomer.Tm <= tm_max) and
                  (oligomer.hairpin.dg > hpn_min) and
                  (oligomer.GC < GC_max) and
                  (oligomer.homodimer.dg > ho_min) and
                  (not oligomer.base_run))

        if not stable:
            continue

        try:
            avg_start = find_oligomer_pos(kmer, genome_sequences, orientation)
        except ValueError:
            continue

        oligomer.start = avg_start
        oligomers.append(oligomer)

    return oligomers


def map_pos_to_oligomer(oligomers, verbose=False):
    pos_oligomer_map = {}
    for oligomer in oligomers:
        mapped_oligomers = pos_oligomer_map.get(oligomer.start, [])
        mapped_oligomers.append(oligomer)
        pos_oligomer_map[oligomer.start] = mapped_oligomers

    return pos_oligomer_map


# TO RELOCATE
def build_pham_gene_map(db_filter, conditionals, phams_in=[], verbose=False):
    gene_col = db_filter.get_column("gene.GeneID")
    pham_col = db_filter.get_column("gene.PhamID")

    phams = db_filter.transpose("gene.PhamID")

    if phams_in:
        phams = list(set(phams_in).intersection(set(phams)))

    pham_gene_map = {}
    for pham in phams:
        pham_conditionals = conditionals + [(pham_col == pham)]

        gene_values = db_filter.build_values(column=gene_col,
                                             where=pham_conditionals)

        if gene_values:
            pham_gene_map[pham] = gene_values

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
