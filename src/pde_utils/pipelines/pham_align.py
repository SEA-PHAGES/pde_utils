"""Pipeline for exporting aligned pham sequences."""
import argparse
import time
from pathlib import Path

from pdm_utils.functions import configfile
from pdm_utils.functions import pipelines_basic

from pde_utils.functions import alignment


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pham_align"


# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def main(unparsed_args_list):
    """Uses parsed_args to run the entirety of the pham align pipeline.

    :param unparsed_args_list: Input a list of command line args.
    :type unparsed_args_list: list[str]
    """
    args = parse_pham_align(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)

    execute_pham_align(alchemist, folder_path=args.folder_path,
                       folder_name=args.folder_name, values=values,
                       filters=args.filters, groups=args.groups,
                       file_type=args.file_type,
                       mat_out=args.distmat_out, tree_out=args.guidetree_out,
                       verbose=args.verbose, dump=args.dump, force=args.force,
                       threads=args.number_threads)


def parse_pham_align(unparsed_args_list):
    DATABASE_HELP = """
        Name of the MySQL database to align phams from
        """

    CONFIG_FILE_HELP = """
        Pham align option that enables use of a config file for sourcing
        credentials
            Follow selection argument with the path to the config file
            specifying MySQL and NCBI credentials.
        """
    DUMP_HELP = """
        Pham align option that dumps exported files directly to the desired
        working directory.
        """
    FORCE_HELP = """
        Pham align option that aggresively creates and overwrites directories.
        """
    VERBOSE_HELP = """
        Pham align option that enables progress print statements
        """
    FOLDER_PATH_HELP = """
        Pham align option to change the path of the directory where the
        exported files are stored.
            Follow selection argument with the path to the desired working
            directory.
        """
    FOLDER_NAME_HELP = """
        Pham align option to change the name of the directory where the
        exported files are stored.
            Follow selection argument with the desired name.
        """
    NUMBER_THREADS_HELP = """
        Pipeline option that allows for multithreading of workflows.
            Follow selection argument with number of threads to be used
        """

    IMPORT_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selectionn argument with path to the csv file containing
            the names of each genome in the first column.
        """
    SINGLE_PHAMS_HELP = """
        Selection input option that imports values from cmd line input.
            Follow seelctiona rgument with space separated phams in the
            database
        """
    WHERE_HELP = """
        Data filtering option that filters data by the inputted expressions.
            Follow selection argument with formatted filter expression:
                {Table}.{Column}={Value}
        """
    GROUP_BY_HELP = """
        Data selection option that sorts data by the inputted columns.
            Follow selection argument with formatted column expressions:
                {Table}.{Column}={Value}
        """

    OUTFILE_TYPE_HELP = """
        Pham align option to change the format of alignment file exported.
            Follow selection argument with a supported alignment file type.
        """
    DISTMAT_OUT_HELP = """
        Pham align option to toggle on pham distance matrix file generation.
        """
    GUIDETREE_OUT_HELP = """
        Pham align option to toggle on pham guidetree file generation.
        """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str, help=DATABASE_HELP)

    parser.add_argument("-c", "--config_file", help=CONFIG_FILE_HELP,
                        type=pipelines_basic.convert_file_path)
    parser.add_argument("-m", "--folder_name", type=str, help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path",  type=Path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)
    parser.add_argument("-d", "--dump", action="store_true", help=DUMP_HELP)
    parser.add_argument("-f", "--force", action="store_true", help=FORCE_HELP)
    parser.add_argument("-np", "--number_threads", type=int,
                        help=NUMBER_THREADS_HELP)

    parser.add_argument("-if", "--import_file", help=IMPORT_FILE_HELP,
                        type=pipelines_basic.convert_file_path, dest="input")
    parser.add_argument("-in", "--import_names", help=SINGLE_PHAMS_HELP,
                        nargs="*", dest="input")
    parser.add_argument("-w", "--where", nargs="?", help=WHERE_HELP,
                        dest="filters")
    parser.add_argument("-g", "--group_by", nargs="*", help=GROUP_BY_HELP,
                        dest="groups")

    parser.add_argument("-ft", "--file_type", type=str,
                        choices=alignment.CLUSTALO_FORMATS,
                        help=OUTFILE_TYPE_HELP)
    parser.add_argument("-mat", "--distmat_out", action="store_true",
                        help=DISTMAT_OUT_HELP)
    parser.add_argument("-tree", "--guidetree_out", action="store_true",
                        help=GUIDETREE_OUT_HELP)

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, input=[], filters="", groups=[],
                        file_type="fasta", number_threads=1)

    args = parser.parse_args(unparsed_args_list[2:])
    return args


def execute_pham_align(alchemist, folder_path=None,
                       folder_name=DEFAULT_FOLDER_NAME, values=None,
                       filters="", groups=[],
                       file_type="fasta", mat_out=False, tree_out=False,
                       threads=1, verbose=False, dump=False, force=False):
    """Executes the entirety of the pham align pipeline.
    :param alchemist: A connected and fully built AlchemyHandler object.
    :type alchemist: AlchemyHandler
    :param folder_path: Path to a valid dir for working dir creation.
    :type folder_path: Path
    :param folder_name: A name for the working directory.
    :type folder_name: str
    :param force: A boolean to toggle aggresive building of directories.
    :type force: bool
    :param values: List of values to filter database results.
    :type values: list[str]
    :param verbose: A boolean value to toggle progress print statements.
    :type verbose: bool
    :param dump: A boolean value to toggle dump in current working dir.
    :type dump: bool
    :param filters: A MySQL formatted WHERE clause string
    :type filters: str
    :param groups: A list of supported MySQL column names to group by.
    :type groups: list[str]
    :param file_type: Format type of sequence alignment file to export.
    :type file_type: str
    :param mat_out: A boolean to toggle distance matrix file generation.
    :type mat_out: bool
    :param tree_out: A boolean to toggle guidetree file generation.
    :type tree_out: bool
    :param threads: Number of processes to spawn during alignment workflow.
    :type threads: int
    """
    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                             values=values, verbose=verbose)
    working_path = pipelines_basic.create_working_path(
                                                    folder_path, folder_name,
                                                    dump=dump, force=force)

    data_cache = {}
    conditionals_map = pipelines_basic.build_groups_map(
                                                db_filter, working_path,
                                                groups=groups, verbose=verbose,
                                                force=force)
    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        db_filter.reset()
        db_filter.values = values

        conditionals = conditionals_map[mapped_path]
        db_filter.values = db_filter.build_values(where=conditionals)

        if db_filter.hits() == 0:
            print(f"No database entries received for '{mapped_path}'")
            continue

        pipelines_basic.create_working_dir(mapped_path, dump=dump, force=force)

    execute_pham_MSA_alignment(alchemist, mapped_path,
                               db_filter.values, data_cache=data_cache,
                               file_type=file_type,
                               mat_out=mat_out, tree_out=tree_out,
                               threads=threads, verbose=verbose)


def execute_pham_MSA_alignment(alchemist, working_dir, phams, data_cache=None,
                               file_type="fasta", mat_out=False,
                               tree_out=False, threads=1, verbose=False):
    if data_cache is None:
        data_cache = {}

    if verbose:
        print("Writing pham amino acid sequences to file...")

    fasta_path_map = alignment.create_pham_fastas(
                                        alchemist.engine, phams, working_dir,
                                        data_cache=data_cache,
                                        threads=threads, verbose=verbose)

    if verbose:
        print("Aligning pham amino acid sequences...")
    alignment.align_fastas(fasta_path_map, override=True,
                           mat_out=mat_out, tree_out=tree_out,
                           file_type=file_type,
                           threads=threads, verbose=verbose)
