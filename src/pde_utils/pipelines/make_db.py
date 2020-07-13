import argparse
import os
import shutil
import shlex
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from subprocess import DEVNULL, PIPE, Popen

from Bio import AlignIO
from pdm_utils.functions import basic
from pdm_utils.functions import fileio
from pdm_utils.functions import pipelines_basic

from pde_utils.classes import clustal
from pde_utils.functions import alignment

#GLOBAL VARIABLES
#-----------------------------------------------------------------------------

DATABASE_TYPES = ["hhsuite", "mmseqs"]

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_db"
DEFAULT_FOLDER_PATH = Path.cwd()

PHAM_FASTA_COLUMNS = ["gene.GeneID", "gene.Translation"]

#MAIN FUNCTIONS
#-----------------------------------------------------------------------------
def main(unparsed_args_list):
    args = parse_make_db(unparsed_args_list)

    alchemist = pipelines_basic.build_alchemist(args.database) 
    values = pipelines_basic.parse_value_input(args.input)

    execute_make_db(alchemist, args.folder_path, args.folder_name,
                    args.db_type, values=values, verbose=args.verbose,
                    filters=args.filters, groups=args.groups, 
                    db_name=args.db_name, threads=args.number_threads)
    
def parse_make_db(unparsed_args_list):
    DATABASE_HELP = """
        Name of the MySQL database to export from.
        """
    DATABASE_TYPE_HELP = """
        Type of the database to make
        """

    VERBOSE_HELP = """
        Export option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Export option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired export directory.
        """
    FOLDER_NAME_HELP = """
        Export option to change the name
        of the directory where the exported files are stored.
            Follow selection argument with the desired name.
        """

    IMPORT_FILE_HELP = """
        Selection input option that imports values from a csv file.
            Follow selection argument with path to the
            csv file containing the names of each genome in the first column.
        """
    IMPORT_NAMES_HELP = """
        Selection input option that imports values from cmd line input.
            Follow selection argument with space separated
            names of genes in the database.
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

    DATABASE_NAME_HELP = """
        Pipeline option that allows for renaming of the created database.
            Follow selection argument with the desired name of the database.
    """
    NUMBER_THREADS_HELP = """
        Pipeline option that allows for multithreading of workflows.
            Follow selection argument with number of threads to be used
        """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("db_type", type=str, choices=DATABASE_TYPES,
                        help=DATABASE_TYPE_HELP)

    parser.add_argument("-m", "--folder_name", type=str, 
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path",
                        type=pipelines_basic.convert_dir_path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)

    parser.add_argument("-dbm", "--database_name", type=str,
                        help=DATABASE_NAME_HELP)
    parser.add_argument("-np", "--number_threads", type=int, nargs="?",
                        help=NUMBER_THREADS_HELP)

    parser.add_argument("-if", "--import_file",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP, dest="input")
    parser.add_argument("-in", "--import_names", nargs="*",
                        help=IMPORT_NAMES_HELP, dest="input")
    parser.add_argument("-f", "--where", nargs="?",
                        help=WHERE_HELP,
                        dest="filters")
    parser.add_argument("-g", "--group_by", nargs="+",
                        help=GROUP_BY_HELP,
                        dest="groups") 

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH,
                        verbose=False, input=[],
                        filters="", groups=[],
                        db_name=None, number_threads=1) 

    parsed_args = parser.parse_args(unparsed_args_list[1:])
    return parsed_args
    
def execute_make_db(alchemist, folder_path, folder_name, db_type,
                    values=None, verbose=False, filters="", groups=[],
                    db_name=None, threads=1):
    if db_name is None:
        db_name = alchemist.database
    
    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                             values=values,
                                             verbose=verbose)

    db_path = folder_path.joinpath(folder_name)
    db_path = basic.make_new_dir(folder_path, db_path, attempt=50)

    data_cache = {}
    conditionals_map = {}
    pipelines_basic.build_groups_map(db_filter, db_path,
                                     conditionals_map,
                                     groups=groups, verbose=verbose)

    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        db_filter.reset()
        db_filter.values = values

        conditionals = conditionals_map[mapped_path]
        db_filter.values = db_filter.build_values(where=conditionals)

        if db_filter.hits() == 0:
            print(f"No database entries received for '{mapped_path}'.")
            continue
  
        if db_type == "hhsuite":
            execute_make_hhsuite_database(alchemist, db_filter.values,
                                          db_path, db_name,
                                          data_cache=data_cache,
                                          threads=threads)

def execute_make_hhsuite_database(alchemist, values, db_dir, db_name, 
                                                              data_cache=None,
                                                              verbose=False,
                                                              threads=1):
    aln_dir = db_dir.joinpath("pham_alignments")
    aln_dir.mkdir()

    fasta_path_map = create_pham_fastas(alchemist, values, aln_dir, 
                                                  data_cache=data_cache,
                                                  verbose=verbose)
    align_pham_fastas(fasta_path_map, threads=threads, 
                                                  verbose=verbose)

    if verbose:
        stdout = PIPE
    else:
        stdout = DEVNULL

    msa_fftuple = create_msa_ffindex(aln_dir, db_dir, db_name, stdout=stdout)
    a3m_fftuple = create_a3m_ffindex(db_dir, db_name, msa_fftuple,
                                                  threads=threads, 
                                                  stdout=stdout)
    hmm_fftuple = create_hmm_ffindex(db_dir, db_name, a3m_fftuple, 
                                                  threads=threads,
                                                  stdout=stdout)
    create_cs219_ffindex(db_dir, db_name, stdout=stdout)

    sorting_file = create_sorting_file(db_dir, db_name, stdout=stdout)
    sort_a3m_ffindex(db_dir, db_name, sorting_file, a3m_fftuple, 
                                                  stdout=stdout)
    sort_hmm_ffindex(db_dir, db_name, sorting_file, hmm_fftuple, 
                                                  stdout=stdout)

    
    msa_fftuple[0].unlink()
    msa_fftuple[1].unlink()
    sorting_file.unlink()

    if not verify_hhsuite_database(db_dir, db_name):
        print(f"Inconsistencies detected in HHsuite database "
              f"at '{db_dir}'.\n  Scrapping database.")
        shutil.rmtree(db_dir)

def create_pham_fastas(alchemist, phams, aln_dir, data_cache=None,
                                                  verbose=False):
    if data_cache is None:
        data_cache = {}

    fasta_path_map = {}
    for pham in phams:
        fasta_path = aln_dir.joinpath(".".join([str(pham), "fasta"]))
        fasta_path_map[pham] = fasta_path
        
        gs_to_ts = data_cache.get(pham) 
        if gs_to_ts is None:
            gs_to_ts = alignment.get_pham_genes(alchemist.engine, pham)
            data_cache[pham] = gs_to_ts

        fileio.write_fasta(gs_to_ts, fasta_path)

    return fasta_path_map

def align_pham_fastas(fasta_path_map, threads=1, verbose=False):
    if verbose:
        verbose = 2
    else:
        verbose = 0

    for pham, fasta_path in fasta_path_map.items():
        alignment.clustalo(fasta_path, fasta_path, outfmt="fasta", 
                                                   threads=threads,
                                                   verbose=verbose)

#HHSUITE DB HELPER FUNCTIONS
#-----------------------------------------------------------------------------
def create_msa_ffindex(fasta_dir, db_dir, db_name, stdout=None):
    if stdout is None:
        stdout = DEVNULL
    
    msa_ffdata = db_dir.joinpath("".join([db_name, "_msa.ffdata"]))
    msa_ffindex = db_dir.joinpath("".join([db_name, "_msa.ffindex"]))
    
    command = f"ffindex_build -s {msa_ffdata} {msa_ffindex} {fasta_dir}"

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (msa_ffdata, msa_ffindex)

def create_a3m_ffindex(db_dir, db_name, msa_fftuple, threads=1, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    a3m_ffdata = db_dir.joinpath("".join([db_name, "_a3m.ffdata"]))
    a3m_ffindex = db_dir.joinpath("".join([db_name, "_a3m.ffindex"]))
  
    os.environ["OMP_NUM_THREADS"] = "1"
    command = (f"mpirun -np {threads} "
               f"ffindex_apply_mpi {msa_fftuple[0]} {msa_fftuple[1]} "
               f"-i {a3m_ffindex} -d {a3m_ffdata} "
                "-- hhconsensus -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0")

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (a3m_ffdata, a3m_ffindex)

def create_hmm_ffindex(db_dir, db_name, a3m_fftuple, threads=1, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    hmm_ffdata = db_dir.joinpath("".join([db_name, "_hmm.ffdata"]))
    hmm_ffindex = db_dir.joinpath("".join([db_name, "_hmm.ffindex"]))
    
    command = (f"mpirun -np {threads} ffindex_apply_mpi "
               f"{a3m_fftuple[0]} {a3m_fftuple[1]} " 
               f"-i {hmm_ffindex} -d {hmm_ffdata} "
                "-- hhmake -i stdin -o stdout -v 0")

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()
    
    return (hmm_ffdata, hmm_ffindex)

def create_cs219_ffindex(db_dir, db_name, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    a3m_ffpath = db_dir.joinpath("".join([db_name, "_a3m"]))
    cs219_ffpath = db_dir.joinpath("".join([db_name, "_cs219"]))
    
    command = ("cstranslate -f -x 0.3 -c 4 -I a3m "
              f"-i {a3m_ffpath} -o {cs219_ffpath}")

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

def create_sorting_file(db_dir, db_name, sorting_file=None, stdout=None):
    if sorting_file is None:
        sorting_file = db_dir.joinpath("sorting.dat")

    cs219_ffindex = db_dir.joinpath("".join([db_name, "_cs219.ffindex"]))
    
    sort_command = f"sort -k3 -n -r {cs219_ffindex} "
    cut_command = "cut -f1"

    split_sort_command = shlex.split(sort_command)
    split_cut_command = shlex.split(cut_command)
    
    sort_process = Popen(args=split_sort_command, stdout=PIPE)
    cut_process = Popen(args=split_cut_command, stdout=PIPE,
                                                  stdin=sort_process.stdout)
    sort_process.stdout.close()

    file_handle = sorting_file.open(mode="w")
    file_handle.write((cut_process.communicate()[0]).decode("utf-8"))
    file_handle.close()
    
    return sorting_file

def sort_a3m_ffindex(db_dir, db_name, sorting_file, a3m_tuple, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    ord_a3m_ffdata = db_dir.joinpath("".join([db_name, "ordered_a3m.ffdata"]))
    ord_a3m_ffindex = db_dir.joinpath("".join([db_name, "ordered_a3m.ffindex"]))

    command = (f"ffindex_order {sorting_file} {a3m_tuple[0]} {a3m_tuple[1]} "
               f"{ord_a3m_ffdata} {ord_a3m_ffindex}")
    
    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    ord_a3m_ffdata.replace(a3m_tuple[0])
    ord_a3m_ffindex.replace(a3m_tuple[1])

def sort_hmm_ffindex(db_dir, db_name, sorting_file, hmm_tuple, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    ord_hmm_ffdata = db_dir.joinpath("".join([db_name, "ordered_hmm.ffdata"]))
    ord_hmm_ffindex = db_dir.joinpath("".join([db_name, "ordered_hmm.ffindex"]))

    command = (f"ffindex_order {sorting_file} {hmm_tuple[0]} {hmm_tuple[1]} "
               f"{ord_hmm_ffdata} {ord_hmm_ffindex}")
    
    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    ord_hmm_ffdata.replace(hmm_tuple[0])
    ord_hmm_ffindex.replace(hmm_tuple[1])

def verify_hhsuite_database(db_dir, db_name, threads=1, stdout=None):
    if stdout == None:
        stdout = DEVNULL

    return True

if __name__ == "__main__":
    main(sys.argv)
