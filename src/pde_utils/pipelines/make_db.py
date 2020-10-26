import argparse
import os
import shutil
import shlex
import sys
import time
from psutil import cpu_count
from pathlib import Path
from subprocess import DEVNULL, PIPE, Popen
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import phameration

from pde_utils.classes import clustal
from pde_utils.functions import alignment
from pde_utils.functions import fileio as pde_fileio

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

DATABASE_TYPES = ["hhsuite", "blastp"]

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_db"
DEFAULT_FOLDER_PATH = Path.cwd()

PHAM_FASTA_COLUMNS = ["gene.GeneID", "gene.Translation"]

# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(unparsed_args_list):
    args = parse_make_db(unparsed_args_list)

    alchemist = pipelines_basic.build_alchemist(args.database)
    values = pipelines_basic.parse_value_input(args.input)

    execute_make_db(alchemist, args.folder_path, args.folder_name,
                    args.db_type, values=values, verbose=args.verbose,
                    filters=args.filters, groups=args.groups,
                    db_name=args.db_name, threads=args.threads,
                    use_mpi=args.use_mpi)


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

    MPI_HELP = """
        Pipeline option that allows for use of ffindex_apply_mpi
        parallelization.
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
    parser.add_argument("-np", "--threads", type=int, nargs="?",
                        help=NUMBER_THREADS_HELP)

    parser.add_argument("-n", "--database_name", type=str,
                        help=DATABASE_NAME_HELP)
    parser.add_argument("-mpi", "--use_mpi", action="store_true",
                        help=MPI_HELP)

    parser.add_argument("-if", "--import_file",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP, dest="input")
    parser.add_argument("-in", "--import_names", nargs="*",
                        help=IMPORT_NAMES_HELP, dest="input")
    parser.add_argument("-w", "--where", nargs="?",
                        help=WHERE_HELP,
                        dest="filters")
    parser.add_argument("-g", "--group_by", nargs="+",
                        help=GROUP_BY_HELP,
                        dest="groups")

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH,
                        verbose=False, input=[],
                        filters="", groups=[],
                        db_name=None, threads=1)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_make_db(alchemist, folder_path, folder_name, db_type,
                    values=None, verbose=False, filters="", groups=[],
                    db_name=None, threads=1, use_mpi=False):
    if db_name is None:
        db_name = alchemist.database

    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                             values=values,
                                             verbose=verbose)

    working_path = pipelines_basic.create_working_path(folder_path,
                                                       folder_name)

    conditionals_map = pipelines_basic.build_groups_map(
                                            db_filter, working_path,
                                            groups=groups, verbose=verbose)

    data_cache = {}
    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        db_filter.reset()
        db_filter.values = values

        conditionals = conditionals_map[mapped_path]
        db_filter.values = db_filter.build_values(where=conditionals)

        if db_filter.hits() == 0:
            print(f"No database entries received for '{mapped_path}'.")
            continue

        pipelines_basic.create_working_dir(mapped_path)

        if db_type == "hhsuite":
            execute_make_hhsuite_database(alchemist, db_filter.values,
                                          mapped_path, db_name,
                                          data_cache=data_cache,
                                          threads=threads, verbose=verbose,
                                          use_mpi=use_mpi)


def execute_make_hhsuite_database(alchemist, values, db_dir, db_name,
                                  data_cache=None, verbose=False, threads=1,
                                  use_mpi=False):
    aln_dir = db_dir.joinpath("pham_alignments")
    aln_dir.mkdir()

    create_pham_alignments(
                    alchemist, values, aln_dir, data_cache=data_cache,
                    threads=threads, verbose=verbose)

    if verbose:
        stdout = PIPE
    else:
        stdout = DEVNULL

    if use_mpi:
        physical_cores = cpu_count(logical=False)

        if physical_cores < threads:
            if verbose:
                print("Designated process count greater than machine's "
                      "number of physical cores...\n"
                      f"STEPPING DOWN TO {physical_cores} PROCESSES")
            threads = physical_cores

    if verbose:
        print("Creating multiple sequence ffindex file...")
    msa_fftuple = create_msa_ffindex(aln_dir, db_dir, db_name, stdout=stdout)

    if verbose:
        print("Condensing multiple sequence ffindex file into a3m file...")
    a3m_fftuple = create_a3m_ffindex(db_dir, db_name, msa_fftuple,
                                     threads=threads, stdout=stdout,
                                     use_mpi=use_mpi)

    if verbose:
        print("Creating Hidden Markov Model ffindex file...")
    hmm_fftuple = create_hmm_ffindex(db_dir, db_name, a3m_fftuple,
                                     threads=threads, stdout=stdout,
                                     use_mpi=use_mpi)

    if verbose:
        print("Creating cs219 ffindex file...")
    create_cs219_ffindex(db_dir, db_name, stdout=stdout)

    if verbose:
        print("Sorting ffindex files for fast database indexing...")
    sorting_file = create_sorting_file(db_dir, db_name, stdout=stdout)
    sort_a3m_ffindex(db_dir, db_name, sorting_file, a3m_fftuple, stdout=stdout)
    sort_hmm_ffindex(db_dir, db_name, sorting_file, hmm_fftuple, stdout=stdout)

    msa_fftuple[0].unlink()
    msa_fftuple[1].unlink()
    sorting_file.unlink()

    if not verify_hhsuite_database(db_dir, db_name):
        print(f"Inconsistencies detected in HHsuite database "
              f"at '{db_dir}'.\n  Scrapping database.")
        shutil.rmtree(db_dir)


def execute_make_mmseqs_database(alchemist, values, db_dir, db_name,
                                 data_cache=None, verbose=False, threads=1):
    pass


def execute_make_blast_protein_database(
                                    alchemist, values, db_dir, db_name,
                                    data_cache=None, verbose=False, threads=1):
    pass

# HELPER FUNCTIONS
# -----------------------------------------------------------------------------


def create_pham_alignments(alchemist, values, aln_dir, data_cache=None,
                           threads=1, verbose=False):
    if data_cache is None:
        data_cache = {}

    if verbose:
        print("Writing pham amino acid sequences to file...")
    fasta_path_map = alignment.create_pham_fastas(
                                        alchemist.engine, values, aln_dir,
                                        data_cache=data_cache, threads=threads,
                                        verbose=verbose)

    if verbose:
        print("Aligning pham amino acid sequences...")
    aln_path_map = alignment.align_fastas(
                                        fasta_path_map, override=True,
                                        threads=threads, verbose=verbose)

    name_path_map = {}
    for pham, path in aln_path_map.items():
        name_path_map[f"pham{pham}"] = path
    if verbose:
        print("Adding names to pham alignments...")
    pde_fileio.name_comment_files(name_path_map, threads=threads,
                                  verbose=verbose)

    return aln_path_map


# HHSUITE DB HELPER FUNCTIONS
# -----------------------------------------------------------------------------
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


def create_a3m_ffindex(db_dir, db_name, msa_fftuple, threads=1, stdout=None,
                       use_mpi=False):
    if stdout is None:
        stdout = DEVNULL

    a3m_ffdata = db_dir.joinpath("".join([db_name, "_a3m.ffdata"]))
    a3m_ffindex = db_dir.joinpath("".join([db_name, "_a3m.ffindex"]))

    command = (f"{msa_fftuple[0]} {msa_fftuple[1]} "
               f"-d {a3m_ffdata} -i {a3m_ffindex} "
               "-- hhconsensus -M 50 -maxres 65535 "
               "-i stdin -oa3m stdout -v 0")

    if use_mpi:
        os.environ["OMP_NUM_THREADS"] = "1"
        command = "".join([(f"mpirun -np {threads} ffindex_apply_mpi "),
                           command])
    else:
        command = "".join(["ffindex_apply ", command])

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (a3m_ffdata, a3m_ffindex)


def create_hmm_ffindex(db_dir, db_name, a3m_fftuple, threads=1, stdout=None,
                       use_mpi=False):
    if stdout is None:
        stdout = DEVNULL

    hmm_ffdata = db_dir.joinpath("".join([db_name, "_hmm.ffdata"]))
    hmm_ffindex = db_dir.joinpath("".join([db_name, "_hmm.ffindex"]))

    command = (f"{a3m_fftuple[0]} {a3m_fftuple[1]} "
               f"-i {hmm_ffindex} -d {hmm_ffdata} "
               "-- hhmake -i stdin -o stdout -v 0")

    if use_mpi:
        os.environ["OMP_NUM_THREADS"] = "1"
        command = "".join([f"mpirun -np {threads} ffindex_apply_mpi ",
                           command])
    else:
        command = "".join(["ffindex_apply ", command])

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
    ord_a3m_ffindex = db_dir.joinpath("".join(
                                            [db_name, "ordered_a3m.ffindex"]))

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
    ord_hmm_ffindex = db_dir.joinpath("".join(
                                            [db_name, "ordered_hmm.ffindex"]))

    command = (f"ffindex_order {sorting_file} {hmm_tuple[0]} {hmm_tuple[1]} "
               f"{ord_hmm_ffdata} {ord_hmm_ffindex}")

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    ord_hmm_ffdata.replace(hmm_tuple[0])
    ord_hmm_ffindex.replace(hmm_tuple[1])


def verify_hhsuite_database(db_dir, db_name, threads=1, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    return True


# MMSEQS DB HELPER FUNCTIONS
# -----------------------------------------------------------------------------

CLUSTER_ARGS = {"tmp_dir": None,
                "identity": None,
                "coverage": None,
                "e_value": None,
                "sens": None,
                "steps": None,
                "threads": None,
                "aln_mode": None,
                "cov_mode": None,
                "clu_mode": None}


def align_centroids(cfasta_path, tmp_dir, threads=1):
    caln_path = tmp_dir.joinpath("centroid.aln")
    cmat_path = tmp_dir.joinpath("centroid.mat")

    alignment.clustalo(cfasta_path, caln_path, mat_out_path=cmat_path,
                       threads=threads)

    caln = clustal.MultipleSequenceAlignment(caln_path)
    caln.parse_alignment()
    cmat = clustal.PercentIdentityMatrix(cmat_path)
    cmat.parse_matrix()

    return caln, cmat


def clust_centroids(cfasta_path, aln_dir):
    db_dir = aln_dir.joinpath("DB_tmp")
    db_dir.mkdir()

    seq_db = db_dir.joinpath("sequenceDB")
    clu_db = db_dir.joinpath("clusterDB")
    psf_db = db_dir.joinpath("seqfileDB")
    pre_nbhds_file = aln_dir.joinpath("pre_nbhds.fasta")

    phameration.mmseqs_createdb(cfasta_path, seq_db)
    phameration.mmseqs_cluster(seq_db, clu_db, CLUSTER_ARGS)
    phameration.mmseqs_createseqfiledb(seq_db, clu_db, psf_db)
    phameration.mmseqs_result2flat(seq_db, seq_db, psf_db,
                                   pre_nbhds_file)
    pre_neighborhoods = phameration.parse_mmseqs_output(
                                        pre_nbhds_file)

    return pre_neighborhoods


def build_centroid_mmseqs_dbs(caln_path, centroid_fasta_map, aln_dir):
    database_dir = aln_dir.joinpath("mmdb")
    database_dir.mkdir()

    total_centroid_db = database_dir.joinpath("centroidDB")
    phameration.mmseqs_createdb(caln_path, total_centroid_db)

    centroid_db_map = {}
    for pham, centroid_fasta in centroid_fasta_map.items():
        centroid_db = database_dir.joinpath(f"{pham}centroidDB")

        phameration.mmseqs_createdb(centroid_fasta, centroid_db)
        centroid_db_map[pham] = centroid_db

    return total_centroid_db, centroid_db_map


if __name__ == "__main__":
    main(sys.argv)
