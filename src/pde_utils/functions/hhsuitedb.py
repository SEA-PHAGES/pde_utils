import os
import shlex
import shutil
from subprocess import DEVNULL, PIPE, Popen

from pde_utils.functions import fileio as pde_fileio


# CREATE DATABASE
# -----------------------------------------------------------------------------
def create_hhsuitedb(aln_dir, db_dir, db_name, threads=1, use_mpi=False,
                     verbose=False, pham_versions=None):
    if verbose:
        stdout = PIPE
    else:
        stdout = DEVNULL

    if verbose:
        print("Creating multiple sequence ffindex file...")
    msa_fftuple = build_msa_ffindex(aln_dir, db_dir, db_name, stdout=stdout)

    if verbose:
        print("Condensing multiple sequence ffindex file into a3m file...")
    a3m_fftuple = convert_a3m_ffindex(db_dir, db_name, msa_fftuple,
                                      threads=threads, stdout=stdout,
                                      use_mpi=use_mpi)

    if verbose:
        print("Creating Hidden Markov Model ffindex file...")
    hmm_fftuple = convert_hmm_ffindex(db_dir, db_name, a3m_fftuple,
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

    db_path = db_dir.joinpath(db_name)
    if pham_versions is not None:
        pde_fileio.create_phamdb_version_file(db_path, pham_versions)

    return db_path


def build_msa_ffindex(fasta_dir, db_dir, db_name, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    msa_ffdata = db_dir.joinpath("".join([db_name, "_msa.ffdata"]))
    msa_ffindex = db_dir.joinpath("".join([db_name, "_msa.ffindex"]))

    command = f"ffindex_build -s {msa_ffdata} {msa_ffindex} {fasta_dir}"

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (msa_ffdata, msa_ffindex)


def build_hmm_ffindex(hmm_dir, db_dir, db_name, stdout=None):
    if stdout is None:
        stdout = DEVNULL

    hmm_ffdata = db_dir.joinpath("".join([db_name, "_hmm.ffdata"]))
    hmm_ffindex = db_dir.joinpath("".join([db_name, "_hmm.ffindex"]))

    command = f"ffindex_build -s {hmm_ffdata} {hmm_ffindex} {hmm_dir}"

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return (hmm_ffdata, hmm_ffindex)


def convert_a3m_ffindex(db_dir, db_name, msa_fftuple, threads=1, stdout=None,
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


def convert_hmm_ffindex(db_dir, db_name, a3m_fftuple, threads=1, stdout=None,
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


def create_cs219_ffindex(db_dir, db_name, threads=1, use_mpi=False,
                         stdout=None):
    if stdout is None:
        stdout = DEVNULL

    a3m_ffpath = db_dir.joinpath("".join([db_name, "_a3m"]))
    cs219_ffpath = db_dir.joinpath("".join([db_name, "_cs219"]))

    command = ("cstranslate -f -x 0.3 -c 4 -I a3m "
               f"-i {a3m_ffpath} -o {cs219_ffpath}")

    if use_mpi:
        command = "".join([f"mpirun -np {threads}", command])

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
