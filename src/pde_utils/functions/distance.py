from pathlib import Path
import shlex
from subprocess import Popen, PIPE

import Levenshtein

# GLOBAL VARIABLES
MASH_TOOL_SETTINGS = {"kmer": 15, "sketch": 25000}


def mash_sketch(seq_path, out_path, kmer=None, sketch=None, verbose=0):
    command = (f"mash sketch {seq_path} -o {out_path}")

    if kmer is not None:
        command = " ".join([command, "-k", str(kmer)])
    if sketch is not None:
        command = " ".join([command, "-s", str(sketch)])

    command = shlex.split(command)
    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    out_path = out_path.with_name("".join([out_path.name, ".msh"]))
    return out_path


def mash_dist(query_seq_path, target_seq_path, outfile=None, kmer=None,
              sketch=None, verbose=0):
    command = (f"mash dist {query_seq_path} {target_seq_path}")

    if kmer is not None:
        command = " ".join([command, "-k", str(kmer)])
    if sketch is not None:
        command = " ".join([command, "-s", str(sketch)])
    if outfile is not None:
        command = " ".join([command, ">", str(outfile)])

    command = shlex.split(command)
    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()

        out = out.decode("utf-8")
        if verbose > 0 and out:
            print(out)
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return out


def calculate_levenshtein(source_str, target_str, identity=False):
    LD = Levenshtein.distance(source_str, target_str)

    source_len = len(source_str)
    target_len = len(target_str)

    if source_len > target_len:
        percent = (LD / source_len) * 100
    else:
        percent = (LD / target_len) * 100

    if identity:
        percent = 100 - percent

    return percent


def calculate_gcs(query_phams, target_phams):
    shared_phams = set(query_phams).intersection(set(target_phams))

    query_spp = len(shared_phams) / len(query_phams)
    target_spp = len(shared_phams) / len(target_phams)

    return (query_spp + target_spp) / 2


def calculate_norm_gcs(query_phams, target_phams):
    pass


def mBed(fasta_file, distmat_out):
    """
    Runs mBed to generate a k-tuple distance matrix and UPGMA algorithm
    guidetree.  Infile is expected to be in Fasta multiple sequence format.
    :param fasta_file: FASTA file containing sequences to be aligned:
    :type fasta_file: Path
    :type fasta_file: str
    :param distmat_out: Desired destination of the mBed distmat
    :type distmat_out: Path
    :type distmat_out: str
    """
    command = (f"""mBed -infile "{fasta_file}" """)

    command = shlex.split(command)

    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()

    default_distmat = Path.cwd().joinpath("distMat.out")

    if isinstance(distmat_out, str):
        distmat_out = Path(distmat_out)

    default_distmat.replace(distmat_out)
    return default_distmat
