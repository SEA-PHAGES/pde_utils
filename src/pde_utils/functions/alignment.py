from pathlib import Path
import shlex
from subprocess import Popen, PIPE

from Bio import AlignIO
from Bio.Emboss import Applications
import Levenshtein
from pdm_utils.functions import fileio as pdm_fileio
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import parallelize

from pde_utils.functions import multithread

# GLOBAL VARIABLES
# ----------------------------------------------------------------------
EMBOSS_TOOLS = ["needle", "water", "stretcher"]

EMBOSS_TOOL_SETTINGS = {"needle": {"gapopen": 10,
                                   "gapextend": 0.5},
                        "stretcher": {"gapopen": 12,
                                      "gapextend": 2},
                        "water": {"gapopen": 10,
                                  "gapextend": 0.5}
                        }

CLUSTALO_FORMATS = ["fasta", "clustal", "clustal", "phylip", "selex",
                    "stockholm", "vienna"]

# PAIRWISE ALIGNMENT TOOLS
# ---------------------------------------------------------------------


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


def pairwise_align(query_seq_path, target_seq_path, outfile_path,
                   tool="needle", gapopen=None, gapextend=None):
    """Aligns two sequences from fasta_files using EMBOSS tools.
    :param query_seq_path: Path to sequence A for pairwise alignment.
    :type query_seq_path: Path
    :param target_seq_path: Path to sequence B for pairwise alignment.
    :type target_seq_path: Path
    :param outfile_path: Path to write the alignment.
    :type outfile_path: Path
    :returns: Returns a Biopython object containing sequence alignment info.
    :rtype: MultipleSequenceAlignment
    """
    settings = EMBOSS_TOOL_SETTINGS.get(tool)
    if settings is None:
        raise ValueError

    if gapopen is None:
        gapopen = settings["gapopen"]
    if gapextend is None:
        gapextend = settings["gapextend"]

    if tool == "needle":
        cline_init = create_needle_cline
    elif tool == "stretcher":
        cline_init = create_stretcher_cline
    elif tool == "water":
        cline_init = create_water_cline

    emboss_cline = cline_init(query_seq_path, target_seq_path, outfile_path,
                              gapopen, gapextend)

    stdout, stderr = emboss_cline()

    alignment = AlignIO.read(outfile_path, "emboss")
    return alignment


def create_needle_cline(query_seq_path, target_seq_path, outfile_path,
                        gapopen, gapextend):
    """Helper function that retrieves the Needle cline command.
    :param query_seq_path: Path to sequence A for pairwise alignment.
    :type query_seq_path: Path
    :param target_seq_path: Path to sequence B for pairwise alignment.
    :type target_seq_path: Path
    :param outfile_path: Path to write the alignment.
    :type outfile_path: Path
    """
    needle_cline = Applications.NeedleCommandline(
                                        asequence=query_seq_path,
                                        bsequence=target_seq_path,
                                        outfile=outfile_path,
                                        gapopen=gapopen, gapextend=gapextend)

    return needle_cline


def create_stretcher_cline(query_seq_path, target_seq_path, outfile_path,
                           gapopen, gapextend):
    """Helper function that retrieves the Needle cline command.
    :param query_seq_path: Path to sequence A for pairwise alignment.
    :type query_seq_path: Path
    :param target_seq_path: Path to sequence B for pairwise alignment.
    :type target_seq_path: Path
    :param outfile_path: Path to write the alignment.
    :type outfile_path: Path
    """
    stretcher_cline = Applications.StretcherCommandline(
                                        asequence=query_seq_path,
                                        bsequence=target_seq_path,
                                        outfile=outfile_path,
                                        gapopen=gapopen, gapextend=gapextend)

    return stretcher_cline


def create_water_cline(query_seq_path, target_seq_path, outfile_path,
                       gapopen, gapextend):
    """Helper function that retrieves the Needle cline command.
    :param query_seq_path: Path to sequence A for pairwise alignment.
    :type query_seq_path: Path
    :param target_seq_path: Path to sequence B for pairwise alignment.
    :type target_seq_path: Path
    :param outfile_path: Path to write the alignment.
    :type outfile_path: Path
    """
    water_cline = Applications.WaterCommandline(
                                        asequence=query_seq_path,
                                        bsequence=target_seq_path,
                                        outfile=outfile_path,
                                        gapopen=gapopen, gapextend=gapextend)

    return water_cline

# MULTIPLE SEQUENCE TOOLS
# ---------------------------------------------------------------------


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
    command = (f"mBed -infile {fasta_file}")

    command = shlex.split(command)

    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()

    default_distmat = Path.cwd().joinpath("distMat.out")

    if isinstance(distmat_out, str):
        distmat_out = Path(distmat_out)

    default_distmat.replace(distmat_out)
    return default_distmat


def clustalo(fasta_file, aln_out_path, mat_out_path=None, tree_out_path=None,
             outfmt="clustal", infmt="fasta", threads=1, verbose=0):
    """
    Runs Clustal Omega to generate a multiple sequence alignment (MSA)
    and percent identity matrix (PIM) for the indicated file. Infile is
    expected to be in FASTA multiple sequence format. MSA will be in
    Clustal format.
    :param fasta_file: FASTA file containing sequences to be aligned
    :type fasta_file: str
    :param aln_out_path: the multiple sequence alignment (MSA) output file
    :type aln_out_path: str
    :type aln_out_path: Path
    :param mat_out_path: The percent identity matrix (PIM) output file
    :type mat_out_path: str
    :type mat_out_path: Path
    :param tree_out_path: The alignment guide tree output file
    :type tree_out_path: str
    :type tree_out_path: Path
    :param outfmt: The file format of the alignment to be exported.
    :type outfmt: str
    :param infmt:  The file format of the sequence file to be read in
    :type infmt: str
    :param threads: number of threads to use for alignment
    :type threads: int
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: Returns the sequence alignment path[0] and distance matrix path[1]
    :rtype: tuple
    """
    # Make sure verbose is in range(0,3)
    if verbose <= 0:
        verbose = 0
    elif verbose > 2:
        verbose = 2

    # Build Clustal Omega command that will produce a clustal-formatted
    # alignment output file and percent identity matrix
    command = f"clustalo -i {fasta_file} -o {aln_out_path} " \
              f" --outfmt={outfmt} --infmt={infmt} "\
              f"--force --output-order=tree-order " \
              f"--threads={threads}"

    if mat_out_path is not None:
        command = " ".join([command, (f"--distmat-out={mat_out_path} "
                                      "--full --percent-id")])
    if tree_out_path is not None:
        command = " ".join([command, (f"--guidetree-out={tree_out_path} ")])

    for _ in range(verbose):
        command += " -v"                # Add verbosity to command
    command = shlex.split(command)      # Convert command to arg list

    # Run the Clustal Omega command as a subprocess
    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    # Return output files so user can find them easily, whether they specified
    # the filenames or not
    return (aln_out_path, mat_out_path)


def hhmake(infile_path, hhm_path, name=None, add_cons=False, seq_lim=None,
           M=50, seq_id=90, verbose=0):
    """Runs HHsuite3 tool hhmake to make a HMM file from a MSA file.

    :param infile_path: FASTA formatted multiple sequence alignment input file.
    :type infile_path: str
    :type infile_path: Path
    :param hhm_path: Path for the HMM file to be exported to.
    :type hhm_path: str
    :type hhm_path: Path
    :param name: Optional naming for the HMM profile.
    :type name: str
    :param add_cons: Option to make the consensus sequence the master sequence.
    :type add_cons: bool
    :param seq_lim: Option to limit the number of query sequences displayed.
    :type seq_lim: int
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: Returns the HHM path
    :rtype: str
    :rtype: Path
    """
    command = f"hhmake -i {infile_path} -o {hhm_path} -v {verbose} -M {M}"

    if name is not None:
        command = " ".join([command, "-name", name])
    if add_cons:
        command = " ".join([command, "-add_cons"])
    if seq_lim is not None:
        command = " ".join([command, "-seq", str(seq_lim)])
    if seq_id is not None:
        command = " ".join([command, "-id", str(seq_id)])

    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return hhm_path


def hhalign(query_path, target_path, hhr_path, verbose=0):
    """Runs HHsuite3 tool hhalign to create a hhr result file from input files.

    :param query_path: HMM profile query input file.
    :type query_path: str
    :type query_path: Path
    :param target_path: HMM profile target input file.
    :type target_path: str
    :type target_path: Path
    :param hhr_path: Path for the result file to be exported to.
    :type hhr_path: str
    :type hhr_path: Path
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: Returns the hhr result path
    :rtype: str
    :rtype: Path
    """
    command = f"hhalign -i {query_path} -t {target_path} -o {hhr_path} " \
              f"-v {verbose}"

    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return hhr_path


def hmmbuild(infile_path, hmm_path, seq_type="dna", verbose=0):
    """Runs HMMER tool hmmbuild to make an HMM file from an MSA file.

    :param infile_path: FASTA formatted multiple sequence alignment input file.
    :type infile_path: str
    :type infile_path: Path
    :param hmm_path: Path for the HMM file to be exported to.
    :type hmm_path: str
    :type hmm_path: Path
    """
    command = "hmmbuild"

    if seq_type.lower() == "dna":
        command = " ".join([command, "--dna"])
    elif seq_type.lower() == "rna":
        command = " ".join([command, "--rna"])
    elif seq_type.lower() == "amino":
        command = " ".join([command, "--amino"])

    command = " ".join([command, str(hmm_path), str(infile_path)])

    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return hmm_path

# FILE I/O TOOLS
# ----------------------------------------------------------------------


def get_pham_genes(engine, phamid):
    """
    Queries the database for the geneids and translations found in the
    indicated pham. Returns a dictionary mapping the pham's geneids to
    their associated translations. All geneid:translation pairs will be
    represented (i.e. no redundant gene filtering is done...).
    :param engine: the Engine allowing access to the database
    :type engine: sqlalchemy Engine
    :param phamid: the pham whose genes are to be returned
    :type phamid: str
    :return: pham_genes
    :rtype: dict
    """
    pham_genes = dict()

    # Query will return the pham's GeneIDs and Translations grouped by genes
    # that share the same sequence
    query = f"SELECT GeneID, Translation FROM gene WHERE PhamID = {phamid} " \
            f"ORDER BY Translation, GeneID ASC"
    query_results = mysqldb_basic.query_dict_list(engine, query)

    for dictionary in query_results:
        geneid = dictionary["GeneID"]
        translation = dictionary["Translation"].decode("utf-8")
        pham_genes[geneid] = translation

    return pham_genes


def create_pham_fastas(engine, phams, aln_dir, data_cache=None, threads=1,
                       verbose=False):
    if data_cache is None:
        data_cache = {}

    work_items = []
    fasta_path_map = {}
    for pham in phams:
        fasta_path = aln_dir.joinpath(".".join([str(pham), "fasta"]))
        fasta_path_map[pham] = fasta_path

        gs_to_ts = data_cache.get(pham)
        if gs_to_ts is None:
            gs_to_ts = get_pham_genes(engine, pham)
            data_cache[pham] = gs_to_ts

        work_items.append((gs_to_ts, fasta_path))

    multithread.multithread(pdm_fileio.write_fasta, work_items, threads,
                            verbose=verbose)

    return fasta_path_map


def align_fastas(fasta_path_map, mat_out=False, tree_out=False,
                 file_type="fasta", mode="clustalo", override=False,
                 outdir=None, threads=1, verbose=False):
    verbose_num = 0

    work_items = []
    aln_path_map = {}
    for pham, fasta_path in fasta_path_map.items():
        if outdir is not None:
            working_dir = outdir
        else:
            working_dir = fasta_path.parent

        fasta_path_name = fasta_path.with_suffix("").name

        if override:
            aln_path = working_dir.joinpath(".".join([fasta_path_name,
                                                      "fasta"]))
        else:
            aln_path = working_dir.joinpath(".".join([fasta_path_name,
                                                      "aln"]))

        aln_path_map[pham] = aln_path

        mat_path = None
        if mat_out:
            mat_path = working_dir.joinpath(".".join([fasta_path_name,
                                                      "mat"]))

        tree_path = None
        if tree_out:
            tree_path = working_dir.joinpath(".".join([fasta_path_name,
                                                       "tree"]))

        work_items.append((fasta_path, aln_path, mat_path, tree_path,
                          file_type, "fasta", 1, verbose_num))

    if mode == "clustalo":
        aln_driver = clustalo
    else:
        raise NotImplementedError("Alignment program not supported.")

    parallelize.parallelize(work_items, threads, aln_driver, verbose=verbose)

    return aln_path_map


def create_hmms(aln_path_map, name=False, outdir=None, M=50, seq_id=90,
                add_cons=False, seq_lim=None, threads=1, verbose=False):
    verbose_num = 0

    work_items = []
    hmm_path_map = {}
    for pham, aln_path in aln_path_map.items():

        if outdir is not None:
            hmm_path_name = aln_path.with_suffix(".hmm").name
            hmm_path = outdir.joinpath(hmm_path_name)
        else:
            hmm_path = aln_path.with_suffix(".hmm")

        hmm_path_map[pham] = hmm_path

        hmm_name = None
        if name:
            hmm_name = str(pham)

        work_items.append((aln_path, hmm_path, hmm_name,
                           add_cons, seq_lim, M, seq_id, verbose_num))

    parallelize.parallelize(work_items, threads, hhmake, verbose=verbose)

    return hmm_path_map
