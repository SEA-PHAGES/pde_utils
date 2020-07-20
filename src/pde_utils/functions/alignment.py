import shlex
import textwrap
from subprocess import Popen, PIPE

from Bio import AlignIO
from Bio.Emboss import Applications
from pdm_utils.functions import fileio
from pdm_utils.functions import mysqldb_basic

#GLOBAL VARIABLES
#----------------------------------------------------------------------
EMBOSS_TOOLS = ["needle", "water", "stretcher"]

EMBOSS_TOOL_SETTINGS = {"needle"    : {"gapopen"    : 10,
                                       "gapextend"  : 0.5},
                        "stretcher" : {"gapopen"    : 12,
                                       "gapextend"  : 2},
                        "water"     : {"gapopen"    : 10,
                                       "gapextend"  : 0.5}
                       }

CLUSTALO_FORMATS = ["fasta", "clustal", "clustal", "phylip", "selex",
                    "stockholm", "vienna"]

#PAIRWISE ALIGNMENT TOOLS
#---------------------------------------------------------------------
def pairwise_align(query_seq_path, target_seq_path, outfile_path, 
                                                        tool="needle",
                                                        gapopen=None,
                                                        gapextend=None):
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

#MULTIPLE SEQUENCE TOOLS
#---------------------------------------------------------------------
def clustalo(fasta_file, aln_out_path, mat_out_path=None, outfmt="clustal",
                                                     infmt="fasta",
                                                     threads=1, verbose=0):
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
    :param mat_out_path: the percent identity matrix (PIM) output file
    :type mat_out_path: str
    :type mat_out_path: Path
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

    if not mat_out_path is None:
        command = " ".join([command, (f"--distmat-out={mat_out_path} "
                                       "--full --percent-id")])

    for x in range(verbose):
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

def hhmake(infile_path, hmm_path, name=None, add_cons=False, seq_lim=None,
                                                           M=50, verbose=0):
    """Runs HHsuite3 tool hhmake to make a HMM file from a MSA file.

    :param infile_path: FASTA formatted multiple sequence alignment input file.
    :type infile_path: str
    :type infile_path: Path
    :param hmm_path: Path for the HMM file to be exported to.
    :type hmm_path: str
    :type hmm_path: Path
    :param name: Optional naming for the HMM profile.
    :type name: str
    :param add_cons: Option to make the consensus sequence the master sequence.
    :type add_cons: bool
    :param seq_lim: Option to limit the number of query sequences displayed.
    :type seq_lim: int
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: Returns the HMM path
    :rtype: str
    :rtype: Path
    """
    command = f"hhmake -i {infile_path} -o {hmm_path} -v {verbose} -M {M}"

    if not name is None:
        command = " ".join([command, "-name {name}"])
    if add_cons:
        command = " ".join([command, "-add_cons"])
    if seq_lim:
        command = " ".join([command, "-seq {seq_lim}"])
    
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return hmm_path

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

#FILE I/O TOOLS
#----------------------------------------------------------------------
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

def create_pham_fasta(engine, phams, aln_dir, data_cache=None,
                                                  verbose=False):
    """
    """
    if data_cache is None:
        data_cache = {}

    fasta_path_map = {}
    for pham in phams:
        fasta_path = aln_dir.joinpath(".".join([str(pham), "fasta"]))
        fasta_path_map[pham] = fasta_path
        
        gs_to_ts = data_cache.get(pham) 
        if gs_to_ts is None:
            gs_to_ts = get_pham_genes(engine, pham)
            data_cache[pham] = gs_to_ts

        fileio.write_fasta(gs_to_ts, fasta_path)

    return fasta_path_map

def align_pham_fastas(fasta_path_map, mat_out=False, threads=1, verbose=False):
    if verbose:
        verbose = 2
    else:
        verbose = 0

    mat_path_map = {}
    for pham, fasta_path in fasta_path_map.items():
        mat_path = None
        if mat_out:
            mat_path = fasta_path.with_name(".".join([str(pham), "mat"]))
            mat_path_map[pham] = mat_path

        clustalo(fasta_path, fasta_path, outfmt="fasta", mat_out_path=mat_path,
                                                   threads=threads,
                                                   verbose=verbose)

    return mat_path_map

def create_pham_hmms(fasta_path_map, name=False, M=50,
                     add_cons=False, seq_lim=None, verbose=False):
    if verbose:
        verbose = 2
    else:
        verbose = 0

    hmm_path_map = {}
    for pham, fasta_path in fasta_path_map.items():

        hmm_path = fasta_path.with_name(".".join([str(pham), "hmm"]))
        hmm_path_map[pham] = hmm_path

        hmm_name = None
        if name:
            hmm_name = str(pham)

        hhmake(fasta_path, hmm_path, name=hmm_name, add_cons=add_cons,
                M=M, seq_lim=seq_lim, verbose=verbose)
   
    return hmm_path_map


