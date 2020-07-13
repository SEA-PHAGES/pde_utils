import shlex
import textwrap
from subprocess import Popen, PIPE

#Imports to be removed
from pathlib import Path
from pdm_utils.classes.alchemyhandler import AlchemyHandler
from pdm_utils.classes import clustal

from Bio import AlignIO
from Bio.Emboss import Applications

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
    :type phamid: intbiopython pairwise distance
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
    :param aln_out: the multiple sequence alignment (MSA) output file
    :type aln_out: str
    :param mat_out: the percent identity matrix (PIM) output file
    :type mat_out: str
    :param outfmt: The file format of the alignment to be exported.
    :type outfmt: str
    :param threads: number of threads to use for alignment
    :type threads: int
    :param verbose: verbosity level (0-2)
    :type verbose: int
    :return: [aln_out, mat_out]
    :rtype: list
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

    if mat_out_path:
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
                                                           verbose=0):
    command = f"hhmake -i {infile_path} -o {hmm_path} -v {verbose}"

    if not name is None:
        if not isinstance(name, str):
            raise TypeError
        command = " ".join([command, name])
    if add_cons:
        command = " ".join([command, "-add_cons"])
    if seq_lim:
        if not isinstance(seq_lim, int):
            raise TypeError
        command = " ".join([command, "-seq {seq_lim}"])

    
    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return hmm_path

def hhalign(query_path, target_path, hhr_path, add_cons=True, verbose=0):
    command = f"hhalign -i {query_path} -t {target_path} -o {hhr_path} " \
              f"-v {verbose}"

    if add_cons:
        command = " ".join([command, "-add_cons"])

    with Popen(args=shlex.split(command), stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

    return hhr_path

#################
# Example Usage #
#################
def not_main():
    # Set up MySQL handler
    alchemist = AlchemyHandler(database="Actinobacteriophage")
    alchemist.connect(pipeline=True)
    engine = alchemist.engine

    # Choose a pham, set up infile name, and construct FASTA
    pham = 11176
    infile = "/tmp/pham_11176.fasta"
    genes = get_pham_genes(engine, pham)
    write_genes_to_fasta(genes, infile)

    # Run clustal and expand the return value to output filenames
    out_aln, out_mat = clustalo(infile)

    # Parse the identity matrix, get the centroid name and average
    # distance from centroid to each gene in the pham
    mat = PercentIdentityMatrix(out_mat)
    mat.parse_matrix()
    centroid = mat.get_centroid()
    neighbors = mat.get_nearest_neighbors(centroid, 50)
    neighbors = "\n".join(neighbors)
    average_dist = round(100 - mat.get_average_identity(centroid), 3)

    # Parse the multiple sequence alignment and get centroid sequence w/o gaps
    aln = MultipleSequenceAlignment(out_aln)
    aln.parse_alignment()
    centroid_seq = aln.get_sequence(centroid, gaps=False)

    print("=======================================")
    print(f"Centroid gene: {centroid}")
    print(f"Centroid sequence: {centroid_seq[:16]}...")
    print(f"Average distance from centroid: {average_dist}%")
    print(f"{neighbors}")
    print("=======================================")

