import argparse
import sys
import time

from pdm_utils.functions import fileio as pdm_fileio
from pdm_utils.functions import querying
from pdm_utils.functions import mysqldb_basic
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import phameration
from psutil import cpu_count

from pde_utils.classes import clustal
from pde_utils.functions import alignment
from pde_utils.functions import blastdb
from pde_utils.functions import fileio as pde_fileio
from pde_utils.functions import hhsuitedb


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------

DATABASE_TYPES = ["hhsuite", "blast"]

DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_db"

PHAM_FASTA_COLUMNS = ["gene.GeneID", "gene.Translation"]

# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(unparsed_args_list):
    args = parse_make_db(unparsed_args_list)

    alchemist = pipelines_basic.build_alchemist(args.database)
    values = pipelines_basic.parse_value_input(args.input)

    execute_make_db(alchemist, args.db_type, folder_path=args.folder_path,
                    folder_name=args.folder_name, values=values,
                    verbose=args.verbose, filters=args.filters,
                    groups=args.groups, db_name=args.database_name,
                    threads=args.threads, use_mpi=args.use_mpi)


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

    subparsers = parser.add_subparsers(dest="db_type", help=DATABASE_TYPE_HELP)

    hhsuite_parser = subparsers.add_parser("hhsuite")
    blast_parser = subparsers.add_parser("blast")

    hhsuite_parser.add_argument("-np", "--threads", type=int, nargs="?",
                                help=NUMBER_THREADS_HELP)
    hhsuite_parser.add_argument("-mpi", "--use_mpi", action="store_true",
                                help=MPI_HELP)

    blast_parser

    for subparser in [hhsuite_parser, blast_parser]:
        subparser.add_argument("-m", "--folder_name", type=str,
                               help=FOLDER_NAME_HELP)
        subparser.add_argument("-o", "--folder_path",
                               type=pipelines_basic.convert_dir_path,
                               help=FOLDER_PATH_HELP)
        subparser.add_argument("-v", "--verbose", action="store_true",
                               help=VERBOSE_HELP)

        subparser.add_argument("-n", "--database_name", type=str,
                               help=DATABASE_NAME_HELP)

        subparser.add_argument("-if", "--import_file",
                               type=pipelines_basic.convert_file_path,
                               help=IMPORT_FILE_HELP, dest="input")
        subparser.add_argument("-in", "--import_names", nargs="*",
                               help=IMPORT_NAMES_HELP, dest="input")
        subparser.add_argument("-w", "--where", nargs="?",
                               help=WHERE_HELP,
                               dest="filters")
        subparser.add_argument("-g", "--group_by", nargs="+",
                               help=GROUP_BY_HELP, dest="groups")

        subparser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                               folder_path=None, verbose=False, input=[],
                               filters="", groups=[],
                               database_name=None, threads=1, use_mpi=False)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_make_db(alchemist, db_type, values=None, folder_path=None,
                    folder_name=DEFAULT_FOLDER_NAME, verbose=False, filters="",
                    groups=[], db_name=None, threads=1, use_mpi=False,
                    mol_type=None, hash_index=False, parse_seqids=True,
                    gi_mask=False, mask_data=None, mask_id=None, logfile=None,
                    tax_id=None, tax_id_map=None):
    if db_name is None:
        db_name = alchemist.database

    if verbose:
        print("Retrieving database version...")
    db_version = mysqldb_basic.get_first_row_data(alchemist.engine, "version")

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
                                          mapped_path, db_name, db_version,
                                          data_cache=data_cache,
                                          threads=threads, verbose=verbose,
                                          use_mpi=use_mpi)
        elif db_type == "blast":
            execute_make_blast_database(
                    alchemist, db_filter.values, mapped_path, db_name,
                    db_version, data_cache=data_cache, verbose=verbose,
                    hash_index=False, parse_seqids=True, gi_mask=False,
                    mask_data=None, mask_id=None, logfile=None, tax_id=None,
                    tax_id_map=None)


def execute_make_hhsuite_database(alchemist, values, db_dir, db_name,
                                  db_version,
                                  data_cache=None, verbose=False, threads=1,
                                  use_mpi=False):
    aln_dir = db_dir.joinpath("pham_alignments")
    aln_dir.mkdir()

    create_pham_alignments(
                    alchemist, values, aln_dir, data_cache=data_cache,
                    threads=threads, verbose=verbose)

    if use_mpi:
        physical_cores = cpu_count(logical=False)

        if physical_cores < threads:
            if verbose:
                print("Designated process count greater than machine's "
                      "number of physical cores...\n"
                      f"STEPPING DOWN TO {physical_cores} PROCESSES")
            threads = physical_cores

    hhsuitedb.create_hhsuitedb(aln_dir, db_dir, db_name, cores=threads,
                               verbose=verbose, versions=db_version)


def execute_make_blast_database(alchemist, values, db_dir, db_name, db_version,
                                mol_type="prot", data_cache={}, verbose=False,
                                hash_index=False, parse_seqids=True,
                                gi_mask=False, mask_data=None, mask_id=None,
                                logfile=None, tax_id=None, tax_id_map=None):
    fasta_path = create_genes_fasta(alchemist, values, db_dir, db_name,
                                    mol_type=mol_type, verbose=verbose,
                                    data_cache=data_cache)

    blastdb.create_blastdb(fasta_path, db_dir, db_name, mol_type=mol_type,
                           hash_index=hash_index, parse_seqids=parse_seqids,
                           gi_mask=gi_mask, mask_data=mask_data,
                           mask_id=mask_id, logfile=logfile, tax_id=tax_id,
                           tax_id_map=tax_id_map, pham_versions=db_version)


def execute_make_mmseqs_database(alchemist, values, db_dir, db_name,
                                 data_cache=None, verbose=False, threads=1):
    pass


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def create_pham_alignments(alchemist, values, aln_dir, data_cache=None,
                           threads=1, verbose=False):
    if data_cache is None:
        data_cache = {}

    if verbose:
        print("Retrieving pham gene translations...")
    pham_ts_to_id = alignment.get_pham_gene_translations(alchemist, values)

    if verbose:
        print("Writing/aligning pham amino acid sequences to file...")

    aln_path_map = alignment.create_pham_alns(
                               alchemist, aln_dir, pham_ts_to_id,
                               cores=threads, verbose=verbose)

    pham_annotation_map = alignment.get_pham_gene_annotations(
                                                         alchemist, values)

    path_name_map = {}
    for pham, path in aln_path_map.items():
        annotation = pham_annotation_map.get(pham, "")
        if annotation == "":
            annotation = "hypothetical protein"

        path_name_map[path] = f"{annotation} [pham {pham}]"

    if verbose:
        print("Adding names to pham alignments...")
    pde_fileio.name_comment_files(path_name_map, threads=threads,
                                  verbose=verbose)

    return aln_path_map


def create_genes_fasta(alchemist, values, fasta_dir, db_name,
                       mol_type="prot", verbose=False, data_cache=None):
    if data_cache is None:
        data_cache = {}

    fasta_path = fasta_dir.joinpath(".".join([db_name, "fasta"]))

    gs_to_ts = map_translations(alchemist, values)

    pdm_fileio.write_fasta(gs_to_ts, fasta_path)

    return fasta_path


def map_translations(alchemist, pham_ids):
    gene = alchemist.metadata.tables["gene"]

    pham_id = gene.c.PhamID
    gene_id = gene.c.GeneID
    translation = gene.c.Translation

    query = querying.build_select(alchemist.graph, [gene_id, translation])
    results = querying.execute(alchemist.engine, query, in_column=pham_id,
                               values=pham_ids)

    gs_to_ts = {}
    for result in results:
        gs_to_ts[result["GeneID"]] = result["Translation"].decode("utf-8")

    return gs_to_ts


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
