import argparse
import sys
import time
from decimal import Decimal
from pathlib import Path

from networkx import MultiDiGraph, Graph
from pdm_utils.functions import basic
from pdm_utils.functions import fileio as pdm_fileio
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import phameration

from pde_utils.classes import clustal
from pde_utils.functions import alignment
from pde_utils.functions import fileio as pde_fileio
from pde_utils.functions import search

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pan"
DEFAULT_FOLDER_PATH = Path.cwd()

NETWORKX_FILE_TYPES = ["csv", "gexf", "gpickle", "graphml", "al", "nl",
                       "json", "cyjs", "yaml", "pajek", "shp"]
EDGE_WEIGHTS = ["DB", "Distance"]

# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(unparsed_args_list):
    args = parse_build_pan(unparsed_args_list)

    alchemist = pipelines_basic.build_alchemist(args.database)
    values = pipelines_basic.parse_value_input(args.input)

    execute_build_pan(alchemist, args.hhsuite_database, args.folder_path,
                      args.folder_name, values=values, verbose=args.verbose,
                      filters=args.filters, groups=args.groups,
                      threads=args.number_threads, M=args.min_percent_gaps,
                      aI=args.avg_identity, mI=args.min_identity,
                      B=args.DB_stiffness, file_format=args.file_format)


def parse_build_pan(unparsed_args_list):
    DATABASE_HELP = """
        Name of the MySQL database to build from.
        """
    HHSUITE_DATABASE_HELP = """
        Name of the HHSUITE database to build from.
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

    NUMBER_THREADS_HELP = """
        Pipeline option that allows for multithreading of workflows.
            Follow selection argument with number of threads to be used
        """
    DB_STIFFNESS_HELP = """
        PAN building option that controls the cutoff for the DB separation
        value required to make a connection.
            Follow seelction argument with the DB stiffness parameter [0-1].
        """
    AVG_IDENTITY_HELP = """
        PAN building option that controls the initial centroid identity cutoff
        between phams during the building of neighborhoods.
            Follow selection argument with the minimum average identity
            [0-100].
        """
    MIN_IDENTITY_HELP = """
        PAN building option that controls the final centroid to target pham
        identity cutoff between phams during the building of neighborhoods.
            Follow selection argument with the minimum identity [0-100].
        """
    MATCH_STATE_CUTOFF_HELP = """
        PAN building option that controls the maximum percentage of gaps
        at a particular column that defines a match state.
            Follow selection argument with the minimum percent gaps [0-100].
        """
    FILE_FORMAT_HELP = """
        PAN export option that selects the output file type.
            Follow selection argument with a valid graph file type.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str,
                        help=DATABASE_HELP)
    parser.add_argument("hhsuite_database", type=convert_hmm_db_path,
                        help=HHSUITE_DATABASE_HELP)

    parser.add_argument("-m", "--folder_name", type=str,
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path",
                        type=pipelines_basic.convert_dir_path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)

    parser.add_argument("-B", "--DB_stiffness", type=Decimal,
                        help=DB_STIFFNESS_HELP)
    parser.add_argument("-aI", "--avg_identity", type=Decimal,
                        help=AVG_IDENTITY_HELP)
    parser.add_argument("-mI", "--min_identity", type=Decimal,
                        help=MIN_IDENTITY_HELP)
    parser.add_argument("-M", "--min_percent_gaps", type=int,
                        help=MATCH_STATE_CUTOFF_HELP)

    parser.add_argument("-t", "--file_format", type=str,
                        choices=NETWORKX_FILE_TYPES,
                        help=FILE_FORMAT_HELP)
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
    parser.add_argument("-np", "--number_threads", type=int, nargs="?",
                        help=NUMBER_THREADS_HELP)

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH,
                        verbose=False, input=[],
                        filters="", groups=[],
                        db_name=None, number_threads=1,
                        DB_stiffness=0.5, avg_identity=20,
                        min_identity=30,
                        min_percent_gaps=50, file_format="graphml")

    parsed_args = parser.parse_args(unparsed_args_list[1:])
    return parsed_args


def execute_build_pan(alchemist, hhdb_path, folder_path, folder_name,
                      values=None, verbose=False, filters="", groups=[],
                      threads=1, M=50, aI=20, mI=30, B=0.5, file_format="csv"):
    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                             values=values)

    pan_path = folder_path.joinpath(folder_name)
    pan_path = basic.make_new_dir(folder_path, pan_path, attempt=50)

    conditionals_map = {}
    pipelines_basic.build_groups_map(db_filter, pan_path,
                                     conditionals_map,
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

        base_tmp = mapped_path.joinpath("base_tmp")
        base_tmp.mkdir()
        maps_tuple = create_pham_files(alchemist.engine, db_filter.values,
                                       base_tmp, threads=threads, M=M,
                                       verbose=verbose)
        pan = build_pan_base(alchemist, db_filter.values, maps_tuple,
                             data_cache=data_cache)

        nbhd_tmp = mapped_path.joinpath("nbhd_tmp")
        nbhd_tmp.mkdir()
        build_pan_neighborhoods(alchemist, pan, nbhd_tmp, aI=aI, mI=mI, B=B,
                                threads=threads, verbose=verbose,
                                data_cache=data_cache)

        pde_fileio.write_graph(pan, file_format, mapped_path, folder_name)


def create_pham_files(engine, values, aln_dir, threads=1, M=50, verbose=False,
                      data_cache=None):
    if data_cache is None:
        data_cache = {}

    fasta_path_map = alignment.create_pham_fasta(engine, values, aln_dir,
                                                 data_cache=data_cache,
                                                 threads=threads)
    aln_path_map = alignment.align_pham_fastas(fasta_path_map, mat_out=True,
                                               threads=threads, override=True)
    hmm_path_map = alignment.create_pham_hmms(fasta_path_map, name=True,
                                              M=M, threads=threads)

    mat_path_map = {}
    for pham, aln_path in aln_path_map.items():
        mat_path = aln_path.with_name(".".join([str(pham), "fasta"]))
        mat_path_map[pham] = mat_path

    return (aln_path_map, mat_path_map, hmm_path_map)


def build_pan_base(alchemist, values, base_maps_tuple, data_cache=None):
    if data_cache is None:
        data_cache = {}

    pan = MultiDiGraph()
    for pham in values:
        aln_path = base_maps_tuple[0].get(pham)
        mat_path = base_maps_tuple[1].get(pham)

        aln = clustal.MultipleSequenceAlignment(aln_path, fmt="fasta")
        aln.parse_alignment()

        if not mat_path.is_file():
            centroid, trunc_centroid_seq = aln.longest_gene()
            spread = 0
        else:
            mat = clustal.PercentIdentityMatrix(mat_path)
            mat.parse_matrix()
            centroid = mat.get_centroid()
            spread = 100 - mat.get_average_identity(centroid)

        pham_genes = data_cache.get(pham)
        if pham_genes is None:
            pham_genes = alignment.get_pham_genes(alchemist.engine, pham)

        centroid_seq = pham_genes.get(centroid)

        pan.add_node(pham, Spread=spread, CentroidGene=centroid,
                     CentroidSeq=centroid_seq)

    return pan


def build_pan_neighborhoods(alchemist, pan, tmp_dir, aI=20, mI=30, B=0.5,
                            threads=1, verbose=False, data_cache=None):
    if data_cache is None:
        data_cache = {}

    centroid_seqs = {}
    cfasta_path = tmp_dir.joinpath("centroid.fasta")

    pan_dict = dict(pan.nodes(data=True))
    for pham, pham_data in pan_dict.items():
        centroid_seqs[pham] = pham_data["CentroidSeq"]

    pdm_fileio.write_fasta(centroid_seqs, cfasta_path)

    caln, cmat = align_centroids(cfasta_path, tmp_dir, threads=threads)

    for pham, pham_data in pan_dict.items():
        centroid = pham_data["CentroidGene"]
        centroid_seq = pham_data["CentroidSeq"]
        locals = cmat.get_nearest_neighbors(str(pham), aI)

        neighbors = {}
        for local in locals:
            t_spread = pan_dict[int(local)]["Spread"]
            distance = cmat.get_distance(str(pham), local)

            print(t_spread)
            print(distance)

            DBsep = ((t_spread)/distance)
            print(DBsep)
            if DBsep >= B:
                neighbors.update({int(local): {"DBsep": t_spread,
                                               "distance": distance}})

        for neighbor, edge_data in neighbors.items():
            fasta_path = tmp_dir.joinpath(f"{neighbor}_{pham}cent.fasta")
            aln_path = tmp_dir.joinpath(f"{neighbor}_{pham}cent.aln")
            mat_path = tmp_dir.joinpath(f"{neighbor}_{pham}cent.mat")

            pham_genes = data_cache.get(neighbor)
            if pham_genes is None:
                pham_genes = alignment.get_pham_genes(
                                        alchemist.engine, neighbor)

            pham_genes[centroid] = centroid_seq
            pdm_fileio.write_fasta(pham_genes, fasta_path)
            alignment.clustalo(fasta_path, aln_path, mat_out_path=mat_path,
                               threads=threads)

            mat = clustal.PercentIdentityMatrix(mat_path)
            mat.parse_matrix()
            min_dist_seqs = mat.get_nearest_neighbors(centroid, mI)

            if min_dist_seqs:
                closest_gene = min_dist_seqs[0]
                min_distance = mat.get_distance(centroid, closest_gene)

                pan.add_edge(pham, str(neighbor),
                             MinDistance=min_distance,
                             LinkingGene=closest_gene,
                             Distance=edge_data["distance"],
                             DBSeparation=edge_data["DBsep"])

    return pan


def convert_hmm_db_path(path):
    return Path(path)

# Methods of clustering centroid sequences


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
