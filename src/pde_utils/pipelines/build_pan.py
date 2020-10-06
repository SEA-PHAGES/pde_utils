import argparse
import sys
import time
from decimal import Decimal
from threading import Lock
from pathlib import Path

from Bio import Phylo
from networkx import Graph
from pdm_utils.functions import configfile
from pdm_utils.functions import fileio as pdm_fileio
from pdm_utils.functions import parallelize
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import phameration

from pde_utils.classes import clustal
from pde_utils.classes import pan_models
from pde_utils.functions import alignment
from pde_utils.functions import multithread
from pde_utils.functions import pan_handling
# from pde_utils.functions import fileio as pde_fileio
# from pde_utils.functions import search

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pan"

# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(unparsed_args_list):
    args = parse_build_pan(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)

    execute_build_pan(alchemist, args.hhsuite_database,
                      folder_path=args.folder_path,
                      folder_name=args.folder_name,
                      values=values, verbose=args.verbose,
                      filters=args.filters, groups=args.groups,
                      threads=args.number_threads, M=args.min_percent_gaps,
                      aI=args.avg_identity, mI=args.min_identity,
                      B=args.DB_stiffness)


def parse_build_pan(unparsed_args_list):
    DATABASE_HELP = """
        Name of the MySQL database to build from.
        """
    HHSUITE_DATABASE_HELP = """
        Name of the HHSUITE database to build from.
        """

    VERBOSE_HELP = """
        Build PAN option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Build PAN option to change the path
        of the directory where the exported files are stored.
            Follow selection argument with the path to the
            desired export directory.
        """
    FOLDER_NAME_HELP = """
        Build PAN option to change the name
        of the directory where the exported files are stored.
            Follow selection argument with the desired name.
        """
    CONFIG_FILE_HELP = """
        Build PAN option that enables use of a config file for sourcing
        credentials
            Follow selection argument with the path to the config file
            specifying MySQL and NCBI credentials.
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
    # FILE_FORMAT_HELP = """
    #    PAN export option that selects the output file type.
    #        Follow selection argument with a valid graph file type.
    # """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str,
                        help=DATABASE_HELP)
    parser.add_argument("hhsuite_database", type=Path,
                        help=HHSUITE_DATABASE_HELP)

    parser.add_argument("-m", "--folder_name", type=str,
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path", type=Path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)
    parser.add_argument("-c", "--config_file", help=CONFIG_FILE_HELP,
                        type=pipelines_basic.convert_file_path)

    parser.add_argument("-B", "--DB_stiffness", type=Decimal,
                        help=DB_STIFFNESS_HELP)
    parser.add_argument("-aI", "--avg_identity", type=Decimal,
                        help=AVG_IDENTITY_HELP)
    parser.add_argument("-mI", "--min_identity", type=Decimal,
                        help=MIN_IDENTITY_HELP)
    parser.add_argument("-M", "--min_percent_gaps", type=int,
                        help=MATCH_STATE_CUTOFF_HELP)

    # parser.add_argument("-t", "--file_format", type=str,
    #                    choices=NETWORKX_FILE_TYPES,
    #                    help=FILE_FORMAT_HELP)
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
    parser.add_argument("-np", "--number_threads", type=int, nargs="?",
                        help=NUMBER_THREADS_HELP)

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, verbose=False, input=[],
                        filters="", groups=[],
                        db_name=None, number_threads=1,
                        DB_stiffness=0.5, avg_identity=20, min_identity=35,
                        min_percent_gaps=50)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_build_pan(alchemist, hhdb_path, pan_name=None,
                      folder_path=None, folder_name=DEFAULT_FOLDER_NAME,
                      values=None, verbose=False, filters="", groups=[],
                      threads=1, M=50, aI=20, mI=35, B=0.5):
    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                             values=values)

    working_path = pipelines_basic.create_working_path(
                                            folder_path, folder_name)

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

        if pan_name is None:
            pan_name = folder_name

        pipelines_basic.create_working_dir(mapped_path)
        pan_path = mapped_path.joinpath(".".join([pan_name, "sqlite"]))

        pan_alchemist = pan_handling.build_pan(pan_path)

        pham_data_dir = mapped_path.joinpath("pham_alns")
        pham_data_dir.mkdir()
        data_maps_tuple = create_pham_alns(
                                        alchemist.engine, db_filter.values,
                                        pham_data_dir, threads=threads, M=M,
                                        verbose=verbose, data_cache=data_cache)

        pan_dict = build_pan_nodes(alchemist, pan_alchemist, db_filter.values,
                                   data_maps_tuple, data_cache=data_cache,
                                   threads=threads, verbose=verbose)

        cent_data_dir = mapped_path.joinpath("cent_alns")
        cent_data_dir.mkdir()
        build_pan_neighborhoods(alchemist, pan_alchemist, pan_dict,
                                cent_data_dir, data_maps_tuple,
                                aI=aI, mI=mI, B=B, threads=threads,
                                verbose=verbose)


def build_pan_nodes(alchemist, pan_alchemist, values, data_maps_tuple,
                    data_cache=None, threads=1, verbose=False):
    if data_cache is None:
        data_cache = {}

    pan_dict = {}
    pan_lock = Lock()

    db_lock = Lock()

    work_items = []
    for pham in values:
        aln_path = data_maps_tuple[0].get(pham)
        mat_path = data_maps_tuple[1].get(pham)
        tree_path = data_maps_tuple[2].get(pham)
        hmm_path = data_maps_tuple[3].get(pham)
        pham_genes = data_cache.get(pham)

        work_items.append((alchemist, pan_dict, pham, aln_path, mat_path,
                           tree_path, hmm_path, pham_genes, pan_lock, db_lock))

    multithread.multithread(build_pan_nodes_threadtask, work_items,
                            threads=threads)

    if verbose:
        print("...Writing pham data to PAN...")
    for cluster in pan_dict.values():
        pan_alchemist.session.add(cluster)

    pan_alchemist.session.commit()
    return pan_dict


def build_pan_neighborhoods(alchemist, pan_alchemist, pan_dict, data_dir,
                            data_maps_tuple, aI=20, mI=35, B=0.5,
                            threads=1, verbose=False):
    cent_graph = create_centroid_graph(alchemist.engine, pan_dict, data_dir,
                                       threads=threads, verbose=verbose)

    if verbose:
        print("...Constructing pham neighborhoods...")
    identity_edges = []
    for pham, target_cluster in pan_dict.items():
        for neighbor, edge_weights in cent_graph[pham].items():
            if edge_weights["CentroidIdentity"] < aI:
                continue

            cent_edge = cent_graph[pham][neighbor]

            t_spread = target_cluster.Spread

            local_cluster = pan_dict[neighbor]
            identity = cent_edge["CentroidIdentity"]

            DBsep = ((t_spread)/(100 - identity))
            if DBsep >= B:
                min_identity = cent_edge.get("MinIdentity")
                if min_identity is None:
                    min_identity = estimate_min_identity(
                                       local_cluster, target_cluster, data_dir)

                    cent_edge["MinIdentity"] = min_identity

                if min_identity >= mI:
                    identity_edges.append(pan_models.IdentityEdge(
                                Source=neighbor, Target=pham,
                                DBSeparation=DBsep, CentroidIdentity=identity,
                                MinIdentity=min_identity))

    if verbose:
        print("...Writing neighborhood data to PAN...")
    pan_alchemist.session.add_all(identity_edges)
    pan_alchemist.session.commit()


# HELPER FUNCTIONS
# -----------------------------------------------------------------------------


def create_pham_alns(engine, values, aln_dir, threads=1, M=50, verbose=False,
                     data_cache=None):
    if data_cache is None:
        data_cache = {}

    if verbose:
        print("...Writing pham fastas...")
    fasta_path_map = alignment.create_pham_fastas(engine, values, aln_dir,
                                                  data_cache=data_cache,
                                                  threads=threads,
                                                  verbose=verbose)

    if verbose:
        print("...Aligning pham fastas...")
    aln_path_map = alignment.align_fastas(fasta_path_map,
                                          mat_out=True, tree_out=True,
                                          threads=threads, override=True,
                                          verbose=verbose)

    if verbose:
        print("...Calculating pham HMM profiles...")
    hmm_path_map = alignment.create_hmms(fasta_path_map, name=True,
                                         M=M, threads=threads,
                                         verbose=verbose)

    mat_path_map = {}
    tree_path_map = {}
    for pham, aln_path in aln_path_map.items():
        mat_path = aln_path.with_name(".".join([str(pham), "mat"]))
        mat_path_map[pham] = mat_path

        tree_path = aln_path.with_name(".".join([str(pham), "tree"]))
        tree_path_map[pham] = tree_path

    return (aln_path_map, mat_path_map, tree_path_map, hmm_path_map)


def create_centroid_graph(engine, pan_dict, aln_dir, threads=1, verbose=False):
    work_items = []
    fasta_path_map = {}
    cent_graph = Graph()
    for pham, cluster in pan_dict.items():
        fasta_path = aln_dir.joinpath("".join([str(pham), "_cent.fasta"]))
        fasta_path_map[pham] = fasta_path

        cent_graph.add_node(pham, CentroidID=cluster.CentroidID,
                            CentroidSeq=cluster.CentroidSeq.decode("utf-8"),
                            Spread=cluster.Spread)

        work_items.append(({pham: cluster.CentroidSeq.decode("utf-8")},
                           fasta_path))

    if verbose:
        print("...Writing pham centroid fastas...")
    multithread.multithread(pdm_fileio.write_fasta, work_items, threads)

    work_items = []
    for source_pham, source_fasta_path in fasta_path_map.items():
        for target_pham, target_fasta_path in fasta_path_map.items():
            if source_pham == target_pham:
                continue

            cent_edge = cent_graph[source_pham].get(target_pham)
            if cent_edge is not None:
                continue

            cent_aln_path = aln_dir.joinpath(
                                    "".join([str(source_pham), "_",
                                             str(target_pham), "_cents.aln"]))

            work_items.append((source_fasta_path, target_fasta_path,
                               cent_aln_path, "needle"))
            cent_graph.add_edge(source_pham, target_pham)

    if verbose:
        print("...Aligning pham centroids...")
    cent_alns = parallelize.parallelize(
                                work_items, threads, alignment.pairwise_align,
                                verbose=verbose)

    for cent_aln in cent_alns:
        identity = (cent_aln.annotations["identity"] /
                    cent_aln.get_alignment_length()) * 100

        source = alignment[0]
        target = alignment[0]

        cent_graph[source.id][target.id]["CentroidIdentity"] = identity

    return cent_graph


def build_pan_nodes_threadtask(alchemist, pan_dict, pham, aln_path, mat_path,
                               tree_path, hmm_path, pham_genes,
                               pan_lock, db_lock):
    aln = clustal.MultipleSequenceAlignment(aln_path, fmt="fasta")
    aln.parse_alignment()

    if not mat_path.is_file():
        mat = None
        tree = None
        centroid, trunc_centroid_seq = aln.longest_gene()
        spread = 0
    else:
        mat = clustal.PercentIdentityMatrix(mat_path)
        mat.parse_matrix()
        centroid = mat.get_centroid()
        spread = 100 - mat.get_average_identity(centroid)

        tree = Phylo.read(tree_path, "newick")

    if pham_genes is None:
        db_lock.acquire()
        pham_genes = alignment.get_pham_genes(alchemist.engine, pham)
        db_lock.release()

    centroid_seq = pham_genes.get(centroid)

    pan_lock.acquire()
    cluster = pan_models.Cluster(pham, Spread=spread, CentroidID=centroid,
                                 CentroidSeq=centroid_seq.encode("utf-8"),
                                 MultipleSequenceAlignment=aln,
                                 PercentIdentityMatrix=mat,
                                 AlignmentGuidetree=tree)
    pan_dict[pham] = cluster
    pan_lock.release()


def estimate_min_identity(source_cluster, local_cluster, data_dir):
    if source_cluster.AlignmentGuidetree is not None:
        source_linker = estimate_linker_sequence(
                                    source_cluster, local_cluster.CentroidID,
                                    local_cluster.CentroidSeq.decode("utf-8"),
                                    data_dir)
        source_linker_seq = source_cluster.MultipleSequenceAlignment\
                                          .get_sequence(source_linker)
    else:
        source_linker = source_cluster.CentroidID
        source_linker_seq = source_cluster.CentroidSeq.decode("utf-8")
    source_len = len(source_linker_seq)

    if local_cluster.AlignmentGuidetree is not None:
        target_linker = estimate_linker_sequence(
                                    local_cluster, source_cluster.CentroidID,
                                    source_cluster.CentroidSeq.decode("utf-8"),
                                    data_dir)
        target_linker_seq = local_cluster.MultipleSequenceAlignment\
                                         .get_sequence(target_linker)
    else:
        target_linker = local_cluster.CentroidID
        target_linker_seq = local_cluster.CentroidSeq.decode("utf-8")
    target_len = len(target_linker_seq)

    pairwise_aln_path = data_dir.joinpath("pairwise_linkers.aln")

    source_path = data_dir.joinpath("source_linker.fasta")
    pdm_fileio.write_fasta({source_linker: source_linker_seq}, source_path)

    target_path = data_dir.joinpath("target_linker.fasta")
    pdm_fileio.write_fasta({target_linker: target_linker_seq}, target_path)

    linker_alignment = alignment.pairwise_align(
                                            source_path, target_path,
                                            pairwise_aln_path, tool="water")

    if source_len < target_len:
        linker_pid = linker_alignment.annotations["identity"] / source_len
    else:
        linker_pid = linker_alignment.annotations["identity"] / target_len

    return linker_pid * 100


def estimate_linker_sequence(source_cluster, target_centroid_id,
                             target_centroid_seq, data_dir):
    pairwise_aln_path = data_dir.joinpath("pairwise_linkers.aln")
    guidetree = source_cluster.AlignmentGuidetree

    centroid_path = data_dir.joinpath("centroid.fasta")
    pdm_fileio.write_fasta({target_centroid_id: target_centroid_seq},
                           centroid_path)

    curr_node = guidetree.clade
    while curr_node.is_bifurcating() and len(curr_node.clades) > 1:
        left_child = curr_node.clades[0]
        left_seq_path = data_dir.joinpath("left_linker.fasta")

        right_child = curr_node.clades[1]
        right_seq_path = data_dir.joinpath("right_linker.fasta")

        left_seq_id = get_furthest_sequence(source_cluster, left_child)
        left_seq = source_cluster.MultipleSequenceAlignment.get_sequence(
                                                left_seq_id, gaps=False)
        pdm_fileio.write_fasta({left_seq_id: left_seq}, left_seq_path)

        right_seq_id = get_furthest_sequence(source_cluster, right_child)
        right_seq = source_cluster.MultipleSequenceAlignment.get_sequence(
                                                right_seq_id, gaps=False)
        pdm_fileio.write_fasta({right_seq_id: right_seq}, right_seq_path)

        left_alignment = alignment.pairwise_align(
                                            left_seq_path, centroid_path,
                                            pairwise_aln_path, tool="needle")
        right_alignment = alignment.pairwise_align(
                                            right_seq_path, centroid_path,
                                            pairwise_aln_path, tool="needle")

        left_pid = (left_alignment.annotations["identity"] /
                    left_alignment.get_alignment_length())

        right_pid = (right_alignment.annotations["identity"] /
                     right_alignment.get_alignment_length())

        if left_pid > right_pid:
            curr_node = left_child
        else:
            curr_node = right_child

    linker_id = get_furthest_sequence(source_cluster, curr_node)
    return linker_id


def get_furthest_sequence(cluster, seq_clade):
    seq_ids = []
    for terminal_clade in seq_clade.get_terminals():
        seq_ids.append(terminal_clade.name)

    centroid_id = cluster.CentroidID

    furthest_seq_id = centroid_id
    max_distance = 0
    for seq_id in seq_ids:
        distance = cluster.PercentIdentityMatrix.get_distance(centroid_id,
                                                              seq_id)
        if furthest_seq_id is None:
            furthest_seq_id = seq_id
            max_distance = distance
            continue

        if distance > max_distance:
            furthest_seq_id = seq_id
            max_distance = distance

    return furthest_seq_id


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
