import argparse
from decimal import Decimal
import math
import multiprocessing
import shutil
import sys
import time
from threading import Lock
from pathlib import Path
import pickle
import random

from Bio.Application import ApplicationError
from pdm_utils.functions import basic
from pdm_utils.functions import configfile
from pdm_utils.functions import fileio as pdm_fileio
from pdm_utils.functions import parallelize
from pdm_utils.functions import pipelines_basic

from pde_utils.classes import pan_models
from pde_utils.functions import alignment
from pde_utils.functions import multithread
from pde_utils.functions import pan_handling
from pde_utils.functions import fileio as pde_fileio
from pde_utils.functions import search

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
TEMP_DIR = Path("/tmp/pde_utils_build_pan_cache")
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pan"

PAN_GRAPH_EDGEWEIGHTS = ["CentroidDistance", "DBSeparation", "MinDistance"]
# MAIN FUNCTIONS
# -----------------------------------------------------------------------------


def main(unparsed_args_list):
    args = parse_build_pan(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)

    values = pipelines_basic.parse_value_input(args.input)

    execute_build_pan(alchemist, hhdb_path=args.hhsuite_database,
                      folder_path=args.folder_path,
                      folder_name=args.folder_name,
                      values=values, verbose=args.verbose,
                      filters=args.filters, groups=args.groups,
                      threads=args.number_threads, M=args.min_percent_gaps,
                      aD=args.avg_distance, mD=args.min_distance,
                      B=args.DB_stiffness, PANgraph_out=args.PANgraph_out)


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
    AVG_DISTANCE_HELP = """
        PAN building option that controls the initial centroid distance cutoff
        between phams during the building of neighborhoods.
            Follow selection argument with the minimum average identity
            [0-100].
        """
    MIN_DISTANCE_HELP = """
        PAN building option that controls the final centroid to target pham
        distance cutoff between phams during the building of neighborhoods.
            Follow selection argument with the minimum identity [0-100].
        """
    MATCH_STATE_CUTOFF_HELP = """
        PAN building option that controls the maximum percentage of gaps
        at a particular column that defines a match state.
            Follow selection argument with the minimum percent gaps [0-100].
        """
    PAN_GRAPH_OUT_HELP = """
        PAN export option that selects the output graph dump file type.
            Follow selection argument with a valid graph file type.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str,
                        help=DATABASE_HELP)

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
    parser.add_argument("-aD", "--avg_distance", type=Decimal,
                        help=AVG_DISTANCE_HELP)
    parser.add_argument("-mD", "--min_distance", type=Decimal,
                        help=MIN_DISTANCE_HELP)
    parser.add_argument("-M", "--min_percent_gaps", type=int,
                        help=MATCH_STATE_CUTOFF_HELP)

    parser.add_argument("-hhdb", "--hhsuite_database", type=Path,
                        help=HHSUITE_DATABASE_HELP)
    parser.add_argument("-gout", "--PANgraph_out", type=str,
                        choices=pde_fileio.NETWORKX_FILE_TYPES,
                        help=PAN_GRAPH_OUT_HELP)

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
                        filters="", groups=[], hhsuite_database=None,
                        db_name=None, number_threads=1,
                        DB_stiffness=0.2, avg_distance=75, min_distance=65,
                        min_percent_gaps=50, PANgraph_out=None)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_build_pan(alchemist, hhdb_path=None, pan_name=None,
                      folder_path=None, folder_name=DEFAULT_FOLDER_NAME,
                      values=None, verbose=False, filters="", groups=[],
                      threads=1, M=50, aD=75, mD=65, B=0.2,
                      PANgraph_out=None):
    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                             values=values)

    working_path = pipelines_basic.create_working_path(
                                            folder_path, folder_name)

    conditionals_map = pipelines_basic.build_groups_map(
                                            db_filter, working_path,
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

        if pan_name is None:
            pan_name = folder_name

        pipelines_basic.create_working_dir(mapped_path)
        pan_path = mapped_path.joinpath(".".join([pan_name, "sqlite"]))

        pan_alchemist = pan_handling.build_pan(pan_path)
        pan_alchemist.expire_on_commit = True

        pham_data_dir = mapped_path.joinpath("pham_alns")
        pham_data_dir.mkdir()
        data_maps_tuple = create_pham_alns(
                                        alchemist.engine, db_filter.values,
                                        pham_data_dir, threads=threads, M=M,
                                        verbose=verbose)

        build_pan_nodes(pan_alchemist, db_filter.values, data_maps_tuple,
                        threads=threads, verbose=verbose)

        cent_data_dir = mapped_path.joinpath("cent_alns")
        cent_data_dir.mkdir()
        build_pan_neighborhoods(alchemist, pan_alchemist, db_filter.values,
                                cent_data_dir, data_maps_tuple,
                                aD=aD, mD=mD, B=B,
                                threads=threads, verbose=verbose)

        hmm_data_dir = mapped_path.joinpath("pham_hhrs")
        hmm_data_dir.mkdir()

        if hhdb_path is not None:
            raise NotImplementedError(
                                "Town building is not implemented yet... :(")
            if verbose:
                print("...Calculating pham HMM profiles...")
            hmm_path_map = alignment.create_hmms(data_maps_tuple[0], name=True,
                                                 M=M, threads=threads,
                                                 verbose=verbose)
            build_pan_towns(alchemist, pan_alchemist, hhdb_path,
                            hmm_data_dir, hmm_path_map,
                            threads=threads, verbose=verbose)

        if PANgraph_out is not None:
            pan_graph = pan_handling.to_networkx(pan_alchemist)
            pde_fileio.write_graph(pan_graph, PANgraph_out,
                                   mapped_path, pan_name,
                                   edge_weights=PAN_GRAPH_EDGEWEIGHTS)

        shutil.rmtree(Path(TEMP_DIR))


def build_pan_nodes(pan_alchemist, values, data_maps_tuple, threads=1,
                    verbose=False):
    pan_nodes_dict = {}
    pan_lock = Lock()

    work_items = []
    for pham in values:
        aln_path = data_maps_tuple[0].get(pham)
        mat_path = data_maps_tuple[1].get(pham)
        tree_path = data_maps_tuple[2].get(pham)

        work_items.append((pan_nodes_dict, pham, aln_path, mat_path,
                           tree_path, pan_lock))

    multithread.multithread(build_pan_nodes_threadtask, work_items,
                            threads=threads)

    if verbose:
        print("...Writing pham data to PAN...")

    clusters = list(pan_nodes_dict.values())
    pan_alchemist.session.add_all(clusters)
    pan_alchemist.session.commit()
    pan_alchemist.session.close()


def build_pan_neighborhoods(alchemist, pan_alchemist, values, data_dir,
                            data_maps_tuple, aD=75, mD=65, B=0.2,
                            threads=1, verbose=False):
    matrix_chunks = create_centroid_graph(pan_alchemist, values, data_dir,
                                          threads=threads, verbose=verbose)

    thread_manager = multiprocessing.Manager()
    mD_cache = thread_manager.dict()
    data_cache = thread_manager.dict()
    path_cache = thread_manager.dict()

    temp_dir = Path(TEMP_DIR).joinpath("linker_files")
    temp_dir.mkdir()

    read_work_set = set()
    aln_work_items = []

    if verbose:
        print("...Constructing base for pham neighborhoods...")
    construct_neighborhood_base(pan_alchemist, matrix_chunks,
                                read_work_set, aln_work_items, data_dir,
                                temp_dir, mD_cache, data_cache, path_cache,
                                aD, mD, B)

    if verbose:
        print("...Reloading pham neighborhood cluster data...")
    for cluster in read_work_set:
        path_cache[cluster] = (data_maps_tuple[0].get(int(cluster)),
                               data_maps_tuple[1].get(int(cluster)),
                               data_maps_tuple[2].get(int(cluster)))

    aln_work_chunks = basic.partition_list(aln_work_items,
                                           int(math.sqrt(len(aln_work_items))))
    if verbose:
        print("...Computing pham cluster minimum distances...")
    identity_edge_chunks = parallelize.parallelize(
                                             aln_work_chunks, threads,
                                             build_neighborhood_edge_process,
                                             verbose=verbose)
    identity_edges = []
    for chunk in identity_edge_chunks:
        identity_edges = identity_edges + chunk

    if verbose:
        print("...Writing neighborhood data to PAN...")
    pan_alchemist.session.add_all(identity_edges)
    pan_alchemist.session.commit()
    pan_alchemist.session.close()


def build_pan_towns(alchemist, pan_alchemist, hhdb_path, pan_dict,
                    hmm_data_dir, data_maps_tuple, threads=1, verbose=False):
    work_items = []
    hhr_path_map = {}
    for pham, hmm_path in data_maps_tuple[3].items():
        hhr_path = hmm_data_dir.joinpath(".".join([str(pham), "hhr"]))
        hhr_path_map[pham] = hhr_path

        work_items.append((hmm_path, hhdb_path, hhr_path, None, False,
                           1, 0, 0, 1))

    if verbose:
        print("...Performing iterations of hhblitz to find HMM-HMM "
              "relationships...")
    parallelize.parallelize(work_items, threads, search.hhblits,
                            verbose=verbose)


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

    mat_path_map = {}
    tree_path_map = {}
    for pham, aln_path in aln_path_map.items():
        mat_path = aln_path.with_name(".".join([str(pham), "mat"]))
        mat_path_map[pham] = mat_path

        tree_path = aln_path.with_name(".".join([str(pham), "tree"]))
        tree_path_map[pham] = tree_path

    return (aln_path_map, mat_path_map, tree_path_map)


def create_centroid_graph(pan_alchemist, clusters, aln_dir, threads=1,
                          verbose=False):
    thread_manager = multiprocessing.Manager()
    cluster_data = thread_manager.list()

    cluster_data_dicts = pan_handling.retrieve_cluster_data(
                                                    pan_alchemist, clusters)

    for cluster_data_dict in cluster_data_dicts:
        cluster_data.append((cluster_data_dict["ClusterID"],
                             cluster_data_dict["CentroidSeq"].decode("utf-8")))

    work_items = []
    for i in range(len(cluster_data)):
        work_items.append((i, cluster_data))
    random.shuffle(work_items)

    temp_dir_path = Path(TEMP_DIR)
    if temp_dir_path.is_dir():
        shutil.rmtree(temp_dir_path)

    temp_dir_path.mkdir()

    if verbose:
        print("...Calculating centroid Levenshtein distances...")
    matrix_chunks = parallelize.parallelize(
                                work_items, threads,
                                create_centroid_graph_process, verbose=verbose)

    return matrix_chunks


def create_centroid_graph_process(index, cluster_data):
    len_cluster_data = len(cluster_data)
    if index == len_cluster_data:
        return []

    source = cluster_data[index][0]
    source_str = cluster_data[index][1]

    size = int(math.sqrt(len_cluster_data - index))
    edges = []
    for i in range(((len_cluster_data - index + 1) + size - 1) // size):
        process_cluster_data = cluster_data[
                                       index+1+(i*size):index+((i+1)*size)]

        for target_data in process_cluster_data:
            edges.append(distance_to_graph(source, source_str,
                                           target_data[0], target_data[1]))

    filepath = Path(TEMP_DIR).joinpath(str(source))
    with filepath.open(mode="wb") as filehandle:
        pickle.dump(edges, filehandle)

    return filepath


def construct_neighborhood_base(pan_alchemist, matrix_chunks,
                                read_work_set, aln_work_items, data_dir,
                                temp_dir, mD_cache, data_cache, path_cache,
                                aD, mD, B):
    for distance_row in matrix_chunks:
        with distance_row.open(mode="rb") as filehandle:
            distance_edges = pickle.load(filehandle)

        for edge in distance_edges:
            distance = edge[2]["CentroidDistance"]
            source = edge[0]
            target = edge[1]

            if source == target:
                continue

            if distance > aD:
                continue

            source_node = data_cache.get(source)
            if source_node is None:
                source_node = pan_handling.retrieve_cluster_data(
                                                    pan_alchemist, [source])[0]
                data_cache[source] = source_node
            s_spread = source_node["Spread"]

            target_node = data_cache.get(target)
            if target_node is None:
                target_node = pan_handling.retrieve_cluster_data(
                                                    pan_alchemist, [target])[0]
                data_cache[target] = target_node
            t_spread = target_node["Spread"]

            if distance < 50:
                print(f"......{source} and {target} seem too friendly...")

            fwd_DBsep = (t_spread/distance)
            rvs_DBsep = (s_spread/distance)
            if fwd_DBsep >= B or rvs_DBsep >= B:
                read_work_set.add(source)
                read_work_set.add(target)

                if fwd_DBsep >= B:
                    aln_work_items.append((
                            source, target, distance, fwd_DBsep, mD, temp_dir,
                            data_dir, data_cache, mD_cache, path_cache))
                if rvs_DBsep >= B:
                    aln_work_items.append((
                            target, source, distance, rvs_DBsep, mD, temp_dir,
                            data_dir, data_cache, mD_cache, path_cache))


def distance_to_graph(source_id, source_str, target_id, target_str):
    distance = alignment.calculate_levenshtein(source_str, target_str)

    return (source_id, target_id, {"CentroidDistance": distance})


def build_pan_nodes_threadtask(pan_dict, pham, aln_path, mat_path,
                               tree_path, pan_lock):
    cluster = pan_handling.parse_cluster(pham, MSA_path=aln_path,
                                         PIM_path=mat_path, GT_path=tree_path)
    cluster.parse_centroid()
    cluster.parse_spread()

    pan_lock.acquire()
    pan_dict[pham] = cluster
    pan_lock.release()


def write_centroids_threadtask(temp_dir, source_path, target_path,
                               source_gs_to_ts, target_gs_to_ts):
    temp_dir.mkdir()

    pdm_fileio.write_fasta(source_gs_to_ts, source_path)
    pdm_fileio.write_fasta(target_gs_to_ts, target_path)


def build_neighborhood_edge_process(work_chunk):
    identity_edges = []

    for work_item in work_chunk:
        identity_edge = build_identity_edge(*work_item)
        if identity_edge is not None:
            identity_edges.append(identity_edge)

    return identity_edges


def build_identity_edge(source, target, avg_dis, DBsep, min_mD,
                        temp_dir, data_dir, data_cache, mD_cache,
                        path_cache):
    source_aln_files = path_cache[source]
    source_cluster = pan_handling.parse_cluster(
                            source, data_cache[source], source_aln_files[0],
                            source_aln_files[1], source_aln_files[2])
    source_cluster.parse_centroid()

    target_aln_files = path_cache[target]
    target_cluster = pan_handling.parse_cluster(
                            target, data_cache[target], target_aln_files[0],
                            target_aln_files[1], target_aln_files[2])
    target_cluster.parse_centroid()

    min_id = get_mD_cache_edge(source, target, mD_cache)

    if min_id is None:
        min_dis = estimate_min_distance(source_cluster, target_cluster,
                                        temp_dir, data_dir)

        mD_cache[source] = {target: min_dis}

    if min_dis <= min_mD:
        return pan_models.IdentityEdge(
                                Source=source, Target=target,
                                CentroidDistance=avg_dis,
                                DBSeparation=DBsep,
                                MinDistance=min_dis)


def get_mD_cache_edge(source, target, mD_cache):
    source_node = mD_cache.get(source)
    target_node = mD_cache.get(source)

    if source_node is not None:
        mI = source_node.get(target)
        if mI is not None:
            return mI
    elif target_node is not None:
        mI = target_node.get(target)
        if mI is not None:
            return mI

    return None


def estimate_min_distance(source_cluster, target_cluster, temp_dir, data_dir,
                          align=False):
    if source_cluster.GT is not None:
        source_leaf = estimate_linker_sequence(source_cluster,
                                               target_cluster.centroid_seq_str,
                                               temp_dir)

        source_linker = source_leaf.name
        source_linker_seq = source_leaf.comment["Sequence"]
        source_linker_path = source_leaf.comment["FilePath"]
    else:
        source_linker = source_cluster.CentroidID
        source_linker_seq = source_cluster.centroid_seq_str
        source_linker_path = temp_dir.joinpath(".".join([source_linker,
                                                         "fasta"]))

    if target_cluster.GT is not None:
        target_leaf = estimate_linker_sequence(target_cluster,
                                               source_cluster.centroid_seq_str,
                                               temp_dir)

        target_linker = target_leaf.name
        target_linker_seq = target_leaf.comment["Sequence"]
        target_linker_path = target_leaf.comment["FilePath"]
    else:
        target_linker = target_cluster.CentroidID
        target_linker_seq = target_cluster.centroid_seq_str
        target_linker_path = temp_dir.joinpath(".".join([target_linker,
                                                         "fasta"]))

    if align:
        source_len = len(source_linker_seq)
        target_len = len(target_linker_seq)

        if not source_linker_path.is_file():
            pdm_fileio.write_fasta({source_linker: source_linker_seq},
                                   source_linker_path)

        if not target_linker_path.is_file():
            pdm_fileio.write_fasta({target_linker: target_linker_seq},
                                   target_linker_path)

        pairwise_path = data_dir.joinpath("".join([source_linker, "__",
                                                   target_linker, ".fasta"]))

        error = True
        while error:
            try:
                linker_alignment = alignment.pairwise_align(
                                        source_linker_path, target_linker_path,
                                        pairwise_path, tool="needle")
                error = False
            except ApplicationError:
                time.sleep(0.2)

        if source_len > target_len:
            linker_pid = linker_alignment.annotations["identity"] / source_len
        else:
            linker_pid = linker_alignment.annotations["identity"] / target_len

        linker_pid = linker_pid * 100
    else:
        linker_pid = alignment.calculate_levenshtein(source_linker_seq,
                                                     target_linker_seq)

    return linker_pid


def estimate_linker_sequence(source_cluster, target_centroid_seq, temp_dir):
    guidetree = source_cluster.GT

    curr_node = guidetree.clade
    while curr_node.is_bifurcating() and len(curr_node.clades) > 1:
        left_child = curr_node.clades[0]
        right_child = curr_node.clades[1]

        left_leaf = get_furthest_sequence(source_cluster, left_child)
        right_leaf = get_furthest_sequence(source_cluster, right_child)

        left_pid = left_leaf.comment.get("TargetIdentity")
        if left_pid is None:
            left_pid = alignment.calculate_levenshtein(
                                          left_leaf.comment["Sequence"],
                                          target_centroid_seq, identity=True)
            left_leaf.comment["TargetIdentity"] = left_pid
            pass

        right_pid = right_leaf.comment.get("TargetIdentity")
        if right_pid is None:
            right_pid = alignment.calculate_levenshtein(
                                          right_leaf.comment["Sequence"],
                                          target_centroid_seq, identity=True)
            right_leaf.comment["TargetIdentity"] = right_pid

        if left_pid > right_pid:
            curr_node = left_child
        else:
            curr_node = right_child

    linker = get_furthest_sequence(source_cluster, curr_node)

    if linker.comment.get("FilePath") is None:
        linker_path = temp_dir.joinpath(".".join([linker.name,
                                                  "fasta"]))
        if not linker_path.is_file():
            pdm_fileio.write_fasta({linker.name: linker.comment["Sequence"]},
                                   linker_path)
        linker.comment["FilePath"] = linker_path

    return linker


def get_furthest_sequence(cluster, seq_clade):
    leaves = []
    for terminal_clade in seq_clade.get_terminals():
        leaves.append(terminal_clade)

    centroid_id = cluster.CentroidID

    furthest_leaf = None
    max_distance = 0
    for leaf in leaves:
        seq_id = leaf.name

        if leaf.comment is None:
            leaf.comment = {}

        distance = leaf.comment.get("SourceIdentity")

        if distance is None:
            distance = cluster.PIM.get_distance(centroid_id, seq_id)
            leaf.comment["SourceIdentity"] = distance

        if furthest_leaf is None:
            furthest_leaf = leaf
            max_distance = distance
        elif distance > max_distance:
            furthest_leaf = leaf
            max_distance = distance

    sequence = furthest_leaf.comment.get("Sequence")
    if sequence is None:
        furthest_leaf.comment["Sequence"] = cluster.MSA.get_sequence(
                                                furthest_leaf.name, gaps=False)
    return furthest_leaf


def correct_identity(identity, source_len, target_len):
    if source_len < target_len:
        shortest = source_len
        longest = target_len
    else:
        shortest = target_len
        longest = source_len

    return (identity * shortest) / longest


if __name__ == "__main__":
    main(sys.argv)
