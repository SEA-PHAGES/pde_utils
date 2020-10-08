import argparse
from decimal import Decimal
import multiprocessing
import random
import sys
import time
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
from pde_utils.functions import fileio as pde_fileio
from pde_utils.functions import search

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pan"

PAN_GRAPH_EDGEWEIGHTS = ["CentroidIdentity", "DBSeparation", "MinIdentity"]
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
    PAN_GRAPH_OUT_HELP = """
        PAN export option that selects the output graph dump file type.
            Follow selection argument with a valid graph file type.
    """

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
                        filters="", groups=[],
                        db_name=None, number_threads=1,
                        DB_stiffness=0.5, avg_identity=20, min_identity=35,
                        min_percent_gaps=50, PANgraph_out=None)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_build_pan(alchemist, hhdb_path, pan_name=None,
                      folder_path=None, folder_name=DEFAULT_FOLDER_NAME,
                      values=None, verbose=False, filters="", groups=[],
                      threads=1, M=50, aI=20, mI=35, B=0.5,
                      PANgraph_out=None):
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
        pan_alchemist.session.expire_on_commit = False

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

        hmm_data_dir = mapped_path.joinpath("pham_hhrs")
        hmm_data_dir.mkdir()
        build_pan_towns(alchemist, pan_alchemist, hhdb_path, pan_dict,
                        hmm_data_dir, data_maps_tuple,
                        threads=threads, verbose=verbose)

        if PANgraph_out is not None:
            pan_graph = pan_handling.to_networkx(pan_alchemist)
            pde_fileio.write_graph(pan_graph, PANgraph_out,
                                   mapped_path, pan_name,
                                   edge_weights=PAN_GRAPH_EDGEWEIGHTS)


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

    clusters = list(pan_dict.values())
    pan_alchemist.session.add_all(clusters)
    pan_alchemist.session.commit()
    pan_alchemist.session.expunge_all()
    pan_alchemist.session.close()

    return pan_dict


def build_pan_neighborhoods(alchemist, pan_alchemist, pan_dict, data_dir,
                            data_maps_tuple, aI=20, mI=35, B=0.5,
                            threads=1, verbose=False):
    cent_graph = create_centroid_graph(alchemist.engine, pan_dict, data_dir,
                                       threads=threads, verbose=verbose)

    if verbose:
        print("...Constructing base for pham neighborhoods...")

    thread_manager = multiprocessing.Manager()
    work_items = []
    work_locks = {}
    min_identities = thread_manager.dict()
    identity_edges = []
    for pham, target_cluster in pan_dict.items():
        for neighbor, edge_weights in cent_graph[pham].items():
            if pham == neighbor:
                continue

            cent_edge = cent_graph[pham][neighbor]

            t_spread = target_cluster.Spread

            local_cluster = pan_dict[neighbor]
            identity = correct_identity(cent_edge["CentroidIdentity"],
                                        len(local_cluster.CentroidSeq),
                                        len(target_cluster.CentroidSeq))

            if edge_weights["CentroidIdentity"] < aI:
                continue

            if identity > 80:
                print(f"......{pham} and {neighbor} seem too friendly...")

            DBsep = (t_spread/(100 - identity))
            cent_edge["DBSeparation"] = DBsep
            if DBsep >= B:
                local_lock = work_locks.get(neighbor)
                if local_lock is None:
                    local_lock = thread_manager.Lock()
                    work_locks[neighbor] = local_lock

                target_lock = work_locks.get(pham)
                if target_lock is None:
                    target_lock = thread_manager.Lock()
                    work_locks[pham] = target_lock

                work_items.append((
                                local_cluster, target_cluster, min_identities,
                                identity, DBsep, mI, data_dir,
                                local_lock, target_lock))

    if verbose:
        print("...Computing pham cluster minimum pairwise distances...")
    random.shuffle(work_items)
    identity_edges = parallelize.parallelize(
                            work_items, threads, test_identity_edge_process,
                            verbose=verbose)
    identity_edges = [edge for edge in identity_edges if edge is not None]

    if verbose:
        print("...Writing neighborhood data to PAN...")
    pan_alchemist.session.add_all(identity_edges)
    pan_alchemist.session.commit()
    pan_alchemist.session.expunge_all()
    pan_alchemist.session.close()


def build_pan_towns(alchemist, pan_alchemist, hhdb_path, pan_dict,
                    hmm_data_dir, data_maps_tuple, threads=1, verbose=False):
    work_items = []
    hhr_path_map = {}
    for pham, hmm_path in data_maps_tuple[3].items():
        hhr_path = hmm_data_dir.joinpath(".".join([str(pham), "hhr"]))
        hhr_path_map[pham] = hhr_path

        work_items.append((hmm_path, hhdb_path, hhr_path))

    if verbose:
        print("...Performing iterations of hhblitz to find HMM-HMM "
              "relationships...")
    parallelize.parallelize(work_items, threads, search.hhblitz,
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
    cent_graph = Graph()

    cent_gs_to_ts = {}

    phams = list(pan_dict.keys())
    for pham in phams:
        cluster = pan_dict[pham]
        cent_gs_to_ts[pham] = cluster.CentroidSeq.decode("utf-8")

        cent_graph.add_node(pham, CentroidID=cluster.CentroidID,
                            CentroidSeq=cluster.CentroidSeq.decode("utf-8"),
                            Spread=cluster.Spread)

    if verbose:
        print("...Writing pham cluster centroids to fasta...")
    cent_fasta_path = aln_dir.joinpath("centroids.fasta")
    pdm_fileio.write_fasta(cent_gs_to_ts, cent_fasta_path)

    if verbose:
        print("...mBedding pham centroids for k-tuple distances...")
    distmat_path = aln_dir.joinpath("centroids.mat")
    alignment.mBed(cent_fasta_path, distmat_path)

    cmat = clustal.PercentIdentityMatrix(distmat_path)
    cmat.parse_matrix(file_type="mbed")

    for source in phams:
        for target in phams:
            distance = cmat.get_identity(str(source), str(target)) * 100
            identity = 100 - distance
            cent_graph.add_edge(source, target, CentroidIdentity=identity)

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


def test_identity_edge_process(local_cluster, target_cluster, min_identities,
                               avg_identity, DBsep, mI, data_dir,
                               local_lock, target_lock):
    source = int(local_cluster.ClusterID)
    target = int(target_cluster.ClusterID)

    source_node = min_identities.get(source)
    target_node = min_identities.get(target)

    if source_node is not None:
        min_identity = source_node.get(target)
        if min_identity is None:
            local_lock.acquire()
            target_lock.acquire()

            min_identity = estimate_min_identity(
                                local_cluster, target_cluster, data_dir)

            local_lock.release()
            target_lock.release()

            min_identities[source] = {target: min_identity}
    elif target_node is not None:
        min_identity = target_node.get(target)
        if min_identity is None:
            local_lock.acquire()
            target_lock.acquire()

            min_identity = estimate_min_identity(
                                local_cluster, target_cluster, data_dir)

            local_lock.release()
            target_lock.release()

            min_identities[source] = {target: min_identity}
    else:
        local_lock.acquire()
        target_lock.acquire()

        min_identity = estimate_min_identity(
                                local_cluster, target_cluster, data_dir)

        local_lock.release()
        target_lock.release()

        min_identities[source] = {target: min_identity}

    if min_identity >= mI:
        return pan_models.IdentityEdge(
                                Source=source, Target=target,
                                CentroidIdentity=avg_identity,
                                DBSeparation=DBsep, MinIdentity=min_identity)


def estimate_min_identity(source_cluster, local_cluster, data_dir):
    if source_cluster.AlignmentGuidetree is not None:
        source_leaf = estimate_linker_sequence(
                                    source_cluster, local_cluster.CentroidID,
                                    local_cluster.CentroidSeq.decode("utf-8"),
                                    data_dir)

        source_linker = source_leaf.name
        source_linker_seq = source_leaf.comment["Sequence"]
    else:
        source_linker = source_cluster.CentroidID
        source_linker_seq = source_cluster.CentroidSeq.decode("utf-8")

    source_len = len(source_linker_seq)
    source_path = data_dir.joinpath(".".join([source_linker, "fasta"]))
    if not source_path.is_file():
        pdm_fileio.write_fasta({source_linker: source_linker_seq}, source_path)

    if local_cluster.AlignmentGuidetree is not None:
        target_leaf = estimate_linker_sequence(
                                    local_cluster, source_cluster.CentroidID,
                                    source_cluster.CentroidSeq.decode("utf-8"),
                                    data_dir)

        target_linker = target_leaf.name
        target_linker_seq = target_leaf.comment["Sequence"]
    else:
        target_linker = local_cluster.CentroidID
        target_linker_seq = local_cluster.CentroidSeq.decode("utf-8")

    target_len = len(target_linker_seq)
    target_path = data_dir.joinpath(".".join([target_linker, "fasta"]))
    if not target_path.is_file():
        pdm_fileio.write_fasta({target_linker: target_linker_seq}, target_path)

    pairwise_path = data_dir.joinpath("".join([source_linker, "__",
                                               target_linker, "fasta"]))

    linker_alignment = alignment.pairwise_align(
                                            source_path, target_path,
                                            pairwise_path, tool="needle")

    if source_len > target_len:
        linker_pid = linker_alignment.annotations["identity"] / source_len
    else:
        linker_pid = linker_alignment.annotations["identity"] / target_len

    return linker_pid * 100


def estimate_linker_sequence(source_cluster, target_centroid_id,
                             target_centroid_seq, data_dir):
    guidetree = source_cluster.AlignmentGuidetree

    centroid_path = data_dir.joinpath(".".join([target_centroid_id, "fasta"]))
    if not centroid_path.is_file():
        pdm_fileio.write_fasta({target_centroid_id: target_centroid_seq},
                               centroid_path)

    curr_node = guidetree.clade
    while curr_node.is_bifurcating() and len(curr_node.clades) > 1:
        left_child = curr_node.clades[0]
        right_child = curr_node.clades[1]

        left_leaf = get_furthest_sequence(source_cluster, left_child)
        left_seq_path = data_dir.joinpath(".".join(
                                                [left_leaf.name, "fasta"]))
        if not left_seq_path.is_file():
            pdm_fileio.write_fasta({left_leaf.name:
                                    left_leaf.comment["Sequence"]},
                                   left_seq_path)

        right_leaf = get_furthest_sequence(source_cluster, right_child)
        right_seq_path = data_dir.joinpath(".".join(
                                                 [right_leaf.name, "fasta"]))
        if not right_seq_path.is_file():
            pdm_fileio.write_fasta({right_leaf.name:
                                    right_leaf.comment["Sequence"]},
                                   right_seq_path)

        left_alignment = left_leaf.comment.get("TargetAlignment")
        if left_alignment is None:
            pairwise_path = data_dir.joinpath("".join(
                                                [left_leaf.name, "__",
                                                 target_centroid_id, ".aln"]))
            left_alignment = alignment.pairwise_align(
                                            left_seq_path, centroid_path,
                                            pairwise_path, tool="needle")
            left_leaf.comment["TargetAlignment"] = left_alignment

        right_alignment = right_leaf.comment.get("TargetAlignment")
        if right_alignment is None:
            pairwise_path = data_dir.joinpath("".join(
                                                [right_leaf.name, "__",
                                                 target_centroid_id, ".aln"]))
            right_alignment = alignment.pairwise_align(
                                            right_seq_path, centroid_path,
                                            pairwise_path, tool="needle")
            right_leaf.comment["TargetAlignment"] = right_alignment

        left_pid = (left_alignment.annotations["identity"] /
                    left_alignment.get_alignment_length())

        right_pid = (right_alignment.annotations["identity"] /
                     right_alignment.get_alignment_length())

        if left_pid > right_pid:
            curr_node = left_child
        else:
            curr_node = right_child

    linker = get_furthest_sequence(source_cluster, curr_node)
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
            distance = cluster.PercentIdentityMatrix.get_distance(centroid_id,
                                                                  seq_id)
            leaf.comment["SourceIdentity"] = distance

        if furthest_leaf is None:
            furthest_leaf = leaf
            max_distance = distance
        elif distance > max_distance:
            furthest_leaf = leaf
            max_distance = distance

    sequence = furthest_leaf.comment.get("Sequence")
    if sequence is None:
        furthest_leaf.comment["Sequence"] = cluster.MultipleSequenceAlignment\
                                                   .get_sequence(
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
