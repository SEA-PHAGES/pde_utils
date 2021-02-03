import argparse
import decimal
import math
from pathlib import Path
import shutil
import string
import textwrap
import time

from pdm_utils.functions import (configfile, multithread, parallelize,
                                 fileio as pdm_fileio, pipelines_basic)
from pdm_utils.pipelines.revise import TICKET_HEADER

from pde_utils.functions import (alignment, clustering, seq_distance,
                                 sql_queries)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_cluster_db")
TEMP_DIR = (f"/tmp/{DEFAULT_FOLDER_NAME}_temp")

CLUSTER_DB_SUBPIPELINES = ["analyze", "cluster"]
DEFAULT_SETTINGS = {"kmer": 15, "sketch": 25000, "gcs": 0.35, "ani": 0.7,
                    "gcsmax": 0.90, "animax": 0.95, "gcsS": 0.8, "gcsM": 2,
                    "aniS": 0, "aniM": 1}

CLUSTER_ANALYSIS_HEADER = ["Subject PhageID", "Subject Cluster",
                           "Query PhageID", "Query Cluster", "GCS",
                           "Estimated ANI"]


def main(unparsed_args_list):
    args = parse_cluster_db(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)
    values = pipelines_basic.parse_value_input(args.input)

    execute_cluster_db(alchemist,
                       folder_path=args.folder_path,
                       folder_name=args.folder_name,
                       values=values, verbose=args.verbose,
                       filters=args.filters, groups=args.groups,
                       threads=args.number_threads, kmer=args.kmer_size,
                       sketch=args.sketch_size,
                       gcs=args.gene_content_similarity_min,
                       ani=args.average_nucleotide_identity_min,
                       gcsmax=args.gene_content_similarity_max,
                       animax=args.average_nucleotide_identity_max,
                       gcsS=args.gcsS, gcsM=args.gcsM, aniS=args.aniS,
                       aniM=args.aniM, evaluate=args.dump_evaluation,
                       mat_out=args.distmat_out, subcluster=args.subcluster,
                       cluster_prefix=args.cluster_prefix)


def parse_cluster_db(unparsed_args_list):
    DATABASE_HELP = """
        Name of the MySQL database to perform clustering analysis on.
        """
    VERBOSE_HELP = """
        Cluster DB option that enables progress print statements.
        """
    FOLDER_PATH_HELP = """
        Cluster DB option to change the path of the directory where the
        exported files are stored.
            Follow selection argument with the path to the desired
            export directory.
        """
    FOLDER_NAME_HELP = """
        Cluster DB option to change the name of the directory where
        the exported files are stored.
            Follow selection argument with the desired name.
        """
    CONFIG_FILE_HELP = """
        Cluster DB option that enables use of a config file for sourcing
        credentials and parameters
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

    GENE_CONTENT_SIMILARITY_MIN_HELP = """
        Cluster DB option to change the lower threshold for analyses of
        phage genome gene content similarity.
            Follow selection argument with the desired percentage [0-1]
        """
    AVERAGE_NUCLEOTIDE_IDENTITY_MIN_HELP = """
        Cluster DB option to change the lower threshold for analyses of
        phage genome average nucleotide identity.
            Follow selection argument with the desired percentage [0-1]
        """
    GENE_CONTENT_SIMILARITY_MAX_HELP = """
        Cluster DB option to change the upper threshold for analyses of
        phage genome gene content similarity.
            Follow selection argument with the desired percentage [0-1]
        """
    AVERAGE_NUCLEOTIDE_IDENTITY_MAX_HELP = """
        Cluster DB option to change the upper threshold for analyses of
        phage genome average nucleotide identity.
            Follow selection argument with the desired percentage [0-1]
        """
    KMER_SIZE_HELP = """
        Cluster DB option to change the k-mer size that determines
        the building block of the hashes used to estimate ANI.
            Follow selection argument with the desired k-mer size.
        """
    SKETCH_SIZE_HELP = """
        Cluster DB option to set the minimum amount of non-redundant
        min-hashes used to calculate ANI.
            Follow selection argument with the desired sketch size.
        """
    GENE_CONTENT_SIMILARITY_EPS_MODIFIER_HELP = """
        """
    GENE_CONTENT_SIMILARITY_MINPTS_MODIFIER_HELP = """
        """
    AVERAGE_NUCLEOTIDE_IDENTITY_EPS_MODIFIER_HELP = """
        """
    AVERAGE_NUCLEOTIDE_IDENTITY_MINPTS_MODIFIER_HELP = """
        """

    parser = argparse.ArgumentParser()

    parser.add_argument("database", type=str, help=DATABASE_HELP)
    parser.add_argument("-c", "--config_file", help=CONFIG_FILE_HELP,
                        type=pipelines_basic.convert_file_path)
    parser.add_argument("-m", "--folder_name", type=str,
                        help=FOLDER_NAME_HELP)
    parser.add_argument("-o", "--folder_path", type=Path,
                        help=FOLDER_PATH_HELP)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help=VERBOSE_HELP)

    parser.add_argument("-if", "--import_file",
                        type=pipelines_basic.convert_file_path,
                        help=IMPORT_FILE_HELP, dest="input")
    parser.add_argument("-in", "--import_names", nargs="*",
                        help=IMPORT_NAMES_HELP, dest="input")

    parser.add_argument("-w", "--where", nargs="?",
                        help=WHERE_HELP, dest="filters")
    parser.add_argument("-g", "--group_by", nargs="+",
                        help=GROUP_BY_HELP, dest="groups")
    parser.add_argument("-np", "--number_threads", type=int, nargs="?",
                        help=NUMBER_THREADS_HELP)

    parser.add_argument("-kmer", "--kmer_size", type=int,
                        help=KMER_SIZE_HELP)
    parser.add_argument("-sketch", "--sketch_size", type=int,
                        help=SKETCH_SIZE_HELP)

    parser.add_argument("-gcs", "--gene_content_similarity_min",
                        type=float,
                        help=GENE_CONTENT_SIMILARITY_MIN_HELP)
    parser.add_argument("-ani", "--average_nucleotide_identity_min",
                        type=float,
                        help=AVERAGE_NUCLEOTIDE_IDENTITY_MIN_HELP)
    parser.add_argument("-gcsmax", "--gene_content_similarity_max",
                        type=float,
                        help=GENE_CONTENT_SIMILARITY_MAX_HELP)
    parser.add_argument("-animax", "--average_nucleotide_identity_max",
                        type=float,
                        help=AVERAGE_NUCLEOTIDE_IDENTITY_MAX_HELP)

    parser.add_argument("-gcsS", type=float,
                        help=GENE_CONTENT_SIMILARITY_EPS_MODIFIER_HELP)
    parser.add_argument("-gcsM", type=float,
                        help=GENE_CONTENT_SIMILARITY_MINPTS_MODIFIER_HELP)
    parser.add_argument("-aniS", type=float,
                        help=AVERAGE_NUCLEOTIDE_IDENTITY_EPS_MODIFIER_HELP)
    parser.add_argument("-aniM", type=float,
                        help=AVERAGE_NUCLEOTIDE_IDENTITY_MINPTS_MODIFIER_HELP)

    parser.add_argument("-mat", "--distmat_out", action="store_true")
    parser.add_argument("-eval", "--dump_evaluation", action="store_true")
    parser.add_argument("-cpre", "--cluster_prefix", type=str)
    parser.add_argument("-sub", "--subcluster", action="store_true")

    parser.set_defaults(
                    folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                    config_file=None, verbose=False, input=[],
                    filters="", groups=[], number_threads=1,
                    kmer_size=DEFAULT_SETTINGS["kmer"],
                    sketch_size=DEFAULT_SETTINGS["sketch"],
                    gene_content_similarity_min=DEFAULT_SETTINGS["gcs"],
                    average_nucleotide_identity_min=DEFAULT_SETTINGS["ani"],
                    gene_content_similarity_max=DEFAULT_SETTINGS["gcsmax"],
                    average_nucleotide_identity_max=DEFAULT_SETTINGS["animax"],
                    gcsS=DEFAULT_SETTINGS["gcsS"],
                    gcsM=DEFAULT_SETTINGS["gcsM"],
                    aniS=DEFAULT_SETTINGS["aniS"],
                    aniM=DEFAULT_SETTINGS["aniM"],
                    cluster_prefix=None)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_cluster_db(
                alchemist, folder_path=None,
                folder_name=DEFAULT_FOLDER_NAME, values=None, verbose=None,
                filters="", groups=[], threads=1,
                kmer=DEFAULT_SETTINGS["kmer"],
                sketch=DEFAULT_SETTINGS["sketch"],
                gcs=DEFAULT_SETTINGS["gcs"], ani=DEFAULT_SETTINGS["ani"],
                gcsmax=DEFAULT_SETTINGS["gcsmax"],
                animax=DEFAULT_SETTINGS["animax"],
                gcsS=DEFAULT_SETTINGS["gcsS"], gcsM=DEFAULT_SETTINGS["gcsM"],
                aniS=DEFAULT_SETTINGS["aniS"], aniM=DEFAULT_SETTINGS["aniM"],
                mat_out=False, evaluate=False, subcluster=False,
                cluster_prefix=None):
    db_filter = pipelines_basic.build_filter(alchemist, "phage", filters,
                                             values=values)

    working_path = pipelines_basic.create_working_path(
                                            folder_path, folder_name)
    temp_dir = create_temp_path(TEMP_DIR)
    conditionals_map = pipelines_basic.build_groups_map(
                                            db_filter, working_path,
                                            groups=groups, verbose=verbose)

    values = db_filter.values
    for mapped_path in conditionals_map.keys():
        db_filter.reset()
        db_filter.values = values

        conditionals = conditionals_map[mapped_path]
        db_filter.values = db_filter.build_values(where=conditionals)

        if verbose:
            print("Querying MySQL database for clustering metadata...")
        cluster_metadata = query_cluster_metadata(db_filter)

        gcs_matrix = calculate_gcs_matrix(alchemist, db_filter.values,
                                          verbose=verbose, cores=threads)

        pipelines_basic.create_working_dir(mapped_path)

        if verbose:
            print("Clustering database genomes...")
        cluster_scheme = gcs_cluster(
                                mapped_path, gcs_matrix,
                                cluster_metadata[0], cluster_metadata[1],
                                gcs=gcs, gcsmax=gcsmax, S=gcsS, M=gcsM,
                                evaluate=evaluate, cores=threads,
                                verbose=verbose,
                                cluster_prefix=cluster_prefix)

        if subcluster:
            sketch_path_map = sketch_genomes(db_filter, temp_dir,
                                             verbose=verbose)

            if verbose:
                print("Subclustering database genomes...")
            ani_subcluster(mapped_path, sketch_path_map, cluster_scheme,
                           cluster_metadata[0], cluster_metadata[1],
                           cluster_metadata[2], cores=threads,
                           verbose=verbose, ani=ani, animax=animax,
                           evaluate=evaluate)

            empty = True
            for _ in mapped_path.iterdir():
                empty = False

            if empty:
                shutil.rmtree(mapped_path)


def query_cluster_metadata(db_filter):
    cluster_data = db_filter.retrieve(["phage.Cluster",
                                       "phage.Subcluster"])
    cluster_lookup = {}
    subcluster_lookup = {}
    for phage_id, data_dict in cluster_data.items():
        cluster_lookup[phage_id] = data_dict["Cluster"][0]
        subcluster_lookup[phage_id] = data_dict["Subcluster"][0]

    seqid_cluster_map = db_filter.group("phage.Cluster")

    return (cluster_lookup, seqid_cluster_map,
            subcluster_lookup)


# CLUSTERING HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def calculate_gcs_matrix(alchemist, phage_ids, cores=1, verbose=False):
    if verbose:
        print("Querying MySQL database for genome gene content...")
    phage_gc_nodes = []
    for phage in phage_ids:
        phage_gc_nodes.append(
                    sql_queries.get_distinct_phams_from_organism(
                                                    alchemist, phage))

    if verbose:
        print("Calculating gene content similarity matrix...")
    gcs_matrix = clustering.build_symmetric_matrix(
                                phage_gc_nodes,
                                seq_distance.calculate_gcs, names=phage_ids,
                                cores=cores, verbose=verbose)
    return gcs_matrix


def calculate_ani_matrix(phage_ids, sketch_path_map, cores=1, verbose=False):
    sketch_paths = []
    for phage_id in phage_ids:
        sketch_paths.append(sketch_path_map[phage_id])

    ani_matrix = clustering.build_symmetric_matrix(
                                sketch_paths, calculate_ani,
                                names=phage_ids, cores=cores, verbose=verbose)

    return ani_matrix


def sketch_genomes(db_filter, working_dir, verbose=False, threads=1,
                   kmer=DEFAULT_SETTINGS["kmer"],
                   sketch=DEFAULT_SETTINGS["sketch"]):
    gs_and_ts = db_filter.select(["phage.PhageID", "phage.Sequence"],
                                 return_dict=False)

    fasta_dir = create_temp_path(str(working_dir.joinpath("fasta")))
    fasta_path_map = write_genome_fastas(
                                 gs_and_ts, fasta_dir, verbose=verbose,
                                 threads=threads)

    sketch_dir = create_temp_path(str(working_dir.joinpath("sketches")))
    sketch_path_map = sketch_genome_fastas(
                                      fasta_path_map, sketch_dir,
                                      verbose=verbose, threads=threads,
                                      kmer=kmer, sketch=sketch)

    return sketch_path_map


def cluster_db(matrix, eps, cores=1, verbose=False, is_distance=False,
               emax=0.9, S=1.6, M=2):
    if verbose:
        print("...Greedily defining clusters...")
    greedy_scheme = clustering.dbscan(matrix, eps, 1,
                                      is_distance=is_distance,
                                      return_matrix=True)

    cluster_counter = 0
    second_scheme = dict()
    work_items = []
    for greedy_cluster, submatrix in greedy_scheme.items():
        if greedy_cluster is None or submatrix.size <= 1:
            noise = second_scheme.get(None, list())
            second_scheme[None] = noise + submatrix.labels
            continue

        work_items.append((submatrix, is_distance, eps, emax, S, M))

    if verbose:
        print("...Performing clustering iterations...")
    layered_schemes = parallelize.parallelize(
                            work_items, cores, iter_cluster_process,
                            verbose=verbose)

    iter_scheme = dict()
    work_items = []
    for scheme in layered_schemes:
        for cluster, cluster_members in scheme.items():
            if cluster is None:
                noise = second_scheme.get(None, list())
                iter_scheme[None] = noise + cluster_members
                continue

            cluster_counter += 1
            iter_scheme[cluster_counter] = cluster_members

    return iter_scheme


def iter_cluster_process(submatrix, is_distance, emin, emax, S, M):
    mean = submatrix.get_mean()
    std_dev = submatrix.get_SD()

    standard_submatrix = submatrix.standardize(metric="z_score")
    std_mean = standard_submatrix.get_mean()
    std_median = standard_submatrix.get_median()

    eps = mean + ((std_mean - std_median) * S * std_dev)
    if eps < emin:
        eps = emin
    elif eps > emax:
        eps = emax

    size = submatrix.size
    minpts = int(round(((math.log2(size) - 1) * M), 0))
    if minpts < 1:
        minpts = 1

    dbscan_scheme = clustering.dbscan(submatrix, eps,
                                      minpts, is_distance=is_distance,
                                      return_matrix=True)

    dbscan_centroids = list()
    noise = list()
    for dbscan_cluster, dbscan_matrix in dbscan_scheme.items():
        if dbscan_cluster is None:
            noise += dbscan_matrix.labels
            continue

        dbscan_centroids.append(dbscan_matrix.get_centroid())

    if not dbscan_centroids:
        return {None: noise}

    kmeans_scheme = clustering.lloyds(submatrix, dbscan_centroids,
                                      is_distance=is_distance,
                                      return_matrix=False)

    iter_scheme = dict()
    for kmeans_cluster, kmeans_members in kmeans_scheme.items():
        if len(kmeans_members) <= 1:
            noise += kmeans_members
            continue

        iter_scheme[kmeans_cluster] = kmeans_members

    iter_scheme[None] = noise

    return iter_scheme


def gcs_cluster(working_dir, gcs_matrix, cluster_lookup, cluster_seqid_map,
                cores=1, verbose=False, gcs=DEFAULT_SETTINGS["gcs"],
                gcsmax=DEFAULT_SETTINGS["gcsmax"], S=DEFAULT_SETTINGS["gcsS"],
                M=DEFAULT_SETTINGS["gcsM"], evaluate=False,
                cluster_prefix=None):
    cluster_scheme = cluster_db(gcs_matrix, gcs, emax=gcsmax, S=S, M=M,
                                cores=cores, verbose=verbose,
                                is_distance=False)

    # TEST CODE
    # =========================================================================
    # new_cluster_lookup = {}
    # for cluster, cluster_members in cluster_scheme.items():
    #    for member in cluster_members:
    #        new_cluster_lookup[member] = cluster

    # for cluster in cluster_seqid_map.keys():
    #    print(f"Old cluster {cluster}")

    #    new_cluster_set = set()
    #    for cluster_member in cluster_seqid_map[cluster]:
    #        new_cluster_set.add(new_cluster_lookup[cluster_member])

    #    print(f"\tEnded up in {new_cluster_set}")

    # shutil.rmtree(working_dir)
    # return
    # =========================================================================
    # TEST CODE

    old_clusters = list(cluster_seqid_map.keys())
    cluster_redistributions = get_cluster_redistributions(
                                cluster_scheme, cluster_lookup, old_clusters)

    cluster_scheme = assign_cluster_names(
                                cluster_scheme, cluster_redistributions,
                                verbose=verbose, cluster_prefix=cluster_prefix)

    # TEST CODE
    # =========================================================================
    # new_cluster_lookup = {}
    # for cluster, cluster_members in cluster_scheme.items():
    #    for member in cluster_members:
    #        new_cluster_lookup[member] = cluster

    # stable_old_clusters = set(cluster_seqid_map.keys())
    # new_clusters = set(cluster_scheme.keys())
    # diff_clusters = stable_old_clusters.difference(new_clusters)

    # for cluster in diff_clusters:
    #     print(f"Neglected cluster {cluster}")

    #    new_cluster_set = set()
    #    for cluster_member in cluster_seqid_map[cluster]:
    #        new_cluster_set.add(new_cluster_lookup[cluster_member])

    #    print(f"\tEnded up in {new_cluster_set}")

    # shutil.rmtree(working_dir)
    # return
    # =========================================================================
    # TEST CODE

    scheme_alterations = diff_cluster_schemes(cluster_scheme, cluster_lookup)

    wrote_eval = False
    if evaluate:
        new_matrix_cache = dict()

        if verbose:
            print("Evaluating clustering scheme...")
        scheme_metadata = evaluate_clustering_scheme(
                            gcs_matrix, cluster_scheme,
                            matrix_cache=new_matrix_cache,
                            cores=cores, verbose=verbose)

        old_scheme_metadata = evaluate_clustering_scheme(
                            gcs_matrix, cluster_seqid_map,
                            cores=cores)

        alteration_metadata = evaluate_scheme_alteration(
                            gcs_matrix, cluster_scheme, cluster_seqid_map,
                            scheme_alterations)

        wrote_eval = write_clustering_evaluation(
                                    working_dir, scheme_metadata,
                                    old_scheme_metadata, alteration_metadata,
                                    "gene content similarity")

    wrote_ticket = write_clustering_update_ticket(
                                     working_dir, scheme_alterations,
                                     field="Cluster")

    wrote = wrote_eval or wrote_ticket
    if not wrote:
        if verbose:
            print("No changes made to current clustering scheme.")
        shutil.rmtree(working_dir)

    return cluster_scheme


def ani_subcluster(working_dir, sketch_path_map, cluster_scheme,
                   cluster_lookup, cluster_seqid_map,
                   subcluster_lookup, cores=1, verbose=False,
                   ani=DEFAULT_SETTINGS["ani"],
                   animax=DEFAULT_SETTINGS["animax"], evaluate=False):

    for cluster, cluster_members in cluster_scheme.items():
        if cluster is None:
            continue

        cluster_members_set = set(cluster_members)
        old_cluster_members = set(cluster_seqid_map.get(cluster, list()))
        noncluster_members = list(cluster_members_set.difference(
                                                        old_cluster_members))

        old_cluster_members = list(old_cluster_members.intersection(
                                                    cluster_members_set))
        noncluster_members = list(noncluster_members)

        old_subclusters = set()
        subcluster_seqid_map = {}
        altered_subcluster_lookup = {}
        for member in old_cluster_members:
            subcluster = subcluster_lookup[member]
            altered_subcluster_lookup[member] = subcluster
            old_subclusters.add(subcluster)

            seqids = subcluster_seqid_map.get(subcluster, list())
            seqids.append(member)
            subcluster_seqid_map[subcluster] = seqids

        for nonmember in noncluster_members:
            nonmember_cluster = cluster_lookup[nonmember]
            altered_subcluster_lookup[nonmember] = nonmember_cluster

            seqids = subcluster_seqid_map.get(None, list())
            seqids.append(nonmember)
            subcluster_seqid_map[None] = seqids

        if verbose:
            print(f"Subclustering {cluster}...")

        ani_matrix = calculate_ani_matrix(cluster_members, sketch_path_map,
                                          cores=cores, verbose=verbose)

        subcluster_scheme = cluster_db(ani_matrix, ani, emax=animax,
                                       cores=cores, verbose=verbose,
                                       is_distance=False)

        subcluster_redistributions = get_cluster_redistributions(
                                            subcluster_scheme,
                                            altered_subcluster_lookup,
                                            old_subclusters)

        subcluster_scheme = assign_cluster_names(
                                            subcluster_scheme,
                                            subcluster_redistributions,
                                            verbose=verbose,
                                            subcluster=cluster)

        scheme_alterations = diff_cluster_schemes(subcluster_scheme,
                                                  altered_subcluster_lookup)

        subcluster_dir = working_dir.joinpath(str(cluster))
        pipelines_basic.create_working_dir(subcluster_dir)

        wrote_ticket = write_clustering_update_ticket(
                                subcluster_dir, scheme_alterations,
                                field="Subcluster")

        wrote_eval = False
        if evaluate:
            new_matrix_cache = dict()

            if verbose:
                print(f"...Evaluating cluster {cluster} "
                      "subclustering scheme...")
            scheme_metadata = evaluate_clustering_scheme(
                                ani_matrix, subcluster_scheme,
                                matrix_cache=new_matrix_cache,
                                cores=cores, verbose=verbose)

            old_scheme_metadata = evaluate_clustering_scheme(
                                ani_matrix, subcluster_seqid_map,
                                cores=cores)

            alteration_metadata = evaluate_scheme_alteration(
                                ani_matrix, subcluster_scheme,
                                subcluster_seqid_map, scheme_alterations)

            write_clustering_evaluation(
                                subcluster_dir, scheme_metadata,
                                old_scheme_metadata, alteration_metadata,
                                "average nucleotide identity")

        wrote = wrote_eval or wrote_ticket

        if not wrote:
            shutil.rmtree(subcluster_dir)


def calculate_ani(subject_path, query_path):
    mash_output = alignment.mash_dist(subject_path, query_path)
    ani_data = mash_output.split("\t")

    if len(ani_data) == 5:
        return 1 - float(ani_data[2])


# CLUSTER-NAMING HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def assign_cluster_names(cluster_scheme, cluster_redistributions,
                         verbose=False, cluster_prefix=None, subcluster=None):
    named_scheme = dict()
    named_scheme[None] = cluster_scheme.pop(None, list())

    old_clusters = list(cluster_redistributions.keys())

    remaining = list([x for x in cluster_scheme.keys() if x is not None])
    assigned = list()

    for old_cluster, cluster_redistribution in cluster_redistributions.items():
        if old_cluster is None:
            continue

        old_cluster_inheretors = sorted(cluster_redistribution.items(),
                                        key=lambda x: len(x[1]), reverse=True)

        for old_cluster_inheretor in old_cluster_inheretors:
            num_cluster = old_cluster_inheretor[0]
            if num_cluster is None:
                continue

            if num_cluster in assigned:
                continue

            named_scheme[old_cluster] = cluster_scheme[num_cluster]
            remaining.remove(num_cluster)
            assigned.append(num_cluster)
            break

    for num_cluster in remaining:
        if subcluster is not None:
            if len(cluster_scheme) <= 1:
                named_scheme[None] = cluster_scheme[num_cluster]
                continue

            new_subcluster = gen_new_subcluster(subcluster, old_clusters)
            named_scheme[new_subcluster] = cluster_scheme[num_cluster]

            if verbose:
                print(f"......Created new subcluster '{new_subcluster}'...")

            continue

        new_cluster = gen_new_cluster(old_clusters,
                                      cluster_prefix=cluster_prefix)
        named_scheme[new_cluster] = cluster_scheme[num_cluster]

        old_clusters.append(new_cluster)
        assigned.append(num_cluster)

        if verbose:
            print(f"...Created new cluster '{new_cluster}'...")

    return named_scheme


def gen_new_cluster(old_clusters, cluster_prefix=None):
    if cluster_prefix is None:
        cluster_prefix = ""
    alphabet = string.ascii_uppercase

    chars = 1
    counter = []
    while True:
        filled = True
        if chars > len(counter):
            for _ in range(chars - len(counter)):
                counter.append(0)

        letters = [cluster_prefix]
        for char_num in counter:
            letters.append(alphabet[char_num])

        new_cluster = "".join(letters)
        if new_cluster not in old_clusters:
            return new_cluster

        for i in range(len(counter)):
            if counter[len(counter)-1-i] < (len(alphabet) - 1):
                counter[len(counter)-1-i] += 1
                filled = False
                break

            counter[len(counter)-1-i] = 0

        if filled:
            chars += 1
            if chars >= 5:
                break


def gen_new_subcluster(cluster, old_subclusters):
    for i in range(1, 1000):
        subcluster = "".join([cluster, str(i)])
        if subcluster not in old_subclusters:
            return subcluster


# CLUSTERING EVALUATION HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def evaluate_clustering_scheme(matrix, cluster_scheme, cores=1, verbose=False,
                               matrix_cache=None):
    if matrix_cache is None:
        matrix_cache = dict()

    work_items = []
    for cluster, cluster_members in cluster_scheme.items():
        if cluster is None:
            continue

        cluster_matrix = matrix_cache.get(cluster)

        if cluster_matrix is None:
            cluster_matrix = matrix.get_submatrix_from_labels(
                                            cluster_scheme[cluster])

        work_items.append((cluster, cluster_matrix))

    evaluations = parallelize.parallelize(work_items, cores,
                                          cluster_evaluation_subprocess,
                                          verbose=verbose)

    return {data[0]: data[1] for data in evaluations}


def cluster_evaluation_subprocess(cluster, matrix):
    data = dict()

    data["centroid"] = matrix.get_centroid()
    data["spread"] = matrix.get_average_edge(data["centroid"])
    data["average_value"] = matrix.get_mean()
    data["standard_deviation"] = matrix.get_SD()
    data["num_members"] = matrix.size

    return cluster, data


def evaluate_scheme_alteration(matrix, cluster_scheme, old_cluster_scheme,
                               scheme_alterations):
    scheme_alteration_metadata = dict()
    for cluster, cluster_alteration_data in scheme_alterations.items():

        altered_cluster_metadata = dict()
        for change in cluster_alteration_data:
            old_cluster_data = altered_cluster_metadata.get(
                                                    change["old_cluster"],
                                                    dict())

            alterations = old_cluster_data.get("alterations", list())
            alterations.append(change)
            old_cluster_data["alterations"] = alterations

            members = old_cluster_data.get("old_members", list())
            members.append(change["id"])
            old_cluster_data["old_members"] = members

            altered_cluster_metadata[change["old_cluster"]] = old_cluster_data

        for old_cluster, old_cluster_data in altered_cluster_metadata.items():
            old_cluster_data["percent_merged"] = decimal.Decimal(
                                        len(old_cluster_data["old_members"]) /
                                        len(old_cluster_scheme[old_cluster]))

        scheme_alteration_metadata[cluster] = altered_cluster_metadata

    return scheme_alteration_metadata


def diff_cluster_schemes(cluster_scheme, cluster_lookup):
    scheme_alterations = dict()

    for new_cluster, cluster_members in cluster_scheme.items():
        for member in cluster_members:
            old_cluster = cluster_lookup[member]
            if new_cluster != old_cluster:
                altered_data = scheme_alterations.get(new_cluster, list())
                data = dict()
                data["id"] = member
                data["old_cluster"] = old_cluster
                data["new_cluster"] = new_cluster

                altered_data.append(data)
                scheme_alterations[new_cluster] = altered_data

    return scheme_alterations


def get_cluster_redistributions(cluster_scheme, cluster_lookup, old_clusters):
    cluster_redistributions = dict()

    for cluster_num, cluster_members in cluster_scheme.items():
        for member in cluster_members:
            old_cluster = cluster_lookup.get(member)

            cluster_redistribution = cluster_redistributions.get(
                                                        old_cluster, dict())

            distributed_members = cluster_redistribution.get(
                                                        cluster_num, list())
            distributed_members.append(member)

            cluster_redistribution[cluster_num] = distributed_members
            cluster_redistributions[old_cluster] = cluster_redistribution

    cluster_redistribution_weights = list()
    for old_cluster, cluster_redistribution in cluster_redistributions.items():
        weight = 0
        for new_cluster, cluster_members in cluster_redistribution.items():
            weight += len(cluster_members)

        cluster_redistribution_weights.append((old_cluster, weight))

    cluster_redistribution_weights.sort(key=lambda x: x[1], reverse=True)

    cluster_redistributions = {
                cluster: cluster_redistributions[cluster]
                for cluster, weight in cluster_redistribution_weights}

    return cluster_redistributions


# WRITING FUNCTIONS
# -----------------------------------------------------------------------------
def write_clustering_evaluation(working_dir, scheme_metadata,
                                old_scheme_metadata, alteration_metadata,
                                metric, max_qlines=5):
    full_data_lines = []
    quick_data_lines = []

    total_old_avg = 0
    total_old_clusters = 0
    for old_cluster, old_scheme_data in old_scheme_metadata.items():
        if old_cluster is None:
            continue

        total_old_clusters += 1
        total_old_avg += old_scheme_data["average_value"]

    if total_old_clusters > 0:
        total_old_avg /= total_old_clusters
    else:
        total_old_avg = 0

    total_new_avg = 0
    total_new_clusters = 0
    for new_cluster, new_scheme_data in scheme_metadata.items():
        if new_cluster is None:
            continue

        new_scheme_data = scheme_metadata.get(new_cluster)

        total_new_clusters += 1
        total_new_avg += new_scheme_data["average_value"]

    if total_new_clusters > 0:
        total_new_avg /= total_new_clusters
    else:
        total_new_avg = 0

    for new_cluster, alteration_data in alteration_metadata.items():
        new_scheme_data = scheme_metadata.get(new_cluster, None)
        old_scheme_data = old_scheme_metadata.get(new_cluster, None)

        if new_scheme_data is not None:
            new_avg = new_scheme_data["average_value"]

            full_data_lines.append(f"Cluster {new_cluster} changes:")

            if old_scheme_data is None:
                line = "\n".join(textwrap.wrap("".join([
                        f"> New cluster {new_cluster} has an average of ",
                        "{:.1f} % intercluster ".format(new_avg*100),
                        metric]), width=75))
                full_data_lines.append(textwrap.indent(line, "\t"))
            else:
                old_avg_gcs = old_scheme_data["average_value"]

                line = "\n".join(textwrap.wrap("".join([
                            f"> Reclustered {new_cluster} has an average of ",
                            "{:.1f} % intercluster ".format(new_avg*100),
                            metric]), width=75))
                full_data_lines.append(textwrap.indent(line, "\t"))

                line = "\n".join(textwrap.wrap("".join([
                            f"> Cluster {new_cluster} genomes had an average ",
                            "of {:.1f}% intercluster ".format(old_avg_gcs*100),
                            metric]), width=75))
                full_data_lines.append(textwrap.indent(line, "\t"))
        else:
            full_data_lines.append("Clustering noise:")

        for old_cluster, old_cluster_redist_data in alteration_data.items():
            percent = old_cluster_redist_data["percent_merged"] * 100
            old_members = old_cluster_redist_data["old_members"]
            num_old_members = len(old_members)
            full_quick_lines = num_old_members <= max_qlines

            if old_cluster is None:
                old_cluster = "Singleton"
                pass

            if not full_quick_lines:
                if new_cluster is None:
                    quick_data_lines.append(
                            f"Removed {round(percent, 1)}% of {old_cluster} ")
                elif old_scheme_data is None:
                    quick_data_lines.append(
                                    f"Created {new_cluster} <- "
                                    f"{old_cluster} ({round(percent, 1)}%)")
                else:
                    quick_data_lines.append(
                                f"Merged {new_cluster} <- "
                                f"{old_cluster} ({round(percent, 1)}%)")
            elif new_cluster is None:
                line = "\n".join(textwrap.wrap("".join([
                        f"> {round(percent, 1)}% of {old_cluster} genomes ",
                        "were labelled as noise"]), width=75))
                full_data_lines.append(textwrap.indent(line, "\t"))
            else:
                line = "\n".join(textwrap.wrap("".join([
                        f"> {round(percent, 1)}% of {old_cluster} genomes ",
                        f"were merged into {new_cluster}"]), width=75))
                full_data_lines.append(textwrap.indent(line, "\t"))

            for member in old_cluster_redist_data["old_members"]:
                if full_quick_lines:
                    if new_cluster is None:
                        quick_data_lines.append(
                                        f"Removed {member} [{old_cluster}] ")
                    elif old_scheme_data is None:
                        quick_data_lines.append(
                                    f"Created {new_cluster} <- "
                                    f"{member} [{old_cluster}]")
                    else:
                        quick_data_lines.append(f"Added {new_cluster} <- "
                                                f"{member} [{old_cluster}]")

                if new_cluster is None:
                    line = "\n".join(textwrap.wrap((
                                f"* {member} [{old_cluster}] "
                                f"was removed (noise)"),
                                width=71))
                    full_data_lines.append(textwrap.indent(line, "\t\t"))
                else:
                    line = "\n".join(textwrap.wrap((
                                f"* {member} [{old_cluster}] "
                                f"was added to cluster {new_cluster}"),
                                width=71))
                    full_data_lines.append(textwrap.indent(line, "\t\t"))

        full_data_lines.append("")

    if (not full_data_lines) and (not quick_data_lines):
        return

    filepath = working_dir.joinpath("log.txt")
    with filepath.open(mode="w") as filehandle:
        filehandle.write("pde_utils genome clustering pipeline\n")
        filehandle.write("".join(["=" * 79, "\n"]))
        filehandle.write(
                    f"Previous number of clusters: {total_old_clusters}\n")
        filehandle.write(f"Total number of clusters: {total_new_clusters}\n")
        filehandle.write("".join([
                            f"Old total intercluster {metric} ",
                            "average: {:.1f} %\n".format(total_old_avg*100)]))
        filehandle.write("".join([
                            f"New total intercluster {metric} ",
                            "average: {:.1f} %\n".format(total_new_avg*100)]))
        filehandle.write("\n\n")

        filehandle.write("Quick summary:\n")
        filehandle.write("".join(["=" * 79, "\n"]))
        for line in quick_data_lines:
            filehandle.write("".join([line, "\n"]))
        filehandle.write("\n\n")

        filehandle.write("Full summary:\n")
        filehandle.write("".join(["=" * 79, "\n"]))
        for line in full_data_lines:
            filehandle.write("".join([line, "\n"]))

    return True


def write_clustering_update_ticket(working_dir, scheme_alterations,
                                   field="Cluster", filename=None):
    if filename is None:
        filename = working_dir.with_suffix(".csv").name
    update_dicts = []

    for cluster, diff_data in scheme_alterations.items():
        for data_dict in diff_data:
            update_dict = {}

            if cluster is None:
                cluster = "NULL"

            update_data = ("phage", field, cluster, "PhageID", data_dict["id"])
            for i in range(len(TICKET_HEADER)):
                update_dict[TICKET_HEADER[i]] = update_data[i]

            update_dicts.append(update_dict)

    if not update_dicts:
        return False

    filepath = working_dir.joinpath(filename)
    pdm_fileio.export_data_dict(update_dicts, filepath, TICKET_HEADER,
                                include_headers=True)

    return True


def write_genome_fastas(gs_and_ts, fasta_dir, verbose=False, threads=1):
    if verbose:
        print("Writing genome fasta files...")

    work_items = []
    fasta_path_map = {}
    for seq_id, seq in gs_and_ts:
        seq_path = fasta_dir.joinpath("".join([seq_id, ".fasta"]))
        fasta_path_map[seq_id] = seq_path

        work_items.append(({seq_id: seq.decode("utf-8")}, seq_path))

    multithread.multithread(work_items, threads, pdm_fileio.write_fasta,
                            verbose=verbose)

    return fasta_path_map


def sketch_genome_fastas(fasta_path_map, sketch_dir, verbose=False,
                         threads=1, kmer=DEFAULT_SETTINGS["kmer"],
                         sketch=DEFAULT_SETTINGS["sketch"]):
    if verbose:
        print("Sketching genome fasta files...")

    work_items = []
    sketch_path_map = {}
    for seq_id, fasta_path in fasta_path_map.items():
        sketch_path = sketch_dir.joinpath(f"{seq_id}.msh")
        sketch_path_map[seq_id] = sketch_path

        work_items.append((fasta_path, sketch_path, kmer, sketch))

    parallelize.parallelize(work_items, threads, alignment.mash_sketch,
                            verbose=verbose)

    return sketch_path_map


# MISC FUNCTIONS
# -----------------------------------------------------------------------------
def create_temp_path(temp_path):
    temp_dir = Path(temp_path)

    if temp_dir.is_dir():
        shutil.rmtree(temp_dir)

    temp_dir.mkdir()
    return temp_dir


# DILAPIDATED FUNCTIONS
# -----------------------------------------------------------------------------

def get_average_intercluster_identity(matrix):
    centroid = matrix.get_centroid()

    return matrix.get_average_value(centroid)


def get_intersecting_average_identity(matrix, subject_labels, query_labels):
    subject_indicies = []
    for label in subject_labels:
        subject_indicies.append(matrix.get_index_from_label(label))

    query_indicies = []
    for label in query_labels:
        query_indicies.append(matrix.get_index_from_label(label))

    intersecting_gcs = 0
    num_intersections = 0
    for subject_index in subject_indicies:
        for query_index in query_indicies:
            num_intersections += 1

            intersecting_gcs += matrix.get_cell(subject_index, query_index)

    if num_intersections > 0:
        return (intersecting_gcs / num_intersections)


def write_cluster_analysis(intracluster_edges, working_dir, file_name=None):
    data_dicts = []
    for edge in intracluster_edges:
        data_dict = {}
        for i in range(len(edge)):
            data_dict[CLUSTER_ANALYSIS_HEADER[i]] = edge[i]

        data_dicts.append(data_dict)

    if file_name is None:
        file_name = working_dir.name

    filepath = working_dir.joinpath(file_name).with_suffix(".csv")
    pdm_fileio.export_data_dict(data_dicts, filepath, CLUSTER_ANALYSIS_HEADER,
                                include_headers=True)


def retrieve_intracluster_edges(db_filter, working_dir, node_names, matrix,
                                lookup_dict, gcs=DEFAULT_SETTINGS["gcs"],
                                kmer=DEFAULT_SETTINGS["kmer"],
                                sketch=DEFAULT_SETTINGS["sketch"],
                                threads=1, verbose=False):
    intracluster_edges = []
    node_name_set = set()

    for i in range(matrix.size):
        target_name = matrix.labels[i]
        target_cluster = lookup_dict[target_name]
        for j in range(i, matrix.size):
            query_name = matrix.labels[j]
            query_cluster = lookup_dict[query_name]

            if (target_cluster is not None) and (query_cluster is not None):
                if query_cluster == target_cluster:
                    continue

            pairwise_gcs = matrix.get_cell(i, j)
            if pairwise_gcs >= gcs:
                node_name_set.add(target_name)
                node_name_set.add(query_name)

                intracluster_edges.append((
                                        target_name, str(target_cluster),
                                        query_name, str(query_cluster),
                                        str(round(pairwise_gcs, 3))))

    db_filter.values = list(node_name_set)
    gs_and_ts = db_filter.select(["phage.PhageID", "phage.Sequence"],
                                 return_dict=False)

    fasta_dir = create_temp_path(str(working_dir.joinpath("fasta")))
    fasta_path_map = write_genome_fastas(
                                 gs_and_ts, fasta_dir, verbose=verbose,
                                 threads=threads)

    sketch_dir = create_temp_path(str(working_dir.joinpath("sketches")))
    sketch_path_map = sketch_genome_fastas(
                                      fasta_path_map, sketch_dir,
                                      verbose=verbose, threads=threads,
                                      kmer=kmer, sketch=sketch)

    work_items = []
    for edge in intracluster_edges:
        work_items.append((sketch_path_map[edge[0]], sketch_path_map[edge[2]],
                           edge))

    if verbose:
        print("Calculating phage genome ANI...")
    intracluster_edges = parallelize.parallelize(
                            work_items, threads, calculate_ani_process,
                            verbose=verbose)

    intracluster_edges.sort(reverse=True, key=lambda x: (
                                                float(x[4]) + float(x[5])) / 2)
    return intracluster_edges


def calculate_ani_process(query_path, subject_path, edge_data):
    ani_data = alignment.mash_dist(query_path, subject_path)
    ani_data = ani_data.split("\t")
    if len(ani_data) != 5:
        return

    ani = 1 - float(ani_data[2])
    edge_data = edge_data + (round(ani, 3),)

    return edge_data
