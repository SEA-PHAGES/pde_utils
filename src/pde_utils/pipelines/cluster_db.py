import argparse
from pathlib import Path
import shutil
import string
import textwrap
import time

from pdm_utils.functions import (basic, configfile, multithread, parallelize,
                                 fileio as pdm_fileio, pipelines_basic)
from pdm_utils.pipelines.revise import TICKET_HEADER

from pde_utils.functions import (alignment, clustering, seq_distance,
                                 sql_queries)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_cluster_db")
TEMP_DIR = (f"/tmp/{DEFAULT_FOLDER_NAME}_temp")

CLUSTER_DB_SUBPIPELINES = ["analyze", "cluster"]
DEFAULT_SETTINGS = {"kmer": 15, "sketch": 25000, "gcs": 0.35, "ani": 0.7}

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
                       evaluate=args.dump_evaluation,
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
                        type=float)

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
                    cluster_prefix=None)

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_cluster_db(
                alchemist, pipeline="cluster", folder_path=None,
                folder_name=DEFAULT_FOLDER_NAME, values=None, verbose=None,
                filters="", groups=[], threads=1,
                kmer=DEFAULT_SETTINGS["kmer"],
                sketch=DEFAULT_SETTINGS["sketch"],
                gcs=DEFAULT_SETTINGS["gcs"], ani=DEFAULT_SETTINGS["ani"],
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
            print("Querying MySQL database for genome metadata...")
        cluster_metadata = query_cluster_metadata(db_filter)

        if verbose:
            print("Calculating gene content similarity matrix...")
        gcs_matrix = calculate_gcs_matrix(
                                                alchemist, db_filter.values,
                                                verbose=verbose, cores=threads)

        pipelines_basic.create_working_dir(mapped_path)
        if pipeline == "analyze":
            intracluster_edges = retrieve_intracluster_edges(
                                        db_filter, temp_dir,
                                        db_filter.values, gcs_matrix,
                                        cluster_metadata[0], threads=threads,
                                        verbose=verbose, kmer=kmer,
                                        sketch=sketch, gcs=gcs)

            write_cluster_analysis(intracluster_edges, mapped_path)

        elif pipeline == "cluster":
            if verbose:
                print("Reclustering database genomes...")
            cluster_scheme = gcs_cluster(
                                    mapped_path, gcs_matrix,
                                    cluster_metadata[0], cluster_metadata[1],
                                    gcs=gcs, evaluate=evaluate,
                                    verbose=verbose,
                                    cluster_prefix=cluster_prefix)

            if subcluster:
                sketch_path_map = sketch_genomes(db_filter, temp_dir,
                                                 verbose=verbose)

                if verbose:
                    print("Subclustering database genomes...")
                ani_subcluster(mapped_path, sketch_path_map, cluster_scheme,
                               cluster_metadata[0], cluster_metadata[1],
                               cluster_metadata[2],
                               verbose=verbose, ani=ani, evaluate=evaluate)

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


# MATRIX HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def calculate_gcs_matrix(alchemist, phage_ids, cores=1, verbose=False):
    phage_gc_nodes = []
    for phage in phage_ids:
        phage_gc_nodes.append(
                    sql_queries.get_distinct_phams_from_organism(
                                                    alchemist, phage))

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


def gcs_cluster(working_dir, gcs_matrix, cluster_lookup, cluster_seqid_map,
                verbose=False, gcs=DEFAULT_SETTINGS["gcs"], evaluate=False,
                cluster_prefix=None):
    cluster_scheme = gcs_matrix.get_clusters(gcs)

    old_clusters = list(cluster_seqid_map.keys())
    cluster_scheme, old_cluster_histograms = recluster(
                                            cluster_scheme, cluster_lookup,
                                            old_clusters, verbose=verbose,
                                            cluster_prefix=cluster_prefix)

    altered_cluster_scheme = {}
    for cluster, cluster_members in cluster_scheme.items():
        if ((cluster not in old_clusters or
             len(old_cluster_histograms[cluster]) > 1)
                and cluster is not None):
            altered_cluster_scheme[cluster] = cluster_members

    wrote_eval = False
    if evaluate:
        wrote_eval = evaluate_clustering(
                            working_dir, gcs_matrix, altered_cluster_scheme,
                            old_cluster_histograms, cluster_lookup,
                            cluster_seqid_map, "gene content similarity")

    wrote_ticket = write_reclustering_update_ticket(
                                     working_dir, cluster_lookup,
                                     altered_cluster_scheme, field="Cluster")

    wrote = wrote_eval or wrote_ticket
    if not wrote:
        shutil.rmtree(working_dir)

    return cluster_scheme


def ani_subcluster(working_dir, sketch_path_map, cluster_scheme,
                   cluster_lookup, cluster_seqid_map,
                   subcluster_lookup,
                   threads=1, verbose=False, ani=DEFAULT_SETTINGS["ani"],
                   evaluate=False):

    for cluster, cluster_members in cluster_scheme.items():
        if cluster is None:
            continue
        old_cluster_members = cluster_seqid_map.get(cluster, list())

        noncluster_members = [x for x in cluster_members
                              if x not in old_cluster_members]

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

            seqids = subcluster_seqid_map.get(nonmember_cluster, list())
            seqids.append(nonmember)
            subcluster_seqid_map[nonmember_cluster] = seqids

        if verbose:
            print(f"Subclustering {cluster}...")

        ani_matrix = calculate_ani_matrix(cluster_members, sketch_path_map,
                                          cores=threads, verbose=verbose)

        subcluster_scheme = ani_matrix.get_clusters(ani)
        subcluster_scheme, old_subcluster_histogram = recluster(
                                        subcluster_scheme,
                                        altered_subcluster_lookup,
                                        list(old_subclusters),
                                        subcluster=cluster, verbose=verbose)

        altered_subcluster_scheme = {}
        for subcluster, subcluster_members in subcluster_scheme.items():
            if ((subcluster not in old_subclusters or
                 len(old_subcluster_histogram[subcluster]) > 1)
                    and subcluster is not None):
                altered_subcluster_scheme[subcluster] = subcluster_members

        subcluster_dir = working_dir.joinpath(str(cluster))
        pipelines_basic.create_working_dir(subcluster_dir)

        wrote_ticket = write_reclustering_update_ticket(
                                subcluster_dir, subcluster_lookup,
                                subcluster_scheme, field="Subcluster")

        wrote_eval = False
        if evaluate:
            wrote_eval = evaluate_clustering(
                                subcluster_dir, ani_matrix,
                                altered_subcluster_scheme,
                                old_subcluster_histogram,
                                altered_subcluster_lookup,
                                subcluster_seqid_map,
                                "average nucleotide identity")

        wrote = wrote_eval or wrote_ticket

        if not wrote:
            shutil.rmtree(subcluster_dir)


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
            if counter[i] < (len(alphabet) - 1):
                counter[i] += 1
                filled = False
                break

        if filled:
            chars += 1
            if chars >= 5:
                break


def gen_new_subcluster(cluster, old_subclusters):
    for i in range(1, 1000):
        subcluster = "".join([cluster, str(i)])
        if subcluster not in old_subclusters:
            return subcluster


def calculate_ani(subject_path, query_path):
    mash_output = alignment.mash_dist(subject_path, query_path)
    ani_data = mash_output.split("\t")

    if len(ani_data) == 5:
        return 1 - float(ani_data[2])


def recluster(cluster_scheme, cluster_lookup, old_clusters,
              verbose=False, cluster_prefix=None, subcluster=None):
    old_cluster_histograms = {}
    new_cluster_scheme = {}
    for tmp_cluster_name, cluster_members in cluster_scheme.items():
        old_cluster_histogram = {}
        for cluster_member in cluster_members:
            old_cluster = cluster_lookup[cluster_member]

            num_cluster = old_cluster_histogram.get(old_cluster, 0)
            num_cluster += 1
            old_cluster_histogram[old_cluster] = num_cluster

        old_clusters_histogram = basic.sort_histogram(old_cluster_histogram)
        merging_old_clusters = list(old_clusters_histogram.keys())
        new_cluster = merging_old_clusters[0]

        if len(cluster_members) > 1:
            if new_cluster is None:
                if len(merging_old_clusters) > 1:
                    for old_cluster in merging_old_clusters[1:]:
                        new_cluster = old_cluster

                        if (new_cluster is not None and
                                new_cluster not in new_cluster_scheme.keys()
                                and new_cluster in old_clusters):
                            break

            if ((new_cluster is None) or
                (new_cluster in new_cluster_scheme.keys()) or
                    (new_cluster not in old_clusters)):
                if subcluster is None:
                    new_cluster = gen_new_cluster(
                                            old_clusters,
                                            cluster_prefix=cluster_prefix)
                    if verbose:
                        print(f"...Created new cluster '{new_cluster}'...")

                    old_clusters.append(new_cluster)

                elif len(cluster_scheme) > 1:
                    new_cluster = gen_new_subcluster(subcluster,
                                                     old_clusters)
                    if verbose:
                        print(f"...Created new subcluster '{new_cluster}'...")

                    old_clusters.append(new_cluster)

            if subcluster is not None and len(cluster_scheme) <= 1:
                new_cluster = None

        else:
            new_cluster = None

        temp_cluster_members = new_cluster_scheme.get(new_cluster, list())

        new_cluster_scheme[new_cluster] = (temp_cluster_members +
                                           list(cluster_members))
        old_cluster_histograms[new_cluster] = old_cluster_histogram

    return (new_cluster_scheme, old_cluster_histograms)


def evaluate_clustering(working_dir, matrix, cluster_scheme,
                        old_cluster_histograms, cluster_lookup,
                        cluster_seqid_map, metric):
    full_data_lines = []
    quick_data_lines = []

    for new_cluster, new_cluster_members in cluster_scheme.items():
        old_cluster_histogram = old_cluster_histograms[new_cluster]

        full_data_lines.append(f"Cluster {new_cluster} changes:")

        new_avg_gcs = get_average_intercluster_identity(
                                matrix, list(new_cluster_members))

        old_cluster_members = cluster_seqid_map.get(new_cluster)
        if old_cluster_members is None:
            line = "\n".join(textwrap.wrap("".join([
                    f"> New cluster {new_cluster} has an average of ",
                    "{:.1f} % intercluster ".format(new_avg_gcs*100),
                    metric]), width=75))
            full_data_lines.append(textwrap.indent(line, "\t"))
            full_data_lines.append("")

            quick_data_lines.append(f"Created {new_cluster}")
            continue

        old_avg_gcs = get_average_intercluster_identity(
                                        matrix, old_cluster_members)

        line = "\n".join(textwrap.wrap("".join([
                       f"> Cluster {new_cluster} genomes had an average of ",
                       "{:.1f}% intercluster ".format(old_avg_gcs*100),
                       metric]), width=75))
        full_data_lines.append(textwrap.indent(line, "\t"))

        for old_cluster, num_members in old_cluster_histogram.items():
            if new_cluster == old_cluster:
                continue

            added_cluster_members = cluster_seqid_map[old_cluster]
            merged_cluster_members = set(added_cluster_members).intersection(
                                                    set(new_cluster_members))

            old_merged_avg_gcs = get_average_intercluster_identity(
                                    matrix,
                                    list(cluster_seqid_map[old_cluster]))

            percent = (round(num_members / len(
                             cluster_seqid_map[old_cluster]), 3) * 100)
            if percent < 50 or old_cluster is None:
                if old_cluster is None:
                    old_cluster = "Singleton"

                for merged_cluster_member in merged_cluster_members:
                    avg_int_gcs = get_intersecting_average_identity(
                                        matrix, list(old_cluster_members),
                                        [merged_cluster_member])
                    quick_data_lines.append(
                                    f"Added {new_cluster} <- "
                                    f"{merged_cluster_member} [{old_cluster}]")

                    line = "\n".join(textwrap.wrap((
                          f"* {merged_cluster_member} [{old_cluster}] "
                          f"was added to cluster {new_cluster}"), width=71))
                    full_data_lines.append(textwrap.indent(line, "\t\t"))

                    line = "\n".join(textwrap.wrap("".join([
                          f"- {merged_cluster_member} has an ",
                          "average of {:.1f}% ".format(avg_int_gcs*100),
                          metric,
                          f" with cluster {new_cluster} genomes"]), width=67))
                    full_data_lines.append(textwrap.indent(line, "\t\t\t"))
                continue

            quick_data_lines.append(f"Merged {new_cluster} <- "
                                    f"{old_cluster} ({percent}%)")

            line = "\n".join(textwrap.wrap((
                             f"* {percent}% of cluster {old_cluster} "
                             f"was merged into cluster {new_cluster}"),
                             width=71))
            full_data_lines.append(textwrap.indent(line, "\t\t"))

            line = "\n".join(textwrap.wrap("".join([
                    f"- Cluster {old_cluster} genomes had an average of ",
                    "{:.1f}% intercluster ".format(old_merged_avg_gcs*100),
                    metric]), width=67))
            full_data_lines.append(textwrap.indent(line, "\t\t\t"))

            avg_int_gcs = get_intersecting_average_identity(
                                        matrix, list(old_cluster_members),
                                        list(merged_cluster_members))
            if avg_int_gcs is not None:
                line = "\n".join(textwrap.wrap("".join([
                        f"- Merged cluster {old_cluster} genomes have an ",
                        "average of {:.1f}% ".format(avg_int_gcs*100),
                        f"{metric} with cluster ",
                        f"{new_cluster} genomes"]), width=67))
                full_data_lines.append(textwrap.indent(line, "\t\t\t"))

        line = "\n".join(textwrap.wrap("".join([
                    f"> Reclustered {new_cluster} has an average of ",
                    "{:.1f} % intercluster ".format(new_avg_gcs*100),
                    metric]), width=75))
        full_data_lines.append(textwrap.indent(line, "\t"))
        full_data_lines.append("")

    if (not full_data_lines) and (not quick_data_lines):
        return

    filepath = working_dir.joinpath("log.txt")
    with filepath.open(mode="w") as filehandle:
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


def write_reclustering_update_ticket(working_dir, cluster_lookup,
                                     cluster_scheme, field="Cluster",
                                     filename=None):
    if filename is None:
        filename = working_dir.with_suffix(".csv").name
    update_dicts = []

    for cluster, cluster_members in cluster_scheme.items():
        for member in cluster_members:
            if cluster != cluster_lookup[member]:
                update_dict = {}

                if cluster is None:
                    cluster = "NULL"

                update_data = ["phage", field, cluster, "PhageID", member]
                for i in range(len(TICKET_HEADER)):
                    update_dict[TICKET_HEADER[i]] = update_data[i]

                update_dicts.append(update_dict)

    if not update_dicts:
        return False

    filepath = working_dir.joinpath(filename)
    pdm_fileio.export_data_dict(update_dicts, filepath, TICKET_HEADER,
                                include_headers=True)

    return True


def get_average_intercluster_identity(matrix, labels):
    label_indicies = []
    for label in labels:
        label_indicies.append(matrix.get_index_from_label(label))

    avg_gcs = 0
    num_edges = 0
    for i in range(len(label_indicies)):
        for j in range(i+1, len(label_indicies)):
            num_edges += 1

            avg_gcs += matrix.get_cell(label_indicies[i],
                                       label_indicies[j])
    if num_edges == 0:
        return 1

    return avg_gcs / num_edges


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


def create_temp_path(temp_path):
    temp_dir = Path(temp_path)

    if temp_dir.is_dir():
        shutil.rmtree(temp_dir)

    temp_dir.mkdir()
    return temp_dir


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
