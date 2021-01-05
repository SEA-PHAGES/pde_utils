import argparse
from pathlib import Path
import shutil
import time

from pdm_utils.functions import (configfile, multithread, parallelize,
                                 fileio as pdm_fileio, pipelines_basic)

from pde_utils.functions import (alignment, phage_similarity)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = (f"{time.strftime('%Y%m%d')}_cluster_db")
TEMP_DIR = (f"/tmp/{DEFAULT_FOLDER_NAME}_temp")

CLUSTER_DB_SUBPIPELINES = ["analyze"]
DEFAULT_SETTINGS = {"kmer": 15, "sketch": 25000, "gcs": 0.35}

CLUSTER_ANALYSIS_HEADER = ["Subject PhageID", "Subject Cluster",
                           "Query PhageID", "Query Cluster", "GCS",
                           "Estimated ANI"]


def main(unparsed_args_list):
    args = parse_cluster_db(unparsed_args_list)

    config = configfile.build_complete_config(args.config_file)

    alchemist = pipelines_basic.build_alchemist(args.database, config=config)
    values = pipelines_basic.parse_value_input(args.input)

    execute_cluster_db(alchemist, pipeline=args.pipeline,
                       folder_path=args.folder_path,
                       folder_name=args.folder_name,
                       values=values, verbose=args.verbose,
                       filters=args.filters, groups=args.groups,
                       threads=args.number_threads, kmer=args.kmer_size,
                       sketch=args.sketch_size,
                       gcs=args.gene_content_similarity_min)


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

    subparsers = []
    subparser_gen = parser.add_subparsers(dest="pipeline", required=True)

    analyze_parser = subparser_gen.add_parser("analyze")
    subparsers.append(analyze_parser)

    for subparser in subparsers:
        subparser.add_argument("database", type=str, help=DATABASE_HELP)
        subparser.add_argument("-c", "--config_file", help=CONFIG_FILE_HELP,
                               type=pipelines_basic.convert_file_path)
        subparser.add_argument("-m", "--folder_name", type=str,
                               help=FOLDER_NAME_HELP)
        subparser.add_argument("-o", "--folder_path", type=Path,
                               help=FOLDER_PATH_HELP)
        subparser.add_argument("-v", "--verbose", action="store_true",
                               help=VERBOSE_HELP)

        subparser.add_argument("-if", "--import_file",
                               type=pipelines_basic.convert_file_path,
                               help=IMPORT_FILE_HELP, dest="input")
        subparser.add_argument("-in", "--import_names", nargs="*",
                               help=IMPORT_NAMES_HELP, dest="input")

        subparser.add_argument("-w", "--where", nargs="?",
                               help=WHERE_HELP, dest="filters")
        subparser.add_argument("-g", "--group_by", nargs="+",
                               help=GROUP_BY_HELP, dest="groups")
        subparser.add_argument("-np", "--number_threads", type=int, nargs="?",
                               help=NUMBER_THREADS_HELP)

    analyze_parser.add_argument("-kmer", "--kmer_size", type=int,
                                help=KMER_SIZE_HELP)
    analyze_parser.add_argument("-sketch", "--sketch_size", type=int,
                                help=SKETCH_SIZE_HELP)
    analyze_parser.add_argument("-gcs", "--gene_content_similarity_min",
                                type=float,
                                help=GENE_CONTENT_SIMILARITY_MIN_HELP)

    for subparser in subparsers:
        subparser.set_defaults(
                        folder_name=DEFAULT_FOLDER_NAME, folder_path=None,
                        config_file=None, verbose=False, input=[],
                        filters="", groups=[], number_threads=1,
                        kmer_size=DEFAULT_SETTINGS["kmer"],
                        sketch_size=DEFAULT_SETTINGS["sketch"],
                        gene_content_similarity_min=DEFAULT_SETTINGS["gcs"])

    parsed_args = parser.parse_args(unparsed_args_list[2:])
    return parsed_args


def execute_cluster_db(
                alchemist, pipeline="analyze", folder_path=None,
                folder_name=DEFAULT_FOLDER_NAME, values=None, verbose=None,
                filters="", groups=[], threads=1,
                kmer=DEFAULT_SETTINGS["kmer"],
                sketch=DEFAULT_SETTINGS["sketch"],
                gcs=DEFAULT_SETTINGS["gcs"]):
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
            print("Calculating gene content similarity matrix...")

        gcs_matrix = phage_similarity.calculate_gcs_matrix(
                                                alchemist, db_filter.values,
                                                verbose=verbose, cores=threads)
        if pipeline == "analyze":
            cluster_data = db_filter.retrieve("phage.Cluster")
            cluster_lookup = {}
            for phage_id, data_dict in cluster_data.items():
                cluster_lookup[phage_id] = data_dict["Cluster"][0]

            intracluster_edges = retrieve_intracluster_edges(
                                            db_filter, temp_dir,
                                            db_filter.values, gcs_matrix,
                                            cluster_lookup, threads=threads,
                                            verbose=verbose, kmer=kmer,
                                            sketch=sketch, gcs=gcs)

            pipelines_basic.create_working_dir(mapped_path)
            write_cluster_analysis(intracluster_edges, mapped_path)


def retrieve_intracluster_edges(db_filter, working_dir, node_names, matrix,
                                lookup_dict, gcs=DEFAULT_SETTINGS["gcs"],
                                kmer=DEFAULT_SETTINGS["kmer"],
                                sketch=DEFAULT_SETTINGS["sketch"],
                                threads=1, verbose=False):
    intracluster_edges = []
    node_name_set = set()

    for i in range(len(matrix)):
        target_cluster = lookup_dict[node_names[i]]
        for j in range(1, len(matrix[i])):
            query_cluster = lookup_dict[node_names[i+j]]

            if (target_cluster is not None) and (query_cluster is not None):
                if query_cluster == target_cluster:
                    continue

            if matrix[i][j] >= gcs:
                node_name_set.add(node_names[i])
                node_name_set.add(node_names[i+j])

                intracluster_edges.append((
                                        node_names[i], str(target_cluster),
                                        node_names[i+j], str(query_cluster),
                                        str(round(matrix[i][j], 3))))

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
                            work_items, threads, intracluster_edge_ani_process,
                            verbose=verbose)

    intracluster_edges.sort(reverse=True, key=lambda x: (
                                                float(x[4]) + float(x[5])) / 2)
    return intracluster_edges


def intracluster_edge_ani_process(query_path, subject_path, edge_data):
    ani_data = alignment.mash_dist(query_path, subject_path)
    ani_data = ani_data.split("\t")
    if len(ani_data) != 5:
        return

    ani = 1 - float(ani_data[2])
    edge_data = edge_data + (round(ani, 3),)

    return edge_data


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
        sketch_path = sketch_dir.joinpath(fasta_path.with_suffix("").name)
        sketch_path_map[seq_id] = sketch_path.with_suffix(".msh")

        work_items.append((fasta_path, sketch_path, kmer, sketch))

    parallelize.parallelize(work_items, threads, alignment.mash_sketch,
                            verbose=verbose)

    return sketch_path_map
