import argparse
import os
import shutil
import sys
import time
from decimal import Decimal
from pathlib import Path

import shlex
from subprocess import PIPE, Popen

from networkx import compose
from networkx import MultiDiGraph, Graph
from networkx import readwrite
from pdm_utils.functions import basic
from pdm_utils.functions import fileio
from pdm_utils.functions import pipelines_basic
from pdm_utils.functions import phameration

from pde_utils.classes import clustal
from pde_utils.functions import alignment
from pde_utils.functions import search

#GLOBAL VARIABLES
#-----------------------------------------------------------------------------
DEFAULT_FOLDER_NAME = f"{time.strftime('%Y%m%d')}_pan"
DEFAULT_FOLDER_PATH = Path.cwd()

NETWORKX_FILE_TYPES = ["csv", "gexf", "gpickle", "graphml", "al", "nl",
                       "json", "cyjs", "yaml", "pajek", "shp"]
EDGE_WEIGHTS = ["DB", "Distance"]

#MAIN FUNCTIONS
#-----------------------------------------------------------------------------
def main(unparsed_args_list):
    args = parse_build_pan(unparsed_args_list)

    alchemist = pipelines_basic.build_alchemist(args.database)
    values = pipelines_basic.parse_value_input(args.input)
    
    execute_build_pan(alchemist, args.hhsuite_database, args.folder_path,
                      args.folder_name, values=values, verbose=args.verbose,
                      filters=args.filters, groups=args.groups,
                      threads=args.number_threads, M=args.min_percent_gaps,
                      aD=args.min_distance, B=args.DB_stiffness,
                      file_format=args.file_format)

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
        PAN building option that controls the cutoff for the DB separation value
        required to make a connection.
            Follow seelction argument with the DB stiffness parameter [0-1].
        """
    MIN_DISTANCE_HELP = """
        PAN building option that controls the initial centroid distance cutoff 
        between phams during the building of neighborhoods.
            Follow selection argument with the minimum average distance [0-1].
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
    
    parser.add_argument("-np", "--number_threads", type=int, nargs="?",
                        help=NUMBER_THREADS_HELP)
    parser.add_argument("-B", "--DB_stiffness", type=Decimal,
                        help=DB_STIFFNESS_HELP)
    parser.add_argument("-aD", "--min_distance", type=Decimal,
                        help=MIN_DISTANCE_HELP)
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

    parser.set_defaults(folder_name=DEFAULT_FOLDER_NAME,
                        folder_path=DEFAULT_FOLDER_PATH,
                        verbose=False, input=[],
                        filters="", groups=[],
                        db_name=None, number_threads=1,
                        DB_stiffness=0.5, min_distance=0.15,
                        min_percent_gaps=50, file_format="graphml")  
    
    parsed_args = parser.parse_args(unparsed_args_list[1:])
    return parsed_args

def execute_build_pan(alchemist, hhdb_path, folder_path, folder_name, 
                      values=None, verbose=False, filters="", groups=[],
                      threads=1, M=50, aD=20, B=0.5, file_format="csv"):
    db_filter = pipelines_basic.build_filter(alchemist, "pham", filters,
                                                           values=values)

    pan_path = folder_path.joinpath(folder_name)
    pan_path = basic.make_new_dir(folder_path, pan_path, attempt=50)

    conditionals_map = {}
    pipelines_basic.build_groups_map(db_filter, pan_path,
                                     conditionals_map,
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

        aln_dir = mapped_path.joinpath("pham_alignments")
        aln_dir.mkdir()

        maps_tuple = build_pan_base(alchemist.engine, db_filter.values, aln_dir, 
                                    threads=threads, M=M, verbose=verbose)
        aln_obj_dict, mat_obj_dict = read_pan_base(db_filter.values, maps_tuple) 
        caln  = write_centroid_fastas(aln_dir, aln_obj_dict, mat_obj_dict,
                                      threads=threads, verbose=verbose)

        cmat = align_centroids(caln, aln_dir, threads=threads)
        sys.exit(1)
 
        pan = build_pan_neighborhoods(alchemist, db_filter.values,
                                      mat_obj_dict, cmat, aD=aD, B=B,
                                      verbose=verbose)
        write_graph(pan, file_format, mapped_path, folder_name)
        
def build_pan_base(engine, values, aln_dir, threads=1, M=50, verbose=False):
    fasta_path_map = alignment.create_pham_fastas(engine, values, aln_dir)
    mat_path_map = alignment.align_pham_fastas(fasta_path_map, mat_out=True, 
                                               threads=threads, verbose=verbose)
    hmm_path_map = alignment.create_pham_hmms(fasta_path_map, name=True,
                                              M=M, verbose=verbose)

    return (fasta_path_map, mat_path_map, hmm_path_map)

def read_pan_base(values, base_maps_tuple):
    aln_obj_dict = {}
    mat_obj_dict = {}
    for pham in values:        
        aln_path = base_maps_tuple[0].get(pham)
        mat_path = base_maps_tuple[1].get(pham)
     
        aln = clustal.MultipleSequenceAlignment(aln_path, fmt="fasta")
        aln.parse_alignment()

        if not mat_path.is_file():
            mat = None
        else:
            mat = clustal.PercentIdentityMatrix(mat_path)
            mat.parse_matrix()

        aln_obj_dict[pham] = aln
        mat_obj_dict[pham] = mat
    
    return aln_obj_dict, mat_obj_dict

def write_centroid_fastas(aln_dir, aln_obj_dict, mat_obj_dict, threads=1,
                                                                 verbose=False):
    caln_path = aln_dir.joinpath("centroids.fasta")

    centroid_seqs = {} 
    for pham in aln_obj_dict.keys():
        aln = aln_obj_dict.get(pham)
        mat = mat_obj_dict.get(pham)

        if mat is None:
            centroid, centroid_seq = aln.longest_gene()
        else:
            centroid = mat.get_centroid()
            centroid_seq = aln.get_sequence(centroid, gaps=False)
            
        centroid_seqs[pham] = centroid_seq

        #lone_centroid_path = aln_dir.joinpath(f"{pham}_centroid.fasta")
        #fileio.write_fasta({pham : centroid_seq}, lone_centroid_path)
    
    fileio.write_fasta(centroid_seqs, caln_path)

    return caln_path

def build_pan_neighborhoods(alchemist, phams, mat_obj_dict, 
                                              cmat, verbose=False,
                                              aD=20, B=0.5): 
    pan = MultiDiGraph()

    spread_dict = {}
    for pham in phams:
        mat = mat_obj_dict.get(pham)
        spread = 0
        if not mat is None:
            centroid = mat.get_centroid()
            spread = 100 - mat.get_average_identity(centroid)

        spread_dict[pham] = spread

        pan.add_node(pham, Spread=spread)

    for pham in phams:  
        neighbors = cmat.get_nearest_neighbors(str(pham), aD)
        centroid_nodes = cmat.node_names
        for neighbor in neighbors: 
            q_spread = spread_dict[pham]
            t_spread = spread_dict[(int(neighbor))]

            query_idx = centroid_nodes.index(str(pham))
            target_idx = centroid_nodes.index(neighbor)
            row = cmat.get_row(query_idx)
            pid = row[target_idx]
            dist = 100 - pid

            DB_sep = round((t_spread)/dist, 3)
            if DB_sep >= B:
                pan.add_edge(pham, int(neighbor), Distance=dist, 
                                                  DB=DB_sep)
    
    return pan

def write_graph(graph, file_format, export_path, file_name):
    file_path = export_path.joinpath(f"{file_name}.{file_format}")

    if file_format == "csv":
        readwrite.edgelist.write_edgelist(graph, file_path, delimiter=",",
                                          data=EDGE_WEIGHTS) 
    elif file_format == "gexf":
        readwrite.gexf.write_gexf(graph, file_path)
    elif file_format == "gml":
        readwrite.gml.write_gml(graph, file_path) 
    elif file_format == "gpickle":
        readwrite.gpickle.write_gpickle(graph, file_path)
    elif file_format == "graphml":
        readwrite.graphml.write_graphml(graph, file_path)
    elif file_format == "json" or file_format == "cyjs":
        if file_format == "json":
            json_data = readwrite.json_graph.node_link_data(graph)
        else:
            json_data = readwrite.json_graph.cytoscape_data(graph)
       
        file_path.touch()
        file_handle = file_path.open(mode="w")
        json.dump(json_data, file_handle)
        file_handle.close()
    elif file_format == "yaml":
        readwrite.nx_yaml.write_yaml(graph, file_path)
    elif file_format == "pajek":
        readwrite.pajek.write_pajek(graph, file_path)
    elif file_format == "shp":
        readwrite.nx_shp.write_shp(graph, export_path)
    elif file_format == "al":
        readwrite.adjlist.write_adjlist(graph, file_path)
    elif file_format == "nl":
        file_path.touch()
        file_handle = file_path.open(mode="w")
        for node in list(graph.nodes()):
            file_handle.write(f"{node}\n")
        file_handle.close()
    else:
        raise ValueError("Graph output format is not recognized.")

def convert_hmm_db_path(path):
    return Path(path)



#Methods of clustering centroid sequences

def align_centroids(caln_path, aln_dir, threads=1):
    cmat_path = aln_dir.joinpath("centroids.mat")

    alignment.clustalo(caln_path, caln_path, mat_out_path=cmat_path,
                       outfmt="fasta", threads=threads)
    
    return cmat_path

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


CLUSTER_ARGS = {"tmp_dir" : None,
                "identity" : None,
                "coverage" : None, 
                "e_value" : None,
                "sens" : None,
                "steps" : None,
                "threads" : None,
                "aln_mode"  : None,
                "cov_mode" : None,
                "clu_mode" : None}


def clust_centroids(caln_path, aln_dir):
    db_dir = aln_dir.joinpath("DB_tmp")
    
    seq_db = db_dir.joinpath("sequenceDB")
    clu_db = db_dir.joinpath("clusterDB")
    psf_db = db_dir.joinpath("seqfileDB")
    pre_neighborhoods_file = db_dir.joinpath("pre_neighborhoods.fasta")
    
    phameration.mmseqs_createdb(caln_path, seq_db)
    phameration.mmseqs_cluster(seq_db, clu_db, CLUSTER_ARGS)
    phameration.mmseqs_createseqfiledb(seq_db, clu_db, psf_db)
    phameration.mmseqs_result2flat(seq_db, seq_db, psf_db, 
                                        pre_neighborhoods_file) 
    pre_neighborhoods = phameration.parse_mmseqs_output(
                                        pre_neighborhoods_file)

    return pre_neighborhoods

if __name__ == "__main__":
    main(sys.argv)
