from math import (sqrt, floor)

from Bio.Phylo import BaseTree
from numpy import zeros
from pdm_utils.functions import (basic, parallelize)

from pde_utils.classes.matricies import SymmetricMatrix


# MATRIX-HANDLING FUNCTIONS
# -----------------------------------------------------------------------------
def build_symmetric_matrix(nodes, distance_function, names=None,
                           cores=1, verbose=False):
    work_items = []

    row_indicies = [i for i in range(len(nodes))]
    chunk_size = int(floor(sqrt(len(nodes))))
    for i in row_indicies:
        subject = nodes[i]
        if len(nodes) - 1 == i:
            work_items.append((distance_function, subject, [], i, 0))
        else:
            query_node_chunks = basic.partition_list(nodes[i+1:], chunk_size)
            for j in range(len(query_node_chunks)):
                work_items.append((distance_function,
                                   subject, query_node_chunks[j], i, j))

    matrix_data = parallelize.parallelize(work_items, cores,
                                          build_matrix_process)

    matrix_data.sort(key=lambda x: (x[1], x[2]))

    if names is None:
        names = row_indicies
    else:
        if len(names) != len(row_indicies):
            names = row_indicies
    matrix = SymmetricMatrix(names)

    for data in matrix_data:
        for i in range(len(data[0])):
            col = (data[2] * chunk_size) + (data[1] + i + 1)
            matrix.fill_cell(data[1], col, data[0][i])

    matrix.fill_diagonal(1)
    return matrix


def build_matrix_process(distance_function, subject, queries, row,
                         col_chunk):
    row_data = zeros(len(queries))

    for i in range(len(queries)):
        row_data[i] = distance_function(subject, queries[i])

    return (row_data, row, col_chunk)


# CLUSTERING ALGORITHMS
# -----------------------------------------------------------------------------
def greedy(matrix, threshold, return_matrix=False):
    """ Clusters by interpretting the threshold as a lower bound for cluster
    inclusion.

    :param threshold: the lower bound for cluster inclusion
    """
    clusters = {}

    # Step 1: build an adjacency matrix
    adj_mat = list()
    for i in range(matrix.size):
        for j in range(i + 1, matrix.size):
            if matrix.matrix[i][j] >= threshold:
                adj_mat.append({matrix.labels[i], matrix.labels[j]})
    # Step 2: build connected components from adjacency matrix
    while len(adj_mat) > 0:
        anchor = adj_mat.pop(0)
        cleanup = list()
        for i, edge in enumerate(adj_mat):
            if anchor.intersection(edge) > set():
                anchor = anchor.union(edge)
                cleanup.append(i)
        # Check whether this component now overlaps an existing one
        merged = False
        for key, component in clusters.items():
            if component.intersection(anchor) > set():
                clusters[key] = component.union(anchor)
                merged = True
                break
        if not merged:
            clusters[len(clusters) + 1] = anchor
        for i in reversed(cleanup):
            adj_mat.pop(i)
    # Step 3: add any individual nodes that are not connected to others
    used = set()
    for key, value in clusters.items():
        used = used.union(value)
    missing = set(matrix.labels).difference(used)
    for node in missing:
        clusters[len(clusters) + 1] = {node}
    # Step 4: If desired, create submatricies for each set of clustered labels
    if return_matrix:
        cluster_matricies = {}
        for cluster, cluster_members in clusters.items():
            cluster_matricies[cluster] = matrix.get_submatrix_from_labels(
                                                        cluster_members)

    return clusters


def upgma(matrix, iterations, metric="DB"):
    tree = create_upgma_tree(matrix)

    curr_clades = [x for x in tree.root.clades]

    clustering_iters_map = dict()
    for _ in range(iterations):
        clade_matricies = []
        clustering_scheme = {}
        for i in range(len(curr_clades)):
            clade = curr_clades[i]
            if clade.matrix is None:
                clade.matrix = clade_to_matrix(clade, matrix)

            clade_matricies.append(clade.matrix)
            clustering_scheme[i+1] = clade.matrix.labels

        if metric == "DB":
            metric = calculate_DB_index(matrix, clade_matricies)

        clustering_iters_map[len(clade_matricies)] = (metric,
                                                      clustering_scheme)

        if not split_weakest_clade(curr_clades):
            break

    num_clusters, cluster_data = min(clustering_iters_map.items(),
                                     key=lambda x: x[1][0])

    return cluster_data[1]


# CLUSTERING METRIC FUNCTIONS
# -----------------------------------------------------------------------------
def calculate_DB_index(matrix, submatricies):
    centroid_spread_map = {}
    for submatrix in submatricies:
        centroid = submatrix.get_centroid()
        centroid_spread_map[centroid] = submatrix.get_average_value(centroid)

    centroid_matrix = matrix.get_submatrix_from_labels(list(
                                                centroid_spread_map.keys()))

    centroid_adj_map = centroid_matrix.create_adjacency_map()

    db_index = 0
    for centroid, spread in centroid_spread_map.items():
        nearest_centroid, dist = min(centroid_adj_map[centroid],
                                     key=lambda x: x[1])

        db_sep = (float(spread + centroid_spread_map[nearest_centroid]) /
                  float(dist))
        db_index += db_sep

    return db_index / len(centroid_spread_map)


# UPGMA HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def create_upgma_tree(matrix):
    adj_map = matrix.create_adjacency_map()
    closest_pairs = create_closest_pairs(adj_map)
    clade_map = create_clade_map(adj_map)

    for i in range(matrix.size-2):
        source, pair_edge = min(closest_pairs.items(), key=lambda x: x[1])
        merge_closest_edge(adj_map, clade_map, closest_pairs,
                           (source, pair_edge[0]), pair_edge[1])

    unmerged_clusters = list(clade_map.keys())
    unmerged_clades = list(clade_map.values())
    root = BaseTree.Clade(branch_length=(
                                adj_map[unmerged_clusters[0]][
                                                unmerged_clusters[1]]),
                          clades=unmerged_clades)
    root.matrix = matrix
    tree = BaseTree.Tree(root=root, rooted=False)

    return tree


def create_closest_pairs(adj_map):
    closest_pairs = dict()
    temp_adj_map = dict()
    for cluster, adj_list in adj_map.items():
        adj_dict = dict()

        closest = adj_list[0][0]
        closest_value = adj_list[0][1]
        for link, value in adj_list:
            adj_dict[link] = value

            if value < closest_value:
                closest = link
                closest_value = value

        closest_pairs[cluster] = (closest, value)
        temp_adj_map[cluster] = adj_dict

    for cluster, adj_dict in temp_adj_map.items():
        adj_map[cluster] = adj_dict

    return closest_pairs


def create_clade_map(adj_map):
    clade_map = dict.fromkeys(adj_map.keys())

    for cluster in adj_map.keys():
        clade = BaseTree.Clade(name=str(cluster))
        clade.matrix = None
        clade_map[cluster] = clade

    return clade_map


def merge_closest_edge(adj_map, clade_map, closest_pairs, closest_pair,
                       pair_value):
    source = closest_pair[0]
    merging = closest_pair[1]
    if closest_pair[1] < source:
        source = closest_pair[1]
        merging = closest_pair[0]

    source_dict = dict()

    for cluster, adj_dict in adj_map.items():
        if cluster in closest_pair:
            continue

        adj_dict[source] = (adj_dict[source] + adj_dict[merging]) / 2
        adj_dict.pop(merging)

        source_dict[cluster] = adj_dict[source]

        ccp = (cluster, closest_pairs[cluster][0])
        if source in ccp or merging in ccp:
            closest_pairs[cluster] = min(adj_dict.items(), key=lambda x: x[1])

    adj_map.pop(merging)
    adj_map[source] = source_dict

    closest_pairs.pop(merging)
    closest_pairs[source] = min(source_dict.items(), key=lambda x: x[1])

    merging_clade = clade_map.pop(merging)
    clade_map[source] = BaseTree.Clade(branch_length=pair_value,
                                       clades=[clade_map[source],
                                               merging_clade])
    clade_map[source].matrix = None


# DISTANCE-BASED TREE FUNCTIONS
# -----------------------------------------------------------------------------
def clade_to_matrix(clade, test_matrix):
    leaves = clade.get_terminals()
    leaf_names = [clade.name for clade in leaves]

    return test_matrix.get_submatrix_from_labels(leaf_names)


def split_weakest_clade(clades):
    weakest_clade = None
    for clade in clades:
        if clade.is_terminal():
            continue

        if weakest_clade is None:
            weakest_clade = clade

        if clade.branch_length > weakest_clade.branch_length:
            weakest_clade = clade

    if weakest_clade is None:
        return False

    clades.remove(weakest_clade)

    for subclade in weakest_clade.clades:
        clades.append(subclade)

    return True
