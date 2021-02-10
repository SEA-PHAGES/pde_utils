import queue
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
                                          build_matrix_process,
                                          verbose=verbose)

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
def greedy(matrix, threshold, is_distance=True, return_matrix=False):
    """ Clusters by interpretting the threshold as a lower bound for cluster
    inclusion.

    :param threshold: the lower bound for cluster inclusion
    """
    clusters = {}

    # Step 1: build an adjacency matrix
    adj_mat = list()
    for i in range(matrix.size):
        for j in range(i + 1, matrix.size):
            if is_distance:
                if matrix.matrix[i][j] <= threshold:
                    adj_mat.append({matrix.labels[i], matrix.labels[j]})
            else:
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
                                                        list(cluster_members))
        clusters = cluster_matricies

    return clusters


def upgma(matrix, iterations, metric="DB", is_distance=True,
          return_matrix=False):
    tree = create_upgma_tree(matrix, is_distance=is_distance)

    curr_clades = [x for x in tree.root.clades]

    clustering_scheme_map = dict()
    for i in range(iterations):
        clade_matricies = []
        clustering_scheme = {}
        for j in range(len(curr_clades)):
            clade = curr_clades[j]
            if clade.matrix is None:
                clade.matrix = clade_to_matrix(clade, matrix)

            clustering_scheme[j+1] = clade.matrix.labels

            if clade.matrix.size <= 1:
                continue
            clade_matricies.append(clade.matrix)

        if metric == "DB":
            metric_val = calculate_DB_index(matrix, clade_matricies,
                                            is_distance=is_distance)

        if metric in ["DB"]:
            if metric_val == 0:
                continue

        outliers = [x[1][0] for x in clustering_scheme.items()
                    if len(x[1]) <= 1]
        clustering_scheme[len(clustering_scheme) + 1] = outliers

        clustering_scheme_map[i] = (metric_val, clustering_scheme,
                                    len(curr_clades))

        if not split_weakest_clade(curr_clades, is_distance=is_distance):
            break

    if metric in ["DB"]:
        num_clusters, cluster_data = min(clustering_scheme_map.items(),
                                         key=lambda x: x[1][0])
    else:
        raise

    clustering_scheme = cluster_data[1]
    if return_matrix:
        cluster_matricies = {}
        for cluster, members in clustering_scheme.items():
            cluster_matricies[cluster] = matrix.get_submatrix_from_labels(
                                                                    members)

        clustering_scheme = cluster_matricies

    return clustering_scheme


def upgma_merge(matrix, submatricies, metric="DB",
                is_distance=True, return_matrix=False):
    clade_centroid_map = {}
    for submatrix in submatricies:
        subclade = create_upgma_tree(submatrix, is_distance=is_distance).root
        clade_centroid_map[submatrix.get_centroid()] = subclade

    centroid_matrix = matrix.get_submatrix_from_labels(list(
                                                clade_centroid_map.keys()))
    centroid_tree = create_upgma_tree(centroid_matrix, is_distance=is_distance)
    centroid_root = centroid_tree.root
    for leaf in centroid_root.get_terminals():
        leaf.clades = clade_centroid_map[leaf.name].clades
        leaf.branch_length = clade_centroid_map[leaf.name].branch_length

    curr_clades = [x for x in centroid_root.clades]

    clustering_scheme_map = dict()
    for i in range(len(submatricies) - 1):
        clade_matricies = []
        clustering_scheme = {}
        for j in range(len(curr_clades)):
            clade = curr_clades[j]
            if clade.matrix is None:
                clade.matrix = clade_to_matrix(clade, matrix)

            if clade.matrix.size <= 1:
                continue

            clade_matricies.append(clade.matrix)
            clustering_scheme[j+1] = clade.matrix.labels

        if metric == "DB":
            metric_val = calculate_DB_index(matrix, clade_matricies,
                                            is_distance=is_distance)
            print(metric_val)

        if metric in ["DB"]:
            if metric_val == 0:
                continue

        outliers = [x[1][0] for x in clustering_scheme.items()
                    if len(x[1]) <= 1]
        clustering_scheme[len(clustering_scheme) + 1] = outliers

        clustering_scheme_map[i] = (metric_val, clustering_scheme,
                                    len(curr_clades))

        if not split_weakest_clade(curr_clades, is_distance=is_distance):
            break

    if metric in ["DB"]:
        num_clusters, cluster_data = min(clustering_scheme_map.items(),
                                         key=lambda x: x[1][0])
    else:
        raise

    clustering_scheme = cluster_data[1]
    if return_matrix:
        cluster_matricies = {}
        for cluster, members in clustering_scheme.items():
            cluster_matricies[cluster] = matrix.get_submatrix_from_labels(
                                                                    members)

        clustering_scheme = cluster_matricies

    return clustering_scheme


def dbscan(matrix, eps, minpts, is_distance=True, return_matrix=False):
    counter = 0
    cluster_lookup = dict()
    for p_label in matrix.labels:
        if cluster_lookup.get(p_label, -1) != -1:
            continue

        neighbors = matrix.get_nearest_neighbors(p_label, eps,
                                                 is_distance=is_distance)
        if len(neighbors) < minpts:
            cluster_lookup[p_label] = None
            continue

        counter += 1
        cluster_lookup[p_label] = counter

        seeds = queue.Queue()

        for q_label in neighbors:
            seeds.put(q_label)

        while not seeds.empty():
            q_label = seeds.get()

            neighbor_cluster = cluster_lookup.get(q_label, -1)
            if neighbor_cluster is None:
                cluster_lookup[q_label] = counter

            if neighbor_cluster != -1:
                continue

            cluster_lookup[q_label] = counter
            neighbors = matrix.get_nearest_neighbors(q_label, eps,
                                                     is_distance=is_distance)
            if len(neighbors) >= minpts:
                for neighbor in neighbors:
                    seeds.put(neighbor)

    cluster_scheme = dict()
    for label, cluster in cluster_lookup.items():
        cluster_members = cluster_scheme.get(cluster, list())
        cluster_members.append(label)

        cluster_scheme[cluster] = cluster_members

    if return_matrix:
        cluster_matricies = {}
        for cluster, cluster_members in cluster_scheme.items():
            cluster_matricies[cluster] = matrix.get_submatrix_from_labels(
                                                        list(cluster_members))
        cluster_scheme = cluster_matricies

    return cluster_scheme


def lloyds(matrix, centroids, eps=None, is_distance=True, return_matrix=False):
    cent_indicies = set()
    for centroid in centroids:
        cent_indicies.add(matrix.get_index_from_label(centroid))

    while True:
        clustering_scheme = dict()
        cent_lookup = dict()
        for i, cent_index in enumerate(cent_indicies):
            clustering_scheme[(i+1)] = [matrix.get_label_from_index(
                                                                cent_index)]
            cent_lookup[cent_index] = (i + 1)

        for i in range(matrix.size):
            if i in cent_indicies:
                continue

            if is_distance:
                closest_func = min
            else:
                closest_func = max

            closest_centroid = closest_func([cent for cent in cent_indicies],
                                            key=lambda x: matrix.matrix[i][x])

            cluster = cent_lookup[closest_centroid]
            if eps is not None:
                if matrix.matrix[i][closest_centroid] < eps:
                    noise = clustering_scheme.get(None, list())
                    noise.append(matrix.get_label_from_index(i))
                    clustering_scheme[None] = noise
                    continue

            clustering_scheme[cluster].append(matrix.get_label_from_index(i))

        scheme_matricies = dict()
        new_cent_indicies = set()
        for cluster, cluster_members in clustering_scheme.items():
            if cluster is None:
                continue

            cluster_members.sort()
            submatrix = matrix.get_submatrix_from_labels(cluster_members)
            scheme_matricies[cluster] = submatrix

            new_cent_indicies.add(matrix.get_index_from_label(
                                                    submatrix.get_centroid()))

        diff_indicies = new_cent_indicies.difference(cent_indicies)
        if not diff_indicies:
            break

        cent_indicies = new_cent_indicies

    if return_matrix:
        return scheme_matricies

    return clustering_scheme


# CLUSTERING METRIC FUNCTIONS
# -----------------------------------------------------------------------------
def calculate_DB_index(matrix, submatricies, is_distance=True):
    centroid_spread_map = {}
    for submatrix in submatricies:
        centroid = submatrix.get_centroid()
        centroid_spread_map[centroid] = submatrix.get_average_edge(centroid)

    centroid_matrix = matrix.get_submatrix_from_labels(list(
                                                centroid_spread_map.keys()))

    centroid_adj_map = centroid_matrix.create_adjacency_map()

    db_index = 0
    for centroid, spread in centroid_spread_map.items():
        adj_map = centroid_adj_map[centroid]
        if not adj_map:
            continue

        if is_distance:
            nearest_centroid, dist = min(adj_map, key=lambda x: x[1])
        else:
            nearest_centroid, dist = max(adj_map, key=lambda x: x[1])

        if dist == 0:
            db_sep = 1
        else:
            db_sep = (float(spread + centroid_spread_map[nearest_centroid]) /
                      float(dist))
        db_index += db_sep

    db_index /= len(centroid_spread_map)
    return db_index


def calculate_silhouette_coeffecient(matrix, submatricies, is_distance=True):
    cent_index_submatrix_map = {}
    for submatrix in submatricies:
        centroid = submatrix.get_centroid()
        centroid_index = matrix.get_index_from_label(centroid)
        cent_index_submatrix_map[centroid_index] = submatrix

    centroid_matrix = matrix.get_submatrix_from_indicies(
                                        list(cent_index_submatrix_map.keys()))

    mean_value_A_maps = dict()
    for centroid_index, query_submatrix in cent_index_submatrix_map.items():
        mean_value_A_map = dict()
        for member in query_submatrix.labels():
            a = query_submatrix.get_average_edge()
            mean_value_A_map[member] = a

        mean_value_A_maps[centroid_index] = mean_value_A_map

    SC = 0
    s_counter = 0
    s_coeffecient_maps = dict()
    for centroid in centroid_matrix.labels:
        s_coeffecient_map = dict()

        cent_index = matrix.get_index_from_label(centroid)
        centroid_distances = centroid_matrix.get_row(cent_index,
                                                     exclude_diagonal=True)

        i_matrix = cent_index_submatrix_map[cent_index]

        nearest_centroid = centroid_matrix.get_label_from_index(
                                    min(range(len(centroid_distances)),
                                        key=lambda x: centroid_distances[x]))
        nearest_cent_index = matrix.get_label_from_index(nearest_centroid)

        j_matrix = cent_index_submatrix_map[nearest_cent_index]
        j_matrix_indicies = list()
        for j_member in j_matrix.labels:
            j_matrix_indicies.append(matrix.get_index_from_label(j_member))

        mean_value_A_map = mean_value_A_maps[cent_index]
        for i_member in i_matrix:
            i_index = matrix.get_index_from_label(i_member)

            b = 0
            for j_index in j_matrix_indicies:
                b += matrix.matrix[i_index][j_index]

            b /= len(j_matrix_indicies)
            a = mean_value_A_map[i_member]

            s_divisor = max([b, a])
            if is_distance:
                s_dividend = b - a
            else:
                s_dividend = a - b

            s = s_dividend / s_divisor
            s_coeffecient_map[i_member] = s
            s_counter += 1
            SC += s

        s_coeffecient_maps[centroid] = s_coeffecient_map

    SC /= s_counter
    return SC, s_coeffecient_maps


def calculate_error_sum_of_squares(matrix, submatricies, is_distance=True,
                                   max_value=1):
    if len(submatricies) < 1:
        return 0

    total_SSE = 0
    for submatrix in submatricies:
        centroid = submatrix.get_centroid()

        SSE = 0
        row = submatrix.get_row(centroid, exclude_diagonal=True)

        if len(row) > 1:
            for value in row:
                if is_distance:
                    SSE += value ** 2
                else:
                    SSE += (max_value - value) ** 2

        total_SSE += SSE

    total_SSE /= len(submatricies)
    return total_SSE
    pass


# UPGMA HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def create_upgma_tree(matrix, is_distance=True):
    adj_map = matrix.create_adjacency_map()
    closest_pairs = create_closest_pairs(adj_map, is_distance=is_distance)
    clade_map = create_clade_map(adj_map)

    for i in range(matrix.size-2):
        if is_distance:
            source, pair_edge = min(closest_pairs.items(), key=lambda x: x[1])
        else:
            source, pair_edge = max(closest_pairs.items(), key=lambda x: x[1])
        merge_closest_edge(adj_map, clade_map, closest_pairs,
                           (source, pair_edge[0]), pair_edge[1],
                           is_distance=is_distance)

    unmerged_clusters = list(clade_map.keys())
    unmerged_clades = list(clade_map.values())
    if len(unmerged_clusters) > 1:
        branch_length = adj_map[unmerged_clusters[0]][unmerged_clusters[1]]

        root = BaseTree.Clade(branch_length=branch_length,
                              clades=unmerged_clades)
    else:
        root = unmerged_clades[0]

    root.matrix = matrix
    tree = BaseTree.Tree(root=root, rooted=False)

    return tree


def create_closest_pairs(adj_map, is_distance=True):
    closest_pairs = dict()
    temp_adj_map = dict()
    for cluster, adj_list in adj_map.items():
        if not adj_list:
            continue

        adj_dict = dict()

        closest = adj_list[0][0]
        closest_value = adj_list[0][1]
        for link, value in adj_list:
            adj_dict[link] = value

            if is_distance:
                if value < closest_value:
                    closest = link
                    closest_value = value
            else:
                if value > closest_value:
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
                       pair_value, is_distance=True):
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
            if is_distance:
                closest_pairs[cluster] = min(adj_dict.items(),
                                             key=lambda x: x[1])
            else:
                closest_pairs[cluster] = max(adj_dict.items(),
                                             key=lambda x: x[1])

    adj_map.pop(merging)
    adj_map[source] = source_dict

    closest_pairs.pop(merging)
    if is_distance:
        closest_pairs[source] = min(source_dict.items(), key=lambda x: x[1])
    else:
        closest_pairs[source] = max(source_dict.items(), key=lambda x: x[1])

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


def split_weakest_clade(clades, is_distance=True):
    weakest_clade = None
    for clade in clades:
        if clade.is_terminal():
            continue

        if weakest_clade is None:
            weakest_clade = clade

        if is_distance:
            if clade.branch_length > weakest_clade.branch_length:
                weakest_clade = clade
        else:
            if clade.branch_length < weakest_clade.branch_length:
                weakest_clade = clade

    if weakest_clade is None:
        return False

    clades.remove(weakest_clade)

    for subclade in weakest_clade.clades:
        clades.append(subclade)

    return True
