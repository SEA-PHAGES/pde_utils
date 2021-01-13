from math import (sqrt, floor)

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
