from sqlalchemy import select

from pdm_utils.functions import (querying, parallelize)

from pde_utils.classes import (matrix)
from pde_utils.functions import alignment


def get_distinct_phams_from_organism(alchemist, organism_id):
    gene_obj = alchemist.metadata.tables["gene"]

    phageid_obj = gene_obj.c.PhageID
    phamid_obj = gene_obj.c.PhamID

    phams_query = select([phamid_obj]).where(
                                        phageid_obj == organism_id).distinct()
    phams = querying.first_column(alchemist.engine, phams_query)
    return phams


def calculate_gcs(query_phams, target_phams):
    shared_phams = set(query_phams).intersection(set(target_phams))

    query_spp = len(shared_phams) / len(query_phams)
    target_spp = len(shared_phams) / len(target_phams)

    return (query_spp + target_spp) / 2


def calculate_gcs_matrix(alchemist, organism_ids, cores=1, verbose=False):
    organisms_gc = []
    for organism_id in organism_ids:
        organisms_gc.append(
                    get_distinct_phams_from_organism(alchemist, organism_id))

    work_items = []
    for i in range(len(organisms_gc)):
        work_items.append((organisms_gc[i:]))

    matricies = parallelize.parallelize(
                            work_items, cores, calculate_gcs_matrix_process,
                            verbose=verbose)

    matricies.sort(key=lambda x: len(x), reverse=True)
    return matricies


def calculate_gcs_symmetric_matrix(alchemist, organism_ids, cores=1,
                                   verbose=False):
    half_matrix = calculate_gcs_matrix(alchemist, organism_ids, cores=cores,
                                       verbose=verbose)

    return half_to_symmetric_matrix(half_matrix, organism_ids)


def half_to_symmetric_matrix(half_matrix, labels):
    symmetric_matrix = matrix.SymmetricMatrix(labels)

    for i in range(len(half_matrix)):
        for j in range(1, len(half_matrix[i])):
            symmetric_matrix.fill_cell(i, i + j, half_matrix[i][j])

    return symmetric_matrix


def calculate_gcs_matrix_process(gc_data):
    if not gc_data:
        return []

    matrix_row = []
    target_gc = gc_data[0]
    for query_gc in gc_data:
        matrix_row.append(calculate_gcs(query_gc, target_gc))

    return matrix_row


def calculate_ani_matrix(sketch_paths, cores=1, verbose=False):
    work_items = []
    for i in range(len(sketch_paths)):
        work_items.append((sketch_paths[i:]))

    matricies = parallelize.parallelize(
                        work_items, cores, calculate_ani_matrix_process,
                        verbose=verbose)
    matricies.sort(key=lambda x: len(x), reverse=True)
    return matricies


def calculate_ani_symmetric_matrix(sketch_paths, cores=1, verbose=False):
    half_matrix = calculate_ani_matrix(sketch_paths, cores=1, verbose=verbose)

    return half_to_symmetric_matrix(
                        half_matrix,
                        [path.with_suffix("").name for path in sketch_paths])


def calculate_ani_matrix_process(sketch_paths):
    if not sketch_paths:
        return []

    matrix_row = []
    target_path = sketch_paths[0]
    for query_path in sketch_paths:
        ani_data = alignment.mash_dist(query_path, target_path)
        ani_data = ani_data.split("\t")

        if len(ani_data) == 5:
            ani = 1 - float(ani_data[2])
            matrix_row.append(ani)

    return matrix_row
