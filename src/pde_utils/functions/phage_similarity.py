from sqlalchemy import select

from pdm_utils.functions import (querying, parallelize)


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


def calculate_gcs_matrix_process(gc_data):
    if not gc_data:
        return []

    matrix_row = []
    target_gc = gc_data[0]
    for query_gc in gc_data:
        matrix_row.append(calculate_gcs(query_gc, target_gc))

    return matrix_row
