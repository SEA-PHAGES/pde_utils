from sqlalchemy import select

from pdm_utils.functions import querying


def get_distinct_phams_from_organism(alchemist, organism_id):
    gene_obj = alchemist.metadata.tables["gene"]

    phageid_obj = gene_obj.c.PhageID
    phamid_obj = gene_obj.c.PhamID

    phams_query = select([phamid_obj]).where(
                                        phageid_obj == organism_id).distinct()
    phams = querying.first_column(alchemist.engine, phams_query)
    return phams
