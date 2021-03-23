from Bio import Phylo
from networkx import MultiDiGraph
from pdm_utils import AlchemyHandler
from pdm_utils.functions import querying
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from pde_utils.classes import clustal
from pde_utils.classes.pan_models import (Base, Cluster)


# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
PAN_GRAPH_EDGEWEIGHTS = ["CentroidDistance", "DBSeparation", "MinDistance",
                         "Probability", "Expect", "HMMAlnCols",
                         "HMMAlnQueryStart", "HMMAlnQueryEnd",
                         "HMMAlnHitStart", "HMMAlnHitEnd"]


def build_pan_engine(pan_path):
    engine_string = f"sqlite+pysqlite:///{pan_path}"
    engine = create_engine(engine_string)

    return engine


def build_pan_session(engine):
    session_maker_obj = sessionmaker(bind=engine)
    session = session_maker_obj()
    return session


def build_pan(pan_path):
    engine = build_pan_engine(pan_path)

    metadata = Base.metadata
    metadata.create_all(engine)

    session = build_pan_session(engine)

    alchemist = AlchemyHandler()
    alchemist._engine = engine
    alchemist._metadata = metadata
    alchemist._session = session

    alchemist.connected = True
    alchemist.has_credentials = True
    alchemist.connected_database = True
    alchemist.has_database = True

    return alchemist


def to_networkx(alchemist):
    pan_graph = MultiDiGraph()

    pan_clusters = alchemist.session.query(Cluster).all()
    for cluster in pan_clusters:
        pan_graph.add_node(cluster.ClusterID, Spread=cluster.Spread,
                           CentroidID=cluster.CentroidID,
                           CentroidSeq=cluster.CentroidSeq.decode("utf-8"))

    for cluster in pan_clusters:
        for identity_edge in cluster.NeighborhoodEdges:
            pan_graph.add_edge(identity_edge.Source, identity_edge.Target,
                               DBSeparation=identity_edge.DBSeparation,
                               CentroidDistance=identity_edge.CentroidDistance,
                               MinDistance=identity_edge.MinDistance,
                               Probability="", Expect="", HMMAlnCols="",
                               HMMAlnQueryStart="", HMMAlnQueryEnd="",
                               HMMAlnHitStart="", HMMAlnHitEnd="")

        for hmm_edge in cluster.TownEdges:
            pan_graph.add_edge(hmm_edge.Source, hmm_edge.Target,
                               DBSeparation="", CentroidDistance="",
                               MinDistance="",
                               Probability=hmm_edge.Probability,
                               Expect=hmm_edge.Expect,
                               HMMAlnCols=hmm_edge.AlignedCols,
                               HMMAlnQueryStart=hmm_edge.SourceStart,
                               HMMAlnQueryEnd=hmm_edge.SourceEnd,
                               HMMAlnHitStart=hmm_edge.TargetStart,
                               HMMAlnHitEnd=hmm_edge.TargetEnd)

    return pan_graph


def retrieve_cluster_data(pan_alchemist, cluster_ids):
    cluster_table = Cluster.__table__

    query = querying.build_select(
                        pan_alchemist.graph,
                        [cluster_table.c.Spread, cluster_table.c.CentroidID,
                         cluster_table.c.CentroidSeq,
                         cluster_table.c.ClusterID])
    results = querying.execute(pan_alchemist.engine, query,
                               in_column=cluster_table.c.ClusterID,
                               values=cluster_ids)

    return results


def parse_cluster(cluster_id, data_dict=None, MSA_path=None, PIM_path=None,
                  GT_path=None):
    aln = None
    if MSA_path is not None:
        if MSA_path.is_file():
            aln = clustal.MultipleSequenceAlignment(MSA_path, fmt="fasta")
            aln.parse_alignment()

    mat = None
    tree = None
    if PIM_path is not None:
        if PIM_path.is_file():
            mat = clustal.PercentIdentityMatrix(PIM_path)
            mat.parse_matrix()

    if GT_path is not None:
        if GT_path.is_file():
            tree = Phylo.read(str(GT_path), "newick")

    centroid_id = None
    centroid_seq = None
    spread = None
    if data_dict is not None:
        centroid_id = data_dict.get("CentroidID")
        centroid_seq = data_dict.get("CentroidSeq")
        spread = data_dict.get("Spread")

    return Cluster(ClusterID=cluster_id, CentroidID=centroid_id,
                   CentroidSeq=centroid_seq, Spread=spread,
                   MSA=aln, PIM=mat, GT=tree)
