from networkx import MultiDiGraph
from pdm_utils import AlchemyHandler
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker

from pde_utils.classes import pan_models


def build_pan_engine(pan_path):
    engine_string = f"sqlite+pysqlite:///{pan_path}"
    engine = create_engine(engine_string)

    return engine


def build_pan_metadata(engine=None):
    metadata = MetaData(bind=engine)

    pan_models.map_pan_models(metadata)

    return metadata


def build_pan_session(engine):
    session_maker_obj = sessionmaker(bind=engine)
    session = session_maker_obj()
    return session


def build_pan(pan_path):
    engine = build_pan_engine(pan_path)
    metadata = build_pan_metadata(engine=engine)

    metadata.create_all()

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

    pan_clusters = alchemist.session.query(pan_models.Cluster).all()
    for cluster in pan_clusters:
        pan_graph.add_node(cluster.ClusterID, Spread=cluster.Spread,
                           CentroidID=cluster.CentroidID,
                           CentroidSeq=cluster.CentroidSeq.decode("utf-8"))

    for cluster in pan_clusters:
        for identity_edge in cluster.NeighborhoodEdges:
            pan_graph.add_edge(identity_edge.Source, identity_edge.Target,
                               DBSeparation=identity_edge.DBSeparation,
                               CentroidIdentity=identity_edge.CentroidIdentity,
                               MinIdentity=identity_edge.MinIdentity)

    return pan_graph
