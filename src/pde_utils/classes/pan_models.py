from sqlalchemy import (Column, ForeignKey, Integer, String, Table)
from sqlalchemy.types import (BLOB, Float)
from sqlalchemy.orm import (mapper, relationship)


class Cluster:
    def __init__(self, ClusterID, Spread=None, CentroidID=None,
                 CentroidSeq=None, PercentIdentityMatrix=None,
                 MultipleSequenceAlignment=None, AlignmentGuidetree=None):
        self.ClusterID = ClusterID
        self.Spread = Spread
        self.CentroidID = CentroidID
        self.CentroidSeq = CentroidSeq

        self.NeighborhoodEdges = []

        self.PercentIdentityMatrix = PercentIdentityMatrix
        self.MultipleSequenceAlignment = MultipleSequenceAlignment
        self.AlignmentGuidetree = AlignmentGuidetree


class IdentityEdge:
    def __init__(self, ID=None, Source=None, Target=None, DBSeparation=None,
                 CentroidIdentity=None, MinIdentity=None):
        self.ID = ID
        self.Source = Source
        self.Target = Target
        self.DBSeparation = DBSeparation
        self.CentroidIdentity = CentroidIdentity
        self.MinIdentity = MinIdentity

        self.Source_node = None
        self.Target_node = None


def map_pan_models(metadata):
    cluster = Table("cluster", metadata,
                    Column("ClusterID", String(15), primary_key=True),
                    Column("Spread", Float),
                    Column("CentroidID", String(35)),
                    Column("CentroidSeq", BLOB))

    identity_edge = Table("identity_edge", metadata,
                          Column("ID", Integer, primary_key=True),
                          Column("Source", String(15),
                                 ForeignKey("cluster.ClusterID"),
                                 nullable=False),
                          Column("Target", String(15),
                                 ForeignKey("cluster.ClusterID"),
                                 nullable=False),
                          Column("DBSeparation", Float),
                          Column("CentroidIdentity", Float),
                          Column("MinIdentity", Float),
                          sqlite_autoincrement=True)

    neighborhood_edges_join = (cluster.c.ClusterID == identity_edge.c.Source)
    cluster_properties = {"ClusterID": cluster.c.ClusterID,
                          "Spread": cluster.c.Spread,
                          "CentroidID": cluster.c.CentroidID,
                          "CentroidSeq": cluster.c.CentroidSeq,
                          "NeighborhoodEdges": relationship(
                                          IdentityEdge,
                                          back_populates="SourceNode",
                                          primaryjoin=neighborhood_edges_join)}

    source_node_join = (identity_edge.c.Source == cluster.c.ClusterID)
    target_node_join = (identity_edge.c.Target == cluster.c.ClusterID)
    identity_edge_properties = {
                        "ID": identity_edge.c.ID,
                        "Source": identity_edge.c.Source,
                        "Target": identity_edge.c.Target,
                        "DBSeparation": identity_edge.c.DBSeparation,
                        "CentroidIdentity": identity_edge.c.CentroidIdentity,
                        "MinIdentity": identity_edge.c.MinIdentity,
                        "SourceNode": relationship(
                                          Cluster,
                                          back_populates="NeighborhoodEdges",
                                          primaryjoin=source_node_join),
                        "TargetNode": relationship(
                                          Cluster,
                                          primaryjoin=target_node_join)}

    mapper(Cluster, cluster, properties=cluster_properties)
    mapper(IdentityEdge, identity_edge, properties=identity_edge_properties)
