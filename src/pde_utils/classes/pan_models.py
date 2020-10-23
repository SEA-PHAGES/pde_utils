from Bio.Phylo import BaseTree
from sqlalchemy import (Column, ForeignKey, Integer, String)
from sqlalchemy.types import (BLOB, Float)
from sqlalchemy.ext.declarative import declarative_base
# from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from pde_utils.classes import clustal

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
Base = declarative_base()

NEIGHBORHOOD_EDGES_JOIN = "cluster.c.ClusterID == identity_edge.c.Source"
SOURCE_NODE_JOIN = "identity_edge.c.Source == cluster.c.ClusterID"
TARGET_NODE_JOIN = "identity_edge.c.Target == cluster.c.ClusterID"


class Cluster(Base):
    """Class to hold data about a sequence cluster."""
    __tablename__ = "cluster"

    ClusterID = Column("ClusterID", String(15), primary_key=True)
    Spread = Column("Spread", Float)
    CentroidID = Column("CentroidID", String(35))
    CentroidSeq = Column("CentroidSeq", BLOB)

    NeighborhoodEdges = relationship("IdentityEdge",
                                     back_populates="SourceNode",
                                     primaryjoin=NEIGHBORHOOD_EDGES_JOIN)

    def __init__(self, ClusterID=None, Spread=None, CentroidID=None,
                 CentroidSeq=None, PIM=None, MSA=None, GT=None):
        self.ClusterID = ClusterID
        self.Spread = Spread
        self.CentroidID = CentroidID
        self.CentroidSeq = CentroidSeq

        if CentroidSeq is not None:
            self._centroid_seq_str = CentroidSeq.decode("utf-8")
        else:
            self._centroid_seq_str = None

        self._centroid_seq_str = None

        if PIM is not None:
            if not isinstance(PIM, clustal.PercentIdentityMatrix):
                raise TypeError(
                        "Cluster object expects a PercentIdentityMatrix "
                        f"object, not {type(PIM)}.")

        if MSA is not None:
            if not isinstance(MSA, clustal.MultipleSequenceAlignment):
                raise TypeError(
                        "Cluster object expects a MultipleSequenceAlignment "
                        f"object, not {type(MSA)}.")

        if GT is not None:
            if not isinstance(GT, BaseTree.TreeElement):
                raise TypeError(
                        "Cluster object expects a Biopython TreeElement "
                        f"object, not {type(GT)}.")

        self.PIM = PIM
        self.MSA = MSA
        self.GT = GT

    @property
    def centroid_seq_str(self):
        if self._centroid_seq_str is None:
            if isinstance(self.CentroidSeq, bytes):
                self._centroid_seq_str = self.CentroidSeq.decode("utf-8")

        seq_str = self._centroid_seq_str
        return seq_str

    @centroid_seq_str.setter
    def centroid_seq_str(self, seq_str):
        if isinstance(seq_str, str):
            self.CentroidSeq = seq_str.encode("utf-8")

        self._centroid_seq_str = seq_str

    def parse_centroid(self):
        if self.MSA is None:
            raise AttributeError(
                        "Cluster object requires a MultipleSequenceAlignment "
                        "to recover a centroid from")
        self.MSA.check_initialization("cluster.parse_centroid")

        if self.PIM is None:
            centroid_id, centroid_seq = self.MSA.longest_gene()
        else:
            self.PIM.check_initialization("cluster.parse_centroid")

            centroid_id = self.PIM.get_centroid()
            centroid_seq = self.MSA.get_sequence(centroid_id, gaps=False)

        self.CentroidID = centroid_id
        self.CentroidSeq = centroid_seq.encode("utf-8")
        self.centroid_seq_str = centroid_seq

    def parse_spread(self):
        if self.PIM is None:
            self.Spread = 0
            return

        if self.CentroidID is None:
            self.parse_centroid_id()

        self.PIM.check_initialization("cluster.parse_spread")
        self.Spread = 100 - self.PIM.get_average_identity(self.CentroidID)


class IdentityEdge(Base):
    """Class to hold data about the identity relationship between clusters."""
    __tablename__ = "identity_edge"

    ID = Column("ID", Integer, primary_key=True)
    Source = Column("Source", String(15), ForeignKey("cluster.ClusterID"),
                    nullable=False)
    Target = Column("Target", String(15), ForeignKey("cluster.ClusterID"),
                    nullable=False)
    DBSeparation = Column("DBSeparation", Float)
    CentroidDistance = Column("CentroidDistance", Float)
    MinDistance = Column("MinDistance", Float)

    SourceNode = relationship("Cluster", back_populates="NeighborhoodEdges",
                              primaryjoin=SOURCE_NODE_JOIN)
    TargetNode = relationship("Cluster", primaryjoin=TARGET_NODE_JOIN)

    def __init__(self, ID=None, Source=None, Target=None, DBSeparation=None,
                 CentroidDistance=None, MinDistance=None):
        self.ID = ID
        self.Source = Source
        self.Target = Target
        self.DBSeparation = DBSeparation
        self.CentroidDistance = CentroidDistance
        self.MinDistance = MinDistance

        self.Source_node = None
        self.Target_node = None
