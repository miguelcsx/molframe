from typing import final
from numpy import float32, int64
from numpy.typing import NDArray

@final
class NodeLevel:
    Atoms: NodeLevel
    Residues: NodeLevel

@final
class EdgeDirection:
    Undirected: EdgeDirection
    Symmetric: EdgeDirection
    Directed: EdgeDirection

@final
class NodeFeature:
    Element: NodeFeature
    FormalCharge: NodeFeature
    PartialCharge: NodeFeature
    BFactor: NodeFeature
    Occupancy: NodeFeature
    AtomCount: NodeFeature
    PositionX: NodeFeature
    PositionY: NodeFeature
    PositionZ: NodeFeature

@final
class EdgeFeature:
    Distance: EdgeFeature
    BondOrder: EdgeFeature

@final
class SpatialBackend:
    BruteForce: SpatialBackend
    CellList: SpatialBackend
    KdTree: SpatialBackend
    NeighborList: SpatialBackend
    Auto: SpatialBackend

@final
class EdgeKind:
    @staticmethod
    def bonds() -> EdgeKind: ...
    @staticmethod
    def contacts(cutoff: float) -> EdgeKind: ...
    @staticmethod
    def radius(cutoff: float) -> EdgeKind: ...
    @staticmethod
    def k_nearest(neighbors: int) -> EdgeKind: ...

@final
class MissingFeaturePolicy:
    @staticmethod
    def error() -> MissingFeaturePolicy: ...
    @staticmethod
    def fill(value: float) -> MissingFeaturePolicy: ...

@final
class GraphOptions:
    def __init__(
        self,
        nodes: NodeLevel,
        edges: EdgeKind,
        *,
        direction: EdgeDirection = EdgeDirection.Symmetric,
        features: tuple[list[NodeFeature], list[EdgeFeature]] = ([], []),
        missing: MissingFeaturePolicy | None = None,
        backend: SpatialBackend = SpatialBackend.Auto,
        periodic: bool = False,
    ) -> None: ...

@final
class Graph:
    @property
    def node_count(self) -> int: ...
    @property
    def edge_count(self) -> int: ...
    @property
    def edge_index(self) -> NDArray[int64]: ...
    @property
    def node_features(self) -> NDArray[float32]: ...
    @property
    def edge_features(self) -> NDArray[float32]: ...
    @property
    def node_feature_schema(self) -> list[NodeFeature]: ...
    @property
    def edge_feature_schema(self) -> list[EdgeFeature]: ...
    cost: object
    def to_pyg(self) -> object: ...
    def to_dgl(self) -> object: ...
