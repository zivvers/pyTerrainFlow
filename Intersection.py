from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Tuple
import itertools

class FeatureType(Enum):
    VERTEX = auto()
    EDGE = auto()

class Graticule(Enum):
    DIAG = 0 
    LON = 1
    LAT = 2

@dataclass
class Intersection:

    #point: Tuple[float, float, float]
    point: Tuple[float, float]

    mesh_feature: FeatureType
    clip_feature: FeatureType

    mesh_edge_index: Optional[int] = None
    mesh_edge_order: Optional[int] = None 
    mesh_vertex_index: Optional[int] = None

    letter : Optional = None

    graticule_type: Optional[Graticule] = None
    grat_pix0: Optional[Tuple[int, int]] = None
    grat_pix1: Optional[Tuple[int, int]] = None

    # For clip polygon features we always
    # know that CW is ClipVerts[clip_poly_id(n, m), clip_poly_id(n+1, m)]
    # so we index the Clip Edge/Vertex Intersection the same?
    clip_poly_id: Optional[int] = None
    clip_poly_order: Optional[int] = None

    mesh_t: Optional[int]  = None
    clip_t: Optional[int]  = None

    # wtf Useful for filtering internal decomposition edges.
    #is_union_boundary: bool = True

    @property
    def kind(self):
        return (
            self.mesh_feature.feature_type,
            self.clip_feature.feature_type,
        )

    @property
    def is_mesh_vertex_on_clip_edge(self):
        return self.kind == (
            FeatureType.VERTEX,
            FeatureType.EDGE,
        )

    @property
    def is_mesh_edge_on_clip_vertex(self):
        return self.kind == (
            FeatureType.EDGE,
            FeatureType.VERTEX,
        )

    @property
    def is_vertex_vertex(self):
        return self.kind == (
            FeatureType.VERTEX,
            FeatureType.VERTEX,
        )

    @property
    def is_edge_edge(self):
        return self.kind == (
            FeatureType.EDGE,
            FeatureType.EDGE,
        )