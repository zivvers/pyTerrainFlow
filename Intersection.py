from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Tuple
import itertools

class FeatureType(Enum):
    VERTEX = auto()
    EDGE = auto()

class GraticuleType(Enum):
    DIAG = 0 
    LON = 1
    LAT = 2

@dataclass
class Graticule:
    grat_type: GraticuleType
    e0: Tuple[int, int]
    e1: Tuple[int, int] 
    id: Optional[int] = None

    def calc_id(self, grat_type, _e0, _e1):

        x0,y0 = _e0
        x1,y1 = _e1

        print(f"pix1: {_e0}, pix2: {_e1}")
        
        match grat_type:
            case GraticuleType.DIAG:
                assert x0-y0 == x1-y1
                return x0 - y0
            case GraticuleType.LON:
                assert x0 == x1
                return x0
            case GraticuleType.LAT:
                assert y0 == y1
                return y0


    def __post_init__(self):
        self.id = self.calc_id(self.grat_type, self.e0, self.e1)





@dataclass
class Intersection:

    #
    # Mandatory params
    point: Tuple[float, float]
    mesh_feature: FeatureType
    clip_feature: FeatureType

    # Everything else
    orig_point: Optional[Tuple[float, float]] = None
    mesh_edge_index: Optional[int] = None
    mesh_edge_order: Optional[int] = None 
    mesh_vertex_index: Optional[int] = None

    letter : Optional[str] = None

    graticule: Optional[Graticule] = None
    grat_pix0: Optional[Tuple[int, int]] = None
    grat_pix1: Optional[Tuple[int, int]] = None

    # For clip polygon features we always
    # know that CW is ClipVerts[clip_poly_id(n, m), clip_poly_id(n+1, m)]
    # so we index the Clip Edge/Vertex Intersection the same?
    clip_poly_id: Optional[int] = None
    clip_poly_order: Optional[int] = None

    mesh_t: Optional[int]  = None
    clip_t: Optional[int]  = None

    degenerate: Optional[bool] = None
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