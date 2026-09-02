    
from config import load_config

import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Polygon
from matplotlib.collections import LineCollection
from pyglm import glm
from collections import defaultdict
from affine import Affine
from MeshComponent import MeshComponent, Triangle

from typing import Optional, Tuple
import math

from MeshComponentCreator import MeshComponentCreator
from dataclasses import dataclass, field
from Intersection import *
from enum import Enum



def raster_extent(transform, width, height, col_offset=0, row_offset=0):
    
    corners = [
        # corners in pixel space, shifted by offsets
        (col_offset,         row_offset),
        (col_offset + width, row_offset),
        (col_offset,         row_offset + height),
        (col_offset + width, row_offset + height),
    ]

    xs = []
    ys = []

    for col, row in corners:
        x, y = transform * (col, row)
        xs.append(x)
        ys.append(y)

    left   = min(xs)
    right  = max(xs)
    bottom = min(ys)
    top    = max(ys)

    return {
        "left": left,
        "right": right,
        "bottom": bottom,
        "top": top
    }


class MeshCut:

    DEFAULT_ROAD_WIDTH = 3 # meters

    MAX_INTERSECT_TRI_EDGE = 6

    #
    # we should be using modified transform here so 0,0 corresponds
    # to the inset view
    #
    def __init__(self, buffer, index, tf, bnds):

        self.buffer_array = buffer
        self.index_array = index # triangle index for raster
        self.transform = tf
        self.bounds = bnds   

        #
        #
        # fields for our triangle edge index / intersection index
        self.num_lat = None
        self.num_lon = None
        self.num_diag = None
        self.num_tri_edges = None
        self.tri_edge_index = None
        self.tri_edge_intersect_index = None

        self.num_cols = self.buffer_array.shape[0]
        self.num_rows = self.buffer_array.shape[1]

        self.num_triangles = (self.num_cols - 1)*(self.num_rows - 1)*2       

        # diagonal identifiers!
        self.min_diag = -1*self.num_rows;
        self.max_diag = self.num_cols;
    
        self.reverse_tri_index = []
        self.create_reverse_index()

        self.create_edge_index()

        self.all_mesh_edges = None

        self.all_mesh_verts = None

        self.all_clip_edges = None

        self.vert_clip_index = None

    def cross_2D(self, a, b ):
        return (a.x * b.y) - (a.y * b.x)

    def convert_to_point_3D(self, _pix_x, _pix_y):

        point_y = _pix_y*self.transform.e + self.transform.f ;
        point_x = _pix_x*self.transform.a + self.transform.c ;
        z = self.buffer_array[_pix_x,_pix_y,2].item()

        return point_x, point_y, z

    def convert_to_point(self, _pix_x, _pix_y):

        point_y = _pix_y*self.transform.e + self.transform.f ;
        point_x = _pix_x*self.transform.a + self.transform.c ;
        
        return ( point_x, point_y )
    
    def populate_clip_indices(self, tri_edges_dict):


        for tri_index, inner_dict in tri_edges_dict.items():
            # 
            # need to make sure not repeating edges
            pass

    def lat_edge_index(self, x, y):
        # (x, y) -> (x + 1, y)
        return int( y * (self.num_cols - 1) + x )


    def lon_edge_index(self, x, y):
        # (x, y) -> (x, y + 1)
        return int(self.num_lat + y * self.num_cols + x)


    def diag_edge_index(self, x, y):
        # (x, y) -> (x + 1, y + 1)

        return int(
            self.num_lat
            + self.num_lon
            + y * (self.num_cols - 1)
            + x
        )


    def triangle_edge_indices(self, tri_index):
        num_cell_cols = self.num_cols - 1

        cell_index = tri_index // 2
        x = cell_index % num_cell_cols
        y = cell_index // num_cell_cols

        top = self.lat_edge_index(x, y)
        bottom = self.lat_edge_index(x, y + 1)

        left = self.lon_edge_index(x, y)
        right = self.lon_edge_index(x + 1, y)

        diagonal = self.diag_edge_index(x, y)

        if tri_index % 2 == 0:
            # upper-right
            return top, right, diagonal
        else:
            # lower-left
            return diagonal, bottom, left 
        

    def create_edge_index( self ):

        self.num_lat = self.num_rows * (self.num_cols - 1) # horizontal
        self.num_lon = (self.num_rows - 1) * self.num_cols # vertical
        self.num_diag = (self.num_rows - 1) * (self.num_cols - 1) #diagonal

        self.num_tri_edges = self.num_lat + self.num_lon + self.num_diag

        print(f"creating total of {self.num_tri_edges} edges!")

        self.tri_edge_index = []

        # Latitude Edges (x, y) , (x + 1, y)
        for y in range(self.num_rows):
            for x in range(self.num_cols - 1):
                v0 = y * self.num_cols + x
                v1 = y * self.num_cols + x + 1
                self.tri_edge_index.append((v0, v1))

        print(f"post lat edge index size: {len(self.tri_edge_index)}")

        # Longitude Edges (x, y) , (x, y + 1)
        for y in range(self.num_rows - 1):
            for x in range(self.num_cols):
                v0 = y * self.num_cols + x
                v1 = (y+1) * self.num_cols + x
                self.tri_edge_index.append((v0, v1))

        print(f"post lon edge index size: {len(self.tri_edge_index)}")

        # Diagonal edges (x, y) , (x + 1, y + 1)
        for y in range(self.num_rows - 1):
            for x in range(self.num_cols - 1):
                v0 = y * self.num_cols + x
                v1 = (y+1) * self.num_cols + (x + 1)
                self.tri_edge_index.append((v0, v1))

        print(f"post diag edge index size: {len(self.tri_edge_index)}")

        assert len(self.tri_edge_index) == self.num_tri_edges

        #self.tri_edge_intersect_index = np.full((self.num_tri_edges, self.MAX_INTERSECT_TRI_EDGE),-1, dtype=np.int32)


    def create_reverse_index(self):

        # index_l = []
        # for y in range(1,self.num_rows):
        #     for x in range(1,self.num_cols):

        #         currIndx = x + y * self.num_cols
        #         prevRowSameColIndx = x + (y - 1) * self.num_cols
        #         prevRowBackColIndx = x - 1 + (y - 1) * self.num_cols
        #         prevIndx = currIndx - 1

        #         index_l.append((prevRowBackColIndx, currIndx, prevRowSameColIndx))
        #         index_l.append((prevRowBackColIndx, prevIndx, currIndx))
        
        self.reverse_tri_index = []
        #  4--5--x 
        #  |\4|\ |  
        #  |3\|5\|
        #  3--v--0
        #  |\2|\0|
        #  | \|1\|
        #  x--2--1
        for y in range(self.num_rows):
            for x in range(self.num_cols):
                currReverse = [-1,-1,-1,-1,-1,-1]
                if (x < self.num_cols - 1 and y < self.num_rows - 1):
                    cell = x + y  * (self.num_cols - 1);
                    currReverse[0] = 2 * cell;
                    currReverse[1] = 2 * cell + 1;


                if (x > 0 and y < self.num_rows - 1):
                    cell = (x-1) + y * (self.num_cols - 1);
                    currReverse[2] = 2 * cell;
 

                if (x > 0 and y > 0 ):
                    cell = (x - 1) + (y - 1) * (self.num_cols - 1);
                    currReverse[3] = 2 * cell + 1 ;
                    currReverse[4] = 2 * cell ;

                if (x < self.num_cols - 1 and y > 0):
                    cell = x  + (y - 1) * (self.num_cols - 1) ;
                    currReverse[5] = 2 * cell + 1 ;

                self.reverse_tri_index.append(tuple(currReverse))

        # VOID
        #return index_l, reverse_l


    def get_type_edge(self, edge_index):

        if edge_index < self.num_lat:
            return "lat"
        elif edge_index < self.num_lat + self.num_lon:
            return "lon"
        else: 
            return "diag"

    #
    # for a given unique edge get the 2 corresponding triangles
    #
    def get_tris_edge_index(self, edge_index):

        v0,v1 = self.tri_edge_index[edge_index]
        edge_typ = self.get_type_edge( edge_index )

        pixel_x , pixel_y = v0 % self.num_cols, v0 // self.num_cols

        tri_indices = self.reverse_tri_index[v0]
        #  4--5--x 
        #  |\4|\ |  
        #  |3\|5\|
        #  3--v--0
        #  |\2|\0|
        #  | \|1\|
        #  x--2--1
        match edge_typ:
            case "lon":
                if (pixel_x == self.num_cols - 1):
                    tris = [2]
                elif pixel_x == 0:
                    tris = [1]
                else:
                    tris = [2,1]

            case "lat":
                if (pixel_y == self.num_rows - 1):
                    tris = [5]
                elif (pixel_y == 0):
                    tris =[0]
                else:
                    tris = [5,0]
            case "diag":
                tris = [0,1]

        candidate_tris = [ tri_indices[_indx] for _indx in tris ]
        return candidate_tris


    #
    # takes the place of explicit triangle index for our purposes
    # in python
    #
    def get_triangle_vertex_indices(self, tri_idx):
        num_cells_x = self.num_cols - 1
        num_cells_y = self.num_rows - 1
        num_triangles = 2 * num_cells_x * num_cells_y

        if not 0 <= tri_idx < num_triangles:
            raise IndexError(f"Triangle index {tri_idx} outside [0, {num_triangles - 1}]")

        cell_idx = tri_idx // 2
        tri_in_cell = tri_idx % 2

        cell_x = cell_idx % num_cells_x
        cell_y = cell_idx // num_cells_x

        # v00 = cell_y * self.num_cols + cell_x
        # v10 = v00 + 1
        # v01 = v00 + self.num_cols
        # v11 = v01 + 1

        v00 = (cell_x, cell_y)
        v10 = (cell_x+1, cell_y)
        v01 = (cell_x, cell_y+1) 
        v11 = (cell_x+1, cell_y+1) 
        if tri_in_cell == 0:
            return v00, v10, v11
        else:
            return v00, v11, v01


    #
    # 
    #
    def order_cw(self, polygon_verts):
        center_x = sum(x for x, _, _ in polygon_verts) / len(polygon_verts)
        center_y = sum(y for _, y, _ in polygon_verts) / len(polygon_verts)

        return sorted( polygon_verts,
            key=lambda p: math.atan2(p[1] - center_y, p[0] - center_x) ,
            reverse=True
        )


    #
    # returns returns [ Graticule(), ... ]
    # of full latitude, longitude, diagonals
    def get_inner_edges(self, polygon_verts):

        x,y = zip(*polygon_verts)
        min_x = np.min(x).item()
        max_x = np.max(x).item()
        min_y = np.min(y).item()
        max_y = np.max(y).item()

        #
        # convert to pixels! notice we're using the modified transform
        #
        min_y_pix = (min_y - self.transform.f ) / self.transform.e;
        max_y_pix = (max_y - self.transform.f ) / self.transform.e;

        largest_y_pix = max(max_y_pix, min_y_pix)
        smallest_y_pix = min(min_y_pix, max_y_pix)

        min_inner_x = math.floor( (min_x - ( self.transform.c ) ) / self.transform.a);
        min_inner_y = math.floor( smallest_y_pix );
        max_inner_x = math.ceil( (max_x - ( self.transform.c ) ) / self.transform.a);
        max_inner_y = math.ceil( largest_y_pix );

        print(f"min_y_pix: {min_y_pix}, max_y_pix: {max_y_pix}")
        print(f"min_inner_y: {min_inner_y},  max_inner_y: { max_inner_y} ")

        x = np.arange(min_inner_x, max_inner_x+1)
        y = np.arange(min_inner_y, max_inner_y+1)

        print(f"let's get all longitude lines: {x}")
        print(f"let's get all latitude lines: {y}")

        mesh_edges = []
        
        for x_indx in x:

            pix1 = glm.vec2( x_indx ,   0 );
            pix2 = glm.vec2( x_indx ,  self.num_rows-1 );
        
            e1 = self.convert_to_point(*pix1);
            e2 = self.convert_to_point(*pix2);

            gr = Lon(e1, e2, x_indx)

            print(f"pixels! for vertical edges: {pix1}, {pix2}")
        
            mesh_edges.append(gr)
        
        for y_indx in y:

            pix1 = glm.vec2( 0 ,   y_indx );
            pix2 = glm.vec2( self.num_cols-1 ,  y_indx );
        
            e1 = self.convert_to_point(*pix1);
            e2 = self.convert_to_point(*pix2);
        
            gr = Lat(e1, e2, y_indx)

            mesh_edges.append( gr )
        
            ###p2 = glm.vec2( 0 ,   self.num_cols )

        #for y in range(max(0, min_inner_y - 1), max_inner_y + 1):

        min_d = min_inner_x - max_inner_y
        max_d = max_inner_x - min_inner_y

        max_x = self.num_rows-1
        max_y = self.num_cols - 1

        for d in range(max(min_d, self.min_diag), min(max_d + 1, self.max_diag)):
    
            x1 = max(0, d)
            x2 = min(max_x, d + max_y)

            #if x1 > x2:
            #    return None  

            y1 = x1 - d
            y2 = x2 - d

            print(f"d: {d}; x1, y1: {(x1,y1)}; x2, y2 {(x2, y2)}")

            e1 = self.convert_to_point(x1, y1);
            e2 = self.convert_to_point(x2, y2);

            gr = Diag(e1, e2, d)
        
            mesh_edges.append(gr)

        return mesh_edges

    #
    #
    # returns [ Graticule(), ... ]
    # of latitude, longitude, diagonals sliced within polygon
    # polygon should be CCW list of 2D tuples
    #
    def slice_edges( self, inner_edges, polygon_verts, poly_id ):

        eps= 1e-7

        n = len(polygon_verts)
        # 
        # 
        new_edges = []

        #
        # for each lat, lon, diag edge
        #
        for mesh_edge in inner_edges: 

            print(mesh_edge)

            # these will be modified as 
            # we clip the mesh lines
            p0 = mesh_edge.p0 #glm.vec2(mesh_edge[0])
            p1 = mesh_edge.p1

            #
            # Lat, Lon, Diag
            typ = type(mesh_edge)

            print(f"MESH EDGE [{mesh_edge.p0}, {mesh_edge.p1}] of TYPE: {typ}")

            #
            # represents mesh edge
            #
            _dir = p1 - p0;

            keep_mesh_line = True

            for edge_i in range(n):
                e_1 = glm.vec2( polygon_verts[edge_i] )
                e_2 = glm.vec2( polygon_verts[(edge_i+1) % n] )

                edge = e_2 - e_1

                start_side = self.cross_2D(edge, p0 - e_1)
                dir_side = self.cross_2D(edge, _dir)

                # they're paralell!!
                if abs(dir_side) <= eps:
                    # need to determine if paralell 
                    # is to right or left
                    # need a point on 
                    sidde = self.cross_2D(edge, p0 - e_1)
                    
                    if sidde < -eps:
                        # eXclude entirely
                        keep_mesh_line = False
                        break;
                    
                    elif sidde > eps:
                        continue;


                t_hit = -start_side / dir_side
                
                t_min = 0.0
                t_max = 1.0
                if dir_side > 0:

                    if t_hit > t_min:
                        t_min = t_hit
                        min_clip_i = edge_i
                    #t_min = max(t_min, t_hit)
                else:
                    if t_hit < t_max:
                        t_max = t_hit
                        max_clip_i = edge_i              
                    
                    t_max = min(t_max, t_hit)

                if t_min > t_max + eps:
                    pass

                if t_min > t_max + eps:
                    keep_mesh_line = False
                    break;

                clipped_start = p0 + t_min * _dir ;
                clipped_end = p0 + t_max * _dir ;

                p0 = clipped_start
                p0_clip_index = 0

                p1 = clipped_end
                print(f"mesh lon/lat/diag originally {p0} {p1}, \n\t cut to {clipped_start}, {clipped_end}")

                _dir = p1 - p0;

            if keep_mesh_line:

                mesh_edge.p0 = p0
                mesh_edge.p1 = p1

                mesh_edge.p0_intersect_id = min_clip_i;
                mesh_edge.p1_intersect_id = max_clip_i;
                mesh_edge.poly_id = poly_id # which edge^ of which polygon?

                new_edges.append( mesh_edge )
                #new_edges.append( [tuple(clipped_start), tuple(clipped_end)] )

        return new_edges

    #
    # for 2 tuples are their x, y the same?
    #
    #
    def same_point(self, p0, p1):

       return math.dist(p0[:2] , p1[:2]) < 10**-6


    #
    # 
    #
    def get_edge_indices(self, tri_index):

        vert_indices = self.index_array[tri_index]

        if tri_index % 2 == 0:
            # *---*
            #  \  |
            #   \ | 
            #    \|  
            #     *
            diag_index = self.diag_edge_index(vert_indices[0] % self.num_cols, vert_indices[0] // self.num_cols)
            lat_index = self.lat_edge_index(vert_indices[0] % self.num_cols, vert_indices[0] // self.num_cols)
            lon_index = self.lon_edge_index(vert_indices[2] % self.num_cols, vert_indices[2] // self.num_cols)

            diag_ccw_order, lat_ccw_order, lon_ccw_order = True, False, False 
        
        else:
            #
            # *\
            # | \
            # |  \
            # *---*
            # e0_index = vert_indices[0]
            # e1_index = vert_indices[2]
            # e2_index = vert_indices[1]

            #
            # return {"edge_id": <int>, "order_ccw": <Bool>}
            #

            diag_index = self.diag_edge_index(vert_indices[0] % self.num_cols, vert_indices[0] // self.num_cols)
            lat_index = self.lat_edge_index(vert_indices[1] % self.num_cols, vert_indices[1] // self.num_cols)
            lon_index = self.lon_edge_index(vert_indices[0] % self.num_cols, vert_indices[0] // self.num_cols)

            diag_ccw_order, lat_ccw_order, lon_ccw_order = False, True, True

        return [{"index":diag_index, "ccw_order": diag_ccw_order }, \
                {"index":lon_index, "ccw_order": lon_ccw_order },\
                {"index":lat_index, "ccw_order": lat_ccw_order } \
                ]

    #
    # returns 2D point!
    # conversion for the original index to new edge way
    #
    def get_edges(self, tri_index):

        print(f"tri_index: {tri_index}")
        vert_indices = self.index_array[tri_index]

        # currIndx = c + r * dem_patch_num_verts;
        # prevRowSameColIndx = c + (r - 1) * dem_patch_num_verts;
        # prevRowBackColIndx = c - 1 + (r - 1) * dem_patch_num_verts;
        # prevIndx = currIndx - 1;
        # index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])
        # index_array.append([prevRowBackColIndx, prevIndx, currIndx])

        if tri_index % 2 == 0:
            # *---*
            #  \  |
            #   \ | 
            #    \|  
            #     *

            e0_index = vert_indices[0]
            e1_index = vert_indices[1]
            e2_index = vert_indices[2]
        else:
            #
            # *\
            # | \
            # |  \
            # *---*
            e0_index = vert_indices[0]
            e1_index = vert_indices[2]
            e2_index = vert_indices[1]
        
        e0_cr = e0_index % self.num_cols, e0_index // self.num_cols
        e1_cr = e1_index % self.num_cols, e1_index // self.num_cols
        e2_cr = e2_index % self.num_cols, e2_index // self.num_cols

        e0 = self.get_point(*e0_cr)[:2]
        e1 = self.get_point(*e1_cr)[:2]
        e2 = self.get_point(*e2_cr)[:2]

        if tri_index % 2 == 0:
            return [(e0, e1), (e2, e1), (e0, e2) ]
        else:
            return [(e0, e2), (e0, e1), (e2, e1) ]

    def bary_x_y(self, _x, _y, 
                v1,v2,v3 # all tuples!
                ):
        vec1 = glm.vec2( v2[0] - v1[0], v2[1] - v1[1]  )
        vec2 = glm.vec2( v3[0] - v1[0], v3[1] - v1[1] )
        vec3 = glm.vec2( _x - v1[0], _y - v1[1])

        dot00 = glm.dot(vec1, vec1);
        dot01 = glm.dot(vec1, vec2);
        dot02 = glm.dot(vec1, vec3);
        dot11 = glm.dot(vec2, vec2);
        dot12 = glm.dot(vec2, vec3);

        invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        _t = (dot11 * dot02 - dot01 * dot12) * invDenom;
        _u = (dot00 * dot12 - dot01 * dot02) * invDenom;
        _s = 1 - _t - _u;

        return _s, _t, _u #_t >= 0 and _u >= 0 and _s >= 0;

    #
    # for list of list of 2 tuples order 
    # assuming no overlapping !!
    #
    def order_clips(self, *point_pairs):
        def pixel_key(point_pair):
            x_m, y_m = point_pair[0][:2]
            pixel_col, pixel_row = (~self.transform) * ( x_m, y_m )

        return sorted(point_pairs, key=pixel_key)

    def get_2D_verts(self, tri_index):
        verts = self.index_array[tri_index]
        verts_c_r = [(i % self.num_cols, i // self.num_cols) for i in verts]

        return [ (self.buffer_array[c,r,0],self.buffer_array[c,r,1]) for c,r in verts_c_r ]

    #
    # return Triangle Index of point
    #
    def get_tri_index_point(self, _x_m, _y_m ):

        col_f, row_f = (~self.transform) * ( _x_m, _y_m )

        pixel_x, pixel_y = math.floor(col_f), math.floor(row_f)
        
        vert_indx = pixel_y*self.num_cols + pixel_x

        tri_indices = self.reverse_tri_index[vert_indx]

        #  4--5--x 
        #  |\4|\ |  
        #  |3\|5\|
        #  3--v--0
        #  |\2|\0|
        #  | \|1\|
        #  x--2--1

        # bottom right
        # *--*
        # |\ |
        # | \|
        # *--*
        if (pixel_x == self.num_cols - 1 and pixel_y == self.num_rows - 1):

            tri_index_0 = 4
            tri_index_1 = None

        # at right border but not bottom
        elif pixel_x == self.num_cols - 1:
            tri_index_0 = 2
            tri_index_1 = None
            #tri_index = self.reverse_tri_index[vert_indx]

        # at bottom border but not right
        elif pixel_y == self.num_rows - 1:
            tri_index_0 = 5
            tri_index_1 = None  

        else:
            tri_index_0 = 0
            tri_index_1 = 1


        if tri_index_1 == None:
            return tri_indices[tri_index_0]

        else: 
            index_tri_0 = tri_indices[tri_index_0]
            index_tri_1 = tri_indices[tri_index_1]
        
        tri_0_verts = self.get_2D_verts(index_tri_0)
        tri_1_verts = self.get_2D_verts(index_tri_1)

        stu_0 = self.bary_x_y(_x_m, _y_m, *tri_0_verts)
        stu_1 = self.bary_x_y(_x_m, _y_m, *tri_1_verts)

        if all(round(e,10) >= 0 for e in stu_0):
            return index_tri_0
        elif all(round(e,10) >= 0 for e in stu_1):
            return index_tri_1
        else:
            raise ValueError("neither triangle contained point!")

    #
    # returns bool for whether point is on segment
    # and t for how far along
    def point_on_segment(self, p: glm.vec2, p0: glm.vec2, p1: glm.vec2, eps=1e-9):

        p = glm.vec2(p)
        p0 = glm.vec2(p0)
        p1 = glm.vec2(p1)

        seg = p1 - p0
        rel = p - p0

        # Must lie on the infinite line
        if abs(self.cross_2D(seg, rel)) > eps:
            return False, 0.0

        # Parameter along p0 -> p1
        seg_len_sq = glm.dot(seg, seg)

        assert(seg_len_sq > 0)

        u = glm.dot(rel, seg) / seg_len_sq

        # Must also lie within the segment
        return -eps <= u <= 1.0 + eps , u


    #
    # for an intersection of a Graticule 
    #
    def get_mesh_edge_inter(self,
                                graticule_type: Graticule, 
                                inters : tuple[float, float] ):

        inter_pix_x = math.floor( (inters[0] - self.transform.c ) / self.transform.a );
        inter_pix_y = math.floor( (inters[1] - self.transform.f ) / self.transform.e );

        print(f"floored pixel: {(inter_pix_x, inter_pix_y)}")

        print(f"post border containment: {(inter_pix_x, inter_pix_y)}")

        match graticule_type:
            case Graticule.LON:
                
                if (inter_pix_y >= self.num_rows):
                    raise ValueError(f"off grid pixel: {(inter_pix_x, inter_pix_y)}")

                # this will happen for border intersection!
                elif (inter_pix_y == self.num_rows-1):
                    inter_pix_y = inter_pix_y-1

                edge_index = self.lon_edge_index(inter_pix_x, inter_pix_y)
            case Graticule.LAT:

                if (inter_pix_x >= self.num_cols):
                    raise ValueError(f"off grid pixel: {(inter_pix_x, inter_pix_y)}")

                # this will happen for border intersection!
                elif (inter_pix_x == self.num_cols-1):
                    inter_pix_x = inter_pix_x-1

                edge_index = self.lat_edge_index(inter_pix_x, inter_pix_y)

            case Graticule.DIAG:
                
                if (inter_pix_x == self.num_cols-1) or (inter_pix_y == self.num_rows-1):
                    inter_pix_x = inter_pix_x-1
                    inter_pix_y= inter_pix_y-1

                edge_index = self.diag_edge_index(inter_pix_x, inter_pix_y)

        #
        # check for degenerate case!
        #
        edge = self.get_edge_points(edge_index)
        on_edge, _ = self.point_on_segment(inters, edge[0], edge[1])

        if not on_edge:
            print(f"Degenerate edge intersection !!")
            vert_index = inter_pix_x + inter_pix_y*self.num_cols
            othr_edges = self.get_edges_vert( vert_index )

            for o_e in othr_edges: # iterating over tuples

                if o_e is not None:
                    o_edge_indx = o_e[0]
                    o_edge = self.get_edge_points(o_edge_indx)
                    on_edge, _ = self.point_on_segment(inters, o_edge[0], o_edge[1])
                    if on_edge:
                        return o_edge_indx
            else:
                raise ValueError("Can't find right edge!")

        else:
            return edge_index

    #
    # returns [] for no intersection and 
    # [ Intersection:
    #       mesh_feature: Edge or Vertex
    #       clip_feature: Edge or Vertex
    #       point       ]
    # p0, p1 represent clip points (CCW! 3D!)
    # e0, e0 represent mesh points
    #
    # this all assumes that e0 != e1 
    def get_intersection(self, 
                            graticule_type: Graticule
                            , p0 : tuple[float, float, float] , p1 :tuple[float, float, float] 
                             , e0 : tuple[float, float, float] , e1 : tuple[float, float, float]
                             , poly_id
                             , poly_num_points):
        eps = 1e-9

        # for now
        #e0_z = e0[2]
        #e1_z = e1[2]

        p0_glm = glm.vec2(p0)
        p1_glm = glm.vec2(p1)
        e0_glm = glm.vec2(e0)
        e1_glm = glm.vec2(e1)

        mesh_dir = e1_glm - e0_glm;
        clip_dir = p1_glm - p0_glm;

        denominator = self.cross_2D(clip_dir, mesh_dir)
        print(f"clip graticule determinant: {denominator}")

        rel_dir = e0_glm - p0_glm

        #
        # 
        # u is parameter for how along mesh edge we are but notice
        # we're slicing entire GRATICULES with this method so we 
        # need to find intersection first

        #
        # paralell!
        if abs(denominator) <= eps:

            if abs(self.cross_2D(mesh_dir, rel_dir)) > eps: # never touch paralell
                return None
            else:
                #
                # on same line (inf intersections)
                # but we only consider the clip point p0!
                # notice due to 32 bit glm vec we need to change
                # this back to a tuple instead of using original p0
                intersection = tuple( p0_glm )   

                edge_index = self.get_mesh_edge_inter(
                                graticule_type, 
                                intersection )

                print(f"intersection: {tuple(intersection)}, grat: {graticule_type} results in edge: {edge_index}")
                print(f"meshy edge: {e0} , {e1}")
                edge_tup = self.tri_edge_index[edge_index]

                edge0_index, edge1_index = edge_tup

                edge0_vert = self.get_point( edge0_index )
                edge1_vert = self.get_point( edge1_index )

                edge0_glm = glm.vec2(edge0_vert)
                edge1_glm = glm.vec2(edge1_vert)

                sub_dir = edge1_glm - edge0_glm
                denom = glm.dot(sub_dir, sub_dir)

                edge_t = glm.dot(p0_glm - edge0_glm, sub_dir) / denom
                clip_t = 0.0
        else:

            clip_t = self.cross_2D(rel_dir, mesh_dir) / denominator
            mesh_t = self.cross_2D(rel_dir, clip_dir) / denominator

            #
            # notice no epsilon for upper bound!
            clip_edge_hit = 0 <= clip_t < 1.0

            graticule_hit = -eps <= mesh_t <= 1.0 + eps

            print(f"intersection clip hit: {clip_t} \n graticule hit: { mesh_t }")


            if clip_edge_hit and graticule_hit:
                intersection = tuple(p0_glm + clip_t * clip_dir)

                edge_index = self.get_mesh_edge_inter(
                                    graticule_type, 
                                    intersection )

                print(f"intersection: {tuple(intersection)}, grat: {graticule_type} results in edge: {edge_index}")
                print(f"meshy edge: {e0} , {e1}")
                edge_tup = self.tri_edge_index[edge_index]

                edge0_index, edge1_index = edge_tup

                edge0_vert = self.get_point( edge0_index )
                edge1_vert = self.get_point( edge1_index )

                edge0_glm = glm.vec2(edge0_vert)
                edge1_glm = glm.vec2(edge1_vert)

                sub_dir = edge1_glm - edge0_glm
                denom = glm.dot(sub_dir, sub_dir)

                edge_t = glm.dot(intersection - edge0_glm, sub_dir) / denom

            else:
                return None

        new_eps = np.spacing( np.float32(intersection[0]) )

        if edge_t <= new_eps:
            mesh_vertex_index = edge0_index
            mesh_feature = FeatureType.VERTEX
        elif edge_t >= 1.0 - new_eps:
            mesh_vertex_index = edge1_index
            mesh_feature = FeatureType.VERTEX
        else:
            mesh_vertex_index = None
            mesh_feature = FeatureType.EDGE

        if clip_t <= new_eps:
            clip_feature = FeatureType.VERTEX
        elif clip_t >= 1.0 - new_eps:
            clip_feature = FeatureType.VERTEX
            poly_id = (poly_id[0], (poly_id[1]+1)%poly_num_points)
        else:
            clip_feature = FeatureType.EDGE            

        return Intersection( point=intersection,
                                mesh_feature=mesh_feature,
                                clip_feature=clip_feature,
                                mesh_edge_index=edge_index,
                                mesh_vertex_index=mesh_vertex_index,
                                graticule_type=graticule_type,
                                mesh_t=edge_t,
                                clip_t=clip_t, 
                                clip_poly_id=poly_id )


    #
    # 2D interpolated coords for (1.5, 1.0) pixel por ejemplo
    def interpolate_xy_coords(self, pix : Tuple[float, float]):

        x,y = (pix[0] * self.transform.a + self.transform.c), (pix[1] * self.transform.e) + self.transform.f
        return x,y

    #
    # one edge index might be wrong lol due to degeneracy lol
    #
    def get_vert_inters_edges(self, point, edge_i, edge_j):

        v0, v1 = self.tri_edge_index[edge_i]
        v2, v3 = self.tri_edge_index[edge_j]

        verts = [v0, v1, v2, v3]

        point_glm = glm.vec2( point )
        min_dist = math.inf
         
        for vi in verts:
            curr_vert = self.get_point(vi)[:2]
            if min_dist > glm.distance( glm.vec2(curr_vert), point_glm ):
                min_dist = glm.distance( glm.vec2(curr_vert), point_glm )
                min_vert = vi

        print(f"resturning vert of index {min_dist}")
        return min_vert 
    
    def get_mesh_edges_vertex(self, edge_i, edge_j):
    
        v0, v1 = self.tri_edge_index[edge_i]
        v2, v3 = self.tri_edge_index[edge_j]

        shared = set((v0, v1)) & set((v2, v3))
        if len(shared) != 1:
            raise ValueError("Can't get common vertex two edges")

        return next(iter(shared))


    # method that returns a list of intersections defined as
    #  >> see Intersection class <<
    # does degenerecy de-duplication!
    #
    #
    def clip_mesh_edges(self, clip_e0, clip_e1, polygon_id, num_clip_points):

        inters_list = []

        _p0 = clip_e0
        _p1 = clip_e1

        #  num_clip_points = len( clip_poly_point[polygon_id[0]] )

        #min_x = np.min(clip_e0[0], clip_e1[0]).item()
        #max_x = np.max(clip_e0[0], clip_e1[0]).item()
        #min_y = np.min(clip_e0[1], clip_e1[1]).item()
        #max_y = np.max(clip_e0[1], clip_e1[1]).item()     

        e0_x_pix = (clip_e0[0] - self.transform.c ) / self.transform.a;
        e1_x_pix = (clip_e1[0] - self.transform.c ) / self.transform.a;

        e0_y_pix = (clip_e0[1] - self.transform.f ) / self.transform.e;
        e1_y_pix = (clip_e1[1] - self.transform.f ) / self.transform.e;

        max_x_pix = max( e0_x_pix, e1_x_pix )
        min_x_pix = min( e0_x_pix, e1_x_pix )

        max_y_pix = max( e0_y_pix, e1_y_pix )
        min_y_pix = min( e0_y_pix, e1_y_pix )

        min_inner_x = math.floor( min_x_pix );
        min_inner_y = math.floor( min_y_pix );
        
        max_inner_x = math.ceil( max_x_pix );
        max_inner_y = math.ceil( max_y_pix );

    
        x = np.arange(min_inner_x, max_inner_x + 1)
        y = np.arange(min_inner_y, max_inner_y + 1)

        # 
        # vertical (lon) lines
        for x_indx in x:

            pix1 = ( int(x_indx) ,  0 );
            pix2 = ( int(x_indx) , self.num_rows-1 );
        
            _e0 = self.convert_to_point(*pix1);
            _e1 = self.convert_to_point(*pix2);

            #
            # here we use _e0, _e1
            _inter = self.get_intersection(Graticule.LON, _p0, _p1\
                                           , _e0, _e1
                                           , polygon_id, num_clip_points)

            #_inter.clip_poly_id = polygon_id
            if _inter is not None:
                _inter.grat_pix0 = pix1
                _inter.grat_pix1 = pix2
                inters_list.append( _inter )
            
        # 
        # horizontal (lat) lines
        for y_indx in y:

            pix1 = ( 0 ,  int( y_indx ) );
            pix2 = ( self.num_cols-1 , int(y_indx) );
        
            _e0 = self.convert_to_point(*pix1);
            _e1 = self.convert_to_point(*pix2);
        
            _inter = self.get_intersection( Graticule.LAT, _p0, _p1\
                                            , _e0, _e1
                                            , polygon_id, num_clip_points )

            #_inter.clip_poly_id = polygon_id

            if _inter is not None:
                _inter.grat_pix0 = pix1
                _inter.grat_pix1 = pix2
                inters_list.append( _inter )

        min_d = min_inner_x - max_inner_y
        max_d = max_inner_x - min_inner_y

        #max_x = self.num_rows - 1
        #max_y = self.num_cols - 1
        
        # SWiTCH

        max_x = self.num_cols - 1
        max_y = self.num_rows - 1

        d_start = max(min_d, self.min_diag)
        d_end = min(max_d, self.max_diag)

        for d in range(d_start, d_end+1):
    
            x1 = max(0, d)
            x2 = min(max_x, d + max_y)

            if x1 == x2:
                continue

            print(f"x1: {x1}, x2: {x2}")

            y1 = x1 - d
            y2 = x2 - d

            print(f"diag pixels: {(x1,y1)}, {(x2,y2)}")

            _e0 = self.convert_to_point(x1, y1);
            _e1 = self.convert_to_point(x2, y2);

            _inter = self.get_intersection(Graticule.DIAG, _p0, _p1\
                                           , _e0, _e1
                                           , polygon_id, num_clip_points )

            if _inter is not None:
                _inter.grat_pix0 = (int(x1),int(y1))
                _inter.grat_pix1 = (int(x2),int(y2))
                inters_list.append( _inter )

        #
        # de-generacy check
        #
        # I DON'T THINK WE NEED ANYMORE
        # groups = defaultdict(list)
        # for _inter in inters_list:
        #     if _inter.mesh_feature == FeatureType.EDGE:
        #         groups[_inter.point].append(_inter)

        # non_distinct_point = [
        #     group for group in groups.values()
        #     if len(group) > 1
        # ]

        # if len(non_distinct_point) > 0:
        #     print(f"EDGES WITH THE SAME INTER POINT")
        #     print(non_distinct_point)
        #     non_distincts = [nd[0].point for nd in non_distinct_point]
        #     print(f"non-distincts: {non_distincts}")
        #     inters_list = [_inter for _inter in inters_list if _inter.point not in non_distincts]

        # for inter_group in non_distinct_point:

        #     changeme = inter_group[0]
        #     # 
        #     # change this to a vertex of intersection of two edges
        #     #
        #     edgeA = inter_group[0].mesh_edge_index
        #     edgeB = inter_group[1].mesh_edge_index

        #     inters = inter_group[0].point
        #     print(f"find edge intersection of {edgeA}, {edgeB}")

        #     vert_indx = self.get_vert_inters_edges(inters, edgeA, edgeB)
        #     false_vert_indx = self.get_mesh_edges_vertex(edgeA, edgeB)

        #     print(f"behold the shared vertex: {vert_indx} vs. false one: {false_vert_indx}")

        #     vert_point = self.get_point(vert_indx)[:2]

        #     changeme.point = vert_point
        #     changeme.mesh_vertex_index = vert_indx
        #     changeme.mesh_feature = FeatureType.VERTEX
        #     inters_list.append(changeme)


        return inters_list

    #
    #
    # returns each edge and the notion of how
    # whether we're coming from first or second vertex
    def get_edges_vert(self, vert_index):

        c,r = vert_index % self.num_cols , vert_index // self.num_cols

        topRight = c == self.num_cols - 1 and r == 0
        topLeft = c == 0 and r == 0
        bottomRight = c == self.num_cols - 1 and r == self.num_rows - 1
        bottomLeft = c == 0 and r == self.num_rows - 1

        #
        # *---*---*
        # |\  |\  |
        # | 4 5 \ |
        # |  \|  \|
        # *-3-C-0-*
        # |\  |\  |
        # | \ 2 1 |
        # |  \|  \|
        # *---*---*

        if topRight:
            edges = [
                None , # 0
                None , # 1
                ( self.lon_edge_index(c, r), Graticule.LON, True) ,  # 2
                ( self.lat_edge_index(c-1, r), Graticule.LAT, False) , # 3
                None ,
                None
            ]

        elif topLeft:
            edges = [
                (self.lat_edge_index(c, r), Graticule.LAT, True) ,  # 0
                (self.diag_edge_index(c, r), Graticule.DIAG, True) , # 1
                (self.lon_edge_index(c, r), Graticule.LON, True),   # 2
                None ,
                None ,
                None          
            ]
        elif bottomRight:
            edges = [
                None ,                                                  # 0
                None ,                                                  # 1
                None ,                                                  # 2
                (self.lat_edge_index(c-1, r), Graticule.LAT, False ) ,  # 3
                (self.diag_edge_index(c-1, r-1), Graticule.DIAG, False),# 4
                (self.lon_edge_index(c, r-1), Graticule.LON, False )    #5
            ]

        elif bottomLeft:
            edges = [
                (self.lat_edge_index(c, r), Graticule.LAT, True) ,  # 0
                None ,
                None ,
                None ,
                None ,
                (self.lon_edge_index(c, r-1), Graticule.LON, False) # 5
            ]

        elif c == self.num_cols - 1: #right column   
            edges = [
                None ,
                None ,
                (self.lat_edge_index(c-1, r), Graticule.LAT, False) ,    # 2
                (self.lon_edge_index(c, r-1), Graticule.LON , False) ,   # 3
                (self.diag_edge_index(c-1, r-1), Graticule.DIAG, False), # 4
                (self.lon_edge_index(c, r), Graticule.LON, True )        # 5
            ] 

        elif (c == 0): #left column
            edges = [
                (self.lat_edge_index(c, r), Graticule.LAT, True)  ,   # 0
                (self.diag_edge_index(c, r), Graticule.DIAG, True),   # 1
                (self.lon_edge_index(c, r), Graticule.LON, True)  ,   # 2
                None ,                                                # 3
                None ,                                                # 4
                (self.lon_edge_index(c, r-1), Graticule.LON, False)   # 5  
            ]    

        elif (r == 0): #top
            edges = [
                (self.lat_edge_index(c, r), Graticule.LAT, True),     # 0
                (self.diag_edge_index(c, r), Graticule.DIAG, True),   # 1
                (self.lon_edge_index(c, r), Graticule.LON, True),     # 2          
                (self.lat_edge_index(c-1, r), Graticule.LAT, False),  # 3
                None ,
                None
            ]
        elif (r == self.num_rows - 1): # bottom

            edges = [
                (self.lat_edge_index(c, r), Graticule.LAT, True),       # 0 
                None ,
                None ,         
                (self.lat_edge_index(c-1, r), Graticule.LAT, False),    # 3
                (self.diag_edge_index(c-1, r-1), Graticule.DIAG, False),# 4
                (self.lon_edge_index(c, r-1), Graticule.LON, False)     # 5  
            ]

        else: # interior

            edges = [
                (self.lat_edge_index(c, r), Graticule.LAT , True),      # 0
                (self.diag_edge_index(c, r), Graticule.DIAG, True),     # 1  
                (self.lon_edge_index(c, r), Graticule.LON, True),       # 2       
                (self.lat_edge_index(c-1, r), Graticule.LAT, False),    # 3
                (self.diag_edge_index(c-1, r-1), Graticule.DIAG, False),# 4
                (self.lon_edge_index(c, r-1), Graticule.LON, False)     # 5 
            ]

        return edges

    #
    # void method that modifies tri_edges_dict
    # 
    def clip_tris( self, clippd_edges, polygon_verts, tri_edges_dict ):
        #
        # holding a tuple representing edge cut
        # like (0, 0.35) or like 
        #
        # But these edges are different from the triangle indexed edges!
        #      2
        #   *----*
        #   * \  |
        # 1 *  0 | 1
        #   *   \|
        #   *----*
        #     2

        for clippd_edge in clippd_edges:

            p0_pix_x = math.floor( ( clippd_edge.p0.x - ( self.transform.c ) ) / self.transform.a);
            p0_pix_y = math.floor( ( clippd_edge.p0.y - self.transform.f ) / self.transform.e );

            p1_pix_x = math.floor( ( clippd_edge.p1.x - ( self.transform.c ) ) / self.transform.a);
            p1_pix_y = math.floor( ( clippd_edge.p1.y - self.transform.f ) / self.transform.e );

            #
            # based on type of Graticule we must iterate over potential edges on given edge thing
            #
            points_list = []
            pixel_points = []

            match clippd_edge:

                case Lat():
                    assert( p0_pix_y == p1_pix_y)
                    pixel_points = [(col, p0_pix_y) for col in range(min(p0_pix_x, p1_pix_x), max(p0_pix_x, p1_pix_x)+ 1)]
                
                    edge_index = 2

                    start_edge_index = self.lat_edge_index(*pixel_points[0])
                    end_edge_index = self.lat_edge_index(*pixel_points[-1])

                case Lon():

                    assert( p0_pix_x == p1_pix_x)
                    pixel_points = [(p0_pix_x, row) for row in range(min(p0_pix_y, p1_pix_y), max(p0_pix_y, p1_pix_y)+ 1)]
                
                    edge_index = 1

                    start_edge_index = self.lon_edge_index(*pixel_points[0])
                    end_edge_index = self.lon_edge_index(*pixel_points[-1])
                case Diag():

                    d0 = p0_pix_x - p0_pix_y
                    d1 = p1_pix_x - p1_pix_y

                    edge_index = 0

                    # this is better?
                    diagonal_id = round((d0 + d1) / 2.0)
                    assert math.isclose(d0, d1, abs_tol=1e-7)

                    for _col in range(min(p0_pix_x, p1_pix_x), max(p0_pix_x, p1_pix_x)+1):
                        _row = _col - diagonal_id

                        # _row = 
                        pixel_points.append( (_col, _row) )

                    start_edge_index = self.diag_edge_index(*pixel_points[0])
                    end_edge_index = self.diag_edge_index(*pixel_points[-1])

                    print(f"For DIAG {diagonal_id} p0x: {p0_pix_x}, p1x: {p1_pix_x}, pixels: {pixel_points}")


            inter_1 = { "edge_index" : start_edge_index ,   }
            inter_2 = { "edge_index" : end_edge_index   }

            max_iter = len(pixel_points) - 1;
            for _i , pixel in enumerate( pixel_points ):

                is_first = _i == 0;
                is_last = _i == max_iter;
            
                vert_indx = pixel[1] * self.num_cols + pixel[0]

                potential_tri_indices = self.reverse_tri_index[vert_indx]

                curr_point = self.get_point(*pixel)

                match clippd_edge:
                    #  4--5--x 
                    #  |\4|\ |  
                    #  |3\|5\|
                    #  3--v--0
                    #  |\2|\0|
                    #  | \|1\|
                    #  x--2--1
                    case Lon():
                        #s = "LON"
                        #tri_index_0 = potential_tri_indices[1]
                        #tri_index_1 = potential_tri_indices[2]
                        next_pixel = ( pixel[0], pixel[1]+1 )

                    case Lat():
                        #s = "LAT"
                        #tri_index_0 = potential_tri_indices[5]
                        #tri_index_1 = potential_tri_indices[0]
                        next_pixel = ( pixel[0]+1, pixel[1] )

                    case Diag():
                        #s = "DIAG"
                        #tri_index_0 = potential_tri_indices[0]
                        #tri_index_1 = potential_tri_indices[1]
                        next_pixel = ( pixel[0]+1, pixel[1]+1 )

                next_point = self.get_point(*next_pixel)

                next_point_2D = next_point[:2]
                curr_point_2D = curr_point[:2]

                curr_vert = self.get_vert_index(*pixel);
                next_vert = self.get_vert_index(*next_pixel);

                _p0 = tuple(clippd_edge.p0)
                _p1 = tuple(clippd_edge.p1)

                #
                # clip part is middle
                # *-------xxxxx-----* 
                if is_first and is_last:

                    self.vert_clip_dict[curr_vert] = self.same_point(_p0, curr_point_2D)
                    self.vert_clip_dict[next_vert] = self.same_point( _p1, next_point_2D)

                    #clip_edge_i = [ [(clippd_edge.poly_id, clippd_edge.p0_intersect_id),(clippd_edge.poly_id, clippd_edge.p1_intersect_id)] ]
                #
                # clip part is end
                # *-------xxxxx* 
                elif is_first:

                    self.vert_clip_dict[curr_vert] = self.same_point(_p0, curr_point_2D)
                    self.vert_clip_dict[next_vert] = True ;

                    #clip_verts = [ self.same_point( _p0, curr_point_2D ), True ]
                    #clip_edge_i = [ [ (clippd_edge.poly_id, clippd_edge.p0_intersect_id), ( None, None) ]  ]

                #
                # clip part is start
                # *-------xxxxx* 
                elif is_last:

                    self.vert_clip_dict[curr_vert] = True; 
                    self.vert_clip_dict[next_vert] = self.same_point( _p1, next_point_2D ) ;

                # clip part is WHOLE THING
                # *xxxxx* 
                else: 

                    self.vert_clip_dict[curr_vert] = True; 
                    self.vert_clip_dict[next_vert] = True ;

                #tri_edges_dict[tri_index_0]["edges"][edge_index]["bool"] = clip_verts

    #
    # for test situations this is helpful to be able to run 
    # multiple times "reseting"
    #
    def perform_clipping(self, poly_clip_list):

        self.all_mesh_edges = defaultdict( list )

        self.all_mesh_verts = defaultdict( list )

        # only have vertex/edge separate indices for mesh!
        # self.all_clip_edges = defaultdict( lambda: defaultdict(self.make_both_edge_dict) )
        self.all_clip_edges = defaultdict( list )

        self.all_clip_verts = defaultdict( list )

        #
        # now a "poly" is a single clipping entity
        for poly_id, poly in enumerate( poly_clip_list ): #[ ccw_quad_points ] ):

            poly_num_edges = len(poly)

            print(f"poly: {poly}")

            for _j in range(len(poly)):

                clip_p0 = poly[_j]
                clip_p1 = poly[(_j+1) % poly_num_edges]

                poly_edge_inters = self.clip_mesh_edges( clip_p0, clip_p1, (0,_j), poly_num_edges )

                self.all_clip_edges[poly_id , _j] = poly_edge_inters

                for _inter in poly_edge_inters:

                    #
                    # MESH
                    #
                    if _inter.mesh_feature == FeatureType.VERTEX:

                        if _inter.mesh_vertex_index not in self.all_mesh_verts:
                            self.all_mesh_verts[_inter.mesh_vertex_index].append(_inter)
                    else:
                        assert _inter.mesh_feature == FeatureType.EDGE, "Unknown feature type detected"
                        self.all_mesh_edges[_inter.mesh_edge_index].append(_inter)

                    #
                    # CLIP
                    #
                    if _inter.clip_feature == FeatureType.VERTEX:
                        if _inter.clip_poly_id not in self.all_clip_verts:
                            self.all_clip_verts[_inter.clip_poly_id].append(_inter)

        #
        # with single clipping polygon shouldn't have multiple
        # intersections per mesh tri vertex
        for vert_i, val in self.all_mesh_verts.items():
            assert len(val) == 1 , f"Multiple intersections found vertex {vert_i}"

        #
        # now we modify in place the edge intersections
        for tri_edge_i, intersections in self.all_mesh_edges.items():

            self.all_mesh_edges[tri_edge_i].sort( key=lambda _inte: _inte.mesh_t )

            # Do we need this?
            for __i,inter in enumerate( self.all_mesh_edges[tri_edge_i] ):

                self.all_mesh_edges[tri_edge_i][__i].mesh_edge_order = __i;





    #
    # for a given polygon of n verts (type glm vec2) get the mesh 
    # vertices within it
    #
    def get_inner_verts(self, polygon_verts):

        x,y = zip(*polygon_verts)
        min_x = np.min(x).item()
        max_x = np.max(x).item()
        min_y = np.min(y).item()
        max_y = np.max(y).item()

        #
        # convert to pixels!
        #
        min_inner_x = math.ceil( (min_x - ( self.transform.c ) ) / self.transform.a);
        min_inner_y_cand = math.ceil( (min_y - self.transform.f ) / self.transform.e );
        max_inner_x = math.floor( (max_x - ( self.transform.c ) ) / self.transform.a);
        max_inner_y_cand = math.floor( (max_y - self.transform.f ) / self.transform.e );

        min_inner_y = min(min_inner_y_cand, max_inner_y_cand)
        max_inner_y = max(min_inner_y_cand, max_inner_y_cand)

        print(f"min_inner_x: {min_inner_x}, max_inner_x: {max_inner_x}")
        print(f"min_inner_y: {min_inner_y}, max_inner_y: {max_inner_y}")

        x = np.arange(min_inner_x, max_inner_x+1)
        y = np.arange(min_inner_y, max_inner_y+1)

        X, Y = np.meshgrid(x, y)
        pixels = np.vstack([X.ravel(), Y.ravel()]).T

        n = len(polygon_verts)

        inside_points = []

        print(f"pixels: {pixels}")
        for pix in pixels:

            # glm vec2
            point = self.convert_to_point(*pix)

            pix_in_poly = True

            print(f"PIXEL: {pix}")
            for edge_i in range(n):
                e_1 = polygon_verts[edge_i]
                e_2 = polygon_verts[(edge_i+1)%n]

                cross_prod = self.cross_2D(e_2 - e_1, point - e_1 )
                print(f"\tcross product: {cross_prod}")

                #
                # I don't know why this test is "inside point should be negative"
                #
                if cross_prod > 0: # should use epsilon?
                    pix_in_poly = False 
                    #break;

            if pix_in_poly:
                inside_points.append({ "point": point, "pixel": pix })

        return inside_points
    

    def get_vert(self, tri_index, vert_index, edge_dict):

        if tri_index % 2 == 0:
            conversion_dict = {
                0: (0, 0),
                1: (1,1),
                2: (2,1)
            }
        else:
             conversion_dict = {
                0: (0, 0),
                1: (0,1),
                2: (2,0)
            }   

        edge_index, _order = conversion_dict[vert_index] 

        return edge_dict["edges"][edge_index]["orig"][_order]



    #
    # what is the index of the edge of clip poly that intersects this
    # triangle edge
    #
    def clip_index_intercept(self, edge_index, _order, edge_dict):

        return edge_dict["edges"][edge_index]["clip_index"][_order][0]
    


    def get_edge_points_between(self, curr_clip_intersection
                                    , _edge_dict
                                    , non_clip_points):
        
        curr_edge_id = curr_clip_intersection["edge_id"]
        curr_edge_or = curr_clip_intersection["edge_order"]

        new_tri_edge = curr_clip_intersection["tri_edge"]

        inter_none = _edge_dict.get(curr_edge_id, None)

        if inter_none is None:
            raise ValueError("No intersections for edge")

        num_inters = len( inter_none )
        edge_local_order = curr_clip_intersection["edge_local_order"]
            
        new_edge_ccw = tri_edge_list[edge_local_order]["ccw_order"]
        
        if new_edge_ccw:

                if curr_clip_intersection["edge_order"] <= num_inters:
                    order_local_index = curr_clip_intersection["edge_order"] + 1 
                    edge_local_order = edge_local_order    
                    found_next_intersection = True

                else:
                    # get next edge!
                    edge_local_order = (edge_local_order + 1)%3
                    order_local_index = 0 
                    found_next_intersection = False
        
        else:
            if curr_clip_intersection["edge_order"] > 0:
                edge_local_order = edge_local_order 
                order_local_index = curr_clip_intersection["edge_order"] - 1
                   
                found_next_intersection = True

            else:
                # get next edge!
                edge_local_order = (edge_local_order + 1)%3
                order_local_index = 0 
                found_next_intersection = False

        if found_next_intersection:
            edge_id = tri_edge_list[edge_local_order]
            return all_tri_edges[edge_id][order_local_index]

        else: 

            while(not found_next_intersection):

                edge_id = tri_edge_list[edge_local_order]
                new_edge_ccw = tri_edge_list[edge_local_order]["ccw_order"]
        

                # # Add appropriate edge vertex
                vertex_rel_edge = 0 if new_edge_ccw else 1
                edge_tup = self.tri_edge_index[edge_id]
                ccw_vertex_index = edge_tup[vertex_rel_edge]

                ccw_vert_cr = ccw_vertex_index % self.num_cols, ccw_vertex_index // self.num_cols
                ccw_vert = self.get_point(*ccw_vert_cr)[:2]

                non_clip_points.append( ccw_vert )

                new_edge_inters = all_tri_edges.get(edge_id, None)

                if new_edge_inters is None:
                    edge_local_order = (edge_local_order+1)%3

                else:

                    if new_edge_ccw:
                        return new_edge_inters[0]

                    else:
                        return new_edge_inters[-1]


    #
    # technical debt builds so quickly...
    # returns new intersection with edge
    # appends to points!
    #
    def get_clip_points_between(self, clip_poly_list
                                    , curr_clip_intersection
                                    , _clip_dict
                                    , poly_points):

        curr_poly_id = curr_clip_intersection["poly_id"]
        curr_poly_order = curr_clip_intersection["poly_order"]

        next_inter_found = False

        print(f"starting at poly: {curr_poly_id}, #: {curr_poly_order}")

        #
        # n verts, n edges
        #
        num_poly_edges = len( clip_poly_list[curr_poly_id[0]] )

        while(not next_inter_found):

            num_inters = len( _clip_dict[curr_poly_id[0]][curr_poly_id[1]] )

            print(f"# intersections clip edge: {num_inters} vs. clip order: {curr_poly_order}")

            if curr_poly_order < num_inters - 1:
                curr_poly_id = curr_poly_id # keep same
                curr_poly_order+=1
                return _clip_dict[curr_poly_id[0]][curr_poly_id[1]][curr_poly_order]

            else:
                curr_poly_order = 0 # reset

                # new_point = clip_order_d[curr_poly_id][1]
                
                curr_poly_edge = (curr_poly_id[1] + 1) % num_poly_edges
                #curr_poly_id = clip_order_d[curr_poly_id][0]

                #
                # we never change the polygon rite?
                curr_poly_id = (curr_poly_id[0] , curr_poly_edge)
                
                new_point = clip_poly_list[curr_poly_id[0]][curr_poly_id[1]]

                poly_points.append( {"point": new_point, "type": "intersection"} )

                print(f"POLY ID: {curr_poly_id}")

                # does this have any intersections?
                if curr_poly_id[0] in _clip_dict and curr_poly_id[1] in _clip_dict[ curr_poly_id[0] ] \
                    and len(_clip_dict[ curr_poly_id[0] ][curr_poly_id[1]]) > 0:

                    return _clip_dict[curr_poly_id[0]][curr_poly_id[1]][curr_poly_order]
                else: 
                    continue




    #
    # for list of points we have edges we want to get intermediary points for 
    # e.g. 1 -> 2 give us point 2
    #      3 -> 0 give us point 0
    def get_points_between(self, start_edge_tup, end_edge_tup, clip_pol_d):

        poly_id, poly_edge_id = start_edge_tup
        end_poly_id, end_poly_edge_id = end_edge_tup

        if poly_id == end_poly_id:
            return []

        #print(f"poly_id: {poly_id}, poly_edge_id: {poly_edge_id}")
        #print( clip_pol_d[poly_id] )
        return [ clip_pol_d[poly_id]["points"][poly_edge_id] ]
        

        # n = len(points)

        # if n == 0:
        #     return []

        # #
        # # 
        # if start_edge_i == end_edge_i:
        #     return []

        # if not 0 <= start_edge_i < n:
        #     raise IndexError(f"Invalid start edge index: {start_edge_i}")

        # if not 0 <= end_edge_i < n:
        #     raise IndexError(f"Invalid end edge index: {end_edge_i}")

        # between = []
        # current_point_i = (start_edge_i + 1) % n

        # while True:
        #     between.append( tuple(points[current_point_i]) )

        #     if current_point_i == end_edge_i:
        #         break

        #     current_point_i = (current_point_i + 1) % n

        # return between

    #
    # does the clip polygon intersect this edge?
    #
    def clip_intersects_edge(self, edge_index, _order, edge_dict):

        if len( edge_dict["edges"][edge_index]['clip'] ) > 0:

            orig_v = edge_dict["edges"][edge_index]["orig"][_order]
            clip_p = edge_dict["edges"][edge_index]["clip"][0][_order]

            if self.same_point(orig_v, clip_p):
                return None 
            else:
                return clip_p   
        
        else: 
            return None


    #
    # for ordering purposes we need to know if we're heading INTO the clip or if 
    #   we're heading OUT OF the clip
    #
    #
    def heading_into_clip(self, edge_local_index, tri_edge_list, curr_intersection, clip_poly_list):

        print(tri_edge_list)

        edge_verts = self.tri_edge_index[tri_edge_list[edge_local_index]["index"]]
        next_vertex_rel_edge = 1 if tri_edge_list[edge_local_index]["ccw_order"] else 0

        end_vert_index = edge_verts[next_vertex_rel_edge]
        end_CR = end_vert_index % self.num_cols, end_vert_index // self.num_cols
        
        poly_id, edge_id = curr_intersection["poly_id"]

        print(f"poly: {poly_id}, edge: {edge_id}")

        num_edges_clip_poly = len(clip_poly_list[poly_id])

        poly_p1 = clip_poly_list[poly_id][edge_id]

        poly_p2 = clip_poly_list[poly_id][(edge_id+1)%num_edges_clip_poly]

        print(f"poly end: {poly_p2}, start: {poly_p1}")

        # switched so we get CCW!
        #poly_edge_dir = glm.vec2(poly_p2) - glm.vec2(poly_p1)
        poly_edge_dir = glm.vec2(poly_p1) - glm.vec2(poly_p2)

        end_vert = glm.vec2( mesh.get_point(*end_CR) )

        print(f"end vert: { tuple(mesh.get_point(*end_CR)) }, intersection: {tuple(curr_intersection["intersection"])}")

        _dir = end_vert -  glm.vec2( curr_intersection["intersection"] )

        _cross = self.cross_2D(poly_edge_dir, _dir)

        return _cross > 0


    def get_next_edge_intersection(self, curr_intersection, base_intersection, tri_edge_list, all_tri_edges, poly_points ):

        edge_local_order = curr_intersection["edge_local_order"]

        print(tri_edge_list)
        print(curr_intersection["tri_edge"] )
        assert tri_edge_list[edge_local_order]["index"] == curr_intersection["tri_edge"]

        next_found = False

        # can we just go to next intersection?
        edge_ccw = tri_edge_list[edge_local_order]["ccw_order"]

        tri_edge = curr_intersection["tri_edge"]

        assert tri_edge in all_tri_edges

        num_inters = len( all_tri_edges[ tri_edge ] )

        if edge_ccw:
            if curr_intersection["edge_order"] < num_inters-1:
                new_edge_intersection_order = curr_intersection["edge_order"] + 1 
                return all_tri_edges[tri_edge][new_edge_intersection_order]

        else: # go backwards
            if curr_intersection["edge_order"] > 0:
                new_edge_intersection_order = curr_intersection["edge_order"] - 1
                return all_tri_edges[tri_edge][new_edge_intersection_order]

        #
        #
        #    
        while(not next_found):

            # next edge
            edge_local_order = (edge_local_order+1)%3

            tri_edge = tri_edge_list[edge_local_order]["index"]
            edge_ccw = tri_edge_list[edge_local_order]["ccw_order"]

            print(f"edge: {tri_edge}")

            #
            # add relevant triangle vertex
            #
            vertex_rel_edge = 0 if edge_ccw else -1
            edge_tup = self.tri_edge_index[tri_edge]
            
            ccw_vertex_index = edge_tup[vertex_rel_edge]

            ccw_vert_cr = ccw_vertex_index % self.num_cols, ccw_vertex_index // self.num_cols
            ccw_vert = self.get_point(*ccw_vert_cr)[:2]

            poly_points.append({"point": ccw_vert, "type": "vert", "pix": ccw_vert_cr } )

            if tri_edge in all_tri_edges:
                num_inters = len( all_tri_edges[ tri_edge ] )
            else:
                num_inters = 0

            if num_inters > 0:

                # start from front
                if edge_ccw:
                    new_edge_intersection_order = 0

                # start from back
                else:
                    new_edge_intersection_order = num_inters - 1

                return all_tri_edges[tri_edge][new_edge_intersection_order]

            else:
                continue # try next edge

    def verify_edges(self, tri_index):

        tri_edge_list = mesh.get_edge_indices(tri_index)

        edge_0_tup = self.tri_edge_index[ tri_edge_list[0]["index"] ]
        edge_1_tup = self.tri_edge_index[ tri_edge_list[1]["index"] ]
        edge_2_tup = self.tri_edge_index[ tri_edge_list[2]["index"] ]
        
        if tri_index % 2 == 0: #even
            # *---*
            #  \  |
            #   \ | 
            #    \|  
            #     *
            assert edge_0_tup[0] == edge_2_tup[0]
            assert edge_0_tup[1] == edge_1_tup[1]
            assert edge_1_tup[0] == edge_2_tup[1]
        else:
            # *
            # |\
            # | \
            # |  \
            # *---*
            assert edge_0_tup[0] == edge_1_tup[0]
            assert edge_0_tup[1] == edge_2_tup[1]
            assert edge_1_tup[1] == edge_2_tup[0]

    #
    # clip_poly_list needs to be list of lists of points in CCW order 
    #
    def get_poly(self, tri_index, tri_edge_list, base_intersection, tri_edge_processed_dict, all_tri_edges,
                              all_clip_edges, clip_poly_list):
        

        curr_intersection = base_intersection

        # 
        # 0, 1, or 2?
        curr_which_edge = curr_intersection[ "edge_local_order" ]

        headin_ta_clip = self.heading_into_clip( curr_which_edge, tri_edge_list, curr_intersection, clip_poly_list  )

        print(f"heading into the clip? : {headin_ta_clip}")

        if not headin_ta_clip:
            return None
        
        curr_edge_order = curr_intersection["edge_order"]

        if (tri_edge_processed_dict[tri_index][curr_which_edge][curr_edge_order]):

            # it's already processed!
            return None;
    
        else:

            tri_edge_processed_dict[tri_index][curr_which_edge][curr_edge_order] = True

        poly_points = []

        curr_edge_uniq_id = curr_intersection["tri_edge"]

        at_end = False 

        print(f"tri index (!!!): {tri_index}")

        while(not at_end):

            if curr_edge_order == 0:
                vertex_rel_edge = 0 if tri_edge_list[curr_which_edge]["ccw_order"] else -1
                edge_tup = self.tri_edge_index[curr_edge_uniq_id]

                ccw_vertex_index = edge_tup[vertex_rel_edge]

                ccw_vertex_index
                ccw_vert_cr = ccw_vertex_index % self.num_cols, ccw_vertex_index // self.num_cols
                ccw_vert = self.get_point(*ccw_vert_cr)[:2]

                # poly_points.append( ccw_vert )
                
                # is the first one is clipped then just go to 
                # next intersection, this allows us to always
                # first ride the clip 
                #
                # by previously checking for the current edge going into 
                #  we know it's not clipped

            print(f"CURR INTERSECTION: {curr_intersection}")

            poly_points.append({"point": curr_intersection["intersection"] \
                                , "type": "intersection"})

            print(f"size poly points before hitting clip edge: {len(poly_points)}")

            curr_intersection = self.get_clip_points_between(clip_poly_list, curr_intersection, all_clip_edges, poly_points)

            print(f"POST CLIP INTERS: {curr_intersection}")

            poly_points.append({"point": curr_intersection["intersection"]
                                    , "type": "intersection"})

            print(f"size poly points before hitting clip edge: {len(poly_points)}")

            curr_intersection = self.get_next_edge_intersection( curr_intersection, base_intersection, tri_edge_list, all_tri_edges, poly_points )

            if curr_intersection == base_intersection:
                at_end = True

        return poly_points

    #
    # this function returns the unclipped shape at the vertex
    #
    #
    def get_all_non_clip_poly(self, tri_index, edge_dict, clip_poly_dict): 

        if tri_index % 2 == 0:
            #
            # each represents edge of corresponding 
            # vert index and order ( 0 or 1 currently )
            edge_vert_map_dict = {\
                0: [(2,0),(0,0), (1,1), (2,1)],\
                1: [(0,1),(1,1), (2,1), (0,0)],\
                2: [(1,0),(2,1), (0,0), (1,1)]       
            }
            #vert_order_list = [0,2,1]
        else:
            edge_vert_map_dict = {\
                0 :[(0,0),(1,0),(2,0),(0,1)],\
                1: [(1,1),(2,0),(0,1),(1,0)],\
                2: [(2,1),(0,1),(1,0),(2,0)] 
            } 
            #vert_order_list = [0,1,2]


        #
        # invariantly Even or Odd, we know triangle V3rt
        #   0 : 0 of 0 edge
        #   1: 1 of 1 edge
        #   2: 1 of 2 edge
        #
        # Do we store these as 2D or 3D
        #
        verts = [ edge_dict["edges"][0]["orig"][0][:2]
                , edge_dict["edges"][1]["orig"][1][:2]
                , edge_dict["edges"][2]["orig"][1][:2] ]

        polygons = []    
        for vert_i in range(3):

            # CCW order
            index1, _order1 = edge_vert_map_dict[vert_i][0] # curr
            index2, _order2 = edge_vert_map_dict[vert_i][1] # next
            index3, _order3 = edge_vert_map_dict[vert_i][2] # opposiTe
            index4, _order4 = edge_vert_map_dict[vert_i][3] # same

            print(f"index1, _order1: {index1, _order1}")
            print(f"index2, _order2: {index2, _order2}")
            print(f"index3, _order3: {index3, _order3}")
            print(f"TRI Index: {tri_index},  V INDEX: {vert_i} ")
            # print(edge_vert_map_dict)

            #
            # is vertex clipped?
            if edge_dict["edges"][index1]["bool"][_order1]:

                # other thing says to clip too
                assert(edge_dict["edges"][index2]["bool"][_order2])
                continue; # next vert!

            else: # not clipping this vertex!
                points = []

                # does have clipping intersection?
                curr_edge_intersect_p = self.clip_intersects_edge( index1, _order1 ,edge_dict)
                next_edge_intersect_p  = self.clip_intersects_edge( index2, _order2,edge_dict)
                opposite_edge_intersect_p  = self.clip_intersects_edge( index3, _order3, edge_dict)
                same_edge_intersect_p = self.clip_intersects_edge( index4, _order4, edge_dict)
                
                if curr_edge_intersect_p is None:
                    continue;
                else:
                    # first edge intersection
                    points.append(curr_edge_intersect_p)

                    # now vertex itself
                    points.append(verts[vert_i])
                    curr_edge_clip_intercept_tup = self.clip_index_intercept( index1, _order1 , edge_dict )

                    if next_edge_intersect_p is not None:
                        #
                        # what is the clip_intercept?
                        # 
                        next_edge_clip_intercept_tup = self.clip_index_intercept( index2, _order2 , edge_dict )
                        
                        print(f"\tCASE 1: tri: {tri_index}, vert: {vert_i}, curr_edge_clip_ {curr_edge_clip_intercept_tup}, next_edge_kl: {next_edge_clip_intercept_tup}")
                        
                        points.append(next_edge_intersect_p)

                        #print(f"case1 Adding {len(poly_verts[curr_edge_clip_intercept_i:next_edge_clip_intercept_i])} clip points")
                        
                        intermed_clip = self.get_points_between(curr_edge_clip_intercept_tup, next_edge_clip_intercept_tup, clip_poly_dict)

                        #
                        # notice it's backwards!
                        #
                        points.extend( intermed_clip[::-1] )
                        polygons.append(points)

                        continue; # next iteration!!!
                    
                    #
                    # next edge is intercepted?!
                    #
                    elif opposite_edge_intersect_p is not None:
                        pass;pass;pass;
                        pass;

                        points.append(verts[(vert_i+1)%3])
                        points.append( opposite_edge_intersect_p )

                        oppo_edge_clip_intercept_i = self.clip_index_intercept( index3, _order3 , edge_dict )
                        
                        intermed_clip = self.get_points_between(curr_edge_clip_intercept_tup, oppo_edge_clip_intercept_i, clip_poly_dict)

                        print(f"number of intermediary clip points: {len(intermed_clip)}")
                        points.extend( intermed_clip[::-1] )
                        polygons.append(points)

                        print(f"\tCASE 2: tri: {tri_index}, vert: {vert_i}, curr_edge_clip_ {curr_edge_clip_intercept_tup}, opposite_edge_kl: {oppo_edge_clip_intercept_i}")
                        continue;
                        #print(f"case2 Adding {len(poly_verts[curr_edge_clip_intercept_i:next_edge_clip_intercept_i])} clip points")
                      
                    #
                    # there cannot be a single interception!
                    # this must be the case
                    #
                    else:
                        assert(same_edge_intersect_p is not None)
                        same_edge_clip_intercept_i = self.clip_index_intercept( index4, _order4 , edge_dict ) 
                        print(f"\tCASE 3: tri: {tri_index}, vert: {vert_i}, curr_edge_clip_ {curr_edge_clip_intercept_tup}, SAME_edge_kl: {same_edge_clip_intercept_i}")
                        
                        intermed_clip = self.get_points_between(curr_edge_clip_intercept_tup, same_edge_clip_intercept_i, clip_poly_dict)

                        # not backwards!
                        points.extend( intermed_clip)

                        points.append(same_edge_intersect_p)

                        polygons.append(points) 
                        break; 
        
        return polygons


    def get_neighbors(self, c, r):

        #
        #  0---1---2
        #  |       |
        #  3---v---4
        #  |       |
        #  5---6---7
        #
        index0 = (c-1, r-1)
        index1 = (c, r-1)
        index2 = (c+1, r-1)
        index3 = (c-1, r)
        index4 = (c+1, r)
        index5 = (c-1, r+1)
        index6 = (c, r+1)
        index7 = (c+1, r+1)

        topRight = c == self.num_cols - 1 and r == 0
        topLeft = c == 0 and r == 0
        bottomRight = c == self.num_cols - 1 and r == self.num_rows - 1
        bottomLeft = c == 0 and r == self.num_rows - 1

        neighbors = []

        if topRight:
            neighbors.append(index3);
            neighbors.append(index5);
            neighbors.append(index6);

        elif (topLeft):
            neighbors.append(index6);
            neighbors.append(index7);
            neighbors.append(index4);

        elif (bottomRight):
            neighbors.append(index3);
            neighbors.append(index0);
            neighbors.append(index1);

        elif (bottomLeft):
            neighbors.append(index1);
            neighbors.append(index2);
            neighbors.append(index4);

        # right most col
        elif c == self.num_cols - 1:

            neighbors.append(index1);
            neighbors.append(index0);
            neighbors.append(index3);
            neighbors.append(index5);
            neighbors.append(index6);
        
        # left most col
        elif (c == 0):
            neighbors.append(index1);
            neighbors.append(index2);
            neighbors.append(index4);
            neighbors.append(index7);
            neighbors.append(index6);
        
        # bottom row
        elif r == self.num_rows - 1:
            neighbors.append(index3);
            neighbors.append(index0);
            neighbors.append(index1);
            neighbors.append(index2);
            neighbors.append(index4);
        
        #  top row
        elif r == 0:
            neighbors.append(index3);
            neighbors.append(index5);
            neighbors.append(index6);
            neighbors.append(index7);
            neighbors.append(index4);
        # interior
        else:
            neighbors.append(index0);
            neighbors.append(index1);
            neighbors.append(index2);
            neighbors.append(index3);
            neighbors.append(index4);
            neighbors.append(index5);
            neighbors.append(index6);
            neighbors.append(index7);

        return neighbors;

    def get_vert_index(self, pix_X, pix_Y):
        return int( pix_Y * self.num_cols + pix_X )


    def get_point(self, indx1, indx2=None):

        if indx2 is None:
            c =  indx1 % self.num_cols  
            r = indx1 // self.num_cols

        else:
            c = int(indx1)
            r = int(indx2)

        _x = self.buffer_array[c, r, 0].item()
        _y = self.buffer_array[c, r, 1].item()
        elev = self.buffer_array[c, r, 2].item()

        return (_x, _y, elev)
    

    def combine_tri_dict(self, tri_d1, tri_d2, poly_dict):

        combined = tri_d1.copy()

        for tri_index in tri_d1.keys() | tri_d2.keys():

            if tri_index not in tri_d1:
                combined[tri_index] = tri_d2[tri_index]
                combine_both = False
            elif tri_index in tri_d1 and tri_index in tri_d2:

                combine_both = True

            edges = mesh.get_edges(tri_index)

            #
            # we are manually adding orig unclipped edges
            # for tris that otherwise have clipping
            for i in range(3):

                # i for edge ONLY
                #

                if combined[tri_index]["edges"][i]["orig"] is None and tri_d2[tri_index]["edges"][i]["orig"] is None:
                    combined[tri_index]["edges"][i]["orig"] = edges[i]
                    combined[tri_index]["edges"][i]["bool"] = [False, False]

                if combine_both:
                    
                    d1_clip = combined[tri_index]["edges"][i]["clip"]
                    d2_clip = tri_d2[tri_index]["edges"][i]["clip"]
                    
                    combined_clip = d1_clip + d2_clip

                    d1_clip_index = combined[tri_index]["edges"][i]["clip_index"]
                    d2_clip_index = tri_d2[tri_index]["edges"][i]["clip_index"]
                    
                    combined_clip_index = d1_clip_index + d2_clip_index

                    print(f"combo clip: {combined_clip}")
                    print(f"combo clip index: {combined_clip_index}")
                    order = sorted(range(len(combined_clip)), key=lambda _k: self.clip_order(combined_clip[_k]))

                    ordered_clip = [combined_clip[_i] for _i in order]
                    ordered_clip_index_tups = [combined_clip_index[_i] for _i in order]

                    filt_ordered_clip_index_tups = [];
                    filt_ordered_clip = [];
                    itr = 0

                    #
                    # combining  elem[i][0] with elem[i+1][1] [[(_,_), (_,_)],[(_,_), (_,_)],[(_,_), (_,_)]]
                    #  based on conditionally "forbidden" edge indices of 
                    #
                    while itr < len(ordered_clip_index_tups):
                        

                        curr_poly_id, curr_edge_i = ordered_clip_index_tups[i][1]
                        # combine conditiones
                        if itr + 1 < len(ordered_clip_index_tups) and \
                            curr_edge_i in poly_dict[curr_poly_id]["ignore"]:

                            curr_clip_index_tups = ordered_clip_index_tups[itr]
                            next_clip_index_tups = ordered_clip_index_tups[itr+1]
                            print(f"curr_poly_id: {curr_poly_id}, curr_edge_i: {curr_edge_i}")
                            filt_ordered_clip_index_tups.append([curr_clip_index_tups[0],next_clip_index_tups[1]])

                            curr_clip = ordered_clip[itr]
                            next_clip = ordered_clip[itr+1]

                            filt_ordered_clip.append([curr_clip[0],next_clip[1]])

                            itr += 2  # both pairs were consumed
                        else:
                            filt_ordered_clip_index_tups.append(ordered_clip_index_tups[itr])
                            filt_ordered_clip.append( ordered_clip[itr] )
                            itr += 1
                    
                    #forbidden_clip = [i for i,qw in enumerate(ordered_clip_index) if qw[0] == 0 or  qw[0] == 0 ]

                    combined[tri_index]["edges"][i]["clip_index"] = filt_ordered_clip_index_tups
                    combined[tri_index]["edges"][i]["clip"] = filt_ordered_clip

                    #
                    # if one poly clips then vert is clipped
                    #
                    other_bool = tri_d2[tri_index]["edges"][i]["bool"] if tri_d2[tri_index]["edges"][i]["orig"] is not None else [False,False]
                    
                    combined[tri_index]["edges"][i]["bool"] = [
                        combined[tri_index]["edges"][i]["bool"][0] | other_bool[0],
                        combined[tri_index]["edges"][i]["bool"][1] | other_bool[1]
                    ]


        return combined
    
    #
    # since we're assuming the edge clips are non-overlapping
    # we can sory be the first point in each!
    #
    def clip_order(self, points_pair):

        print(f"points_pair: {points_pair}")
        return (~self.transform) * points_pair[0]

    #
    #  don't feed it a border vertex! !
    # "higher" means increasing pixel
    def get_next_edge_grat_dir(self
                                , vert_index
                                , grat_type
                                # False : higher to lower pixel
                                # True: lower to higher pixel
                                , direction ):

        edges = self.get_edges_vert( vert_index )

        match grat_type:
            #
            # *---*---*
            # |\  |\  |
            # | 4 5 \ |
            # |  \|  \|
            # *-3-C-0-*
            # |\  |\  |
            # | \ 2 1 |
            # |  \|  \|
            # *---*---*
            case Graticule.LON:
                edge_indx_rel_tri = [2,5]
            case Graticule.LAT:
                edge_indx_rel_tri = [0,3]
            case Graticule.DIAG:
                edge_indx_rel_tri = [1,4]

        # lower to higher pixel
        if (direction):
            i = 0
        else: # higher to lower pixel
            i = 1

        return edges[edge_indx_rel_tri[i]]


    def inside_clip(self, intersection, poly_list, curr_vert):
        curr_poly_index, curr_poly_edge = intersection.clip_poly_id

        poly = poly_list[curr_poly_index]
        n = len(poly)

        curr_clip_point = poly[curr_poly_edge]
        next_clip_point = poly[(curr_poly_edge + 1) % n]
        prev_clip_point = poly[(curr_poly_edge - 1) % n]

        curr_point = self.get_point(curr_vert)[:2]
        curr_vert_glm = glm.vec2(curr_point)

        if intersection.clip_feature == FeatureType.VERTEX:
            # Assuming curr_poly_edge is the index of the clip vertex hit.
            v = glm.vec2(curr_clip_point)
            prev_v = glm.vec2(prev_clip_point)
            next_v = glm.vec2(next_clip_point)

            incoming_edge_dir = v - prev_v      # prev -> vertex
            outgoing_edge_dir = next_v - v      # vertex -> next

            # Direction from the clip vertex back toward curr_vert
            rel = curr_vert_glm - v

            side_incoming = self.cross_2D(incoming_edge_dir, rel)
            side_outgoing = self.cross_2D(outgoing_edge_dir, rel)

            # CW polygon: inside is right side of BOTH adjacent edges.
            inside_incoming = side_incoming <= 0
            inside_outgoing = side_outgoing <= 0

            return inside_incoming and inside_outgoing

        else:
            intersection_glm = glm.vec2(intersection.point)

            print(f"next clip: {next_clip_point}, curr clip: {curr_clip_point}")
            print(f"curr vert: {curr_vert}, inter: {intersection.point}")
            
            clip_dir = glm.vec2(next_clip_point) - glm.vec2(curr_clip_point)
            mesh_dir = curr_vert_glm - intersection_glm

            side = self.cross_2D(clip_dir, mesh_dir)
            print(f"side: {side}")
            return side <= 0



    #
    # returns True for being in clip , F otherwse
    #
    def follow_geodesic_clip(self, direction # False : higher to lower pixel
                                            # True: lower to higher pixel
                             , grat_type
                             , vert_index
                             , poly_list
                            ,  prev_vert = None ):


        #
        # a vert by itself intersected 
        #   by a clip is considered clipped
        #
        if vert_index in self.all_mesh_verts:
            if prev_vert is None:
                return True
            else:
                print(f"checking vertex intersection!")
                intersection = self.all_mesh_verts[vert_index][0]
                print(f"intersection: {intersection}")
                return self.inside_clip(intersection, poly_list, prev_vert)


        _pixel_x,_pixel_y = (int(vert_index % self.num_cols), 
                    int(vert_index // self.num_cols))

        #
        # reached border without being clipped
        #
        if ((_pixel_x == 0 and not direction) 
            or (_pixel_x == self.num_cols - 1 and direction)) and \
            (grat_type == Graticule.LAT or grat_type == Graticule.DIAG):
            return False
        if (( _pixel_y == 0 and not direction) or 
            ( _pixel_y == self.num_rows - 1 and direction)) and \
                (grat_type == Graticule.LON or grat_type == Graticule.DIAG):
            return False

        print(f"getting edge for vert: {vert_index}, dir: {direction}")

        curr_edge_tup = self.get_next_edge_grat_dir( vert_index, grat_type, direction )

        curr_edge_indx = curr_edge_tup[0]

        print(f"curr_edge_tup: {curr_edge_tup}")

        curr_edge_verts = self.tri_edge_index[curr_edge_indx] 

        # False : higher to lower pixel
        # True: lower to higher pixel
        if direction:
            next_vert = curr_edge_verts[1]
            curr_vert = curr_edge_verts[0]
        else:
            next_vert = curr_edge_verts[0]
            curr_vert = curr_edge_verts[1]

        print(f" direction: {direction} current vertex: {vert_index} edge vertices: {curr_edge_verts}")

        assert vert_index == curr_vert

        if curr_edge_indx in self.all_mesh_edges: #or curr_vert in vertices_dict:

            # True: lower to higher pixel
            # False : higher to lower pixel
            if direction:
                intersection = self.all_mesh_edges[curr_edge_indx ][0]
            else:
                intersection = self.all_mesh_edges[curr_edge_indx ][-1]

            return self.inside_clip(intersection, poly_list, curr_vert)

        else:

            return self.follow_geodesic_clip( direction, \
                                             grat_type, \
                                             next_vert, \
                                             poly_list, \
                                             curr_vert )




    #
    #
    #
    def sort_inter_edge(self, _inter):

        edge_index = _inter.tri_edge_index

        point = _inter.intersection

        print(f"edge index: {edge_index}")

        te0,te1 = self.tri_edge_index[edge_index]

        pix0_x, pix0_y = te0 % self.num_cols, te0 // self.num_cols  
        pix1_x, pix1_y = te1 % self.num_cols, te1 // self.num_cols 

        #
        # TODO: why are our points indexed like this? get in line w/ openGL
        #
        _e0x, _e0y = self.buffer_array[pix0_x, pix0_y, 0], self.buffer_array[pix0_x, pix0_y, 1]
        _e1x, _e1y = self.buffer_array[pix1_x, pix1_y, 0], self.buffer_array[pix1_x, pix1_y, 1]

        dx = _e1x - _e0x
        dy = _e1y - _e0y

        length_squared = dx * dx + dy * dy

        if length_squared == 0:
            raise ValueError("0-length edge detected!")

        return ((point[0] - _e0x) * dx + (point[1] - _e0y) * dy) / length_squared

    #
    # to make default dict referenced by
    # tri_edge (of all unique tri edges)
    def make_both_edge_dict(self):
        return { "intersections": [] }       

    def make_edges_dict(self):
        return {
             "edges" : {
                 edge_index: {
                     "orig": None, # should be 2 points (tuples) this is also how we check
                                   # for edge initialization
                     "bool": [], # should be 2 bool (did we clip vert?)
                     "clip": [], # should be list of list of tuples 
                                 # so edge can be clipped > 1
                     "clip_index" : [], # should be 2 indices that represent
                                        # the CCW indices of the clipping 
                                        # polygon
                 }
                 for edge_index in range(3)
             }
        }

    def get_left_right_vectors(self, prev_vec, next_vec):
        #left_rot_matrix = glm.mat2(0, 1, -1, 0)
        #right_rot_matrix = glm.mat2(0, -1, 1, 0)

        #return glm.normalize(left_rot_matrix * vec2), glm.normalize(right_rot_matrix * vec2)
        #return  glm.normalize(glm.vec2(-vec2.y, vec2.x)), glm.normalize(glm.vec2(vec2.y, -vec2.x))

        if prev_vec is None:
            avg_vec = next_vec
        elif next_vec is None:
            avg_vec = prev_vec
        else:
            avg_vec = (next_vec+prev_vec)/2

        #
        # image coordinates go down!
        #
        #return glm.normalize(left_rot_matrix * vec2), glm.normalize(right_rot_matrix * vec2)
        return glm.normalize(glm.vec2(avg_vec.y, -avg_vec.x)), glm.normalize(glm.vec2(-avg_vec.y, avg_vec.x))

    def create_quad_from_points(self, rl_points, rl_index, width=DEFAULT_ROAD_WIDTH):

        quad_4points = []
        for i in range(rl_index[0], rl_index[1]+1):

            currZ = rl_points[i][2]
            currPoint = glm.vec2( rl_points[i] )
    
            if (i == 0):#first
                prevVec = None 
                nextVec = glm.vec2( rl_points[i+1] ) - currPoint

            elif (i == rl_index[1]): # last
                prevVec = currPoint - glm.vec2( rl_points[i-1] )
                nextVec = None

            else:
                prevVec = currPoint - glm.vec2(rl_points[i-1])
                nextVec = glm.vec2( rl_points[i+1] ) - currPoint


            left_vec, right_vec = self.get_left_right_vectors(prevVec, nextVec)

            print(f"left,R vector: {left_vec}, {right_vec}")             
            l_point = currPoint  + (left_vec * width )
            r_point = currPoint  + (right_vec * width )


            l_point = (*l_point, currZ)
            r_point = (*r_point, currZ)
            

            quad_4points.extend([l_point, r_point])

        # I don't think this is right, CW order depends on the direction 
        # of the points!
        #points_CCW =  [   quad_4points[revers_i][1] for revers_i in range(len(quad_4points)-1, -1, -1)  ] \
        #                +  [   quad_4points[forwar_i][0] for forwar_i in range(len(quad_4points))  ]
        #
        index_CW = list( range( len(quad_4points) ) )

        points_CW = self.order_cw(quad_4points)

        return points_CW, index_CW

    #
    # returns tuple of vec2
    def get_edge_points(self, edge_index):

        edge = self.tri_edge_index[  edge_index ]
        vert1_pix = edge[0] % self.num_cols, edge[0] // self.num_cols
        vert2_pix = edge[1] % self.num_cols, edge[1] // self.num_cols

        vert1 = self.get_point(*vert1_pix)[:2]
        vert2 = self.get_point(*vert2_pix)[:2]

        return (glm.vec2(vert1), glm.vec2(vert2))


    # returns vertices, x_center, y_center
    def get_points_centered_scaled(self, scale):

        points = []
        cnt = 0

        for r in range(self.num_rows):
            for c in range(self.num_cols):
     
                p = self.buffer_array[c, r , 0:3]

                # UNIFORM SCALING SO EVERYTHING HAS SAME NORMALS!

                new_z = p[2] * scale # not centering Z

                points.append( ( ( p[0]-self.x_center ) * scale, ( p[1]-self.y_center)* scale, new_z ) )
        return points

    def get_points(self):

        x_max, y_max, z_max = -math.inf, -math.inf, -math.inf
        x_min, y_min, z_min = math.inf, math.inf, math.inf

        points = []

        for r in range(self.num_rows):
            for c in range(self.num_cols):

                p = self.buffer_array[c, r , 0:3]
                points.append( ( p[0].item(), p[1].item(), p[2].item() ) )

                x_max = max(x_max, p[0])
                y_max = max(y_max, p[1])
                z_max = max(z_max, p[2])
                x_min = min(x_min, p[0])
                y_min = min(y_min, p[1])
                z_min = min(z_min, p[2])

        x_center = (x_max + x_min) / 2
        y_center = (y_max + y_min) / 2
        z_center = (z_max + z_min) / 2

        center = (x_center, y_center, z_center)

        return points, center

    #
    #
    #
    def create_face_obj(self, scale, path, tri_clip_dict):

        vertices, center = self.get_points( )

        curr_vert_index = len(vertices)

        faces = []
        for i, tri in enumerate(self.index_array):

            if i in tri_clip_dict: # it's clipped!
                
                polys = tri_clip_dict[i]

                for poly in polys:
                    poly_str = "f"

                    for p in poly:

                        if p["type"] == "vert":
                            pix = p["pix"]
                            indx = int( pix[0] + (pix[1] * self.num_cols) ) +1 # 1-based!
                            poly_str+=" "+str(indx)
                        else:
                            poly_str+=" "+str(curr_vert_index+1)
                            vertices.append(p["point"])
                            curr_vert_index+=1

                    faces.append( poly_str )
                
            else: # not clipped

                faces.append( f"f {tri[0]+1} {tri[1]+1} {tri[2]+1}") # 1 based!

        #
        # center, scale
        vertices = [(scale*(vert[0] - center[0]), scale*(vert[1] - center[1]), scale*vert[2]) for vert in vertices]

        with open(path + ".obj", "w") as f:
            f.write("# OBJ DEM mesh w/ clipped road\n")

            for v in vertices:
                f.write(f"v {v[0]} {v[1]} {v[2]}\n")

            f.write( "\n".join(faces) )


    #
    # Pixel must be within raster bounds, can be at the borders tho
    #
    def get_candidate_tris(self, pixel_x, pixel_y):
    
        if (pixel_x == self.num_cols - 1 and pixel_y == self.num_rows - 1):
            tri1_v1_st = (pixel_x, pixel_y)
            tri1_v2_st = (pixel_x-1, pixel_y-1)
            tri1_v3_st = (pixel_x-1, pixel_y)
            tri1_sts = [tri1_v1_st, tri1_v2_st, tri1_v3_st]
            tri1_coords = [self.get_point(*st) for st in tri1_sts]

            #print(tri1_coords)
            tri1 = Triangle(tri1_coords, tri1_sts)
            return tri1, None

        elif pixel_x == self.num_cols - 1:
            tri1_v1_st = (pixel_x, pixel_y)
            tri1_v2_st = (pixel_x, pixel_y+1)
            tri1_v3_st = (pixel_x-1, pixel_y)

            tri1_sts = [tri1_v1_st, tri1_v2_st, tri1_v3_st]
            tri1_coords = [self.get_point(*st) for st in tri1_sts]

            #print(tri1_coords)
            tri1 = Triangle(tri1_coords, tri1_sts)
            return tri1, None

        # bottom
        elif pixel_y == self.num_rows - 1:
            tri1_v1_st = (pixel_x, pixel_y)
            tri1_v2_st = (pixel_x+1, pixel_y)
            tri1_v3_st = (pixel_x, pixel_y-1)

            tri1_sts = [tri1_v1_st, tri1_v2_st, tri1_v3_st]
            tri1_coords = [self.get_point(*st) for st in tri1_sts]

            #print(tri1_coords)
            tri1 = Triangle(tri1_coords, tri1_sts)
            return tri1, None

        else:
            tri1_v1_st = (pixel_x, pixel_y)
            tri1_v2_st = (pixel_x+1, pixel_y)
            tri1_v3_st = (pixel_x+1, pixel_y+1)
        
            tri1_sts = [tri1_v1_st, tri1_v2_st, tri1_v3_st]
            tri1_coords = [self.get_point(*st) for st in tri1_sts]

            print(tri1_coords)
            tri1 = Triangle(tri1_coords, tri1_sts)
        
            tri2_v1_st = (pixel_x, pixel_y)
            tri2_v2_st = ((pixel_x+1), (pixel_y+1))
            tri2_v3_st = (pixel_x, (pixel_y+1))
        
            tri2_sts = [tri2_v1_st, tri2_v2_st, tri2_v3_st]

            tri2_coords = [self.get_point(*st) for st in tri2_sts]
        
            tri2 = Triangle(tri2_coords, tri2_sts)
        
            return tri1, tri2

    #
    # not floored!
    def get_pixel(self, point):
        _x_m, _y_m = point
        pixel_y = ( _y_m - self.transform.f ) / self.transform.e;
        pixel_x = ( _x_m - ( self.transform.c ) ) / self.transform.a;

        return pixel_x, pixel_y
                
    #
    #
    #
    def interp_elev(self, point):

        print(f"Point: {point}")
        _x_m, _y_m = point
        pixel_x, pixel_y = self.get_pixel(point)
        pixel_x, pixel_y = (math.floor(pixel_x), math.floor(pixel_y))

        cand_tri1, cand_tri2 = self.get_candidate_tris( pixel_x, pixel_y )

        # happens if on border
        if cand_tri2 == None:

            tri = cand_tri1
            stu = tri.bary_x_y(_x_m, _y_m)

        else:
            tri1_stu = cand_tri1.bary_x_y(_x_m, _y_m)
            print(f"TRI1 BARYCENTRIC (curr point): {tri1_stu[0]} {tri1_stu[1]} {tri1_stu[2]}")
            
            tri2_stu = cand_tri2.bary_x_y(_x_m, _y_m)

            #print(f"TRI1 BARYCENTRIC (curr point): {tri1_stu[0]} {tri1_stu[1]} {tri1_stu[2]}")
            print(f"TRI2 BARYCENTRIC: {tri2_stu[0]} {tri2_stu[1]} {tri2_stu[2]}")

            # round to nearest 10 decimal places for barycentric assertion
            assert( all(round(e,10) >= 0 for e in tri1_stu) \
                    or all(round(e,10) >= 0  for e in tri2_stu) )

            stu = tri1_stu if all(e >= 0  for e in tri1_stu) else tri2_stu
            tri = cand_tri1 if all(e >= 0  for e in tri1_stu) else cand_tri2


        interp_z =  stu[0] * tri.point_s[0][2] + stu[1] * tri.point_s[1][2] + stu[2] * tri.point_s[2][2];      
        
        return interp_z

    def is_vert_clipped( self, vert_index, clip_polys ):

        lat_up = self.follow_geodesic_clip( True 
                , Graticule.LAT
                , vert_index
                , clip_polys )
        lat_down = self.follow_geodesic_clip( False 
                , Graticule.LAT
                , vert_index
                , clip_polys )
        lon_up = self.follow_geodesic_clip( True 
                , Graticule.LON
                , vert_index
                , clip_polys )
        lon_down = self.follow_geodesic_clip( False 
                , Graticule.LON
                , vert_index
                , clip_polys )
        diag_up = self.follow_geodesic_clip( True 
                , Graticule.DIAG
                , vert_index
                , clip_polys )
        diag_down = self.follow_geodesic_clip( False 
                , Graticule.DIAG
                , vert_index
                , clip_polys )
        
        return any([ lat_up, lat_down, \
                lon_up, lon_down, \
                diag_up, diag_down ])

    #
    # have to run `perform_clipping` first!
    # i think it does
    def find_clipped_verts( self , clip_polys ):
        self.vert_clip_index = {}

        for row in range(self.num_rows):
            for col in range(self.num_cols):

                vert_index = row * self.num_cols + col

                print(f"processing vert: {vert_index}")

                vert_clipped_bool = self.is_vert_clipped( vert_index, clip_polys )

                self.vert_clip_index[vert_index] = vert_clipped_bool

    def within_bounds(self, point):
        return (
            self.bounds.left <= point[0] <= self.bounds.right
            and self.bounds.bottom <= point[1] <= self.bounds.top
        )


if __name__ == "__main__":

    cfg = load_config()

    TILE_SIZE = 1024

    #raster_file = "dem_EPSG_26910_542354_4944208.tif"
    raster_file = "roof_raster_fake3.tif"
    texture_file = "osip_EPSG_26910_542354_4944208.tif"

    raster_path = cfg.input_dir / raster_file
    texture_path = cfg.input_dir / texture_file

    #offset_x, offset_y = 90, 71

    offset_x, offset_y = 0,0 

    tex_offset_x , tex_offset_y = offset_x*10, offset_y*10

    dem_resolution = 10 #meters per pixel
    tex_resolution = 1 #meters per pixel
    
    #
    # edit this variable!
    # now we're showing the entire raster
    # 
    dem_num_cells = 2

    tex_row_offset = offset_y * 10#scale
    tex_col_offset = offset_x * 10;####scale


    # what is sqrt area of mesh
    dem_patch_size = dem_resolution * dem_num_cells

    # how many verts per side need? 
    #   e.g. 2 verts per side for 10x10 raster
    dem_patch_num_verts = dem_num_cells+1 # per side

    # each texel should be 1x1 meter!
    tex_patch_size = (10 * (dem_num_cells)) #+ 1

    tex_patch_num_verts = (10 * (dem_num_cells)) + 1
    #num_tex_rows = num_tex_cols = 1*10

    buffer_array = np.zeros((dem_patch_num_verts, dem_patch_num_verts, MeshComponentCreator.ENTRIES_PER_BUFFER), dtype=np.float32)

    num_tri = (dem_patch_num_verts - 1) * (dem_patch_num_verts-1) * 2

    num_tex_tri = (tex_patch_num_verts - 1) * (tex_patch_num_verts-1) * 2

    index_array = []

    tri_edge_array = []

    with rasterio.open(raster_path) as src:
        band = src.read(1)
        transf = src.transform

        bounds = src.bounds

        num_rows_total = band.shape[0]
        num_cols_total = band.shape[1]

        band = band[offset_y:offset_y+dem_patch_num_verts,offset_x:offset_x+dem_patch_num_verts]
        h = src.height
        w = src.width

        rows = np.arange( dem_patch_num_verts ) + offset_y
        cols = np.arange( dem_patch_num_verts ) + offset_x

        cols, rows = np.meshgrid(cols, rows)
        xs, ys = rasterio.transform.xy(src.transform, rows, cols, offset='ul')

        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
        x_min, x_max, y_min, y_max = extent  
        
        x_max = xs[-1];
        y_min = ys[-1];

        num_rows_inset = band.shape[0]
        num_cols_inset = band.shape[1]
    
        orig_transf = src.transform 
        
        new_transf =  orig_transf * Affine.translation(
            offset_x,   # column offset
            offset_y    # row offset
        )

        for r in range(num_rows_inset):
            for c in range(num_cols_inset):

                # should be in meters
                x_m, y_m = rasterio.transform.xy(new_transf, r, c, offset='ul')
                elev = band[r, c]

                if elev < 0:
                    elev = 0

                buffer_array[c, r, 0] = x_m #2ND
                buffer_array[c, r, 1] = y_m
                buffer_array[c, r, 2] = elev
                buffer_array[c, r, 3] = 0.0
                buffer_array[c, r, 4] = 0.0 # downslope angle local
                buffer_array[c, r, 5] = 0.0 # downslope angle global
                buffer_array[c, r, 6] = 0.0 # downslope slope
                buffer_array[c, r, 7] = 0.0 # downslope Bool
                buffer_array[c, r, 8] = 0.0 # downslope pixel x
                buffer_array[c, r, 9] = 0.0 # downslope pixel y
                buffer_array[c, r, 10] = 0.0 # downslope deviation x
                buffer_array[c, r, 11] = 0.0 # downslope deviation y

                buffer_array[c, r, 12] = 0.0 # number of upslope
                buffer_array[c, r, 13] = 0.0 # x of adjusted position
                buffer_array[c, r, 14] = 0.0 # y of adjusted position
                buffer_array[c, r, 15] = 0.0 # interp elev of adjusted position
                buffer_array[c, r, 16] = 0.0 # been visited by flow line drawing

                buffer_array[c, r, 17] = 0.0 # x normal
                buffer_array[c, r, 18] = 0.0 # y normal
                buffer_array[c, r, 19] = 0.0 # z normal
                buffer_array[c, r, 20] = 0.0

    tex_buffer_array = np.zeros((tex_patch_num_verts, tex_patch_num_verts, 3), dtype=np.float32)
    tex_index_array = []

            
    with rasterio.open(texture_path) as tex_src:

        tex_band = tex_src.read(1)
        tex_transf = tex_src.transform
        print(f"band shape: {tex_band.shape}")

        tex_band = tex_band[tex_row_offset : tex_row_offset + tex_patch_size ,
            tex_col_offset : tex_col_offset + tex_patch_size ]

        tex_h = tex_src.height
        tex_w = tex_src.width

        _rows = np.arange(tex_patch_num_verts) + tex_row_offset
        _cols = np.arange(tex_patch_num_verts) + tex_col_offset

        tex_cols, tex_rows = np.meshgrid(_cols, _rows)
        tex_xs, tex_ys = rasterio.transform.xy(tex_src.transform, tex_rows, tex_cols, offset='ul')

        xs = np.asarray(xs)
        ys = np.asarray(ys)

        tex_xs = np.asarray(tex_xs)
        tex_ys = np.asarray(tex_ys)


        for r in range(tex_patch_num_verts):
            for c in range(tex_patch_num_verts):

                #tex_x_m, tex_y_m = tex_xs[c], tex_ys[r] #rasterio.transform.xy(tex_transf, r, c)
                #
                tex_x_m = tex_xs[c]
                tex_y_m = tex_ys[r]

                #tex_buffer_array[c, r, 0] = tex_x_m #2ND
                #tex_buffer_array[c, r, 1] = tex_y_m
                #tex_buffer_array[c, r, 2] = elev


                if ( r < dem_patch_num_verts and c < dem_patch_num_verts):

                    #elev = band[r, c]
                    x_m, y_m = xs[c], ys[r] #rasterio.transform.xy(transf, r, c)

                    if (r ==0 and c ==0):
                        print(f"tex ({tex_x_m}, {tex_y_m }), dem: ({x_m}, {y_m })")


                    #buffer_array[c, r, 0] = x_m #2ND
                    #buffer_array[c, r, 1] = y_m
                   #buffer_array[c, r, 2] = elev


                if (r > 0 and c > 0):
                    currIndx = c + r * tex_patch_num_verts;
                    prevRowSameColIndx = c + (r - 1) * tex_patch_num_verts;
                    prevRowBackColIndx = c - 1 + (r - 1) * tex_patch_num_verts;
                    prevIndx = currIndx - 1;
                    tex_index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])
                    tex_index_array.append([prevRowBackColIndx, prevIndx, currIndx])


                    if ( r < dem_patch_num_verts and c < dem_patch_num_verts):
                        currIndx = c + r * dem_patch_num_verts;
                        prevRowSameColIndx = c + (r - 1) * dem_patch_num_verts;
                        prevRowBackColIndx = c - 1 + (r - 1) * dem_patch_num_verts;
                        prevIndx = currIndx - 1;
                        index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])
                        index_array.append([prevRowBackColIndx, prevIndx, currIndx])


                        #tri_edge_array.append([])
                        #tri_edge_array.append([])

        fig, ax = plt.subplots( figsize=(14, 14) )
    
        # ds = gdal.Open( str(raster_path) )
        # width = ds.RasterXSize
        # height = ds.RasterYSize
        # gt = ds.GetGeoTransform()

        # minx = gt[0]
        # maxy = gt[3]
        # maxx = gt[0] + width * gt[1] + height * gt[2]
        # miny = gt[3] + width * gt[4] + height * gt[5]
            
        dem_extent_d = raster_extent(src.transform, dem_num_cells,dem_num_cells, offset_x, offset_y)
        tex_extent_d = raster_extent(tex_src.transform, dem_num_cells*10, dem_num_cells*10, tex_offset_x, tex_offset_y )
        #im = ax.imshow(band, extent=[dem_extent_d["left"], dem_extent_d["right"], dem_extent_d["bottom"], dem_extent_d["top"]], cmap='viridis', origin='upper')

        print(f"index array size: {index_array}")

        tex_triang = mtri.Triangulation(
            tex_xs.ravel(),
            tex_ys.ravel(),
            np.asarray(tex_index_array)
        )
        
        
        mesh = MeshCut(buffer_array, index_array, new_transf, bounds)

        #
        # validation of edges being smaller -> larger
        for edge in mesh.tri_edge_index:

            v0,v1 = edge
            v0_cr = v0 % mesh.num_cols, v0 // mesh.num_cols
            v1_cr = v1 % mesh.num_cols, v1 // mesh.num_cols
        
            v0_m = mesh.get_point(*v0_cr)[:2]
            v1_m = mesh.get_point(*v1_cr)[:2]

            v0_pix_dist = v0_cr[0]**2 + v0_cr[1]**2
            v1_pix_dist = v1_cr[0]**2 + v1_cr[1]**2

            v0_dist = v0_m[0]**2 + v0_m[1]**2
            v1_dist = v1_m[0]**2 + v1_m[1]**2

            print(f"v0_dist: {v0_dist}, v1_dist: {v1_dist}")
            print(f"edge!: {v0}, {v1}");
            
            #assert v0_dist < v1_dist , "hmm..." ;



        # road_test_points = [[ (543270.0+1, 4943470.0-5)
        #                      , (543270.0, 4943470.0)
        #                      , (543270.0, 4943474.0)
        #                      , (543274.0, 4943479.0) ]]

        road_test_points =  [[(543304.0+1, 4943443.0), (543304.0+1, 4943483.0)]]
        #road_test_points = [[(543284.0, 4943438.0), (543284.0, 4943488.0)]]
        road_test_points = [[(543284.0, 4943488.0), (543284.0, 4943448.0)]]

        road_test_points = [ [mesh.interpolate_xy_coords((0.5,0.25)), mesh.interpolate_xy_coords((2,0.25)) ,\
                              mesh.interpolate_xy_coords((2,1.25)), mesh.interpolate_xy_coords((0.5,1.25)) ] ]

        road_test_points = [ [mesh.interpolate_xy_coords((0.5 , 0)), mesh.interpolate_xy_coords((2 , 0)) ,\
                              mesh.interpolate_xy_coords((2 , 1.25)), mesh.interpolate_xy_coords((0.5 , 1.25)) ] ]


        wit_in = [ mesh.within_bounds(p) for p in road_test_points[0] ]

        #road_test_points = [[(),()]]

        #pixel_y = ( _y_m - self.transform.f ) / self.transform.e;
        #pixel_x = ( _x_m - ( self.transform.c ) ) / self.transform.a;

        road_test_points = [ [(p[0], p[1], mesh.interp_elev(p)) for p in road_test_points[0]] ]

        #road_test_index = [[0, 3]]
        road_test_index = [[0, 1]]

        ccw_quad_points, ccw_quad_index = mesh.create_quad_from_points( road_test_points[0], road_test_index[0] )


        '''
        x,y = zip( *clip_poly_dict[0]["points"]+[clip_poly_dict[0]["points"][0]] )#[:2] )

        plt.plot(x, y, zorder=13, color="purple")

        x,y = zip( *clip_poly_dict[1]["points"]+[clip_poly_dict[1]["points"][0]] )#[:2] )

        plt.plot(x, y, zorder=12, color="yellow")

        '''

        ccw_quad_points = road_test_points[0]

        x , y , z = zip( *ccw_quad_points )
        
        plt.plot(x, y, 'o', zorder=20)

        x,y,z = zip( *ccw_quad_points + [ccw_quad_points[0]] )

        plt.plot(x, y, zorder=19, color="green")

        triang = mtri.Triangulation(
            xs.ravel(),
            ys.ravel(),
            np.asarray(index_array)
        )

        plt.triplot(triang, 'go-', label='DEM Mesh', color='black', linewidth=1)

        ed = defaultdict(int)

        print(f"# unique edges: {len(mesh.tri_edge_index)}")

        for edge_i, e in enumerate(mesh.tri_edge_index):

            e1,e2 = mesh.get_edge_points(edge_i)
            e_avg = (e1+e2)/2

            typ = mesh.get_type_edge(edge_i)

            ed[typ]+=1
            plt.text(e_avg.x, e_avg.y, f"{edge_i}", fontsize=12, color="red",  zorder=20)
            print(f"plotted edge: {edge_i}")
        
        for tri_index, tri in enumerate(index_array):
            tot_x = 0.0
            tot_y = 0.0

            for vert_index in tri:
                _x,_y = vert_index % dem_patch_num_verts,vert_index // dem_patch_num_verts
                tot_x += buffer_array[_x, _y,0]
                tot_y += buffer_array[_x, _y,1]

            plt.text(tot_x/3, tot_y/3, f"{tri_index}", fontsize=12, color="yellow",  zorder=20)
        

        plt.show()

        # p1 = [(543270.0, 4943474.0),
        #         (543273.0, 4943476.0),
        #         (543265.0, 4943475.0),
        #         (543266.0, 4943474.0)]
        # p2 = [ (543271.0, 4943469.0)
        #       , (543273.0, 4943476.0)
        #       , (543270.0, 4943474.0)
        #       , (543270.0, 4943470.0)
        #       ]
        # p2 = [ (543270.0, 4943470.0), (543271.0, 4943469.0) , (543273.0, 4943476.0), (543270.0, 4943474.0) ]

        # p2 = [ (543270.0, 4943470.0), (543271.0, 4943469.0)]

        #x,y = zip(*clip_polys[0])

        #plt.plot(x, y, zorder=10, color="pink")


        tri_edge_processed_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: False)))

        # 
        # now we iterate through! again to unzip
        # 

        outside_clip_polys = []

        tri_dict = defaultdict(list)
        test_dict = defaultdict(int)

        mesh.perform_clipping(  [ ccw_quad_points ] )


        #
        # can we shoot GRATICULES?
        vert_index = 4

        in_clip_4_right = mesh.follow_geodesic_clip( True 
                             , Graticule.LAT
                             , vert_index
                             , [ ccw_quad_points ] )

        in_clip_4_left = mesh.follow_geodesic_clip( False 
                        , Graticule.LAT
                        , vert_index
                        , [ ccw_quad_points ] )

        in_clip_4_up = mesh.follow_geodesic_clip( True 
                             , Graticule.LON
                             , vert_index
                             , [ ccw_quad_points ] )

        in_clip_4_down = mesh.follow_geodesic_clip( False 
                        , Graticule.LON
                        , vert_index
                        , [ ccw_quad_points ] )

        in_clip_4_diag_down = mesh.follow_geodesic_clip( False 
                                , Graticule.DIAG
                                , vert_index
                                , [ ccw_quad_points ] )

        in_clip_4_diag_up = mesh.follow_geodesic_clip( True 
                                , Graticule.DIAG
                                , vert_index
                                , [ ccw_quad_points ] )
                        

        '''
        for tri_edge_i, intersections in all_tri_edges.items():

            candidate_tris = mesh.get_tris_edge_index( tri_edge_i )
            print(f"precipitated by edge: { tri_edge_i } ")
            print(f"candidate triangles are: {candidate_tris[0]}, {candidate_tris[1]}")

            #
            # every unique edges has 1 or 2 triangles
            # the CCW order of which is different
            #

            for c_tri_indx in candidate_tris:

                #if c_tri_indx == 20:
                #    continue

                test_dict[c_tri_indx]+=1

                per_tri_clip = [] 

                # indices of edges in CCW order
                tri_edge_list = mesh.get_edge_indices(c_tri_indx)

                for i in range(3):
                    vertex_rel_edge = 0 if tri_edge_list[i]["ccw_order"] else -1

                    uniq_edge_index = tri_edge_list[i]["index"]

                    edge_tup = mesh.tri_edge_index[uniq_edge_index]

                    tri_edge_inters = all_tri_edges.get( uniq_edge_index, None)

                    if tri_edge_inters is None:
                        continue # to next edge

                    print(tri_edge_inters)

                    for inter_index in range( len( tri_edge_inters["intersections"] ) ):

                        print(f"# tri edge intersections: {len(tri_edge_inters["intersections"])}")

                        orig_inter = tri_edge_inters["intersections"][inter_index]

                        print(f"tri: {c_tri_indx}, edge: {uniq_edge_index}, intersection #: {inter_index}")

                        poss_poly = mesh.get_poly(c_tri_indx, tri_edge_list, orig_inter, tri_edge_processed_dict, all_tri_edges,
                              all_clip_edges, [ ccw_quad_points ])

                        print(f"possible poly: {poss_poly}")

                        if poss_poly is None:
                            pass
                        else:
                            per_tri_clip.append(poss_poly)
                            outside_clip_polys.append(poss_poly)


                tri_dict[c_tri_indx].extend(per_tri_clip)
        '''

        #
        #
        # OBJ creation
        #
        scale = 0.0008
        # mesh.create_face_obj(scale, "clipped_road_8-12-26", tri_dict)

        #
        # what are the clipping points in each
        #tri_inner_dict = defaultdict(lambda:[()])

        # set
        plot_set = set()

        all_polygons = []
 
        #te = tri_edges_dict[18]["edges"]
        #edge_0_clip = te[0]["clip"]
        #edge_1_clip = te[1]["clip"] 

        colors = ['pink', 'orange', 'green', 'blue']

        
        for _i,poly in enumerate(outside_clip_polys):
            
            tup_poly = [tuple(point["point"][:2]) for point in poly]   
            poly = Polygon(tup_poly, facecolor=colors[_i%4], edgecolor='black', linewidth=2, zorder=8)  
            ax.add_patch(poly)  
            # break     
        plt.show()
        



        

