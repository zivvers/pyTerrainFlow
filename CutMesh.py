    
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

import math

from MeshComponentCreator import MeshComponentCreator
from Graticule import *



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
    #
    # we should be using modified transform here so 0,0 corresponds
    # to the inset view
    #
    def __init__(self, buffer, index, tf, bnds):

        self.buffer_array = buffer
        self.index_array = index # triangle index for raster
        self.transform = tf
        self.bounds = bnds   

        self.num_cols = self.buffer_array.shape[0]
        self.num_rows = self.buffer_array.shape[1]

        self.num_triangles = (self.num_cols - 1)*(self.num_rows - 1)*2       

        # diagonal identifiers!
        self.min_diag = -1*self.num_rows;
        self.max_diag = self.num_cols;
    
        self.reverse_tri_index = []
        self.create_reverse_index()

    def cross_2D(self, a, b ):
        return (a.x * b.y) - (a.y * b.x)

    def convert_to_point(self, _pix_x, _pix_y):

        point_y = _pix_y*self.transform.e + self.transform.f ;
        pixel_x = _pix_x*self.transform.a + self.transform.c ;
        
        return glm.vec2(pixel_x, point_y)
    
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
    def order_ccw(self, polygon_verts):
        center_x = sum(x for x, y in polygon_verts) / len(polygon_verts)
        center_y = sum(y for x, y in polygon_verts) / len(polygon_verts)

        return sorted( polygon_verts,
            key=lambda p: math.atan2(p[1] - center_y, p[0] - center_x)
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

        ## BIV
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

        # for x in range(max(0, min_inner_x - 1), max_inner_x + 1):

        #     for x in range(max(0, min_inner_x - 1), max_inner_x + 1):


        #         print(f"possible diagonal: y {y} x {x}")
        #         print(f"pixel1 : {tuple(pix1)}, pix2: {tuple(pix2)}")
        #         print(f"diagonal identifier {x-y}\n")
        #         pix1 = glm.vec2(x, y)
        #         pix2 = glm.vec2(x + 1, y + 1)

        #         e1 = self.convert_to_point(*pix1);
        #         e2 = self.convert_to_point(*pix2);
        
        #         #mesh_edges.append([tuple(e1), tuple(e2)])

        #         diag_identifiers.add(x-y)

        # print(diag_identifiers)
        return mesh_edges

        #point_y = _pix_y*self.transform.e + self.transform.f ;
        #pixel_x = _pix_x * self.transform.a + self.transform.c ;

    #
    #
    # returns [ Graticule(), ... ]
    # of latitude, longitude, diagonals sliced within polygon
    # polygon must be CCW!
    def slice_edges( self, inner_edges, polygon_verts ):

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
                e_1 = polygon_verts[edge_i]
                e_2 = polygon_verts[(edge_i+1) % n]

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

                mesh_edge.p0 =p0
                mesh_edge.p1 =p1

                mesh_edge.p0_intersect_id = min_clip_i;
                mesh_edge.p1_intersect_id = max_clip_i;

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
    # conversion for the original index to new edge way
    #
    def get_edges(self, tri_index):

        vert_indices = self.index_array[tri_index]


        if tri_index % 2 == 0:
            e0_index = vert_indices[0]
            e1_index = vert_indices[1]
            e2_index = vert_indices[2]
        else:
            e0_index = vert_indices[0]
            e1_index = vert_indices[2]
            e2_index = vert_indices[1]
        
        e0_cr = e0_index % self.num_cols, e0_index // self.num_cols
        e1_cr = e1_index % self.num_cols, e1_index // self.num_cols
        e2_cr = e2_index % self.num_cols, e2_index // self.num_cols

        e0 = self.get_point(*e0_cr)
        e1 = self.get_point(*e1_cr)
        e2 = self.get_point(*e2_cr)

        if tri_index % 2 == 0:
            return [(e0, e1), (e2, e1), (e0, e2) ]
        else:
            return [(e0, e1), (e0, e2), (e2, e1) ]

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
    # void method that modifies tri_edges_dict
    # 
    def clip_tris( self, clippd_edges, polygon_verts, tri_edges_dict ):
#[(543264.0, 4943478.0), (543274.0, 4943468.0)]
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
                #  4--5--x 
                #  |\4|\ |  
                #  |3\|5\|
                #  3--v--0
                #  |\2|\0|
                #  | \|1\|
                #  x--2--1
                case Lat():
                    assert( p0_pix_y == p1_pix_y)
                    pixel_points = [(col, p0_pix_y) for col in range(min(p0_pix_x, p1_pix_x), max(p0_pix_x, p1_pix_x)+ 1)]
                
                    edge_index = 2
                case Lon():

                    assert( p0_pix_x == p1_pix_x)
                    pixel_points = [(p0_pix_x, row) for row in range(min(p0_pix_y, p1_pix_y), max(p0_pix_y, p1_pix_y)+ 1)]
                
                    edge_index = 1

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

                    print(f"For DIAG {diagonal_id} p0x: {p0_pix_x}, p1x: {p1_pix_x}, pixels: {pixel_points}")

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
                        s = "LON"
                        tri_index_0 = potential_tri_indices[1]
                        tri_index_1 = potential_tri_indices[2]
                        next_pixel = ( pixel[0], pixel[1]+1 )

                    case Lat():
                        s = "LAT"
                        tri_index_0 = potential_tri_indices[5]
                        tri_index_1 = potential_tri_indices[0]
                        next_pixel = ( pixel[0]+1, pixel[1] )

                    case Diag():
                        s = "DIAG"
                        tri_index_0 = potential_tri_indices[0]
                        tri_index_1 = potential_tri_indices[1]
                        next_pixel = ( pixel[0]+1, pixel[1]+1 )

                #
                # currently tri_index_0 and tri_index_1 might be -1 which means the respective
                # triangle doesn't exist (e.g. tri 5 on latitude pixel (0, _))
                #

                print(f"TRI INDEXS: {tri_index_0}, {tri_index_1} FOR {s} pix: {pixel}")
                #tri_index_set.add(tri_index_0)
                #tri_index_set.add(tri_index_1)

                next_point = self.get_point(*next_pixel)

                next_point_2D = next_point[:2]
                curr_point_2D = curr_point[:2]

                #is_curr_vert = abs(clippd_edge.p0-curr_point) < 10**-5
                #is_next_vert = abs(clippd_edge.p1-next_point) < 10**-5

                tri_edges_dict[tri_index_0]["edges"][edge_index]["orig"] = [curr_point_2D, next_point_2D]
                tri_edges_dict[tri_index_1]["edges"][edge_index]["orig"] = [curr_point_2D, next_point_2D]


                _p0 = tuple(clippd_edge.p0)
                _p1 = tuple(clippd_edge.p1)
                #
                # clip part is middle
                # *-------xxxxx-----* 
                if is_first and is_last:
                    
                    tri_edges_dict[tri_index_0]["edges"][edge_index]["clip"].append( ( _p0, _p1 ) );
                    tri_edges_dict[tri_index_1]["edges"][edge_index]["clip"].append( ( _p0, _p1) );
                    
                    clip_verts = [ self.same_point(_p0, curr_point_2D) , math.dist( _p1, next_point_2D) < 1e-6]
                    clip_edge_i = [ clippd_edge.p0_intersect_id,  clippd_edge.p1_intersect_id  ]
                #
                # clip part is end
                # *-------xxxxx* 
                elif is_first:
                    tri_edges_dict[tri_index_0]["edges"][edge_index]["clip"].append( ( _p0, tuple(next_point_2D)) )
                    tri_edges_dict[tri_index_1]["edges"][edge_index]["clip"].append( ( _p0, tuple(next_point_2D)) )

                    clip_verts = [ self.same_point( _p0, curr_point_2D ), True ]
                    clip_edge_i = [ clippd_edge.p0_intersect_id,  None  ]

                #
                # clip part is start
                # *-------xxxxx* 
                elif is_last:
                    tri_edges_dict[tri_index_0]["edges"][edge_index]["clip"].append( (tuple(curr_point_2D), tuple(clippd_edge.p1)) )
                    tri_edges_dict[tri_index_1]["edges"][edge_index]["clip"].append( (tuple(curr_point_2D), tuple(clippd_edge.p1)) )  

                    clip_verts = [ True, self.same_point( _p1, next_point_2D ) ]
                    clip_edge_i = [ None,  clippd_edge.p1_intersect_id  ]
                #
                # clip part is WHOLE THING
                # *xxxxx* 
                else: 
                    tri_edges_dict[tri_index_0]["edges"][edge_index]["clip"].append( (tuple(curr_point_2D), tuple(next_point_2D) ) )
                    tri_edges_dict[tri_index_1]["edges"][edge_index]["clip"].append( (tuple(curr_point_2D), tuple(next_point_2D) ) )  

                    clip_verts = [True, True]
                    clip_edge_i = [ None, None ]

                #tri_edges_dict[tri_index_0]["edges"][edge_index]["bool"] = clip_verts

                print(f"for TRI: {tri_index_0} , pixel: {pixel} isFirst?: {is_first} isLast?: {is_last} ;; we are setting {edge_index} clippd: {clip_verts}")

                tri_edges_dict[tri_index_0]["edges"][edge_index]["bool"] = clip_verts 
                tri_edges_dict[tri_index_1]["edges"][edge_index]["bool"] = clip_verts

                tri_edges_dict[tri_index_0]["edges"][edge_index]["clip_index"] = clip_edge_i
                tri_edges_dict[tri_index_1]["edges"][edge_index]["clip_index"] = clip_edge_i

        #
        # have to place clip vertices
        for clip_p in polygon_verts:
            
            tri_i = self.get_tri_index_point(*clip_p)
            print(f"adding cingle clip point: {clip_p} to {tri_i}")
            tri_edges_dict[tri_i]["points"].append(clip_p)



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


    def get_point(self, pix_X, pix_Y):
        c = int(pix_X)
        r = int(pix_Y)
        lon = self.buffer_array[c, r, 0].item()
        lat = self.buffer_array[c, r, 1].item()
        elev = self.buffer_array[c, r, 2].item()

        return (lon, lat, elev)


    def make_edges_dict(self):
        return {
            "points": [],
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


if __name__ == "__main__":

    cfg = load_config()

    TILE_SIZE = 1024

    raster_file = "dem_EPSG_26910_542354_4944208.tif"
    texture_file = "osip_EPSG_26910_542354_4944208.tif"

    plot_points = [[542781.25,4943800.5], [542774,4943798]]


    raster_path = cfg.input_dir / raster_file
    texture_path = cfg.input_dir / texture_file

    offset_x, offset_y = 90, 71

    tex_offset_x , tex_offset_y = offset_x*10, offset_y*10

    dem_resolution = 10 #meters per pixel
    tex_resolution = 1 #meters per pixel
    
    #
    # edit this variable!
    dem_num_cells = 4

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

        num_rows_total = band.shape[0]
        num_cols_total = band.shape[1]

        band = band[offset_y:offset_y+dem_patch_num_verts,offset_x:offset_x+dem_patch_num_verts]
        h = src.height
        w = src.width

        rows = np.arange( dem_patch_num_verts ) + offset_y
        cols = np.arange( dem_patch_num_verts ) + offset_x

        cols, rows = np.meshgrid(cols, rows)
        xs, ys = rasterio.transform.xy(src.transform, rows, cols,offset='ul')

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


        #extent = [tex_src.bounds.left, tex_src.bounds.right, tex_src.bounds.bottom, tex_src.bounds.top]
        #x_min, x_max, y_min, y_max = extent  
        #x_max = xs[-1];
        #y_min = ys[-1];


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
        im = ax.imshow(tex_band, extent=[dem_extent_d["left"], dem_extent_d["right"], dem_extent_d["bottom"], dem_extent_d["top"]], cmap='viridis', origin='upper')

        print(f"index array size: {index_array}")

        tex_triang = mtri.Triangulation(
            tex_xs.ravel(),
            tex_ys.ravel(),
            np.asarray(tex_index_array)
        )

        #plt.triplot(tex_triang, 'g-', label='Texture Mesh',color='pink', linewidth=0.25)

        bounds = [ buffer_array[0,0,0].item(), buffer_array[0, num_rows_inset-1, 1].item(), 
                    buffer_array[num_cols_inset-1,0 , 0].item(), buffer_array[0, 0, 1].item() ]

        # test_points = [ [542354.0 - 100, 4944198.0],
        #                 [542356.5 - 100, 4944200.9],
        #                 [542359.0 - 100, 4944203.6],
        #                 [542361.5 - 100, 4944205.9],
        #                 [542364.0 - 100, 4944208.0] ] 
        
        # for i,p in enumerate(test_points):
        #     x,y=p
        #     if (bounds[0] <= x <= bounds[2]) and (bounds[1] <= y <= bounds[3]):
        #         print(f"point {i} within bounds!")
        
        # test_index = np.array( [[0, 4]] , dtype=int )

        #promesheus = MeshComponentCreator(buffer_array, transf, bounds)

        # road_lines_comp_3D = promesheus.create_lines( test_points, test_index )

        # test_points = [[543357, 4943473.0],
        #                [543367, 4943463.0],
        #                [543380.0, 4943472.0],
        #                  ]
        
        test_points = [ glm.vec2(543369 - 100, 4943477.0), glm.vec2(543378.0 - 100, 4943472.0)
                       , glm.vec2(543367 - 100, 4943463.0), glm.vec2(543357 - 100, 4943473.0) ]

        mesh = MeshCut(buffer_array, index_array, new_transf, bounds)

        test_points_ccw = mesh.order_ccw(test_points)

        # inner_points = mesh.get_inner_verts(test_points);

        inner_edges = mesh.get_inner_edges(test_points_ccw);

        cut_edges = mesh.slice_edges( inner_edges, test_points_ccw );

        orig_mesh_lines = LineCollection([ (l.p0, l.p1) for l in inner_edges ], colors='crimson', linewidths=1.5)
        
        cut_mesh_lines = LineCollection([ (l.p0, l.p1) for l in cut_edges] , colors='blue', linewidths=2)

        ax.add_collection( orig_mesh_lines )

        ax.add_collection( cut_mesh_lines )
  

        tri_edges_dict = defaultdict(mesh.make_edges_dict)

        #
        # in place modification of tri_edges_dict
        #
        mesh.clip_tris( cut_edges, test_points_ccw, tri_edges_dict )

        #index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])

        #index_array = [val for idx, val in enumerate(index_array) if idx not in tri_dict.keys()]

        for tri_index, tri in enumerate(index_array):
            tot_x = 0.0
            tot_y = 0.0

            for vert_index in tri:
                _x,_y = vert_index % dem_patch_num_verts,vert_index // dem_patch_num_verts
                tot_x += buffer_array[_x, _y,0]
                tot_y += buffer_array[_x, _y,1]

            #print(f"plot at ({tot_x/3}, {tot_y/3})")
            plt.text(tot_x/3, tot_y/3, f"{tri_index}", fontsize=12, color="yellow",  zorder=20)


        triang = mtri.Triangulation(
            xs.ravel(),
            ys.ravel(),
            np.asarray(index_array)
        )

        plt.triplot(triang, 'go-', label='DEM Mesh', color='black', linewidth=1)

        #
        # what are the clipping points in each
        #tri_inner_dict = defaultdict(lambda:[()])

        # set
        plot_set = set()


        for tri_index, inner_dict in tri_edges_dict.items():

            # NNV
            print(f"tri {tri_index}")
            print(f"orig edge: {inner_dict["edges"][0]["orig"]}")
            print(f"clip bool e0: {inner_dict["edges"][0]["bool"]}")
            print(f"clip bool e1: {inner_dict["edges"][1]["bool"]}")
            print(f"clip bool e2: {inner_dict["edges"][2]["bool"]}")
            print(f"clip index e0: {inner_dict["edges"][0]["clip_index"]}")
            print(f"clip index e1: {inner_dict["edges"][1]["clip_index"]}")
            print(f"clip index e2: {inner_dict["edges"][2]["clip_index"]}")
            print(f"clip bool: {inner_dict["edges"][0]["clip"]}")

            edges = mesh.get_edges(tri_index)

            for i in range(3):
                if inner_dict["edges"][i]["orig"] is None:
                    pass
                else:
                    d_e = inner_dict["edges"][i]["orig"]
                    _e = edges[i]
                    assert(mesh.same_point(d_e[0], _e[0]))
                    assert(mesh.same_point(d_e[1], _e[1]))


            keep_vert_0 = True
            keep_vert_1 = True
            keep_vert_2 = True

            if tri_index % 2 == 0:


                if inner_dict["edges"][0]["orig"] is None:
                    if inner_dict["edges"][2]["orig"] is not None:
                        assert(inner_dict["edges"][2]["bool"][0] == False)

                else:
                    keep_vert_0 = not inner_dict["edges"][0]["bool"][0]

                if inner_dict["edges"][0]["orig"] is None:
                    if inner_dict["edges"][1]["orig"] is not None:
                        assert(inner_dict["edges"][1]["bool"][1] == False)

                else:
                    keep_vert_1 = not inner_dict["edges"][0]["bool"][1]


                if inner_dict["edges"][1]["orig"] is None:
                    if inner_dict["edges"][2]["orig"] is not None:
                        assert(inner_dict["edges"][2]["bool"][1] == False)
                else:
                    keep_vert_2 = not inner_dict["edges"][1]["bool"][0]              
                
                # if inner_dict["edges"][0]["orig"] is None and inner_dict["edges"][2]["orig"] is not None:
                #     vert_0_bool0 = inner_dict["edges"][0]["bool"][0]
                #     vert_0_bool1 = inner_dict["edges"][2]["bool"][0]
                #     assert(vert_0_bool0 == vert_0_bool1)
                
                # else:
                # if inner_dict["edges"][0]["orig"] is not None and inner_dict["edges"][1]["orig"] is not None:
                #     vert_1_bool0 = inner_dict["edges"][0]["bool"][1]
                #     vert_1_bool1 = inner_dict["edges"][1]["bool"][1]
                #     assert(vert_1_bool0 == vert_1_bool1)

                if inner_dict["edges"][1]["orig"] is not None and inner_dict["edges"][2]["orig"] is not None:
                    vert_2_bool0 = inner_dict["edges"][1]["bool"][0]
                    vert_2_bool1 = inner_dict["edges"][2]["bool"][1]  
                    assert(vert_2_bool0 == vert_2_bool1)
               
            else:
                # if inner_dict["edges"][0]["orig"] is not None and inner_dict["edges"][1]["orig"] is not None:
                #     vert_0_bool0 = inner_dict["edges"][0]["bool"][0]
                #     vert_0_bool1 = inner_dict["edges"][1]["bool"][0]
                #     assert(vert_0_bool0 == vert_0_bool1)
                if inner_dict["edges"][0]["orig"] is None:
                    if inner_dict["edges"][1]["orig"] is not None:
                        assert(inner_dict["edges"][1]["bool"][0] == False)
                else:
                    keep_vert_0 = not inner_dict["edges"][0]["bool"][0]

                if inner_dict["edges"][0]["orig"] is not None and inner_dict["edges"][2]["orig"] is not None:
                    vert_1_bool0 = inner_dict["edges"][0]["bool"][1]
                    vert_1_bool1 = inner_dict["edges"][2]["bool"][1]
                    assert(vert_1_bool0 == vert_1_bool1)

                if inner_dict["edges"][0]["orig"] is None:
                    if inner_dict["edges"][2]["orig"] is not None:
                        assert(inner_dict["edges"][2]["bool"][2] == False)

                else:
                    keep_vert_1 = not inner_dict["edges"][0]["bool"][1]

                if inner_dict["edges"][1]["orig"] is not None and inner_dict["edges"][2]["orig"] is not None:
                    vert_2_bool0 = inner_dict["edges"][2]["bool"][0]
                    vert_2_bool1 = inner_dict["edges"][1]["bool"][1]  
                    assert(vert_2_bool0 == vert_2_bool1)

                if inner_dict["edges"][1]["orig"] is None:
                    if inner_dict["edges"][2]["orig"] is not None:
                        assert(inner_dict["edges"][2]["bool"][0] == False)
                else:
                    keep_vert_2 = not inner_dict["edges"][1]["bool"][1] 
            
            
            #assert(vert_2_bool0 == vert_2_bool1)
            print(f"TRIANGLE # {tri_index}, keep v0: {keep_vert_0}, v1: {keep_vert_1}, v2: {keep_vert_2}")

        #
        #    vert_indices = mesh.get_triangle_vertex_indices( tri_indx )
        #
        #    verts = [mesh.get_point(x,y)[:2] for x,y in vert_indices ]
        #    print(f"vert indices: {vert_indices}")
        #    print(f"adding triangle: {verts}")
        #    poly = Polygon(verts, facecolor='cyan', edgecolor='black', linewidth=2, zorder=8)
        #
        #    ax.add_patch(poly)


            #math.dist(tri_edges[1][0][:2], tri_edges[0][0][:2]) < .01
            
        # for tri_index, edge_dict in tri_dict.items():
        #     if (tri_index != -1):
        #         pass
        #         #print(f"removing : {tri_index}")
        #         #del 

        #     tri_verts = index_array[tri_index]

        #     pixs = []
        #     coords = []
        #     for vert_index in tri_verts:
        #         _pix_x,_pix_y = vert_index % dem_patch_num_verts, vert_index // dem_patch_num_verts
        #         _x,_y = buffer_array[_pix_x,_pix_y,0].item(),buffer_array[_pix_x,_pix_y,1].item()

        #         pixs.append((_pix_x,_pix_y))
        #         coords.append((_x,_y,0))

        #     tri1 = Triangle(coords, pixs)

        #     # what clip points are in tris
        #     #tri_inner_dict[ tri_index ] = 

        #     #tri1.ccw_verts_to_downward_edge(tri_index, edge_indx)

        #     #
        #     #
        #     # {<tri index>: {0: {"orig" : [(),()], "clip" : [(), ()], 1: ....}}
        #     #
        #     #

        #     tri_edges_dict = {
        #         "clip_points": [],
        #         "edges" : {
        #             edge_index: {
        #                 "orig": None,
        #                 "clip": None,
        #             }
        #             for edge_index in range(3)
        #         }
        #     }

        #     # check all triangles for clip_point
        #     for clip_p in test_points_ccw:
        #         stu = tri1.bary_x_y(*clip_p)

        #         if all(e >= 0  for e in stu):
        #             print(f"Triangle: {tri_index} contains clip poly vert: {clip_p}")
        #             # tri_inner_dict[ tri_index ].append(tuple(clip_p))

        #             plot_set.add(tuple(clip_p))

        #             tri_edges_dict["clip_points"].append(clip_p)

        #     tri_edges = []
        #     clipped_edges = []

        #     for edge_i in range(3):

        #         # translate CCW indexed points into our expected
        #         # "downwards" edge shape
        #         orig_edge = tri1.ccw_verts_to_downward_edge(tri_index, edge_i)
                
        #         tri_edges_dict["edges"][edge_i]["orig"] = orig_edge 

        #         # we have clipped_edge
        #         if edge_i in edge_dict:
        #             clipped_edge = edge_dict[edge_i]

        #             tri_edges_dict["edges"][edge_i]["clip"] = clipped_edge

        #             start_diff = math.dist(clipped_edge[0][:2] , orig_edge[0][:2])
        #             end_diff = math.dist(clipped_edge[1][:2] , orig_edge[1][:2])

        #             assert start_diff > 0.001 or end_diff > 0.001
                    
        #             plot_set.add(clipped_edge[0][:2])
        #             plot_set.add(clipped_edge[1][:2])

        #             # add first vertex
        #             if start_diff:
        #                 plot_set.add( orig_edge[0][:2] )
        #             # add second vertex
        #             if end_diff:
        #                 plot_set.add( orig_edge[1][:2] )

        #         # no clipping edge!
        #         # add vertices
        #         else:
        #             plot_set.add( orig_edge[0][:2] )
        #             plot_set.add( orig_edge[1][:2] )
                
        #         tri_edges.append( orig_edge )
        #         print(f"clip edge: {clipped_edge}\n\t, origiedge: {orig_edge}, ")

        #     all_tri_edges_dict[tri_index] = tri_edges_dict


            #break;


            #assert math.dist(tri_edges[1][0][:2], tri_edges[0][0][:2]) < .01
            

            #
            # VERT 0 is at start of edge 0 and 1
            # if it's clipped by one it should be by both
            #if math.dist(tri_edges[0][0][:2], clipped_edges[0][0][:2]) > .0002:
            #    assert(math.dist(tri_edges[1][0][:2], clipped_edges[1][0][:2]) > .0002)


            #print(f"tri index: {tri_index}, in set? {}")
            #pass

            # now I have the clipped edges and the inner clip points
            # reconstruct the inner shapes
            



        #for tri_indx in tri_indices:
        #
        #    vert_indices = mesh.get_triangle_vertex_indices( tri_indx )
        #
        #    verts = [mesh.get_point(x,y)[:2] for x,y in vert_indices ]
        #    print(f"vert indices: {vert_indices}")
        #    print(f"adding triangle: {verts}")
        #    poly = Polygon(verts, facecolor='cyan', edgecolor='black', linewidth=2, zorder=8)
        #
        #    ax.add_patch(poly)

        #neu_x, neu_y = zip(*list(plot_set))

        
        te = tri_edges_dict[18]["edges"]
        edge_0_clip = te[0]["clip"]
        edge_1_clip = te[1]["clip"] 



        poly = Polygon(verts, facecolor='cyan', edgecolor='black', linewidth=2, zorder=8)
        ax.add_patch(poly)
        #
        # 
        # plt.scatter(neu_x, neu_y, c="pink", zorder= 30)

        x,y = zip(*test_points_ccw)

        plt.plot(x, y, zorder=10)

        # test_points = [ [543354, 4943488-2.5], [543364, 4943488-2.5]
        #                 , [543364, 4943488-7.5], [543354, 4943488-7.5] ]
        # x,y = zip(*test_points)

        # plt.plot(x, y, zorder=10)

        # #
        # # get example clipping tri
        # #
        # tri_pnts = []

        # for i in range(3):
        #     tri_pnts.append( glm.vec2(xs[index_array[1][i]].item(), ys[index_array[1][i]].item()) )


        # _dir = glm.vec2(test_points[2]) - glm.vec2(test_points[3])

        # _curr_pnt = glm.vec2(543364, 4943488-7.5);

        # for i in range(3):
        #     e_pnt1 =  tri_pnts[i]
        #     e_pnt2 =  tri_pnts[(i+1) % 3]

        #     seg = e_pnt2 - e_pnt1

        #     den = cross_2D(_dir, seg)
        #     print(f"edge {i} determin: {den}")

        #     if den == 0:
        #         continue ;

        #     t = cross_2D(e_pnt1 - _curr_pnt, seg) / den
        #     u = cross_2D(e_pnt1 - _curr_pnt, _dir) / den

        #     print(f"t: {t}, u: {u}")

        #seg = e_pnt2 - e_pnt1
        #den = cross_2D(_dir, seg)

        plt.show()

