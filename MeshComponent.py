import rasterio
import numpy as np
from osgeo import gdal
from pyglm import glm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import math
from collections import Counter
from osgeo import ogr

from matplotlib.animation import FuncAnimation
from matplotlib.patches import Polygon

#
# should have separate utilities library!
def cross_2D( a, b):
    return (a.x * b.y) - (a.y * b.x)


#
# MeshComponent is naive to whole index
# Just the Vertex, Edge, Triangle we will make in MeshComponentCreator
#
class MeshComponent:

    def __init__(self, pts=None, sts=None):
        self.valid = False
        
        self.normal = None
        if pts is None:
            self.point_s = None
            self.pixel_s = None

        else:
            self.valid = True
            #self.mesh = mesh #ThreeDMesh outer class #NOMas
            self.pixel_s = [(int(tup[0]), int(tup[1]) ) for tup in sts]
            self.point_s = [(float(tup[0]), float(tup[1]), float(tup[2])) for tup in pts]
            

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        else:
            return True

    # method to over-ride
    def get_edges(self):

        return list()

    # void method to overide
    def calc_norm_axes(self):
        return NotImplemented

        def print_obj(self):
            print("MeshComponent")
 
class Vertex( MeshComponent ):

    def __init__(self, pts, sts):
            
        print(f"MAKING VERTEX, num STS: {len(sts)}")
        super().__init__(pts, sts)

    def __eq__(self, other):

        # if the other one is a vertex
        if isinstance(other, Vertex):
            # vertex single pixel in edge
            #return all(elem in self.pixel_s for elem in other.pixel_s)
            return Counter(self.pixel_s) == Counter(other.pixel_s)
        # if the other one is an Edge
        if isinstance(other, Edge):
            return False
            

class Edge(MeshComponent):
    def __init__(self, pts, sts):
            
        print(f"MAKING EDGE, num STS: {len(sts)}")
        super().__init__(pts, sts)

    def __eq__(self, other):

        # if the other one is a vertex
        if isinstance(other, Vertex):
            # vertex single pixel in edge
            #return all(elem in self.pixel_s for elem in other.pixel_s)
            return False
        # if the other one is an Edge
        if isinstance(other, Edge):
            return Counter(self.pixel_s) == Counter(other.pixel_s)
                

    def calc_norm_axes(self):
        print("calc EDGE normal")





class Triangle(MeshComponent):

    def __init__(self, pts=None, sts=None):
        super().__init__(pts, sts)
        if sts ==None:
            print(f"MAKING TRIANGLE Invalid: {self.valid}")
        else:  
            print(f"MAKING TRIANGLE, num STS: {len(sts)}")
        

        # these will be GLM vec3
        self.local_x_axis = None
        self.local_y_axis = None

    def print_obj(self):
        print(f"TRI {self.pixel_s[0]} ; {self.pixel_s[1]} ; {self.pixel_s[2]}")


    def bary_x_y(self, _x, _y):

        v1,v2,v3 = self.point_s
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

    def calc_norm_axes(self):

        pnt1, pnt2, pnt3 = self.point_s

        seg1 = glm.vec3(pnt2) - glm.vec3(pnt1);
        seg2 = glm.vec3(pnt3) - glm.vec3(pnt1);
        tri_normal = glm.normalize(glm.cross(seg1, seg2));

        # This is sorta cheating, assume up
        vec_from_center = glm.vec3(0,0,1)

        # flip it!
        if ( glm.dot( tri_normal , vec_from_center ) < 0 ):
            tri_normal = -tri_normal

        self.normal = tri_normal

    def get_edges( self ):
        v_1 = self.point_s[0]
        v_2 = self.point_s[1]
        v_3 = self.point_s[2]

        v_1_st = self.pixel_s[0]
        v_2_st = self.pixel_s[1]
        v_3_st = self.pixel_s[2]

        # return edges!
        return [ Edge( [v_1, v_2], [v_1_st, v_2_st] ), Edge( [v_2, v_3], [v_2_st, v_3_st] ) \
                , Edge( [v_3, v_1], [v_3_st, v_1_st] ) ]
    #
    # for vertex on triangle check if the other two
    # vertices are intersected by our dir vec
    #
    def check_vert_intersection(self, vert, _dir):
        v_1 = self.point_s[0]
        v_2 = self.point_s[1]
        v_3 = self.point_s[2]

        v_1_st = self.pixel_s[0]
        v_2_st = self.pixel_s[1]
        v_3_st = self.pixel_s[2]

        if vert.pixel_s[0] == v_1_st:
            other_verts = [ Vertex( [v_2], [v_2_st] ), Vertex( [v_3], [v_3_st] ) ]
        elif vert.pixel_s[0] == v_2_st:
            other_verts = [ Vertex( [v_1], [v_1_st] ), Vertex( [v_3], [v_3_st] ) ]
        elif vert.pixel_s[0] == v_3_st:
            other_verts = [ Vertex( [v_1], [v_1_st] ), Vertex( [v_2], [v_2_st] ) ]
        else:
            other_verts = []

        if other_verts:
            v_curr = glm.vec3(vert.point_s[0]).xy
            v_other_1 = glm.vec3(other_verts[0].point_s[0]).xy
            v_other_2 = glm.vec3(other_verts[1].point_s[0]).xy

            e1 = v_curr - v_other_1
            e2 = v_curr - v_other_2

            colinear1 = abs(cross_2D(_dir, e1))
            colinear2 = abs(cross_2D(_dir, e2))

            forward1 = glm.dot(_dir, e1)
            forward2 = glm.dot(_dir, e2)

            print(f"EDGE COLINEAR 1: {colinear1:.5f}, CO0LINEAR 2: {colinear2:.5f}")


        else:
            #return []
            print(f"wrong vertex for triangle")



class MeshComponentCreator:

    ENTRIES_PER_BUFFER = 21 #!
    ENTRIES_PER_INDEX = 3

    def __init__(self, buffer, tf, bnds):

        self.buffer_array = buffer
        self.transform = tf
        self.bounds = bnds
        self.num_cols = self.buffer_array.shape[0]
        self.num_rows = self.buffer_array.shape[1]

        self.create_index() # right now just for visualization

    # returns vertices, x_center, y_cente
    def get_points_centered_scaled(self, scale):

        points = []
        cnt = 0

        for r in range(self.num_rows):
            for c in range(self.num_cols):
     
                p = promesheus.buffer_array[c, r , 0:3]

                # UNIFORM SCALING SO EVERYTHING HAS SAME NORMALS!

                new_z = p[2] * scale # not centering Z

                points.append( ( ( p[0]-self.x_center ) * scale, ( p[1]-self.y_center)* scale, new_z ) )
        return points

    def create_index(self):

        self.test_index_array = []

        num_tri = (self.num_cols - 1) * (self.num_rows-1) * 2
        self.index_array = np.empty((num_tri, self.ENTRIES_PER_INDEX), dtype=np.int32)

        indx = 0
        for r in range( 1, self.num_rows):
            for c in range(1,  self.num_cols):
    
                # storing vertex points 
                currIndx = c + r * self.num_cols;
                prevRowSameColIndx = c + (r - 1) * self.num_cols;
                prevRowBackColIndx = c - 1 + (r - 1) * self.num_cols;
                prevIndx = currIndx - 1;

                self.index_array[indx, 0] = prevRowBackColIndx
                self.index_array[indx, 1] = currIndx
                self.index_array[indx, 2] = prevRowSameColIndx
                #self.index_array[indx, 3] = -1 # blank

                #self.add_reverse_index( c - 1, r - 1, indx) #prevRowBackColIndx
                #self.add_reverse_index( c, r, indx) #currIndx
                #self.add_reverse_index( c, r - 1 , indx) #prevRowSameColIndx


                self.test_index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])

                indx+=1
                self.index_array[indx, 0] = prevRowBackColIndx
                self.index_array[indx, 1] = prevIndx
                self.index_array[indx, 2] = currIndx
                #self.index_array[indx, 3] = -1 # blank

                #self.add_reverse_index( c - 1, r - 1, indx) #prevRowBackColIndx
                #self.add_reverse_index( c - 1, r , indx) #prevIndx
                #self.add_reverse_index( c , r , indx) # currIndx

                self.test_index_array.append([prevRowBackColIndx, prevIndx, currIndx])

                indx+=1


    def convert_to_pixel(self, _x_m, _y_m):

        pixel_y = ( _y_m - ( self.transform[5] + 0.5*self.transform[4] ) ) / self.transform[4];
        pixel_x = ( _x_m - ( self.transform[2] + 0.5*self.transform[0] ) ) / self.transform[0];
        
        return (math.floor(pixel_x), math.floor(pixel_y))

    #
    # Assumes that params are already ints!
    #
    def get_point(self, pix_X, pix_Y):
        c = int(pix_X)
        r = int(pix_Y)
        lon = self.buffer_array[c, r, 0].item()
        lat = self.buffer_array[c, r, 1].item()
        elev = self.buffer_array[c, r, 2].item()

        return (lon, lat, elev)

    def get_norm_pix(self, pix_X, pix_Y):
        c = pix_X
        r = pix_Y

        return glm.vec3( self.buffer_array[c, r, 17], self.buffer_array[c, r, 18], self.buffer_array[c, r, 19] );


    def get_norm_vert_edge(self, vert_or_edge):

        # pull norms from vertices
        norm_s = [self.get_norm_pix(*st) for st in vert_or_edge.pixel_s]
        total_norm = sum(norm_s, glm.vec3(0))
        return glm.normalize(total_norm)

    def calculate_vertex_normals(self):    
        x_max, y_max, z_max = -math.inf, -math.inf, -math.inf
        x_min, y_min, z_min = math.inf, math.inf, math.inf

        # calc center of all points
        for r in range(self.num_rows):
            for c in range(self.num_cols):
                print(f"x max: {x_max} x min {x_min}, y max {y_max} y min {y_min}")

                p = self.get_point(c, r )
                x_max = max(x_max, p[0])
                y_max = max(y_max, p[1])
                z_max = max(z_max, p[2])
                x_min = min(x_min, p[0])
                y_min = min(y_min, p[1])
                z_min = min(z_min, p[2])

        self.x_center = (x_max + x_min) / 2
        self.y_center = (y_max + y_min) / 2
        self.z_center = (z_max + z_min) / 2

        center = glm.vec3(self.x_center, self.y_center, self.z_center)

        for r in range(self.num_rows):
            for c in range(self.num_cols):

                p = glm.vec3(self.get_point(c, r ))

                #vec_from_center = p - center

                topRight = c == self.num_cols - 1 and r == 0
                topLeft = c == 0 and r == 0
                bottomRight = c == self.num_cols - 1 and r == self.num_rows - 1
                bottomLeft = c == 0 and r == self.num_rows - 1

                triangles = []

                # get all triangles given vertex position
                if topRight:
                    print("top right")
                    triangles = [self.make_tri_index(c,r,2)]

                elif topLeft:
                    print("top left")
                    triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1)]

                elif bottomRight:
                    print("bottom rigght")
                    triangles = [self.make_tri_index(c,r,4), self.make_tri_index(c,r,3)]

                elif bottomLeft:
                    print("bottom left")
                    triangles = [ self.make_tri_index(c,r,5) ]

                elif c == self.num_cols - 1: #right column
                    print("right col")
                    triangles = [self.make_tri_index(c,r,2), self.make_tri_index(c,r,3), self.make_tri_index(c,r,4)]

                elif (c == 0): #left column
                    triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1), self.make_tri_index(c,r,5)]
                    print("left col")
                elif (r == self.num_rows - 1): # bottom row
                    triangles = [self.make_tri_index(c,r,5), self.make_tri_index(c,r,4)]
                    print("bottom row")
                elif (r == 0): # top row
                    triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1), self.make_tri_index(c,r,2)]
                    print("top row")
                else: # interior
                    print("interior")
                    triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1), self.make_tri_index(c,r,2), \
                                  self.make_tri_index(c,r,3), self.make_tri_index(c,r,4), self.make_tri_index(c,r,5)  ]
                
                vertex_norm = glm.vec3(0,0,0)
                for t in triangles:

                    t.calc_norm_axes()
                    assert(t.normal is not None)
                    vertex_norm = vertex_norm + t.normal
                    
                vertex_norm = glm.normalize(vertex_norm)

                # vec from center stuff isn't working
                vec_from_center = glm.vec3(0,0,1)

                # flip it!
                if (glm.dot(vertex_norm , vec_from_center ) < 0):
                    vertex_norm = -vertex_norm

                self.buffer_array[c, r, 17] = vertex_norm.x
                self.buffer_array[c, r, 18] = vertex_norm.y
                self.buffer_array[c, r, 19] = vertex_norm.z



                




    def print_stats(self):
        print(f"mesh of cols: {self.num_cols}, rows: {self.num_rows}")


    def get_tri_vert_indx(self, indx, c, r):
        match indx:
            case 0:
                return (c+1, r)
            case 1:
                return (c+1, r+1)
            case 2:
                return (c, r+1)
            case 3:
                return (c-1, r)
            case 4:
                return (c-1, r-1)
            case 5:
                return (c, r-1)


    #  4--5--x 
    #  |\4|\ |  
    #  |3\|5\|
    #  3--v--0
    #  |\2|\0|
    #  | \|1\|
    #  x--2--1
    def make_tri_index( self, c, r, i) -> Triangle:
            
        match i:
            case 0:
                sts = [(c,r) , self.get_tri_vert_indx(0, c,r) , self.get_tri_vert_indx(1, c,r)]
                print(f"STs: {sts}")
                coords = [self.get_point(*st) for st in sts]
            case 1:
                sts = [(c,r) , self.get_tri_vert_indx(1, c,r) , self.get_tri_vert_indx(2, c,r)]
                coords = [self.get_point(*st) for st in sts]
            case 2:
                sts = [(c,r) , self.get_tri_vert_indx(2, c,r) , self.get_tri_vert_indx(3, c,r)]
                coords = [self.get_point(*st) for st in sts]
            case 3:
                sts = [(c,r) , self.get_tri_vert_indx(3, c,r) , self.get_tri_vert_indx(4, c,r)]
                coords = [self.get_point(*st) for st in sts]
            case 4:
                sts = [(c,r) , self.get_tri_vert_indx(4, c,r) , self.get_tri_vert_indx(5, c,r)]
                coords = [self.get_point(*st) for st in sts]
            case 5:
                sts = [(c,r) , self.get_tri_vert_indx(5, c,r) , self.get_tri_vert_indx(0, c,r)]
                coords = [self.get_point(*st) for st in sts]

        tri = Triangle(coords, sts)
        return tri


    def create_reverse_index(self):
        #num_rows = band.shape[0]
        #num_cols
        index_l = []
        for y in range(1,self.num_rows):
            for x in range(1,self.num_cols):

                currIndx = x + y * self.num_cols
                prevRowSameColIndx = x + (y - 1) * self.num_cols
                prevRowBackColIndx = x - 1 + (y - 1) * self.num_cols
                prevIndx = currIndx - 1

                index_l.append((prevRowBackColIndx, currIndx, prevRowSameColIndx))
                index_l.append((prevRowBackColIndx, prevIndx, currIndx))
        
        reverse_l = []
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

                reverse_l.append(tuple(currReverse))

        return index_l, reverse_l

    #
    # for a given point on a VERTEX or EDGE  and direction dir
    # find the 2 other verts of the triangle that is next
    #
    def get_next_tri( self, vertex_or_edge, _point, _dir ):
        
        # edge
        if isinstance(vertex_or_edge, Edge):
            print("getting next Tri EDGE")
            pnt1 = vertex_or_edge.point_s[0]
            pnt2 = vertex_or_edge.point_s[1]

            pix1 = vertex_or_edge.pixel_s[0]
            pix2 = vertex_or_edge.pixel_s[1]

            e1_x = pix1[0] 
            e1_y = pix1[1]
            e2_x = pix2[0] 
            e2_y = pix2[1]

            pnt1_2D = glm.vec3(pnt1).xy
            pnt2_2D = glm.vec3(pnt2).xy

            print(f"point1: {pnt1_2D} point2: {pnt2_2D}")

            seg = pnt2_2D - pnt1_2D;

            #if seg.y == 0:
            #    perp = glm.vec2(-1,0) * seg
            #else:
            #    perp = glm.vec2(0, -1) * seg

            perp = glm.rotate(glm.vec3(seg, 0.0), glm.radians(90.0), glm.vec3(0, 0, 1)).xy
        
            #
            # A---*
            #  \  E\ 
            #   \ E \
            #    \E  \
            #     *---B
            # 
            if (e1_x == e2_x):
                print("scenario 1")
                vert_b = (e1_x+1, max(e1_y,e2_y)) # right
                vert_a = (e1_x-1, min(e1_y,e2_y)) # left
            #
            #   A
            #   |\
            #   | \
            #   |  \      ^
            #   *EEE*     |
            #    \  |    -y
            #     \ |     |
            #      \|     
            #       B
            elif e1_y == e2_y:
                print("scenario 2")
                assert(e1_x != e2_x)
                vert_b = (max(e1_x,e2_x), e1_y+1) # right
                vert_a = (min(e1_x,e2_x), e1_y-1) # left
            #
            #   *---A
            #   |E  |   ^
            #   | E |   |
            #   |  E|  -y
            #   B---*   |
            else:
                print("scenario 3")
                min_x = min(e1_x, e2_x)
                min_y = min(e1_y, e2_y)

                max_x = max(e1_x, e2_x)
                max_y = max(e1_y, e2_y)

                vert_b = (min_x, max_y)
                vert_a = (max_x, min_y)
                
            print(f"VERT A: {vert_a}, VERT B : {vert_b}")


            vert_b_pnt = glm.vec3(*self.get_point(*vert_b)) if self.pixel_within_bounds(vert_b) else None
            vert_a_pnt = glm.vec3(*self.get_point(*vert_a)) if self.pixel_within_bounds(vert_a) else None

            side_b = None if vert_b_pnt is None else glm.dot(vert_b_pnt.xy - pnt1_2D , perp)

            side_a = None if vert_a_pnt is None else glm.dot(vert_a_pnt.xy - pnt1_2D , perp)
        
            print(f"DIR: {_dir} PERP: {perp}")
            ray_side = glm.dot(_dir, perp)

            if side_b is not None and ray_side * side_b > 0:
                
                pixels = [pix1, pix2, vert_b]
                coords = [self.get_point(*st) for st in pixels]
                tri = Triangle(coords, pixels)

            elif side_a is not None and ray_side * side_a > 0:
                pixels = [pix1, pix2, vert_a] 
                coords = [self.get_point(*st) for st in pixels]
                tri = Triangle(coords, pixels)
            else:

                tri = Triangle() # non-valid

            return tri

        else:
            assert( isinstance(vertex_or_edge, Vertex) )
            print("getting next Tri VERT")
            # pnt1 = vertex_or_edge.point_s[0]
            c, r = vertex_or_edge.pixel_s[0]

            topRight = c == self.num_cols - 1 and r == 0
            topLeft = c == 0 and r == 0
            bottomRight = c == self.num_cols - 1 and r == self.num_rows - 1
            bottomLeft = c == 0 and r == self.num_rows - 1

            triangles = []

            #print(f"index {index0} {index1} {index2} {index3} {index4} {index5}")

            # get all triangles given vertex position
            if topRight:
                print("top rigght")
                triangles = [self.make_tri_index(c,r,2)]

            elif topLeft:
                print("top left")
                triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1)]

            elif bottomRight:
                print("bottom rigght")
                triangles = [self.make_tri_index(c,r,4), self.make_tri_index(c,r,3)]

            elif bottomLeft:
                print("bottom left")
                triangles = [self.make_tri_index(c,r,5)]

            elif c == self.num_cols - 1: #right column
                print("right col")
                triangles = [self.make_tri_index(c,r,2), self.make_tri_index(c,r,3), self.make_tri_index(c,r,4)]

            elif (c == 0): #left column
                print("left col")
                triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1), self.make_tri_index(c,r,5)]
                
            elif (r == 0):  # top row
                print("bottom row")
                triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1), self.make_tri_index(c,r,2)]
            
            elif (r == self.num_rows - 1): # bottom row
                print("top row")
                triangles = [self.make_tri_index(c,r,3), self.make_tri_index(c,r,4), self.make_tri_index(c,r,5)]
                
            else: # interior
                print("interior")
                triangles = [self.make_tri_index(c,r,0), self.make_tri_index(c,r,1), self.make_tri_index(c,r,2), \
                              self.make_tri_index(c,r,3), self.make_tri_index(c,r,4), self.make_tri_index(c,r,5)  ]
                               
            chosen_tri = Triangle() # non-valid
            
            for i, _tri in enumerate(triangles):
                print(f"C R: {c}, {r}")
                print(f"VERTEX POINT: {_tri.point_s[0]}")
                print(f"VERTEX PIXEL: {_tri.pixel_s[0]}")
                print(f"EDGE 1 PIXEL: {_tri.pixel_s[1]}")
                print(f"EDGE 2 PIXEL: {_tri.pixel_s[2]}")
                e1 = glm.normalize(glm.vec3(_tri.point_s[1]).xy - glm.vec3(_tri.point_s[0]).xy)
                e2 = glm.normalize(glm.vec3(_tri.point_s[2]).xy - glm.vec3(_tri.point_s[0]).xy)


                e1_cross = cross_2D( _dir, e1)
                e2_cross = cross_2D( e2, _dir)
                print(f"E1 CROSS: {e1_cross}, E2 Cross: {e2_cross}")
                if e1_cross >= 0 and e2_cross > 0:
                    print(f"tri {i} passed")
                    chosen_tri = _tri
                    #return _tri
            return chosen_tri

    def point_on(self, pnt, edgePnt1, edgePnt2):

        AP = edgePnt1 - pnt
        AB = edgePnt1 - edgePnt2

        cros = cross_2D(AP, AB)

        if (abs(cros) > 1e-10):
            return False

        dot_both = glm.dot(AP, AB)
        dot_self = glm.dot(AB, AB)

        if dot_both < 0:
            return False

        if dot_both > dot_self:
            return False

        return True


    def pixel_within_bounds(self, pix):

        if pix[0] < 0 or pix[0] >= self.num_cols or\
            pix[1] < 0 or pix[1] >= self.num_rows:
            return False
        else:
            return True

    def get_candidate_tris(self, pixel_x, pixel_y):
    
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

    def get_component(self, _x_m, _y_m, _pixel_x, _pixel_y):

        tri1, tri2 = self.get_candidate_tris( _pixel_x, _pixel_y )

        tri1_stu = tri1.bary_x_y(_x_m, _y_m)
        tri2_stu = tri2.bary_x_y(_x_m, _y_m)

        print(f"TRI1 BARYCENTRIC (curr point): {tri1_stu[0]} {tri1_stu[1]} {tri1_stu[2]}")
        print(f"TRI2 BARYCENTRIC: {tri2_stu[0]} {tri2_stu[1]} {tri2_stu[2]}")

        # round to nearest 10 decimal places for barycentric assertion
        assert( all(round(e,10) >= 0 for e in tri1_stu) \
                or all(round(e,10) >= 0  for e in tri2_stu) )

        stu = tri1_stu if all(e >= 0  for e in tri1_stu) else tri2_stu
        tri = tri1 if all(e >= 0  for e in tri1_stu) else tri2
        print(f"CHOSEN TRIANGLE S: {stu[0]} T: {stu[1]} U: {stu[2]}")


        interp_z =  stu[0] * tri.point_s[0][2] + stu[1] * tri.point_s[1][2] + stu[2] * tri.point_s[2][2];      

        # on a vertex?

        tri.print_obj()

        print(f"points {tri.point_s[0]}, {tri.point_s[1]}, {tri.point_s[2]}")

        if (tri.point_s[0][0], tri.point_s[0][1]) == (_x_m, _y_m): #s
            print("returning ")
            point = [tri1.point_s[0]]
            st = [tri.pixel_s[0]]
            return Vertex(point, st), tri.point_s[0]

        elif (tri.point_s[1][0], tri.point_s[1][1]) == (_x_m, _y_m): #t
            point = [tri.point_s[1]]
            st = [tri.pixel_s[1]]
            return Vertex(point, st), tri.point_s[1]

        elif (tri.point_s[2][0], tri.point_s[2][1]) == (_x_m, _y_m): #u
            point = [tri.point_s[2]]
            st = [tri.pixel_s[2]]
            return Vertex(point, st), tri.point_s[2]


        # on an edge? previous was using barycentric coords for this
        # but would miss some edges
       
        pnt = glm.vec2(_x_m, _y_m)
        pnt_3D = (_x_m, _y_m, interp_z)

        if self.point_on( pnt, glm.vec2(tri.point_s[0]), glm.vec2(tri.point_s[1])):
            st = [tri.pixel_s[0], tri.pixel_s[1]]
            points = [tri.point_s[0], tri.point_s[1]]
            return Edge(points, st), pnt_3D

        elif self.point_on( pnt, glm.vec2(tri.point_s[1]), glm.vec2(tri.point_s[2])):
            st = [tri.pixel_s[1], tri.pixel_s[2]]
            points = [tri.point_s[1], tri.point_s[2]]
            return Edge(points, st), pnt_3D

        elif self.point_on( pnt, glm.vec2(tri.point_s[2]), glm.vec2(tri.point_s[0])):
            st = [tri.pixel_s[2], tri.pixel_s[0]]
            points = [tri.point_s[2], tri.point_s[0]]
            
            return Edge(points, st), pnt_3D

        else:
            return tri, pnt_3D

    def check_paralell(self, face, _curr_point, _dir=None):
        chosen_pixel = None
        for i,p in enumerate(face.point_s):
            p_xy = glm.vec2(p)
            
            if (p_xy == _curr_point):
                chosen_pixel = face.pixel_s[i]

        print(f"Chosen pixel: { chosen_pixel }")

        edges = face.get_edges()
        for j, e in enumerate(edges):
            pnt2 = e.point_s[1]
            pnt1 = e.point_s[0]
            e_pnt2 = glm.vec3(pnt2).xy
            e_pnt1 = glm.vec3(pnt1).xy
            seg = e_pnt2 - e_pnt1
            den = cross_2D(_dir, seg)

            print(f"edge: {e.pixel_s}")
            print(f"determinant: {den}")

            t = cross_2D(e_pnt1 - _curr_point, seg) / den
            u = cross_2D(e_pnt1 - _curr_point, _dir) / den
            print(f"u: {u:.7f} t: {t:.7f}")


     
    def check_edges( self, face, _dir, _curr_point, _curr_component ):
        
        if isinstance(_curr_component , Vertex):
            face.check_vert_intersection( _curr_component, _dir)
        
        
        edges = face.get_edges() # tuples

        min_t = math.inf 

        print(f"DIR: {_dir}")

        best_component = None # edge or vertex we hit 
        intersection_point = None

        for i, e in enumerate(edges):
            print(f"checking edge {i} w/ pixels {e.pixel_s}")

            if e == _curr_component:
                print("CURRENT COMPONENT EDGE")
                continue


            pnt2 = e.point_s[1]
            pnt1 = e.point_s[0]
            e_pnt2 = glm.vec3(pnt2).xy
            e_pnt1 = glm.vec3(pnt1).xy
            seg = e_pnt2 - e_pnt1
            den = cross_2D(_dir, seg)

            print(f"determinant: {den}")

            # if they're paralell might be vert to vert
            if den == 0: # paralell
            #    continue 
                if isinstance(_curr_component, Vertex):
                    # get other point
                    vert_pixels = _curr_component.pixel_s[0]
                    if vert_pixels == e.pixel_s[0]:
                        other_pixel = e.pixel_s[1]
                        other_point = e.point_s[1]

                    elif vert_pixels == e.pixel_s[1]:
                        other_pixel = e.pixel_s[0]
                        other_point = e.point_s[0]

                    print("FOUND PARALELL EDGE FOR VERT TO VERT")
                    best_component = Vertex([other_point], [other_pixel])

                    return best_component, other_point

                else:
                    continue

            t = cross_2D(e_pnt1 - _curr_point, seg) / den
            u = cross_2D(e_pnt1 - _curr_point, _dir) / den

            print(f"u: {u:5f} t: {t:5f}")
            if ( t > 0 and 0 <= u <= 1 ) and (t < min_t):
                print(f"resetting best edge to edge {i} w/ pixels {e.pixel_s}")
                min_t = t
                z = pnt1[2] + u * (pnt2[2] - pnt1[2])
                intersection_point_xy = _curr_point + t * _dir
                intersection_point = glm.vec3(intersection_point_xy, z)
                ix,iy,iz = intersection_point
                
                e1x,e1y = e_pnt1

                e2x,e2y = e_pnt2

                print(f"intersecion: {ix,iy} e1: {e1x,e1y} e2: {e2x,e2y}")


                # what is the right epsilon?

                #if glm.distance(intersection_point_xy , e_pnt1) < 0.20 :
                if intersection_point_xy == e_pnt1:
                    print("Edge intersection determined to be VERTEX 1 Intersection")
                    best_component = Vertex([e.point_s[0]], [e.pixel_s[0]])
                #elif glm.distance(intersection_point_xy , e_pnt2) < 0.20:
                elif intersection_point_xy == e_pnt2:
                    print("Edge intersection determined to be VERTEX 2 Intersection")
                    best_component = Vertex([e.point_s[1]], [e.pixel_s[1]])

                # if u == 0:
                #     best_component = Vertex([e.point_s[0]], [e.pixel_s[0]])
                # elif u == 1:
                #     best_component = Vertex([e.point_s[1]], [e.pixel_s[1]])
                else:
                    best_component = e
        return best_component, intersection_point

    #
    # returns curr_component (MeshComponent)
    #           , curr_face (MeshComponent)
    #           , curr_point_3D glm.vec3
    def get_next_point( self, curr_component, curr_point, end_point ):

        _dir = glm.normalize(end_point - curr_point);

        if curr_component is None:

            _pixel_x, _pixel_y = self.convert_to_pixel(*curr_point)
                
            
            if _pixel_x < 0 or _pixel_x >= (self.num_cols - 1) \
                or _pixel_y < 0 or _pixel_y >= (self.num_rows - 1):
                # current point outside of raster bounds!
                
                return None, None, None

                # _pixel_end_x, _pixel_end_y = self.convert_to_pixel(*curr_point)

                # if _pixel_end_x < 0 or _pixel_end_x >= (self.num_cols - 1) \
                #     or _pixel_end_y < 0 or _pixel_end_y >= (self.num_rows - 1):
                #     # end point outside of raster bounds!
                #     return None, None, None

                # else:
                #     # switch start and end points!
                

            curr_component, point = self.get_component(*curr_point, _pixel_x, _pixel_y)

            if isinstance(curr_component, Triangle):
                curr_face = curr_component
                return curr_component, glm.vec3(point), curr_face
            else:
                return curr_component, glm.vec3(point), None
                # curr_face = self.get_next_tri( curr_component, curr_point, _dir )

        elif isinstance(curr_component, Triangle):

            curr_face = curr_component 

        else:
            assert( isinstance(curr_component, Edge) or isinstance(curr_component, Vertex))
            
            # replace edge or vert with current face
            curr_face = self.get_next_tri( curr_component, curr_point, _dir )

        if not curr_face.valid:
            print("invalid face!!")
            return None, None, None
        else:
            print(f"Curr face pixels: {curr_face.pixel_s}")

        print(f"Curr face type: {type(curr_face)}, Curr component type: {type(curr_component)}")
        stu =  curr_face.bary_x_y(*curr_point)
        print(f"current face barycentric: {stu}")

        end_stu = curr_face.bary_x_y(*end_point)

        if ( all(elem >= 0 for elem in end_stu) ):
            print( "reached end face" )

            intersection_z =  end_stu[0] * curr_face.point_s[0][2] + end_stu[1] * curr_face.point_s[1][2] + end_stu[2] * curr_face.point_s[2][2];      
            return  None, glm.vec3(end_point, intersection_z) , curr_face

        best_component, intersection_point = self.check_edges( curr_face,  _dir, glm.vec2(curr_point), curr_component )
        
        if best_component is not None:
            print(f"from checking edges we see {best_component.pixel_s}")
        # 
        return best_component, intersection_point, curr_face

    #
    # using DEM create arrays of the steepest downhill slope &
    # angle relative to pixel per paper "A new method for the determination of flow directions
    #   and upslope areas in grid digital elevation model" (Tarboton, 1997)
    #
    def calculate_flows(self):

        steepest_s = 0

        for row in range(self.num_rows):
            for col in range(self.num_cols):
                
                #if row==0 and col==0:
                #    pass
                #else:
                #    continue

                existing_facets = self.get_avail_facets(col, row)

                steepest_up = 0 # upslope positive
                steepest_down = 0 # downslope negative
                ang_up  = 0
                ang_g_up = 0
                ang_down  = 0
                ang_g_down = 0

                e1_up = None
                e2_up = None
                e1_do_st = None #pixel loc
                e2_do_st = None #pixel loc
                e1_do = None # actual 3d point
                e2_do = None # actual 3d point

                up_bool = False
                down_bool = False
           
                for i in existing_facets:
                    print(f"{col} {row} Facet: {i}")
                    param_dict = self.get_param_index(col, row, i)

                    ac = param_dict["ac"]
                    af = param_dict["af"]

                    print(f"e0 pixel: {param_dict['e0']}")
                    print(f"e1 pixel: {param_dict['e1']}")
                    print(f"e2 pixel: {param_dict['e2']}")
                    e0 = glm.vec3(self.get_point(*param_dict["e0"]))
                    e1 = glm.vec3(self.get_point(*param_dict["e1"]))
                    e2 = glm.vec3(self.get_point(*param_dict["e2"]))

                    #e0 = e0_v.point
                    #e1 = e1_v.point
                    #e2 = e2_v.point

                    d1 = glm.distance(e0 , e1)
                    d2 = glm.distance(e1 , e2)

                    s1_up = (e1.z - e0.z) / d1;
                    s2_up = (e2.z - e1.z) / d2

                    # slope direction
                    r_up = math.atan2(s2_up, s1_up)
                    s_up = math.sqrt( s1_up**2 + s2_up**2 )

                    if r_up < 0: # below d1
                        r_up = 0
                        s_up = s1_up

                    if r_up > math.atan2(d2, d1):
                        r_up = math.atan2(d2, d1)
                        s_up = (e2.z - e0.z) / (d1**2 + d2**2)

                    r_g_up = af * r_up + ac * (math.pi / 2)

                    if steepest_up < s_up and s_up > 0:
                        up_bool = True
                        steepest_up = s_up
                        ang_g_up = r_g_up
                        ang_up = r_up
                        e1_up = param_dict["e1"]
                        e2_up = param_dict["e2"]

                    s1_down = (e0.z - e1.z) / d1;
                    s2_down = (e1.z - e2.z) / d2

                    # slope direction
                    r_down = math.atan2(s2_down, s1_down)
                    s_down = math.sqrt( s1_down**2 + s2_down**2 )

                    if r_down < 0: # below d1
                        r_down = 0
                        s_down = s1_down

                    if r_down > math.atan2(d2, d1):

                        r_down = math.atan2(d2, d1)
                        s_down = (e0.z - e2.z) / (d1**2 + d2**2)

                    r_g_down = af * r_down + ac * (math.pi / 2)

                    if steepest_down < s_down and s_down > 0:

                        #print(f"new best steep facet #{i} slope of {s_down}, local angle of {r_down} global angle of {r_g_down}")
                        down_bool = True
                        steepest_down = s_down
                        steepest_down_ad = param_dict["ad"]
                        ang_g_down = r_g_down
                        ang_down = r_down
                        e1_do_st = param_dict["e1"]
                        e2_do_st = param_dict["e2"]

                        e1_do = e1 
                        e2_do = e2

                #
                if up_bool:
                    #currently not capturing
                    pass
                    #up_s = Slope(steepest_up, ang_up, ang_g_up, e1_up, e2_up)
                    #e0_.upslope = up_s                    

                if down_bool:
                    #do_s = Slope(steepest_down, ang_down, ang_g_down, e1_do, e2_do)
                    #e0_.downslope = do_s

                    _which_e, deviation = self.get_delta(e0, e1_do, e2_do, ang_down, steepest_down_ad)
                    
                    downslope_st = e1_do_st if _which_e == 1 else e2_do_st

                    print(f"col: {col} row: {row} downslope pixel ({downslope_st.x}, {downslope_st.y})")

                    self.buffer_array[col, row, 4] = ang_g_down # downslope angle local
                    self.buffer_array[col, row, 5] = ang_down # downslope angle global
                    self.buffer_array[col, row, 6] = steepest_down # downslope slope
                    self.buffer_array[col, row, 7] = 1 # downslope exists
                    self.buffer_array[col, row, 8] = downslope_st.x # downslope pixel x
                    self.buffer_array[col, row, 9] = downslope_st.y # downslope pixel y
                    self.buffer_array[col, row, 10] = deviation.x # downslope deviation x
                    self.buffer_array[col, row, 11] = deviation.y # downslope deviation y
                
                else:

                    self.buffer_array[col, row, 7] = 0 # NO downslope slope


                if (steepest_down < steepest_s):
                    steepest_s = steepest_down

        # normalize flow magnitude
        #flow_mag = flow_mag / steepest_s

        self.steepest_slope = steepest_s

    # From Wu's paper "Nondispersive Drainage Direction Simulation Based
    #   on Flexible Triangular Facets"
    # Is angle pointing more toward e1 or e2 
    # what is the length difference from either e1/e2?
    #
    def get_delta(self, e0, e1, e2, ang, a_d):

        _e0 = glm.vec2(e0)
        _e1 = glm.vec2(e1)
        _e2 = glm.vec2(e2)

        adj_edge = glm.length( glm.vec2(_e1) - glm.vec2(_e0) );
        opp_edge = glm.length( glm.vec2(_e2) - glm.vec2(_e1) );

        tot_ang = math.atan2(opp_edge, adj_edge) 

        if ang > ( tot_ang/2 ): #e2 case
            deviation = opp_edge - (adj_edge * math.tan(ang)) 
            dest_vert = 2
            chosen_a_d = a_d[1]
        else: #e1 case
            deviation = adj_edge * math.tan(ang)
            dest_vert = 1
            chosen_a_d = a_d[0]

        print(f"destin vert : { dest_vert }")
        print(f"Before multiply deviation : {deviation }")
        return (dest_vert, deviation * chosen_a_d)

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
        elif c == 0:
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

    #
    # for given point get the point that is downslope, ()
    #
    def get_downslope_point(self, c, r):

        if (self.buffer_array[c, r, 7] == 1):
            return (int(self.buffer_array[c, r, 8]), int(self.buffer_array[c, r, 9]))
        else:
            return ()


    #
    # recursive way to get # of upslope pixels
    #
    def get_upslope(self, c , r):

        neighbors = self.get_neighbors(c, r)
        tot_pix = 1 
        for n_c, n_r in neighbors:

            candidate_tup = self.get_downslope_point(n_c, n_r)

            # does it drain into c,r ?
            if candidate_tup and candidate_tup[0] == c and candidate_tup[1] == r:
                
                tot_pix += self.get_upslope(n_c, n_r)

        return tot_pix


    #
    # for every pixel get the number of pixels that
    # drain into it (have to run calculate_flows first)
    #
    def calc_upslope_new_point(self):
        for row in range(self.num_rows):
            for col in range(self.num_cols):

                upslope_amt = self.get_upslope(col, row)

                self.buffer_array[col, row, 12] = upslope_amt

        for row in range(self.num_rows):
            for col in range(self.num_cols):

                neighbors = self.get_neighbors(col, row)

                adjusted_point = glm.vec2( self.get_point(col,row ) )
                deviation_tot = glm.vec2(0,0)
                upslope_tot = 1
                for n_c, n_r in neighbors:

                    candidate_tup = self.get_downslope_point(n_c, n_r)

                    if candidate_tup and candidate_tup[0] == col and candidate_tup[1] == row:
                        
                        deviation = glm.vec2(self.buffer_array[n_c, n_r, 10], self.buffer_array[n_c, n_r, 11])
                        upslope_amt = self.buffer_array[col, row, 12]

                        upslope_tot += upslope_amt
                        deviation_tot += upslope_amt*deviation

                # update
                adjusted_point += (deviation_tot / upslope_tot)

                # tri, interp_pnt = td.get_tri_interp_point( adjusted_point.x, adjusted_point.y )

                self.buffer_array[col, row, 13] = adjusted_point.x
                self.buffer_array[col, row, 14] = adjusted_point.y

                # notice we're not getting interpolated z value, 
                # that will be part of line drawing

    #
    # read in GPKG file
    # filter to TIF bounding box
    # return array of array of point tuples
    # 
    def read_filter_gpkg(self, gpkg_file_path):

        datasource = ogr.Open(gpkg_file_path, 0)

        layer = datasource.GetLayerByName("lines")

        if layer is None:
            print(f"Layer not found!")

        linear_ring = ogr.Geometry(ogr.wkbLinearRing)
        layer.SetSpatialFilterRect(*self.bounds)

        layer.ResetReading()
        linestring_coords = []

        for feature in layer:
            
            geom = feature.GetGeometryRef()
            # seems like all linestrings

            # following not working for our created LineStrings
            # geom.GetGeometryType() == ogr.wkbLineString
            #
            if geom is not None and geom.GetGeometryName() == "LINESTRING":
                coords = []
                for p in range(geom.GetPointCount()):
                    # Get X and Y coordinates (Z is optional)
                    coords.append((geom.GetX(p), geom.GetY(p)))
                linestring_coords.append(coords)

        return linestring_coords
    #
    # recursive function
    #
    def create_flowline(self, c, r):
        # mark as visited
        downslope_tup = self.get_downslope_point(c, r)

        adj_pnt = (self.buffer_array[c, r, 13].item() , self.buffer_array[c, r, 14].item()) #, self.buffer_array[c, r, 15])  


        if self.buffer_array[c, r, 16] == 1: # been visited
            return [adj_pnt]
        elif downslope_tup: # has downslope
            self.buffer_array[c, r, 16] = 1
            return [adj_pnt] + self.create_flowline(downslope_tup[0], downslope_tup[1])
        else:
            return [adj_pnt]


    def get_flowlines(self):

        #self.buffer_array[c, r, 12] = 0.0 # number of upslope
        #self.buffer_array[c, r, 13] = 0.0 # x of adjusted position
        #self.buffer_array[c, r, 14] = 0.0 # y of adjusted position
        #self.buffer_array[c, r, 15] = 0.0 # interp elev of adjusted position
        all_lines = []
        for row in range(self.num_rows):
            for col in range(self.num_cols):

                # only get where there's no upslope and is downslope
                if self.buffer_array[col, row, 12] == 1 and self.buffer_array[col, row, 7] == 1:
                    
                    flow_line = self.create_flowline( col, row )
                    if flow_line:
                        all_lines.append(flow_line)
        return all_lines
    
    #
    # for each of 8 triangular facets get vertex index + conversion factors for 
    # global direction
    #
    def get_param_index(self, c, r , facet):
        #
        # returns e_0, e_1, e_2, a_c, a_f
        # o is nother orientation (1 or -1)
        #
        e_0 = glm.vec2(int(c), int(r))
        match facet:
            case 1:
                #r * self.num_cols)+c;
                e_1 = glm.vec2(int(c+1), int(r)) #(r * self.num_cols)+c+1;
                e_2 = glm.vec2(int(c+1), int(r-1)) #((r - 1) * self.num_cols)+c+1;
                a_f = 1
                a_c = 0
                a_d = (glm.vec2(0, 1), glm.vec2(0, -1) )
            case 2:
                e_1 = glm.vec2(int(c), int(r-1)) #((r - 1) * self.num_cols)+c;
                e_2 = glm.vec2(int(c+1), int(r-1)) #((r - 1) * self.num_cols)+c+1;
                a_f = -1
                a_c = 1
                a_d = (glm.vec2(1, 0), glm.vec2(-1, 0) )
            case 3:
                e_1 = glm.vec2(int(c), int(r-1))  #((r - 1) * self.num_cols)+c;
                e_2 = glm.vec2(int(c-1), int(r-1))  #((r - 1) * self.num_cols)+c-1;
                a_f = 1
                a_c = 1
                a_d = ( glm.vec2(-1, 0), glm.vec2(1, 0) )
            case 4:
                e_1 = glm.vec2(int(c-1), int(r))  #(r * self.num_cols)+c-1;
                e_2 = glm.vec2(int(c-1), int(r-1))  #((r - 1) * self.num_cols)+c-1;
                a_f = -1
                a_c = 2 # same for 4 and 5
                a_d = ( glm.vec2(0, 1), glm.vec2(0, -1) )
            case 5:
                e_1 = glm.vec2(int(c-1), int(r))  #(r * self.num_cols)+c-1;
                e_2 = glm.vec2(int(c-1), int(r+1))  #((r + 1) * self.num_cols)+c-1;
                a_f = 1
                a_c = 2 # same for 4 and 5
                a_d = ( glm.vec2(0, -1), glm.vec2(0, 1) )
            case 6:
                e_1 = glm.vec2(int(c), int(r+1)) #((r + 1) * self.num_cols)+c;
                e_2 = glm.vec2(int(c-1), int(r+1)) #((r + 1) * self.num_cols)+c-1;
                a_f = -1
                a_c = 3
                a_d = ( glm.vec2( -1,0 ), glm.vec2(1, 0) )
            case 7:
                e_1 = glm.vec2(int(c), int(r+1)) #((r + 1) * self.num_cols)+c;
                e_2 = glm.vec2(int(c+1), int(r+1)) #((r + 1) * self.num_cols)+c+1;
                a_f = 1
                a_c = 3
                a_d = ( glm.vec2( 1,0 ), glm.vec2(-1, 0) )
            case 8:
                e_1 = glm.vec2(int(c+1), int(r)) #(r * self.num_cols)+c+1;
                e_2 = glm.vec2(int(c+1), int(r+1)) #((r + 1) * self.num_cols)+c+1;
                a_f = -1
                a_c = 4
                a_d = ( glm.vec2( 0,-1), glm.vec2(0, 1) )

        return {"e0": e_0, "e1": e_1, "e2": e_2, "ac": a_c, "af" : a_f, "ad" : a_d}

    # given that a pixel might be on DEM border
    def get_avail_facets(self, c, r):

        topRight = c == self.num_cols - 1 and r == 0;
        topLeft = c ==  0 and r == 0;
        bottomRight = c == self.num_cols - 1 and r == self.num_rows - 1;
        bottomLeft = c == 0 and r == self.num_rows - 1;

        if topRight:
            facets = [5,6]
        elif topLeft:
            facets = [7,8]
        elif bottomRight:
            facets = [3,4]
        elif bottomLeft:
            facets = [1,2]
        elif c == self.num_cols - 1: 
            #right column
            facets = [3,4,5,6]
        elif (c == 0): 
            #left column
            facets = [1,2,7,8]
        elif (r == self.num_rows - 1): 
            # bottom row
            facets = [1,2,3,4]
        elif (r == 0): 
            # top row
            facets = [5,6,7,8]
        else: # interior
            facets = list(range(1,9))

        return facets


    #
    # with multiple line strings of 2D points create lines
    #
    def create_lines(self, line_strings):

        all_points = []

        for _it, ls in enumerate(line_strings):
                
            section_points = []
            for i in range(len(ls)-1):

                reverse = False

                print(ls[i])
                pointA = glm.vec2(ls[i])
                pointB = glm.vec2(ls[i+1])
                curr_component, intersection, curr_face = self.get_next_point(None, pointA, pointB )

                # this means we're starting outside raster
                if curr_component is None:

                    print("Checking reverse!")
                    reverse = True
                    pointA = glm.vec2(ls[i+1])
                    pointB = glm.vec2(ls[i])

                    curr_component, intersection, curr_face = self.get_next_point(None, pointA, pointB )

                    print(f"reverse intersection: {intersection}")

                    if curr_component is None:
                        print("Reverse was NULL TOO")
                        continue # next iter
                        


                reached_end = False

                if (i == 0 ): # we should have intersection from previous segment

                    # AHHHHH-ful
                    if curr_face is not None:
                        curr_face.calc_norm_axes()

                    if curr_component is not None:
                        if isinstance(curr_component, Triangle):
                            curr_face.calc_norm_axes()
                        else:
                            curr_component.normal = self.get_norm_vert_edge(curr_component)


                    section_points = [ {"point": intersection, "face": curr_face, "component": curr_component} ]

                while(not reached_end):
                    curr_component, intersection, curr_face = self.get_next_point( curr_component , glm.vec3(intersection).xy , pointB )

                    # means we've reached second point
                    # or border
                    if (curr_component is None):
                        reached_end = True

                    else:
                        print(f"\tLS: {_it}, Point Index: {i} More than 1 point!")
                    
                    print(f"intersection: {intersection}")

                    # intersection will be None only if we're at 
                    # a border and are going in direction outside it
                    if (intersection is not None):

                        if curr_face is not None:
                            curr_face.calc_norm_axes()

                        if curr_component is not None:
                            if isinstance(curr_component, Triangle):
                                curr_face.calc_norm_axes()
                            else:
                                curr_component.normal = self.get_norm_vert_edge(curr_component)


                        section_points.append({"point": intersection, "face": curr_face, "component": curr_component})

                    #reached_end = True 

                print(f"line string {i} size: {len(section_points)}")
            
            all_points.append(section_points)

        return all_points


    def write_line_obj(self, _lines, obj_file_name ,scale):
        with open(obj_file_name, "w") as f2:
            f2.write("# OBJ flow lines\n")

            vertex_index = 1  # OBJ is 1-based

            for line in _lines:
                if len(line) < 2:
                    continue  # need at least 2 points

                indices = []

                # write vertices

                for pnt_dict in line:

                    point = pnt_dict["point"]
                    curr_face = pnt_dict["face"]
                    curr_comp = pnt_dict["component"]

                    #{"point": intersection, "face": curr_face, "component": curr_component}
                    if curr_comp is not None:
                        norm = curr_comp.normal
                    else:
                        norm = curr_face.normal

                    assert(norm is not None)

                    print(f"normal: {norm}")
                    normal_shifted_point = point + 5 * norm

                    print(f"original point: { point}")
                    v_new = (normal_shifted_point - glm.vec3(self.x_center, self.y_center, 0)) * scale

                    print(f"adjusted point: {v_new}")

                    #v_new = glm.vec3( ( v.x - _x_center ) * scale, ( v.y - _y_center ) * scale, v.z * scale )

                    f2.write(f"v {v_new.x} {v_new.y} {v_new.z}\n")

                    indices.append(vertex_index)

                    vertex_index += 1

                # write line element
                indices_str = " ".join(map(str, indices))
                f2.write(f"l {indices_str}\n\n")


    def write_obj( self, path , flow_lines, road_lines, scale):

        vertices = self.get_points_centered_scaled( scale )
        faces = self.test_index_array

        with open(path + ".obj", "w") as f:
            f.write("# OBJ DEM mesh\n")

            # write vertices
            for v in vertices:
                f.write(f"v {v[0]} {v[1]} {v[2]}\n")

            # write faces (OBJ is 1-based)
            
            for i, j, k in faces:
                f.write(f"f {i+1} {j+1} {k+1}\n")

        print("writing flow lines to OBJ")
        self.write_line_obj(  flow_lines, path + "_flow_lines.obj", scale )
        
        print("writing road lines to OBJ")
        self.write_line_obj(  road_lines, path + "_road_lines.obj", scale )

 
    def plot_point_comp(self, fig_name, fpath, curr_component, intersection, curr_face ):

        fig, ax = plt.subplots(figsize=(14, 14))
        with rasterio.open( fpath ) as src:

            height = src.shape[0]
            width = src.shape[1]
            
            band = src.read(1)

            cols, rows = np.meshgrid(np.arange(width), np.arange(height))
            xs, ys = rasterio.transform.xy(src.transform, rows, cols)

            extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
            x_min, x_max, y_min, y_max = extent            
            im = ax.imshow(band, extent=extent, cmap='viridis', origin='upper')
            triang = mtri.Triangulation(xs, ys, self.test_index_array)
            #plt.tripcolor(triang, xs, ys, edgecolors='k', cmap='viridis')
            plt.triplot(triang, 'go-', label='Triangular Mesh',color='grey', linewidth=1)



            im = ax.imshow(band, extent=extent, cmap='viridis', origin='upper')

            if isinstance(curr_component , Triangle):

                curr_face = curr_component

            if curr_face is None:
                pass 

            else:
                face_points_2D = [ glm.vec3(elem).xy for elem in curr_face.point_s ]
                triangle = Polygon(face_points_2D , facecolor='orange', edgecolor='orange', linewidth=2)
                ax.add_patch(triangle)


            if curr_component is None:
                pass
            else:
                x, y = self.sep_x_y([ glm.vec3(elem).xy for elem in curr_component.point_s ] )

                ax.plot(x,y, marker='o', markersize=8, color="red", linewidth=2.0, zorder=4)


            p_x , p_y = intersection.xy
            ax.plot([p_x], [p_y], marker='*', markersize=2, color="black", linewidth=1.0, zorder=5)

            left_limit = extent[0]
            right_limit = extent[0] + (extent[1] - extent[0]) / 3.0

            left_limit = x_min
            right_limit = x_min + (x_max - x_min) / 3
            bottom_limit = y_max - (y_max - y_min) / 3
            top_limit = y_max
            ax.set_xlim(left_limit, right_limit)
            ax.set_ylim(bottom_limit, top_limit)

            plt.savefig( fig_name, bbox_inches='tight', dpi=100)

    def sep_x_y(self, list_points):
        if list_points:
            x, y = zip(*list_points)
            return list(x), list(y)
        else:
            return [],[]

    def animation_plot( self, line_comps, save_loc ):

        fig, ax = plt.subplots(figsize=(12, 12))

        num_frames = len(line_comps) - 1 

        def sep_x_y(list_points):
            if list_points:
                x, y = zip(*list_points)
                return list(x), list(y)
            else:
                return [],[]

        with rasterio.open(file_path) as src:
            height = src.shape[0]
            width = src.shape[1]
            
            cols, rows = np.meshgrid(np.arange(width), np.arange(height))
            xs, ys = rasterio.transform.xy(src.transform, rows, cols)

            extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

            im = ax.imshow(band, extent=extent, cmap='viridis', origin='upper')
            triang = mtri.Triangulation(xs, ys, self.test_index_array)
            #plt.tripcolor(triang, xs, ys, edgecolors='k', cmap='viridis')
            plt.triplot(triang, 'go-', label='Triangular Mesh',color='grey', linewidth=1)
            ax.autoscale(False) 
            _init_face_points = [ [0,0] ]
            #_init_face_points_2D = [ glm.vec3(elem).xy for elem in _init_face_points ]

            print(_init_face_points)
            triangle = Polygon(_init_face_points, facecolor='orange', edgecolor='orange', linewidth=2)
            ax.add_patch(triangle)
            ln_comp, = ax.plot([], [], marker='o', markersize=8, color="red", linewidth=2.0)
            ln_road, = ax.plot([], [], 'black')


        def update(frame):
        
            # frame should be an INT

            i = frame #min(len(line_comps)-1, frame)

            road_points = [ glm.vec3(elem["point"]).xy for elem in line_comps[:i] ]
            #print(road_points)

            if line_comps[i]["face"] is not None:
                face_points = line_comps[i]["face"].point_s
                face_points_2D = [ glm.vec3(elem).xy for elem in face_points ]
            else:
                face_points_2D = _init_face_points

            print(f"i: {i}, fram: {frame}")

            if line_comps[i]["component"] is not None:

                if isinstance(line_comps[i]["component"], Vertex):
                    print("Vertex Component!!")
                comp_points = sep_x_y([ glm.vec3(elem).xy for elem in line_comps[i]["component"].point_s ])
            else:
                comp_points = [],[]

            triangle.set_xy( face_points_2D )
            
            ln_road.set_data( sep_x_y( road_points ) )
            ln_comp.set_data(  comp_points  )
            #ln_face.set_data(sep_x_y( face_points_2D ))

            return ln_road,ln_comp, triangle


        ani = FuncAnimation(fig, update, frames=np.linspace(0, num_frames, num_frames+1, dtype=int), blit=True, interval=1000)
        ani.save(save_loc, writer='pillow')
        #plt.show()






def plot(file_path, mc, points):
       with rasterio.open(file_path) as src:
            height = src.shape[0]
            width = src.shape[1]
            
            cols, rows = np.meshgrid(np.arange(width), np.arange(height))
            xs, ys = rasterio.transform.xy(src.transform, rows, cols)
            
            fig, ax = plt.subplots( figsize=(12, 12) )

            extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

            im = ax.imshow(band, extent=extent, cmap='viridis', origin='upper')
            triang = mtri.Triangulation(xs, ys, mc.test_index_array)
            plt.triplot(triang, 'go-', label='Triangular Mesh')
            for i in range(height):
                for j in range(width):

                    index = i * width + j

                    text_label = f"({xs[index]:.1f}, {ys[index]:.1f})"
        
                    #ax.text(xs[index], ys[index], text_label, ha='center', va='center', fontsize=8, color='red')


            _x_points, _y_points = zip(*points)
            plt.plot( _x_points, _y_points, marker='o', color='red' )
            plt.show()
                    #-0.5, pnt.col+0.35


def create_esri_raster(self, file_path, low_pnt=(0,0)):

    x_origin, y_origin = 0, 0
    pixel_width = 45
    pixel_height = -45

    height, width = 3,3 

    data = np.full((height, width), 200).astype(np.float32) # Example: random data
    transform = Affine(pixel_width, 0, x_origin, 0, pixel_height, y_origin)

    data[1, 1] = 0
    data[int(low_pnt[0]), int(low_pnt[1])] = 500
    #data[0, 0] = 200
    #data[0, 1] = 200
    #data[0, 2] = 200
    #data[1, 0] = 200
    #data[1, 1] = 500
    #data[1, 2] = 200
    #data[2, 0] = 0
    #data[2, 1] = 200
    #data[2, 2] = 200
    esri_102009_crs = rasterio.crs.CRS.from_string("ESRI:102009")
    # Define the raster profile
    profile = {
        'driver': 'GTiff',  # GeoTIFF format
        'dtype': data.dtype,
        'height': height,  # rows
        'width': width,  # columns
        'count': 1,  # bands
        'crs': esri_102009_crs,
        'transform': transform  # Affine transform
    }

    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(data, 1)  # Write the data to band 1



if __name__ == "__main__":
    
    #file_path = "Accessibility/data/sphere_11_new_v2_reproject_.tif"
    #file_path = "Accessibility/data/clipped_collawash.tif"

    #file_path = r"../Accessibility/data/clipped_collawash_v2.tif"
    #file_path = r"../Accessibility/data/wu_raster.tif"
    file_path = r"data/clipped_collawash_v3.tif"
    gpkg_file = r"data/ziv_esri_roads_v8.gpkg"

    # assumes a geotiff is in ESRI:102009!
    # VOID method updating self.buffer_array
    #  
    with rasterio.open(file_path) as src:

        band = src.read(1)
        bounds = src.bounds
            
        transf = src.transform

        x_min, y_min = rasterio.transform.xy(transf, 0, 0)
        _min = glm.vec3(x_min, y_min, 0)

        # set global
        num_rows = band.shape[0]
        num_cols = band.shape[1]

        buffer_array = np.zeros((num_cols, num_rows, MeshComponentCreator.ENTRIES_PER_BUFFER), dtype=np.float32)

        #ds = gdal.Open(file_path, gdal.GA_ReadOnly)
        #geotransform = ds.GetGeoTransform()

        for r in range(num_rows):
            for c in range(num_cols):

                x_m, y_m = rasterio.transform.xy(transf, r, c)
                elev = band[r, c]

                if elev < 0:
                    elev = 0

                # should be in meters
                # x_m, y_m = transformer.transform(lon, lat)

                # !!!
                #x_m = x_m - (0.5*geotransform[1])
                #y_m = y_m - (0.5*geotransform[5])

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


        promesheus = MeshComponentCreator(buffer_array, transf, bounds)

        # we pretty much only need vertex normals for moving line off
        promesheus.calculate_vertex_normals()

        promesheus.print_stats()
        
        road_line_strings = promesheus.read_filter_gpkg(gpkg_file)


        #plot(file_path, promesheus, [[buffer_array[10, 6, 0], buffer_array[10, 6, 1]]] )

        promesheus.calculate_flows()

        promesheus.calc_upslope_new_point()

        flow_lines = promesheus.get_flowlines()


        # convert road lines to 3D lines
        road_lines_comp_3D = promesheus.create_lines( road_line_strings )

        # print(f"number road lines: {len(road_lines_3D)}")

        # convert flow lines to 3D lines
        flow_lines_comp_3D = promesheus.create_lines( flow_lines )

        obj_path = r"../Accessibility/data/ziv_collawash_flow_mesh_50_v9"
        scale = 0.008

    
        promesheus.write_obj( obj_path, flow_lines_comp_3D, road_lines_comp_3D , scale )


        # promesheus.animation_plot( road_lines_comp_3D[0], "../Accessibility/data/road_line_animation_january_v4.gif")

 
        #SCRIPT
        print("script")
        points = [ [(-1922879.25, 854802.3125), (-1922326.0, 854249.125)] ]

        pointA = glm.vec2(points[0][0])
        pointB = glm.vec2(points[0][1])

        curr_component, intersection, curr_face = promesheus.get_next_point(None, pointA, pointB )

        #fig_ = "geodesic_plt1.png"
        #promesheus.plot_point_comp( fig_ , file_path, curr_component, intersection, curr_face )
        
        
        curr_component, intersection, curr_face = promesheus.get_next_point( curr_component , glm.vec3(intersection).xy , pointB )
        #fig_ = "geodesic_plt2.png"
        #promesheus.plot_point_comp( fig_ , file_path, curr_component, intersection, curr_face )
        
        
        curr_component, intersection, curr_face = promesheus.get_next_point( curr_component , glm.vec3(intersection).xy , pointB )
        #fig_ = "geodesic_plt3.png"
        #promesheus.plot_point_comp( fig_ , file_path, curr_component, intersection, curr_face )
        
        #fig_ = "geodesic_plt4.png"
        #promesheus.plot_point_comp( fig_ , file_path, curr_component, intersection, curr_face )
        
        curr_component, intersection, curr_face = promesheus.get_next_point( curr_component , glm.vec3(intersection).xy , pointB )

        #fig_ = "geodesic_plt5.png"
        #promesheus.plot_point_comp( fig_ , file_path, curr_component, intersection, curr_face )
        
        curr_component, intersection, curr_face = promesheus.get_next_point( curr_component , glm.vec3(intersection).xy , pointB )

        fig_ = "geodesic_plt6.png"
        promesheus.plot_point_comp( fig_ , file_path, curr_component, intersection, curr_face )
              


        #promesheus.animation_plot( flow_lines_comp_3D[21], "Accessibility/data/flow_line_animation_final.gif")

        #vertices = promesheus.get_points_centered_scaled(scale)

        #faces = promesheus.test_index_array


        # road_ls=[]
        # for i,line in enumerate(road_lines_comp_3D):
        #     d = {"index": i, "line": line}
        #     print(i)
        #     if line:
        #         y_max = max([item["point"].y for item in line])
        #         y_min = min([item["point"].y for item in line])
        #     else:
        #         y_max = 1000000
        #     d["y_max"] = y_max
        #     road_ls.append(d)


        # for i, line in enumerate(road_line_strings):
        #     y_max_road = max([item[1] for item in line])
        #     y_min_road = min([item[1] for item in line])
        #     x_max_road = max([item[0] for item in line])
        #     x_min_road = min([item[0] for item in line])
        # print(f"ROAD min x: {x_min_road}, max x: {x_max_road}, min y: {y_min_road}, max y: {y_max_road}")

        # for i, line in enumerate(flow_lines):
        #     y_max_flow = max([item[1] for item in line])
        #     y_min_flow = min([item[1] for item in line])
        #     x_max_flow = max([item[0] for item in line])
        #     x_min_flow = min([item[0] for item in line])
        # print(f"FLOW min x: {x_min_flow}, max x: {x_max_flow}, min y: {y_min_flow}, max y: {y_max_flow}")

        
       

