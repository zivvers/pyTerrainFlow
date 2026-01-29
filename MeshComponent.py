import numpy as np
from pyglm import glm
from collections import Counter

#
# should have separate utilities library!
def cross_2D( a, b ):
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


#
# this is how we're going to keep track
#
class MeshComponentContainer:


    def __init__(self, curr_face=None , curr_component=None , curr_point=None):
        self.face = curr_face
        self.component = curr_component
        self.point = curr_point

    def is_none( self ):

        if self.face is None and self.component is None :
            return True
        else:
            return False


