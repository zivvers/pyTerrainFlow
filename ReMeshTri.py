import numpy as np
from MeshComponent import MeshComponent, Triangle
from shapely.geometry import Polygon
from shapely.ops import triangulate

class ReMeshTri:

    #
    # triangle: Triangle (to clip)
    # , clip_points (array of 2D tuples)
    #
    #
    def sutherland_hodgman_clip(self, tri, clip_points ):
        
        # convert triangle verts from 3D to 2D
        poly_points = [ pnt[:2] for pnt in tri.point_s ]

        for i in range(len(clip_points)):
            k = (i+1) % len(clip_points) 

            poly_points = self.clip( poly_points, clip_points[i], clip_points[k])

        return poly_points

    # intersection with edge and the second point is added
    def get_intersection(self, pnt1, pnt2, pnt3, pnt4):

        x1,y1 = pnt1
        x2,y2 = pnt2
        x3,y3 = pnt3
        x4,y4 = pnt4

        x_num = (x1*y2 - y1*x2) * (x3-x4) - (x1-x2) * (x3*y4 - y3*x4)
        x_den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4)

        y_num = (x1*y2 - y1*x2) * (y3-y4) - (y1-y2) * (x3*y4 - y3*x4)
        y_den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4)

        return ( x_num/x_den, y_num/y_den )

    # # Function to return y-value of point of intersection of two lines
    # def y_intersect(x1, y1, x2, y2, x3, y3, x4, y4):
    #     num = (x1*y2 - y1*x2) * (y3-y4) - (y1-y2) * (x3*y4 - y3*x4)
    #     den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4)
    #     return num/den

    #
    # triangle_verts: vertices of triangle to clip
    # , edge: current clip edge
    #
    #
    def clip(self, poly_verts, clip_point1, clip_point2):

        poly_size = len(poly_verts)

        clp1_x, clp1_y = clip_point1
        clp2_x, clp2_y = clip_point2
        new_points = []
        #for i, vert in enumerate(tri_verts):
        for i in range(poly_size):

            k = (i+1) % poly_size
            ix, iy = poly_verts[i]
            kx, ky = poly_verts[k]

            # get position of first point with regards to clipper line
            i_pos = (clp2_x - clp1_x) * (iy - clp1_y) - (clp2_y - clp1_y) * (ix - clp1_x)

            # get position of second point with regards to clipper line
            k_pos = (clp2_x - clp1_x) * (ky - clp1_y) - (clp2_y - clp1_y) * (kx - clp1_x)


            # case 1 - both points are outside
            if i_pos >= 0 and k_pos >= 0:
                new_points.append((kx, ky))

            # case 2 - only first point is inside
            elif i_pos < 0 and k_pos >= 0:

                new_pnt = self.get_intersection( clip_point1, clip_point2, (ix, iy), (kx, ky) )

                new_points.append( new_pnt )
                new_points.append( (kx, ky) )

            # case 3 - only second point is inside
            elif i_pos >= 0 and k_pos < 0:
                # Only point of intersection with edge is added
                new_pnt = self.get_intersection( clip_point1, clip_point2, (ix, iy), (kx, ky) )
                new_points.append( new_pnt )

            # case 4 - both points are outside
            else:
                pass

        print(f"new points {new_points}")
        return new_points


if __name__ == "__main__":

    clip_points = [(-1922879.25, 854800.25),
                     (-1922880.25, 854801.25),
                     (-1922878.25, 854803.375),
                     (-1922877.125, 854802.3125)]
    clip_points.reverse()

    pts = [(-1922904.375, 854802.3125, 106.70935821533203),
                (-1922879.25, 854802.3125, 107.40161895751953),
                (-1922879.25, 854777.125, 108.56053161621094)]

    pts_2D = [pt[:2] for pt in pts]
    pts_2D.reverse()

    sts = [(3, 4), (4, 4), (4, 5)]

    target_tri = Triangle(pts, sts)

    rm = ReMeshTri()

    new_points = rm.sutherland_hodgman_clip( target_tri, clip_points )




