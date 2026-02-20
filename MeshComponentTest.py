import unittest
from matplotlib import container
import rasterio
import numpy as np
from pyglm import glm
from MeshComponentCreator import MeshComponentCreator
from MeshComponent import MeshComponentContainer, Vertex, Edge, Triangle
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import sys as sys


# def plot(file_path, mc, point):
#        with rasterio.open(file_path) as src:
#             height = src.shape[0]
#             width = src.shape[1]
            
#             cols, rows = np.meshgrid(np.arange(width), np.arange(height))
#             xs, ys = rasterio.transform.xy(src.transform, rows, cols)
            
#             fig, ax = plt.subplots( figsize=(12, 12) )
#             triang = mtri.Triangulation(xs, ys, mc.test_index_array)

#             for i in range(height):
#                 for j in range(width):

#                     index = i * width + j

#                     text_label = f"({xs[index]:.1f}, {ys[index]:.1f})"
        
#                     ax.text(xs[index], ys[index], text_label, ha='center', va='center', fontsize=8, color='red')

#             _x_t, _y_t= point
#             plt.plot( _x_t, _y_t, marker='o', color='blue' )
                    #-0.5, pnt.col+0.35


def instantiate_from_tif(file_path):

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

        for r in range(num_rows):
            for c in range(num_cols):

                x_m, y_m = rasterio.transform.xy(transf, r, c)
                elev = band[r, c]

                # should be in meters
                # x_m, y_m = transformer.transform(lon, lat)

                buffer_array[c, r, 0] = x_m #2ND
                buffer_array[c, r, 1] = y_m
                buffer_array[c, r, 2] = elev
                buffer_array[c, r, 3] = 0.0
                #self.buffer_array[c, r, 4] = r # Umm do we need pixel loc? #3RD
                #self.buffer_array[c, r, 5] = c
                #self.buffer_array[c, r, 6] = 0.0
                #self.buffer_array[c, r, 7] = 0.0
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
        return promesheus 



#
# integration type tests
#
class MeshComponentTest(unittest.TestCase):

    __test__ = True
    
    @classmethod
    def setUpClass(cls):
        cls.file_path = "data/sphere_11_new_reproj.tif"
        cls.promesheus = instantiate_from_tif( cls.file_path )

        cls.file_path = "data/clipped_collawash_v3.tif"
        cls.promesheus2 = instantiate_from_tif( cls.file_path )        

    def test_get_component(self):

        tri1, tri2 = self.promesheus.get_candidate_tris(0,0)

        point1 = (-1803732.875, 761263.25)
        pixel_x, pixel_y = self.promesheus.convert_to_pixel(*point1)
        _1_container = self.promesheus.get_component( *point1 , pixel_x, pixel_y)

        point2 = (-1799191.375 + 10 , 761263.25)
        pixel_x, pixel_y = self.promesheus.convert_to_pixel(*point2)
        _2_container = self.promesheus.get_component(*point2, pixel_x, pixel_y)

        point3 = (-1840065.125, 779429.375)
        pixel_x, pixel_y = self.promesheus.convert_to_pixel(*point3)
        _3_container = self.promesheus.get_component(*point3, pixel_x, pixel_y)

        _1 = _1_container.component
        _2 = _2_container.component
        _3 = _3_container.component

        assert(_1.pixel_s[0] == (9, 5))

        print(f"first type: {type(_1)} index: {_1.pixel_s[0]}")
        print(f"second type: {type(_2)}")
        print(f"third type: {type(_3)}")

    def test_get_component_2(self):

        point1 = (-1803732.875, 761263.25)
        pixel_x, pixel_y = self.promesheus.convert_to_pixel(*point1)
        _1_container = self.promesheus.get_component( *point1 , pixel_x, pixel_y )

        _1 = _1_container.component

        assert(isinstance(_1, Vertex))


    def test_get_next_triangle(self):

        pointA = glm.vec2(-1803732.875, 761263.25)
        pointB = pointA + glm.vec2(100,150)

        _1_container = self.promesheus.get_next_point(MeshComponentContainer(), pointA, pointB )
        _1 = _1_container.component
        assert(isinstance(_1, Vertex))

        _dir = glm.normalize(pointB - pointA)

        curr_face = self.promesheus.get_next_tri(_1, pointA, _dir )
        stu =  curr_face.bary_x_y(*pointB)
        
        #plot(self.file_path, self.promesheus, pointB)
        
        assert( all(elem > 0 for elem in stu) )
        
       

    # @unittest.skip("the reason is you")

    #
    # at vertex move into next triangle
    #
    def test_get_next_point_1(self):

        pointA = glm.vec2(-1803732.875, 761263.25)
        pointB = pointA + glm.vec2(10,0)

        container  = self.promesheus.get_next_point(MeshComponentContainer(), pointA, pointB )

        self.assertIsInstance(container.component, Vertex)

        container = self.promesheus.get_next_point( container , glm.vec3(container.point).xy , pointB )

        print(f"second take type: {type(container.component)}")
        self.assertIsNone(container.component)

    #
    # at vertex move into next triangle and end at vertex there
    #
    def test_get_next_point_2(self):

        pointA = glm.vec2(-1803732.875, 761263.25)
        pointB = glm.vec2(-1799191.375, 765804.8125) # vertex at 10,4

        container  = self.promesheus.get_next_point(MeshComponentContainer(), pointA, pointB )

        #self.assertIsInstance(
        
        self.assertIsInstance(container.component, Vertex)

        container = self.promesheus.get_next_point( container , container.point.xy , pointB )

        self.assertIsInstance(container.component, Edge)


        container = self.promesheus.get_next_point( container , container.point.xy , pointB )
        intersection = container.point

        self.assertIsNone(container.component)

        self.assertEqual(intersection.x , pointB[0] )
        self.assertEqual(intersection.y , pointB[1] )

    #
    # at vertex move along edge to point on that first edge
    #
    def test_get_next_point_3(self):
        pointA = glm.vec2(-1803732.875, 761263.25)
        pointB = glm.vec2(-1803732.875, 761263.25+2000)

        container = self.promesheus.get_next_point( MeshComponentContainer(), pointA , pointB )

        self.assertIsInstance(container.component,  Vertex)

        intersection = container.point
        container = self.promesheus.get_next_point( container , glm.vec3(intersection).xy , pointB )

        self.assertIsNone(container.component)

    #
    # at vertex move along edge to point on that first edge
    #
    def test_get_next_point_4(self):
        pointA = glm.vec2(-1803732.875, 761263.25)
        pointB = glm.vec2(-1803732.875, 761263.25+2000)

        container = self.promesheus.get_next_point( MeshComponentContainer(), pointA , pointB )

        self.assertIsInstance(container.component,  Vertex)
        intersection = container.point
        container = self.promesheus.get_next_point( container , glm.vec3(intersection).xy , pointB )

        print(f"LAST INTERSECTION: {intersection}")
        self.assertIsNone( container.component )

    #
    # at edge move along edge to point on same edge
    #
    def test_get_next_point_5(self):
        pointA = glm.vec2(-1803732.875, 761263.25+2000)
        pointB = glm.vec2(-1803732.875, 761263.25+4000)

        container = self.promesheus.get_next_point( MeshComponentContainer(), pointA , pointB )

        assert( isinstance(container.component, Edge) )

        intersection = container.point
        container = self.promesheus.get_next_point( container.component , glm.vec3(intersection).xy , pointB )

        print(f"LAST INTERSECTION: {intersection}")
        assert(container.component is None)

    #
    # start at vertex move slightly askew to edge
    #
    def test_get_next_point_5(self):
        pointA = glm.vec2(-1817944.125, 774887.875)
        pointB = glm.vec2(-1817357.375, 779429.375)

        container = self.promesheus.get_next_point( MeshComponentContainer(), pointA , pointB )

        self.assertIsInstance( container.component, Edge )

        container = self.promesheus.get_next_point(  container , glm.vec3(container.point).xy , pointB )

        self.assertIsInstance(container.component, Edge)

        container = self.promesheus.get_next_point( container , glm.vec3(container.point).xy , pointB )

        self.assertIsNone( container.component )


    #
    # test point inside raster -> point outside raster
    #
    def test_get_next_point_outside(self):
        pointA = glm.vec2(-1803732.875, 761263.25+200)
        pointB = glm.vec2(-1803732.875+10000, 761263.25+4000)

        container = self.promesheus.get_next_point( MeshComponentContainer(), pointA , pointB )

        self.assertIsInstance( container.component, Edge )

        container  = self.promesheus.get_next_point( container , glm.vec3(container.point).xy , pointB )

        self.assertIsInstance( container.component, Edge )
        self.assertIsInstance( container.point , glm.vec3 )

        container  = self.promesheus.get_next_point( container , glm.vec3(container.point).xy, pointB )

        self.assertIsInstance( container.component, Edge )
        self.assertIsInstance( container.point , glm.vec3 )

        container  = self.promesheus.get_next_point( container , glm.vec3(container.point).xy , pointB )

        self.assertIsInstance( container.component, Edge )
        self.assertIsInstance( container.point , glm.vec3 )

        container  = self.promesheus.get_next_point( container , glm.vec3(container.point).xy , pointB )

        self.assertIsInstance( container.component, Edge )
        self.assertIsInstance( container.point , glm.vec3 )

        container = self.promesheus.get_next_point( container , glm.vec3(container.point).xy , pointB )

        self.assertIsNone( container.component )
        self.assertIsNone( container.point )

    #
    # test point inside on raster border -> point outside raster
    #
    def test_get_next_point_border(self):
        pointA = glm.vec2(-1794649.75, 764914.8125)
        pointB = glm.vec2(-1803732.875+10000, 761263.25+4000)
        
        container = self.promesheus.get_next_point( MeshComponentContainer(), pointA , pointB )

        self.assertIsNone( container.component )
        self.assertIsNone( container.point )

    #
    # test point inside on raster border -> point outside raster
    #
    def test_get_next_point_border_2(self):
        pointA = glm.vec2(-1794649.75, 764914.8125)
        pointB = glm.vec2(-1793732.875, 765263.25)

        pixels = [(11, 4), (11, 5)]
        points = [self.promesheus.get_point(*elem) for elem in pixels]
        edge = Edge(points, pixels)

        container = MeshComponentContainer(None, edge, pointA)
        
        container = self.promesheus.get_next_point( container, pointA , pointB )

        self.assertIsNone( container.component )
        self.assertIsNone( container.point )

    #
    # 2 point single line string
    def test_get_right_number_intersections_1(self):

        points = [[(-1803732.875, 761463.25), (-1799191.375, 763189.0625)]]

        [conv_points] = self.promesheus.create_lines( points )

        self.assertEqual(len(conv_points), 3)

    #
    # 2 point single line string
    def test_get_right_number_intersections_2(self):

        points = [ [(-1803732.875, 761463.25), (-1799191.375, 763189.0625), (-1793732.875, 765263.25)] ]

        [conv_points] = self.promesheus.create_lines( points )

        self.assertEqual(len(conv_points), 5)


    #
    # 2 point single line string from vertex to another vertex
    # but they do collect some edges along the way
    #   
    def test_direct_geodesics(self):

        points = [ [(-1922879.25, 854802.3125), (-1922326.0, 854249.125)] ]


        pointA = glm.vec2(points[0][0])
        pointB = glm.vec2(points[0][1])

        # 0
        container = self.promesheus2.get_next_point( MeshComponentContainer(), pointA, pointB )

        self.assertIsInstance( container.component, Vertex )
        self.assertIsNotNone( container.face )

        # 1
        container = self.promesheus2.get_next_point( container, glm.vec3(container.point).xy , pointB )
        self.assertIsInstance( container.point, glm.vec3 )
        self.assertIsInstance( container.component, Vertex )
        
        # 2
        container =  self.promesheus2.get_next_point( container , glm.vec3(container.point).xy , pointB )
        self.assertIsInstance( container.point, glm.vec3 )
        self.assertIsInstance( container.component, Vertex )
        
        # 3
        container = self.promesheus2.get_next_point( container , glm.vec3(container.point).xy , pointB )
        self.assertIsInstance( container.point, glm.vec3 )
        
        self.assertIsInstance( container.component, Vertex )
        
        # 4
        # If the line is directly on vertex to vertex then why do we have 
        # to visit some edges
        container = self.promesheus2.get_next_point( container , glm.vec3(container.point).xy , pointB )
        self.assertIsInstance( container.point, glm.vec3 )      
        
        self.assertIsInstance( container.component, Edge )

    #
    # check to make sure we're not repeating any points!
    #
    def test_direct_geodesics_2(self):

        points = [ [(-1803732.875, 761463.25), (-1799191.375, 763189.0625), (-1793732.875, 765263.25)] ]

        [conv_points] = self.promesheus2.create_lines( points )

        for i in range(1,len(conv_points)):
            prev_x, prev_y = conv_points[i-1].xy
            curr_x, curr_y = conv_points[i].xy
        
            self.assertNotEqual( (prev_x, prev_y) , (curr_x, curr_y) )


    def test_left_right_vectors(self):

        vec = glm.vec2(0, 1)

        l_vec, r_vec = self.promesheus2.get_left_right_vectors(vec)

        self.assertEqual( l_vec, glm.vec2(-1,0) )
        self.assertEqual( r_vec, glm.vec2(1,0) )


    def test_get_index_tri(self):
        self.promesheus.create_index()


        for i , indices in enumerate( self.promesheus.test_index_array ):

            indx1, indx2, indx3 = indices

            col1,row1 = indx1 % self.promesheus.num_cols, indx1 // self.promesheus.num_cols
            col2,row2 = indx2 % self.promesheus.num_cols, indx2 // self.promesheus.num_cols
            col3,row3 = indx3 % self.promesheus.num_cols, indx3 // self.promesheus.num_cols

            sts = [(col1,row1), (col2,row2), (col3,row3)]
            t = Triangle([], sts)

            tri_index = self.promesheus.get_tri_index(t)

            self.assertEqual( i, tri_index )


if __name__ == "__main__":
    sys.stdout = sys.__stdout__
    unittest.main(argv=["first-arg-is-ignored"], exit=False)

