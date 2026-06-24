    
from config import load_config

import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from pyglm import glm

import math

from MeshComponentCreator import MeshComponentCreator



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

    def __init__(self, buffer, tf, bnds):

        self.buffer_array = buffer
        self.transform = tf
        self.bounds = bnds
        self.num_cols = self.buffer_array.shape[0]
        self.num_rows = self.buffer_array.shape[1]

        self.num_triangles = (self.num_cols - 1)*(self.num_rows - 1)*2       

    def cross_2D(self, a, b ):
        return (a.x * b.y) - (a.y * b.x)

    def convert_to_point(self, _pix_x, _pix_y):

        point_y = _pix_y*self.transform.e + self.transform.f ;
        pixel_x = _pix_x * self.transform.a + self.transform.c ;
        
        return glm.vec2(pixel_x, point_y)

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

    # for a given pixel get all edges
    def get_edges(self, pix):
        
        neighbor_indices = self.get_neighbors(*pix)






if __name__ == "__main__":
    cfg = load_config()

    TILE_SIZE = 1024

    raster_file = "dem_EPSG_26910_542354_4944208.tif"
    texture_file = "osip_EPSG_26910_542354_4944208.tif"

    plot_points = [[542781.25,4943800.5], [542774,4943798]]


    raster_path = cfg.input_dir / raster_file
    texture_path = cfg.input_dir / texture_file

    offset_x, offset_y = 100, 72
    tex_offset_x , tex_offset_y = offset_x*10, offset_y*10

    dem_resolution = 10 #meters per pixel
    tex_resolution = 1 #meters per pixel
    
    #
    # edit this variable!
    dem_num_cells = 3

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

    with rasterio.open(raster_path) as src:
        band = src.read(1)
        transf = src.transform

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

        num_rows = band.shape[0]
        num_cols = band.shape[1]
    
        transf = src.transform

        for r in range(num_rows):
            for c in range(num_cols):

                # should be in meters
                x_m, y_m = rasterio.transform.xy(transf, r, c, offset='ul')
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

        triang = mtri.Triangulation(
            xs.ravel(),
            ys.ravel(),
            np.asarray(index_array)
        )

        plt.triplot(triang, 'go-', label='DEM Mesh',color='black', linewidth=1)
        plt.triplot(tex_triang, 'g-', label='Texture Mesh',color='pink', linewidth=0.25)

        bounds = [ buffer_array[0,0,0].item(), buffer_array[0, num_rows-1, 1].item(), 
                    buffer_array[num_cols-1,0 , 0].item(), buffer_array[0, 0, 1].item() ]

        test_points = [ [542354.0, 4944198.0],
                        [542356.5, 4944200.9],
                        [542359.0, 4944203.6],
                        [542361.5, 4944205.9],
                        [542364.0, 4944208.0] ] 
        
        for i,p in enumerate(test_points):
            x,y=p
            if (bounds[0] <= x <= bounds[2]) and (bounds[1] <= y <= bounds[3]):
                print(f"point {i} within bounds!")
        
        test_index = np.array( [[0, 4]] , dtype=int )

        promesheus = MeshComponentCreator(buffer_array, transf, bounds)

        # road_lines_comp_3D = promesheus.create_lines( test_points, test_index )

        # test_points = [[543357, 4943473.0],
        #                [543367, 4943463.0],
        #                [543380.0, 4943472.0],
        #                  ]
        
        test_points = [ glm.vec2(543369, 4943477.0), glm.vec2(543378.0, 4943472.0)
                       , glm.vec2(543367, 4943463.0), glm.vec2(543357, 4943473.0) ]


        mesh = MeshCut(buffer_array, transf, bounds)

        test_points_ccw = mesh.order_ccw(test_points)

        inner_points = mesh.get_inner_verts(test_points);

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