
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from osgeo import gdal
from config import load_config


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


def point_inside_extent(x, y, extent):
    return (
        extent["left"] <= x <= extent["right"] and
        extent["bottom"] <= y <= extent["top"]
    )


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

    
    tex_row_offset = offset_y * 10#scale
    tex_col_offset = offset_x * 10;####scale
    dem_patch = 6
    num_rows = num_cols = dem_patch
    tex_patch = (10 * (dem_patch - 1)) + 1




    #num_tex_rows = num_tex_cols = 1*10

    buffer_array = np.zeros((num_rows, num_cols, 3), dtype=np.float32)

    num_tri = (num_cols - 1) * (num_rows-1) * 2

    num_tex_tri = (tex_patch - 1) * (tex_patch-1) * 2

    index_array = []

    with rasterio.open(raster_path) as src:
        band = src.read(1)
        transf = src.transform

        band = band[offset_y:offset_y+dem_patch,offset_x:offset_x+dem_patch]
        h = src.height
        w = src.width

        rows = np.arange( dem_patch ) + offset_y
        cols = np.arange( dem_patch ) + offset_x

        cols, rows = np.meshgrid(cols, rows)
        xs, ys = rasterio.transform.xy(src.transform, rows, cols,offset='ul')

        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
        x_min, x_max, y_min, y_max = extent  
        
        x_max = xs[-1];
        y_min = ys[-1];

    tex_buffer_array = np.zeros((tex_patch, tex_patch, 3), dtype=np.float32)
    tex_index_array = []


    with rasterio.open(texture_path) as tex_src:

        tex_band = tex_src.read(1)
        tex_transf = tex_src.transform
        print(f"band shape: {tex_band.shape}")

        print

        tex_band = tex_band[tex_row_offset : tex_col_offset + tex_patch - 1,
            tex_col_offset : tex_col_offset + tex_patch - 1 ]

        tex_h = tex_src.height
        tex_w = tex_src.width

        _rows = np.arange(tex_patch) + tex_row_offset
        _cols = np.arange(tex_patch) + tex_col_offset

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


        for r in range(tex_patch):
            for c in range(tex_patch):

                #tex_x_m, tex_y_m = tex_xs[c], tex_ys[r] #rasterio.transform.xy(tex_transf, r, c)
                #
                tex_x_m = tex_xs[c]
                tex_y_m = tex_ys[r]

                #tex_buffer_array[c, r, 0] = tex_x_m #2ND
                #tex_buffer_array[c, r, 1] = tex_y_m
                #tex_buffer_array[c, r, 2] = elev


                if ( r < num_rows and c < num_cols):

                    #elev = band[r, c]
                    x_m, y_m = xs[c], ys[r] #rasterio.transform.xy(transf, r, c)

                    if (r ==0 and c ==0):
                        print(f"tex ({tex_x_m}, {tex_y_m }), dem: ({x_m}, {y_m })")


                    #buffer_array[c, r, 0] = x_m #2ND
                    #buffer_array[c, r, 1] = y_m
                   #buffer_array[c, r, 2] = elev


                if (r > 0 and c > 0):
                    currIndx = c + r * tex_patch;
                    prevRowSameColIndx = c + (r - 1) * tex_patch;
                    prevRowBackColIndx = c - 1 + (r - 1) * tex_patch;
                    prevIndx = currIndx - 1;
                    tex_index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])
                    tex_index_array.append([prevRowBackColIndx, prevIndx, currIndx])


                    if ( r < num_rows and c < num_cols):
                        currIndx = c + r * num_cols;
                        prevRowSameColIndx = c + (r - 1) * num_cols;
                        prevRowBackColIndx = c - 1 + (r - 1) * num_cols;
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
        
        dem_extent_d = raster_extent(src.transform, num_cols, num_rows, offset_x, offset_y)
        tex_extent_d = raster_extent(tex_src.transform, tex_patch, tex_patch, offset_x*10, offset_y*10)

        im = ax.imshow(tex_band, extent=[tex_extent_d["left"], tex_extent_d["right"], tex_extent_d["bottom"], tex_extent_d["top"]], cmap='viridis', origin='upper')

        print(f"index array size: {index_array}")

        tex_X, tex_Y = np.meshgrid(tex_xs, tex_ys)
        #triang = mtri.Triangulation(xs, ys, index_array)
        tex_triang = mtri.Triangulation(
            tex_X.ravel(),
            tex_Y.ravel(),
            np.asarray(tex_index_array)
        )
        triang = mtri.Triangulation(
            xs.ravel(),
            ys.ravel(),
            np.asarray(index_array)
        )
        plt.triplot(triang, 'go-', label='DEM Mesh',color='black', linewidth=1)
        plt.triplot(tex_triang, 'g-', label='Texture Mesh',color='pink', linewidth=0.25)

        #px = [p[0] for p in plot_points]
        #py = [p[1] for p in plot_points]

        #[ [542764.00000, 4943808.00000], [542774.0000, 4943798.0000] ,  [542774.00000, 4943808.00000] )


        #new_points = [[542774.00000, 4943808.00000], [542774.0000, 4943798.0000]]
        #new_points = [ [542774.00000, 4943808.00000], [542774.0000, 4943798.0000], [542784.00000, 4943798.00000] ]
        #new_points = [[543224,4943378], [543234,4943380.5]]
        #new_points = [[543224.0, 4943378.0], [543232.0, 4943380.0], [543232.125	, 4943380], [543234.0, 4943380.5]]
        new_points = [ [543358.875,4943464.5], [543359.4375,4943468], [543360.0625,4943472] ]
        #new_points = [[543358.875,4943464.5], [543361.25,4943468], [543362.375,4943469.5]]
        px2 = [p[0] for p in new_points]
        py2 = [p[1] for p in new_points]
        ax.scatter(px2, py2, color="green", s=300, zorder=1)


        #ax.scatter(px, py, color="red", s=80, zorder=10)
        #ax.scatter([543224.00000], [4943378.00000],  color="magenta", s=500, zorder=1)
        plt.show()