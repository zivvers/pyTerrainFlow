
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from osgeo import gdal
from config import load_config


def raster_extent(transform, width, height):
    # corners in pixel space
    corners = [
        (0, 0),
        (width, 0),
        (0, height),
        (width, height),
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

    return { "left": left , "right" : right, "bottom" : bottom, "top": top }


if __name__ == "__main__":

    cfg = load_config()

    TILE_SIZE = 1024

    raster_file = "dem_EPSG_26910_542354_4944208.tif"
    texture_file = "osip_EPSG_26910_542354_4944208.tif"

    raster_path = cfg.input_dir / raster_file
    texture_path = cfg.input_dir / texture_file

    
    dem_patch = 4
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

        band = band[:dem_patch, :dem_patch]
        h = src.height
        w = src.width

        rows = np.arange( dem_patch )
        cols = np.arange( dem_patch )

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

        tex_band = tex_band[:tex_patch-1, :tex_patch-1]

        tex_h = tex_src.height
        tex_w = tex_src.width

        tex_cols, tex_rows = np.meshgrid(np.arange(tex_patch), np.arange(tex_patch))
        tex_xs, tex_ys = rasterio.transform.xy(tex_src.transform, tex_rows, tex_cols, offset='ul')

        #extent = [tex_src.bounds.left, tex_src.bounds.right, tex_src.bounds.bottom, tex_src.bounds.top]
        #x_min, x_max, y_min, y_max = extent  
        
        #x_max = xs[-1];
        #y_min = ys[-1];


        for r in range(tex_patch):
            for c in range(tex_patch):

                tex_x_m, tex_y_m = tex_xs[c], tex_ys[r] #rasterio.transform.xy(tex_transf, r, c)
                #

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
                   # buffer_array[c, r, 2] = elev


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
        im = ax.imshow(tex_band, extent=[x_min, x_max, y_min, y_max], cmap='viridis', origin='upper')
    
        # ds = gdal.Open( str(raster_path) )
        # width = ds.RasterXSize
        # height = ds.RasterYSize
        # gt = ds.GetGeoTransform()

        # minx = gt[0]
        # maxy = gt[3]
        # maxx = gt[0] + width * gt[1] + height * gt[2]
        # miny = gt[3] + width * gt[4] + height * gt[5]
        
        dem_extent_d = raster_extent(src.transform, num_cols, num_rows)
        tex_extent_d = raster_extent(tex_src.transform, tex_patch, tex_patch)

        print(f"index array size: {index_array}")

        triang = mtri.Triangulation(xs, ys, index_array)
        tex_triang = mtri.Triangulation(tex_xs, tex_ys, tex_index_array)

        plt.triplot(triang, 'go-', label='DEM Mesh',color='black', linewidth=1)

        plt.triplot(tex_triang, 'g-', label='Texture Mesh',color='pink', linewidth=0.25)

        plt.show()