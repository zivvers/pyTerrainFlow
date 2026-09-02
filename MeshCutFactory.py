from MeshCut import MeshCut
import rasterio
import rasterio
import numpy as np
import matplotlib.tri as mtri
from affine import Affine


class MeshCutFactory:

    def get_mesh(self):
        return self.mesh

    def get_triang(self):
        return self.triang

    #
    # I guess there's an update method to update mesh?
    #
    def __init__(self, tif_file, dem_num_cells, resolution=10):
        self.mesh = None
        self.triang = None # Matplotlib styLE
        self.tif_file = None
        self.dem_num_cells = None
        self.resolution = None

        self.update(tif_file, dem_num_cells, resolution)

    def update(self, tif_file, dem_num_cells, resolution=10):
        self.mesh, self.triang = self.build_mesh(tif_file, dem_num_cells, resolution)

    def build_mesh(self, tif_file, dem_num_cells, resolution=10):

        # not parametizing for now
        offset_x,offset_y = 0, 0

        dem_patch_num_verts = dem_num_cells+1
                
        buffer_array = np.zeros((dem_patch_num_verts, dem_patch_num_verts, 3), dtype=np.float32)
        
        num_tri = (dem_patch_num_verts - 1) * (dem_patch_num_verts-1) * 2

        with rasterio.open(tif_file) as src:
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

            num_rows_inset = band.shape[0]
            num_cols_inset = band.shape[1]

            orig_transf = src.transform 
            
            new_transf =  orig_transf * Affine.translation(
                offset_x,   # column offset
                offset_y    # row offset
            )

            index_array = []

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

                    if (r > 0 and c > 0) and ( r < dem_patch_num_verts and c < dem_patch_num_verts):
                        currIndx = c + r * dem_patch_num_verts;
                        prevRowSameColIndx = c + (r - 1) * dem_patch_num_verts;
                        prevRowBackColIndx = c - 1 + (r - 1) * dem_patch_num_verts;
                        prevIndx = currIndx - 1;
                        index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])
                        index_array.append([prevRowBackColIndx, prevIndx, currIndx])

                rows = np.arange( dem_patch_num_verts ) + offset_y
                cols = np.arange( dem_patch_num_verts ) + offset_x

                cols, rows = np.meshgrid(cols, rows)
                xs, ys = rasterio.transform.xy(src.transform, rows, cols, offset='ul')


        triang = mtri.Triangulation(
            xs.ravel(),
            ys.ravel(),
            np.asarray(index_array)
        )

        return MeshCut(buffer_array, index_array, new_transf, bounds), triang
