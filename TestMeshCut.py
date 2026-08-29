from CutMesh import MeshCut
from config import load_config
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys as sys
from Intersection import *

from shapely.geometry import Polygon

from affine import Affine
import unittest, string


def add_letter(inters):
    for i,_inter in enumerate(inters):
        _inter.letter = string.ascii_uppercase[i]

def plt_triang_base(triang , mesh):

    plt.triplot(triang, 'go-', label='DEM Mesh', color='black', linewidth=1)

    # EDGE LABELS
    for edge_i, e in enumerate(mesh.tri_edge_index):
        e1,e2 = mesh.get_edge_points(edge_i)
        e_avg = (e1+e2)/2
        plt.text(e_avg.x, e_avg.y, f"{edge_i}",\
                    fontsize=11, style="italic", zorder=20)

    for tri_index, tri in enumerate(mesh.index_array):
        tot_x = 0.0
        tot_y = 0.0

        for vert_index in tri:
            #_x,_y = vert_index % dem_patch_num_verts,vert_index // dem_patch_num_verts
            tot_x += mesh.get_point(vert_index)[0]
            tot_y += mesh.get_point(vert_index)[1]

        plt.text(tot_x/3, tot_y/3, f"{tri_index}"\
                 , fontsize=12, weight="bold",  zorder=20)
    
    return plt

#
# format our intersection to show on the
# plot
def format_inter(inter: Intersection) -> str:

    mesh_label = inter.mesh_edge_index if inter.mesh_feature == FeatureType.EDGE else inter.mesh_vertex_index

    return rf'''From Graticule ${inter.graticule_type.name} between pixels {inter.grat_pix0}, {inter.grat_pix1}: 
            Intersection ${inter.letter}$: mesh {inter.mesh_feature.name} id {mesh_label} $t_m$: {inter.mesh_t : .3f} , clip {inter.clip_feature.name} id {inter.clip_poly_id[1]} $t_c$: {inter.clip_t : .3f}'''
    

'''
Intersection(point=(543264.0, 4943493.0), mesh_feature=<FeatureType.EDGE: 2>, clip_feature=<FeatureType.VERTEX: 1>
    , mesh_edge_index=7, mesh_edge_order=None, mesh_vertex_index=None
    , graticule_type=<Graticule.LON: 1>, clip_poly_id=None, clip_poly_order=None
    , mesh_t=0.5, clip_t=0.0), 
'''

#
# the PNTs have to be 2D
#
def plot_points(triang , mesh, pnts, inters, plot_path):

    fig, ax = plt.subplots()
    p = plt_triang_base(triang, mesh )

    #
    # clip points
    x,y = zip( *pnts )
    p.plot(x, y, zorder=10, color="green")

    frmt_inters = [format_inter(_inter) for _inter in inters]

    inter_points = [ _inter.point for _inter in inters ]

    inter_x,inter_y = zip( *inter_points )

    if len(inter_points) >= 5:
        bottom_perc = 0.45
    elif len(inter_points) >= 3:
        bottom_perc = 0.35
    else:
        bottom_perc = 0.3
    #bottom_perc = 0.35 if len(inter_points) >= 3 else 0.3
    
    bottom_start = 0.05 if len(inter_points) > 1 else 0.125

    p.plot(inter_x,inter_y, 'ro', zorder=12)

    p.figtext(
        0.45,                                       # x-coordinate (centered)
        bottom_start,                                      # y-coordinate (near the absolute bottom)
        '\n'.join(frmt_inters),  # LaTeX text (must use 'r' for raw string)
        ha="center",
        color="#1f42b4",
        fontsize=10
    )

    p.tight_layout()
    p.subplots_adjust(bottom=bottom_perc)

    if plot_path is None:
        p.show()
    else:
        p.savefig(plot_path, dpi=300, bbox_inches='tight')

#
#
#
def plot_points_arrows( triang , mesh, clip_pnts, dir, plot_path ):

    fig, ax = plt.subplots()
    p = plt_triang_base(triang, mesh )

    #
    # clip points
    x,y = zip( *clip_pnts )
    p.plot(x, y, zorder=10, color="green")

    ax.annotate("", 
            xy=dir[1],      # Arrow tip destination
            xytext=dir[0],  # Arrow tail starting point
            arrowprops=dict(arrowstyle="->", color="red", lw=2))

    if plot_path is None:
        p.show()
    else:
        p.savefig(plot_path, dpi=300, bbox_inches='tight')


def plot_clip_inters( triang , mesh, clip_verts, inters, remaining_clip, plot_path ):

    CLIPPED_VERT_COLOR = "#D55E00" 
    COL_EDGE_INTER_COLOR = "#0072B2"
    COL_CLIP_ONLY_COLOR = "#CC79A7"

    fig, ax = plt.subplots()
    p = plt_triang_base(triang, mesh )

    x,y = zip( *clip_verts )
    p.plot(x, y, "o", zorder=10, color=CLIPPED_VERT_COLOR, label="Clipped Mesh Verts")

    x,y = zip( *inters )
    p.plot(x, y, "o", zorder=10, color=COL_EDGE_INTER_COLOR, label="Mesh Edge Intersections")

    x,y = zip( *remaining_clip )
    p.plot(x, y, "o", zorder=10, color=COL_CLIP_ONLY_COLOR, label="Remaining Clip Verts")

    p.legend( loc="lower center" )
    p.savefig(plot_path, dpi=300, bbox_inches='tight')


class TestMeshCut(unittest.TestCase):

    __test__ = True
    
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config()

        raster_file = "roof_raster_fake3.tif"

        raster_path = cls.cfg.input_dir / raster_file

        #
        # edit these 2 variables
        offset_x, offset_y = 0,0 
        dem_num_cells = 2

        dem_resolution = 10

        dem_patch_size = dem_resolution * dem_num_cells

        dem_patch_num_verts = dem_num_cells+1
        
        buffer_array = np.zeros((dem_patch_num_verts, dem_patch_num_verts, 3), dtype=np.float32)

        num_tri = (dem_patch_num_verts - 1) * (dem_patch_num_verts-1) * 2

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


            cls.mesh = MeshCut(buffer_array, index_array, new_transf, bounds)

            rows = np.arange( dem_patch_num_verts ) + offset_y
            cols = np.arange( dem_patch_num_verts ) + offset_x

            cols, rows = np.meshgrid(cols, rows)
            xs, ys = rasterio.transform.xy(src.transform, rows, cols, offset='ul')

            # 9 PM
            firstCoord = (buffer_array[1,0, 0].item(), ((buffer_array[1,1, 1] + buffer_array[1,0, 1])/2).item() )
            # noon
            secondCoord = (((buffer_array[1,0,0]+buffer_array[2,0, 0])/2).item(), buffer_array[1,0,1].item() )
            # 3 PM
            thirdCoord = (buffer_array[2,0,0].item() , ((buffer_array[1,1, 1] + buffer_array[1,0, 1])/2).item()  )

            fourthCoord = (((buffer_array[1,0, 0]+buffer_array[2,0, 0])/2).item(), buffer_array[1,1, 1].item() )

            cls.diamond_clip = [firstCoord , secondCoord, thirdCoord, fourthCoord ]
            cls.square_clip = [ cls.mesh.get_point(1,0)[:2], cls.mesh.get_point(2,0)[:2],
                                cls.mesh.get_point(2,1)[:2], cls.mesh.get_point(1,1)[:2] ]

            cls.bigger_rect = [ cls.mesh.interpolate_xy_coords((0.5,0)), cls.mesh.interpolate_xy_coords((2,0)) ,\
                              cls.mesh.interpolate_xy_coords((2,1)), cls.mesh.interpolate_xy_coords((0.5,1)) ]

            cls.bigger_rect2 = [ [cls.mesh.interpolate_xy_coords((0.5, 0))[:2], cls.mesh.interpolate_xy_coords((2, 0))[:2] ,\
                              cls.mesh.interpolate_xy_coords((2, 1.25))[:2], cls.mesh.interpolate_xy_coords((0.5, 1.25))[:2] ] ]


            cls.triang = mtri.Triangulation(
                xs.ravel(),
                ys.ravel(),
                np.asarray(index_array)
            )


    def get_point_categories(self , poly_clip ):

        clp_verts = []
        for vert, clip_bool in self.mesh.vert_clip_index.items():
            if clip_bool:
                clp_verts.append( self.mesh.get_point(vert)[:2] )

        inters = []
        for edge_indx, edge_dict in self.mesh.all_tri_edges.items():

            for _inter in edge_dict:
                inters.append( _inter.point[:2] ) 

        remaining_clip_verts = list( range( len( poly_clip[0] ) ) )

        for clip_edge, clip_dict in self.mesh.all_clip_edges.items():

            for _inter in clip_dict:
                clip_poly_edge_id = _inter.clip_poly_id[1]

                if _inter.clip_feature == FeatureType.VERTEX:
                    if clip_poly_edge_id in remaining_clip_verts:
                        remaining_clip_verts.remove(clip_poly_edge_id)

        remaining_clip = []
        for i in remaining_clip_verts: # make into points
           remaining_clip.append( poly_clip[0][i][:2] ) 

        return clp_verts, inters, remaining_clip

    #
    # returns [Intersection] ,
    #           [tup2D, tup2D]
    #
    def fetch_intersections_clip_i( self, clip_points, i ):

        curr_i = i % 4 # ASSUMING 4 clip points
        next_i = (i+1) % 4
        clip_p0 = clip_points[curr_i]
        clip_p1 = clip_points[next_i]

        clip_edge_inters = self.mesh.clip_mesh_edges(clip_p0 , clip_p1, (0,i))         
        return clip_edge_inters, [ clip_p0 , clip_p1 ]

    def test_interpolate_xy(self):

        interp_coords = self.mesh.interpolate_xy_coords((0.5,0))

        value = (543254.0+5, 4943498.0)

        self.assertAlmostEqual(interp_coords, value)


    def test_create_graticule_plots(self):

        fig, ax = plt.subplots()

        p = plt_triang_base(self.triang 
                    , self.mesh)

        #
        # diag
        diag1 = [self.mesh.get_point(1,0)[:2], self.mesh.get_point(2,1)[:2]]
        diag2 = [self.mesh.get_point(0,0)[:2], self.mesh.get_point(2,2)[:2]]
        diag3 = [self.mesh.get_point(0,1)[:2], self.mesh.get_point(1,2)[:2]]

        print(f"diags: {diag1}, {diag2}, {diag3}")

        x,y = zip(*diag1)
        p.plot(x, y, color='pink', linewidth=2)
        x,y = zip(*diag2)
        p.plot(x, y, color='pink', linewidth=2) 
        x,y = zip(*diag3)
        p.plot(x, y, color='pink', linewidth=2) 

        file_path =  self.cfg.output_dir / f"diagonal_graticule.png"
        p.savefig(file_path, dpi=300, bbox_inches='tight') 
        p.close()
        plt.close('all')
        p = plt_triang_base(self.triang 
                    , self.mesh)
        
        hor1 = [self.mesh.get_point(0,0)[:2], self.mesh.get_point(2,0)[:2]]
        hor2 = [self.mesh.get_point(0,1)[:2], self.mesh.get_point(2,1)[:2]]
        hor3 = [self.mesh.get_point(0,2)[:2], self.mesh.get_point(2,2)[:2]]

        x,y = zip(*hor1)
        p.plot(x, y, color='pink', linewidth=2)
        x,y = zip(*hor2)
        p.plot(x, y, color='pink', linewidth=2) 
        x,y = zip(*hor3)
        p.plot(x, y, color='pink', linewidth=2) 

        file_path =  self.cfg.output_dir / f"latitude_graticule.png"
        p.savefig(file_path, dpi=300, bbox_inches='tight')
        p.close()
        plt.close('all')
        p = plt_triang_base(self.triang 
                    , self.mesh)
        
        ver1 = [self.mesh.get_point(0,0)[:2], self.mesh.get_point(0,2)[:2]]
        ver2 = [self.mesh.get_point(1,0)[:2], self.mesh.get_point(1,2)[:2]]
        ver3 = [self.mesh.get_point(2,0)[:2], self.mesh.get_point(2,2)[:2]]

        x,y = zip(*ver1)
        p.plot(x, y, color='pink', linewidth=2)
        x,y = zip(*ver2)
        p.plot(x, y, color='pink', linewidth=2) 
        x,y = zip(*ver3)
        p.plot(x, y, color='pink', linewidth=2) 

        file_path =  self.cfg.output_dir / f"longtiude_graticule.png"
        p.savefig(file_path, dpi=300, bbox_inches='tight')
        p.close()
        plt.close('all')

    #
    #
    #
    def test_get_edges_vert( self ):

        req_grats = [Graticule.LAT, Graticule.DIAG, Graticule.LON, Graticule.LAT, Graticule.DIAG, Graticule.LON]

        i = 4
        edges = self.mesh.get_edges_vert(i)

        for e_i, edge in enumerate(edges):
            if edge is not None:
                self.assertEqual( edge[1], req_grats[e_i], f"edge {i} is wrong!" ) 


    #
    #
    # next tests make sure from a point can get edge and
    # know which point in edge is current vert
    def test_get_edges_vert_interior( self ):
        # interior
        vert_index = 4
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    def test_get_edges_vert_top_left( self ):

        vert_index = 0
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    def test_get_edges_vert_top_right( self ):

        vert_index = 2
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    def test_get_edges_vert_bottom_right( self ):

        vert_index = 8
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    def test_get_edges_vert_bottom_left( self ):
        vert_index = 6
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    #
    # bottom row, not far right/left
    def test_get_edges_vert_bottom( self ):
        # interior
        vert_index = 7
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    # top row, not far right/left
    def test_get_edges_vert_top( self ):
        # interior
        vert_index = 2
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    # right col, not top/bottom
    def test_get_edges_vert_right( self ):
        # interior
        vert_index = 5
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)

    # left col, not top/bottom
    def test_get_edges_vert_left( self ):
        # interior
        vert_index = 4
        edges = self.mesh.get_edges_vert(vert_index)

        for e_i, edge_tup in enumerate(edges):
            if edge_tup is not None:
                i = 0 if edge_tup[2] else 1
                self.assertEqual(self.mesh.tri_edge_index[ edge_tup[0] ][i], vert_index)


    def test_glancing(self):

        pix1 = ( 1 ,  0 )
        pix2 = ( 1 , self.mesh.num_rows-1 )
    
        _e0 = self.mesh.convert_to_point(*pix1);
        _e1 = self.mesh.convert_to_point(*pix2);
            # starting on the mesh edge
        glancing_inter1 = self.mesh.get_intersection(Graticule.LON, self.diamond_clip[0], self.diamond_clip[1], _e0, _e1)

        assert glancing_inter1 is not None 

    def test_diamond_clip0(self):

        i = 0
        clip_edge_inters,clip_edge = self.fetch_intersections_clip_i( self.diamond_clip, i )      

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        self.assertEqual(clip_edge_inters[0].mesh_edge_index , 7)
        self.assertEqual(clip_edge_inters[1].mesh_edge_index , 13)  

        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"diamond_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    def test_diamond_clip1(self):

        i = 1
        clip_edge_inters , clip_edge = self.fetch_intersections_clip_i( self.diamond_clip, i )        

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        self.assertEqual( clip_edge_inters[0].mesh_edge_index , 1 )

        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"diamond_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    def test_diamond_clip2(self):

        i = 2
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.diamond_clip, i )        

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        self.assertEqual( clip_edge_inters[0].mesh_edge_index , 8 ) 
        self.assertEqual( clip_edge_inters[1].mesh_edge_index , 13 )
        add_letter(clip_edge_inters)

        file_path = self.cfg.output_dir / f"diamond_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    def test_diamond_clip3(self):

        i = 3

        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.diamond_clip, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        self.assertEqual( clip_edge_inters[0].mesh_edge_index , 3)

        add_letter(clip_edge_inters)

        file_path = self.cfg.output_dir / f"diamond_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    #
    # square clip time
    #
    def test_square_clip0(self):
        i = 0

        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.square_clip, i )        

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        self.assertEqual( clip_edge_inters[0].mesh_vertex_index , 1)

        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"square_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    def test_square_clip1(self):
        i = 1
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.square_clip, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        self.assertEqual( clip_edge_inters[0].mesh_vertex_index , 2 )

        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"square_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path )

    def test_square_clip2(self):
        i = 2
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.square_clip, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        
        self.assertEqual(clip_edge_inters[0].mesh_vertex_index, 5 )

        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"square_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path )

    def test_square_clip3(self):

        i = 3
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.square_clip, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        
        self.assertEqual(clip_edge_inters[0].mesh_vertex_index, 4 )
        
        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"square_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    #
    #
    #
    def test_bigger_clip0(self):
        i = 0
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.bigger_rect, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        
        #self.assertEqual(clip_edge_inters[0].mesh_vertex_index, 4 )
        self.assertTrue( any(clp.mesh_edge_index == 0  for clp in clip_edge_inters) )
        self.assertTrue( any(clp.mesh_vertex_index == 1  for clp in clip_edge_inters) )
        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"bigger_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    def test_bigger_clip1(self):
        i = 1
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.bigger_rect, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        
        #self.assertEqual(clip_edge_inters[0].mesh_vertex_index, 4 )

        self.assertEqual(  len(clip_edge_inters) , 2 )
        self.assertTrue( all(clp.mesh_vertex_index == 2  for clp in clip_edge_inters) )
        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"bigger_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 

    def test_bigger_clip2(self):
        i = 2
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.bigger_rect, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        
        #self.assertEqual(clip_edge_inters[0].mesh_vertex_index, 4 )

        self.assertEqual(  len(clip_edge_inters) , 5 )
        self.assertTrue( any(clp.mesh_vertex_index == 5  for clp in clip_edge_inters) )
        self.assertTrue( any(clp.mesh_vertex_index == 4  for clp in clip_edge_inters) )

        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"bigger_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path ) 


    def test_bigger_clip3(self):

        i = 3
        clip_edge_inters, clip_edge = self.fetch_intersections_clip_i( self.bigger_rect, i )          

        print(f"test {i}, # intersections {len(clip_edge_inters)}")

        print(clip_edge_inters)
        #self.assertEqual(clip_edge_inters[0].mesh_vertex_index, 4 )

        self.assertEqual(  len(clip_edge_inters) , 2 )
        self.assertTrue( any(clp.mesh_edge_index == 2  for clp in clip_edge_inters) )
        self.assertTrue( any(clp.mesh_edge_index == 12 for clp in clip_edge_inters) )
        
        add_letter(clip_edge_inters)

        file_path =  self.cfg.output_dir / f"bigger_clip{i}.png"
        
        plot_points( self.triang 
                    , self.mesh
                    , clip_edge
                    , clip_edge_inters
                    , file_path )

    #
    # plot arrows coming from points so we can show 
    # clipped vertex process
    #
    def test_central_vert_inside(self):

        vert_index = 4

        self.mesh.perform_clipping( self.bigger_rect2 )

        in_clip_4_right = self.mesh.follow_geodesic_clip( True 
                        , Graticule.LAT
                        , vert_index
                        , self.bigger_rect2 )

        self.assertTrue( in_clip_4_right )

        #
        # only use for plotting!
        full_poly_clip = self.bigger_rect2[0] + [ self.bigger_rect2[0][0] ]

        #
        # LAT right
        arrow_dir = [ self.mesh.interpolate_xy_coords((1,1))[:2], \
                      self.mesh.interpolate_xy_coords((1.8,1))[:2] ]

        file_path =  self.cfg.output_dir / f"inside_clip_test_vert{vert_index}_lat_right.png"

        plot_points_arrows( self.triang , self.mesh, full_poly_clip, \
                            arrow_dir, file_path )

        # LAT left
        in_clip_4_left = self.mesh.follow_geodesic_clip( False 
                        , Graticule.LAT
                        , vert_index
                        , self.bigger_rect2 )

        self.assertTrue( in_clip_4_left )

        arrow_dir = [ self.mesh.interpolate_xy_coords((1,1))[:2], \
                      self.mesh.interpolate_xy_coords((0.6,1))[:2] ]

        file_path =  self.cfg.output_dir / f"inside_clip_test_vert{vert_index}_lat_left.png"

        plot_points_arrows( self.triang , self.mesh, full_poly_clip, \
                            arrow_dir, file_path )

        # LON up
        in_clip_4_up = self.mesh.follow_geodesic_clip( False 
                        , Graticule.LON
                        , vert_index
                        , self.bigger_rect2 )

        self.assertTrue( in_clip_4_up )

        arrow_dir = [ self.mesh.interpolate_xy_coords((1,1))[:2], \
                      self.mesh.interpolate_xy_coords((1,0.2))[:2] ]

        file_path =  self.cfg.output_dir / f"inside_clip_test_vert{vert_index}_lon_up.png"

        plot_points_arrows( self.triang , self.mesh, full_poly_clip, \
                            arrow_dir, file_path )

        # LON down
        in_clip_4_down = self.mesh.follow_geodesic_clip( True 
                        , Graticule.LON
                        , vert_index
                        , self.bigger_rect2 )

        self.assertTrue( in_clip_4_down )

        arrow_dir = [ self.mesh.interpolate_xy_coords((1,1))[:2], \
                      self.mesh.interpolate_xy_coords((1,1.20))[:2] ]

        file_path =  self.cfg.output_dir / f"inside_clip_test_vert{vert_index}_lon_down.png"

        plot_points_arrows( self.triang , self.mesh, full_poly_clip, \
                            arrow_dir, file_path )

        # DIAG up
        in_clip_4_diag_up = self.mesh.follow_geodesic_clip( False 
                        , Graticule.DIAG
                        , vert_index
                        , self.bigger_rect2 )

        self.assertTrue( in_clip_4_diag_up )

        arrow_dir = [ self.mesh.interpolate_xy_coords((1,1))[:2], \
                      self.mesh.interpolate_xy_coords((0.6,0.6))[:2] ]

        file_path =  self.cfg.output_dir / f"inside_clip_test_vert{vert_index}_diag_up.png"

        plot_points_arrows( self.triang , self.mesh, full_poly_clip, \
                            arrow_dir, file_path )

        # DIAG down
        in_clip_4_diag_down = self.mesh.follow_geodesic_clip( False 
                        , Graticule.DIAG
                        , vert_index
                        , self.bigger_rect2 )

        self.assertTrue( in_clip_4_diag_down )

        arrow_dir = [ self.mesh.interpolate_xy_coords((1,1))[:2], \
                      self.mesh.interpolate_xy_coords((1.2,1.2))[:2] ]

        file_path =  self.cfg.output_dir / f"inside_clip_test_vert{vert_index}_diag_down.png"

        plot_points_arrows( self.triang , self.mesh, full_poly_clip, \
                            arrow_dir, file_path )

        self.mesh.find_clipped_verts( self.bigger_rect2 )

        self.assertTrue(self.mesh.vert_clip_index[4])
        self.assertTrue(self.mesh.vert_clip_index[1])
        self.assertTrue(self.mesh.vert_clip_index[2])
        self.assertTrue(self.mesh.vert_clip_index[5])

        self.assertFalse(self.mesh.vert_clip_index[0])
        self.assertFalse(self.mesh.vert_clip_index[3])
        self.assertFalse(self.mesh.vert_clip_index[6])
        self.assertFalse(self.mesh.vert_clip_index[7])
        self.assertFalse(self.mesh.vert_clip_index[8])

        plot_path = self.cfg.output_dir / "vertex_category.png"

        clp_verts, inters, remaining_clip = \
            self.get_point_categories( self.bigger_rect2 )

        print(f"remaining clip: {remaining_clip}")

        #
        # should be 1 "remaining" clip point
        self.assertTrue( len(remaining_clip) > 0 )

        plot_clip_inters( self.triang , self.mesh , clp_verts \
                         , inters, remaining_clip, plot_path )

if __name__ == "__main__":
    sys.stdout = sys.__stdout__
    unittest.main(argv=["first-arg-is-ignored"], exit=False)