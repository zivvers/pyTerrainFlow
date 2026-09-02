from MeshCut import MeshCut
from MeshCutFactory import *
from pyglm import glm
from config import load_config
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys as sys
from Intersection import *
import math
from shapely.geometry import Polygon

from affine import Affine
import unittest, string
from dataclasses import dataclass

from matplotlib.animation import FuncAnimation

from adjustText import adjust_text


@dataclass(frozen=True)
class Line:
    p0: tuple[float, float]
    p1: tuple[float, float]
    pix0: tuple[int, int]
    pix1: tuple[int, int]


def plot_points_anim(triang , mesh, pnts, inters, clip_points , plot_path, edge_limiter=None, tri_bool=True, extent=None):

    #p = plt_triang_clip(triang, mesh, clip_points+[clip_points[0]] )
    #fig = plt.figure(figsize=(8, 10), dpi=120)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=120)

    # gs = fig.add_gridspec(
    #     2,
    #     1,
    #     height_ratios=[3, 1],
    #     hspace=0.08,
    # )

    ax.triplot(
        triang,
        color="0.65",
        linewidth=0.8,
        linestyle="-",
        marker="o",
        markersize=2.5,
        label="Tri Mesh",
        zorder=1,
    )

    # for edge_i, e in enumerate(mesh.tri_edge_index):
    #     e1,e2 = mesh.get_edge_points(edge_i)
    #     e_avg = (e1+e2)/2
    #     ax.text(e_avg.x, e_avg.y, f"{edge_i}",\
    #                 fontsize=11, style="italic", zorder=20)

    all_clip_points = clip_points+[clip_points[0]]
    clip_xs,clip_ys = zip(*all_clip_points)

    ax.plot(
        clip_xs,
        clip_ys,
        color="0.15",
        linewidth=1.5,
        linestyle="--",
        zorder=2,
        label="Clip Poly",
    )

    #
    # clip points
    x,y = zip( *pnts )
    ax.plot(x, y, zorder=10, color="green")

    frmt_inters = [format_inter(_inter) for _inter in inters]

    inter_points = [ _inter.point for _inter in inters ]

    inter_x,inter_y = zip( *inter_points )

    texts = []

    # -------------------------
    # Edge labels
    # -------------------------

    for edge_i, e in enumerate(mesh.tri_edge_index):
        e1, e2 = mesh.get_edge_points(edge_i)
        e_avg = (e1 + e2) / 2

        if edge_limiter is not None and edge_i not in edge_limiter:
            continue;

        txt = ax.text(
            e_avg.x,
            e_avg.y,
            f"{edge_i}",
            fontsize=11,
            fontstyle="italic",
            zorder=20,
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.10",
                facecolor="white",
                edgecolor="none",
                alpha=0.65,
            ),
        )

        texts.append(txt)

    # -------------------------
    # Triangle labels
    # -------------------------
    if tri_bool:
        for tri_index, tri in enumerate(mesh.index_array):
            tot_x = 0.0
            tot_y = 0.0

            for vert_index in tri:
                p = mesh.get_point(vert_index)
                tot_x += p[0]
                tot_y += p[1]

            tri_x = tot_x / 3.0
            tri_y = tot_y / 3.0

            txt = ax.text(
                tri_x,
                tri_y,
                f"{tri_index}",
                fontsize=12,
                fontweight="bold",
                zorder=21,
                ha="center",
                va="center",
                bbox=dict(
                    boxstyle="round,pad=0.12",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.75,
                ),
            )

            texts.append(txt)


    # -------------------------
    # Intersection labels too
    # -------------------------
    for i in range(len(inter_x)):
        ax.scatter(
            inter_x[i],
            inter_y[i],
            color="red",
            s=100,
            zorder=25,
        )

        txt = ax.text(
            inter_x[i],
            inter_y[i],
            f"{inters[i].letter}",
            fontsize=11,
            fontweight="bold",
            zorder=30,
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.15",
                facecolor="white",
                edgecolor="none",
                alpha=0.85,
            ),
        )

        texts.append(txt)

    #
    #
    #
    start= mesh.interpolate_xy_coords((1, 1))
    start_width = np.spacing( np.float32(start[0]) ) / 5000

    multipliers = 10 ** np.arange(1, 10)

    width_series = multipliers * start_width

    extent = [ start[0] - width_series[0]\
              , start[0] + width_series[0]\
              , start[1] - 1.5
              , start[1] + 1.5 ]

    if extent is not None:  
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
    
    #ax.set_aspect("equal", adjustable="box")
    ax.set_aspect("auto")


    adjust_text(
        texts,
        ax=ax,
        avoid_self=True,
        ensure_inside_axes=True,
        expand=(1.10, 1.20),
        force_text=(0.25, 0.35),
        force_static=(0.25, 0.35),
        force_pull=(0.02, 0.02),
        max_move=(12, 12),
        iter_lim=200,
    )


    
    #p.subplots_adjust(bottom=bottom_perc)

    def update(frame_i):
        width = width_series[frame_i % len(width_series)]

        ax.set_xlim(
            start[0] - width,
            start[0] + width,
        )

        # ax.set_title(f"x window center = {xc:.3f}")

        # Return empty tuple when not using blit
        return ()
    
    num_frames = 150
    anim = FuncAnimation(
        fig,
        update,
        frames=num_frames,
        interval=1000,   # milliseconds per frame
        blit=False,    # use False when changing axis limits
    )

    plt.show()



def add_letter(inters):
    for i,_inter in enumerate(inters):
        _inter.letter = string.ascii_uppercase[i]


def plt_triang_diag_base(triang , mesh, all_clip_points):
    plt.triplot(
        triang,
        color="0.65",
        linewidth=0.8,
        linestyle="-",
        marker="o",
        markersize=2.5,
        label="Tri Mesh",
        zorder=1,
    )

    clip_xs,clip_ys = zip(*all_clip_points)

    plt.plot(
        clip_xs,
        clip_ys,
        color="0.15",
        linewidth=1.5,
        linestyle="--",
        zorder=2,
        label="Clip Poly",
    )

    pos_adj = 0.05

    max_d = max(mesh.num_rows, mesh.num_cols)

    # 0th 
    loc = mesh.interpolate_xy_coords((0,0))
    plt.text(loc[0],
                loc[1],
                f"{0}",
                fontsize=12,
                style="italic",
                ha="right",
                va="bottom",
                zorder=20,
            )

    label_offset = 0.00
    for i in range(1, max_d-1):

        loc = mesh.interpolate_xy_coords((-label_offset, i))
        d = -i
        plt.text(
            loc[0],
            loc[1],
            f"{d}",
            fontsize=12,
            style="italic",
            ha="right",
            va="center",
            zorder=20,
        )

        d = i 
        loc = mesh.interpolate_xy_coords((i, -label_offset))
        plt.text(
            loc[0],
            loc[1],
            f"{d}",
            fontsize=12,
            style="italic",
            ha="center",
            va="bottom",
            zorder=20,
        )

    return plt

def plt_triang_clip(triang , mesh, all_clip_points, tri_label_bool=False):

    clip_xs,clip_ys = zip(*all_clip_points)

    p = plt_triang_base(triang , mesh, tri_label_bool)

    p.plot(
        clip_xs,
        clip_ys,
        color="0.15",
        linewidth=1.5,
        linestyle="--",
        zorder=2,
        label="Clip Poly",
    )
    return p

def plt_triang_base(triang , mesh, tri_label_bool=False):

    plt.triplot(
        triang,
        color="0.65",
        linewidth=0.8,
        linestyle="-",
        marker="o",
        markersize=2.5,
        label="Tri Mesh",
        zorder=1,
    )

    # EDGE LABELS
    # for edge_i, e in enumerate(mesh.tri_edge_index):
    #     e1,e2 = mesh.get_edge_points(edge_i)
    #     e_avg = (e1+e2)/2
    #     plt.text(e_avg.x, e_avg.y, f"{edge_i}",\
    #                 fontsize=11, style="italic", zorder=20)

    # if tri_label_bool:
    #     for tri_index, tri in enumerate(mesh.index_array):
    #         tot_x = 0.0
    #         tot_y = 0.0

    #         for vert_index in tri:
    #             #_x,_y = vert_index % dem_patch_num_verts,vert_index // dem_patch_num_verts
    #             tot_x += mesh.get_point(vert_index)[0]
    #             tot_y += mesh.get_point(vert_index)[1]

    #         plt.text(tot_x/3, tot_y/3, f"{tri_index}"\
    #                 , fontsize=12, weight="bold",  zorder=20)
    
    return plt

#
# format our intersection to show on the
# plot
def format_inter(inter: Intersection) -> str:

    mesh_label = inter.mesh_edge_index if inter.mesh_feature == FeatureType.EDGE else inter.mesh_vertex_index

    return rf'''Intersection {inter.letter} from {inter.graticule_type.name} b/w pixels {inter.grat_pix0}, {inter.grat_pix1}:
      $M$ {{type: {inter.mesh_feature.name}, id: {mesh_label}, $t_m$: {inter.mesh_t : .2f}}}, $C$ {{type: {inter.clip_feature.name}, id: {inter.clip_poly_id[1]}, $t_c$: {inter.clip_t : .2f}}}'''
    

'''
Intersection(point=(543264.0, 4943493.0), mesh_feature=<FeatureType.EDGE: 2>, clip_feature=<FeatureType.VERTEX: 1>
    , mesh_edge_index=7, mesh_edge_order=None, mesh_vertex_index=None
    , graticule_type=<Graticule.LON: 1>, clip_poly_id=None, clip_poly_order=None
    , mesh_t=0.5, clip_t=0.0), 
'''

#
# the PNTs have to be 2D
#
def plot_points(triang , mesh, pnts, inters, clip_points , plot_path, edge_limiter=None, tri_bool=True, extent=None):


    #p = plt_triang_clip(triang, mesh, clip_points+[clip_points[0]] )
    fig = plt.figure(figsize=(8, 10), dpi=120)

    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3, 1],
        hspace=0.08,
    )

    ax = fig.add_subplot(gs[0])
    ax_text = fig.add_subplot(gs[1])
    ax.triplot(
        triang,
        color="0.65",
        linewidth=0.8,
        linestyle="-",
        marker="o",
        markersize=2.5,
        label="Tri Mesh",
        zorder=1,
    )

    # for edge_i, e in enumerate(mesh.tri_edge_index):
    #     e1,e2 = mesh.get_edge_points(edge_i)
    #     e_avg = (e1+e2)/2
    #     ax.text(e_avg.x, e_avg.y, f"{edge_i}",\
    #                 fontsize=11, style="italic", zorder=20)

    all_clip_points = clip_points+[clip_points[0]]
    clip_xs,clip_ys = zip(*all_clip_points)

    ax.plot(
        clip_xs,
        clip_ys,
        color="0.15",
        linewidth=1.5,
        linestyle="--",
        zorder=2,
        label="Clip Poly",
    )

    #
    # clip points
    x,y = zip( *pnts )
    ax.plot(x, y, zorder=10, color="green")

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

    # offsets = [
    #     (10, 10),    # upper right
    #     (-10, 10),   # upper left
    #     (10, -10),   # lower right
    #     (-10, -10),  # lower left
    #     (16, 0),     # right
    #     (-16, 0),    # left
    #     (0, 16),     # above
    #     (0, -16),    # below
    # ]

    # for i in range(len(inter_x)):
    #     dx, dy = offsets[i % len(offsets)]

    #     ha = "left" if dx > 0 else "right" if dx < 0 else "center"
    #     va = "bottom" if dy > 0 else "top" if dy < 0 else "center"
    #     ax.scatter(inter_x[i],inter_y[i], color='red', s=100)
    #     ax.annotate(
    #         f'{inters[i].letter}',  
    #         xy=( inter_x[i],inter_y[i] ) ,
    #         xytext=(dx, dy) ,
    #         textcoords="offset points" ,       
    #         ha=ha,
    #         va=va,
    #     )

    # p.figtext(
    #     0.45,                    # x-coordinate (centered)
    #     bottom_start,            # y-coordinate (near the absolute bottom)
    #     '\n'.join(frmt_inters),  # LaTeX text (must use 'r' for raw string)
    #     ha="center",
    #     color="#1f42b4",
    #     fontsize=10
    # )


    texts = []

    # -------------------------
    # Edge labels
    # -------------------------

    for edge_i, e in enumerate(mesh.tri_edge_index):
        e1, e2 = mesh.get_edge_points(edge_i)
        e_avg = (e1 + e2) / 2

        if edge_limiter is not None and edge_i not in edge_limiter:
            continue;

        txt = ax.text(
            e_avg.x,
            e_avg.y,
            f"{edge_i}",
            fontsize=11,
            fontstyle="italic",
            zorder=20,
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.10",
                facecolor="white",
                edgecolor="none",
                alpha=0.65,
            ),
        )

        texts.append(txt)

    # -------------------------
    # Triangle labels
    # -------------------------
    if tri_bool:
        for tri_index, tri in enumerate(mesh.index_array):
            tot_x = 0.0
            tot_y = 0.0

            for vert_index in tri:
                p = mesh.get_point(vert_index)
                tot_x += p[0]
                tot_y += p[1]

            tri_x = tot_x / 3.0
            tri_y = tot_y / 3.0

            txt = ax.text(
                tri_x,
                tri_y,
                f"{tri_index}",
                fontsize=12,
                fontweight="bold",
                zorder=21,
                ha="center",
                va="center",
                bbox=dict(
                    boxstyle="round,pad=0.12",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.75,
                ),
            )

            texts.append(txt)


    # -------------------------
    # Intersection labels too
    # -------------------------
    for i in range(len(inter_x)):
        ax.scatter(
            inter_x[i],
            inter_y[i],
            color="red",
            s=100,
            zorder=25,
        )

        txt = ax.text(
            inter_x[i],
            inter_y[i],
            f"{inters[i].letter}",
            fontsize=11,
            fontweight="bold",
            zorder=30,
            ha="center",
            va="center",
            bbox=dict(
                boxstyle="round,pad=0.15",
                facecolor="white",
                edgecolor="none",
                alpha=0.85,
            ),
        )

        texts.append(txt)

    if extent is not None:  
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
    
    #ax.set_aspect("equal", adjustable="box")
    ax.set_aspect("auto")


    adjust_text(
        texts,
        ax=ax,
        avoid_self=True,
        ensure_inside_axes=True,
        expand=(1.10, 1.20),
        force_text=(0.25, 0.35),
        force_static=(0.25, 0.35),
        force_pull=(0.02, 0.02),
        max_move=(12, 12),
        iter_lim=200,
    )

    ax_text.axis("off")
    ax_text.text(
        0.02,
        0.95,
        '\n'.join(frmt_inters),
        transform=ax_text.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        family="sans-serif",
    )
    
    #p.subplots_adjust(bottom=bottom_perc)

    plt.savefig(plot_path, dpi=150,facecolor="white")

#
#
#
def plot_points_arrows( triang , mesh, clip_pnts, dir, plot_path ):

    fig, ax = plt.subplots()
    p = plt_triang_clip(triang, mesh, clip_pnts )

    #
    # clip points
    #x,y = zip( *clip_pnts )
    #p.plot(x, y, zorder=10, color="green")

    ax.annotate("", 
            xy=dir[1],      # tip destination
            xytext=dir[0],  # arrow starting point
            arrowprops=dict(arrowstyle="->", color="red", lw=2))

    if plot_path is None:
        p.show()
    else:
        p.savefig(plot_path, dpi=300, bbox_inches='tight')


def plot_clip_inters( triang , mesh, clip_verts, inters, remaining_clip, full_clip, plot_path ):

    CLIPPED_VERT_COLOR = "#D55E00" 
    COL_EDGE_INTER_COLOR = "#0072B2"
    COL_CLIP_ONLY_COLOR = "#CC79A7"

    fig, ax = plt.subplots()
    p = plt_triang_clip(triang, mesh, full_clip )

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

        cls.mesh_cut_facto = MeshCutFactory(raster_path, 2)
        #
        # edit these 2 variables
        # offset_x, offset_y = 0,0 
        # dem_num_cells = 2

        # dem_resolution = 10

        # dem_patch_size = dem_resolution * dem_num_cells

        # dem_patch_num_verts = dem_num_cells+1
        
        # buffer_array = np.zeros((dem_patch_num_verts, dem_patch_num_verts, 3), dtype=np.float32)

        # num_tri = (dem_patch_num_verts - 1) * (dem_patch_num_verts-1) * 2

        # with rasterio.open(raster_path) as src:
        #     band = src.read(1)
        #     transf = src.transform

        #     bounds = src.bounds

        #     num_rows_total = band.shape[0]
        #     num_cols_total = band.shape[1]

        #     band = band[offset_y:offset_y+dem_patch_num_verts,offset_x:offset_x+dem_patch_num_verts]
        #     h = src.height
        #     w = src.width

        #     rows = np.arange( dem_patch_num_verts ) + offset_y
        #     cols = np.arange( dem_patch_num_verts ) + offset_x

        #     num_rows_inset = band.shape[0]
        #     num_cols_inset = band.shape[1]

        #     orig_transf = src.transform 
            
        #     new_transf =  orig_transf * Affine.translation(
        #         offset_x,   # column offset
        #         offset_y    # row offset
        #     )

        #     index_array = []

        #     for r in range(num_rows_inset):
        #         for c in range(num_cols_inset):

        #             # should be in meters
        #             x_m, y_m = rasterio.transform.xy(new_transf, r, c, offset='ul')
        #             elev = band[r, c]

        #             if elev < 0:
        #                 elev = 0

        #             buffer_array[c, r, 0] = x_m #2ND
        #             buffer_array[c, r, 1] = y_m
        #             buffer_array[c, r, 2] = elev

        #             if (r > 0 and c > 0) and ( r < dem_patch_num_verts and c < dem_patch_num_verts):
        #                 currIndx = c + r * dem_patch_num_verts;
        #                 prevRowSameColIndx = c + (r - 1) * dem_patch_num_verts;
        #                 prevRowBackColIndx = c - 1 + (r - 1) * dem_patch_num_verts;
        #                 prevIndx = currIndx - 1;
        #                 index_array.append([prevRowBackColIndx, currIndx, prevRowSameColIndx])
        #                 index_array.append([prevRowBackColIndx, prevIndx, currIndx])


        cls.mesh = cls.mesh_cut_facto.get_mesh()
        cls.triang = cls.mesh_cut_facto.get_triang()
        cls.mesh_cut_facto.update(raster_path, 3)

        cls.second_mesh = cls.mesh_cut_facto.get_mesh()
        cls.second_triang = cls.mesh_cut_facto.get_triang()

        cls.mesh_cut_facto.update(raster_path, 2)

        cls.eps_mesh = cls.mesh_cut_facto.get_mesh()
        cls.eps_triang = cls.mesh_cut_facto.get_triang()

        # 9 PM
        # firstCoord = (buffer_array[1,0, 0].item(), ((buffer_array[1,1, 1] + buffer_array[1,0, 1])/2).item() )
        # # noon
        # secondCoord = (((buffer_array[1,0,0]+buffer_array[2,0, 0])/2).item(), buffer_array[1,0,1].item() )
        # # 3 PM
        # thirdCoord = (buffer_array[2,0,0].item() , ((buffer_array[1,1, 1] + buffer_array[1,0, 1])/2).item()  )

        # fourthCoord = (((buffer_array[1,0, 0]+buffer_array[2,0, 0])/2).item(), buffer_array[1,1, 1].item() )

        # cls.diamond_clip = [firstCoord , secondCoord, thirdCoord, fourthCoord ]

        cls.diamond_clip = [
            (cls.mesh.get_point(1,0)[0], cls.mesh.interpolate_xy_coords((0,0.5))[1]) ,\
            (cls.mesh.interpolate_xy_coords((1.5,0))[0], cls.mesh.get_point(0,0)[1]) ,\
            (cls.mesh.get_point(2,0)[0], cls.mesh.interpolate_xy_coords((0,0.5))[1]) ,\
            (cls.mesh.interpolate_xy_coords((1.5,0))[0], cls.mesh.get_point(1,1)[1])
        ]

        cls.square_clip = [ cls.mesh.get_point(1,0)[:2], cls.mesh.get_point(2,0)[:2],
                            cls.mesh.get_point(2,1)[:2], cls.mesh.get_point(1,1)[:2] ]

        cls.bigger_rect = [ cls.mesh.interpolate_xy_coords((0.75,0)), cls.mesh.interpolate_xy_coords((2,0)) ,\
                            cls.mesh.interpolate_xy_coords((2,1)), cls.mesh.interpolate_xy_coords((0.75,1)) ]

        cls.bigger_rect2 = [ [cls.mesh.interpolate_xy_coords((0.3, 0))[:2], cls.mesh.interpolate_xy_coords((2, 0))[:2] ,\
                            cls.mesh.interpolate_xy_coords((2, 1.25))[:2], cls.mesh.interpolate_xy_coords((0.3, 1.25))[:2] ] ]





    def get_point_categories(self , poly_clip ):

        clp_verts = []
        for vert, clip_bool in self.mesh.vert_clip_index.items():
            if clip_bool:
                clp_verts.append( self.mesh.get_point(vert)[:2] )

        inters = []
        for edge_indx, edge_dict in self.mesh.all_mesh_edges.items():

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
    #
    #
    def get_diagonals_searched(self, clip_points):

        poly_num_edges = len(clip_points)

        max_x_all, max_y_all = -math.inf, -math.inf
        min_x_all, min_y_all = math.inf, math.inf
        for _j in range( len(clip_points) ):
            clip_p0 = clip_points[_j]
            clip_p1 = clip_points[(_j+1) % poly_num_edges]

            p0_x_pix = (clip_p0[0] - self.mesh.transform.c ) / self.mesh.transform.a;
            p1_x_pix = (clip_p1[0] - self.mesh.transform.c ) / self.mesh.transform.a;

            p0_y_pix = (clip_p0[1] - self.mesh.transform.f ) / self.mesh.transform.e;
            p1_y_pix = (clip_p1[1] - self.mesh.transform.f ) / self.mesh.transform.e;

            max_x_pix = max( p0_x_pix, p1_x_pix )
            min_x_pix = min( p0_x_pix, p1_x_pix )

            max_y_pix = max( p0_y_pix, p1_y_pix )
            min_y_pix = min( p0_y_pix, p1_y_pix )

            max_x_all = max(max_x_pix, max_x_all)
            max_y_all = max(max_y_pix, max_y_all)
            min_x_all = min(min_x_pix, min_x_all)
            min_y_all = min(min_y_pix, min_y_all)


        max_x_all = math.ceil(max_x_all)
        max_y_all = math.ceil(max_y_all)
        min_x_all = math.floor(min_x_all)
        min_y_all = math.floor(min_y_all)
        
        min_d = min_x_all - max_y_all
        max_d = max_x_all - min_y_all

        diag_l = []
        for d in range(max(min_d, self.mesh.min_diag), min(max_d + 1, self.mesh.max_diag)):

            x1 = max(0, d)
            x2 = min(self.mesh.num_cols - 1, d + self.mesh.num_rows - 1)

            y1 = x1 - d
            y2 = x2 - d

            _pix0 = x1, y1
            _pix1 = x2, y2
            
            _p0 = self.mesh.convert_to_point(*_pix0)[:2];
            _p1 = self.mesh.convert_to_point(*_pix1)[:2];

            diag_l.append( Line(_p0, _p1, _pix0, _pix1) )

        return diag_l

    #
    # returns [Intersection] ,
    #           [tup2D, tup2D]
    #
    def fetch_intersections_clip_i( self, clip_points, i ):

        curr_i = i % 4 # ASSUMING 4 clip points
        next_i = (i+1) % 4
        clip_p0 = clip_points[curr_i]
        clip_p1 = clip_points[next_i]

        clip_edge_inters = self.mesh.clip_mesh_edges(clip_p0 , clip_p1, (0,i), 4)         
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
        glancing_inter1 = self.mesh.get_intersection(Graticule.LON, self.diamond_clip[0], self.diamond_clip[1] \
                                                        , _e0, _e1
                                                        , (0,1), 4)

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
                    , self.diamond_clip
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
                    , self.diamond_clip
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
                    , self.diamond_clip
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
                    , self.diamond_clip
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
                    , self.square_clip
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
                    , self.square_clip
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
                    , self.square_clip
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
                    , self.square_clip
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
                    , self.bigger_rect
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
                    , self.bigger_rect
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
                    , self.bigger_rect
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
                    , self.bigger_rect
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
                      self.mesh.interpolate_xy_coords((0.4,1))[:2] ]

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
                         , inters, remaining_clip \
                         , self.bigger_rect2[0] + [self.bigger_rect2[0][0]],  plot_path )


    def test_diagonal_graticules(self):

        all_diags = self.get_diagonals_searched( self.bigger_rect2[0] )

        all_clips = self.bigger_rect2[0] + [self.bigger_rect2[0][0]]

        p = plt_triang_diag_base(self.second_triang , self.second_mesh, all_clips )

        for line in all_diags:
            x,y = zip(* [line.p0, line.p1] )
            p.plot(x, y, color='pink', linewidth=2)

        self.assertCountEqual(all_diags, set(all_diags))
        plot_path = self.cfg.output_dir / "diagonals_rectangle.png"
        p.savefig(plot_path, dpi=300, bbox_inches='tight')

        p.close()
        #
        # now for square
        all_diags = self.get_diagonals_searched( self.square_clip )

        all_clips = self.square_clip + [ self.square_clip[0] ]

        p = plt_triang_diag_base(self.second_triang , self.second_mesh, all_clips )
        self.assertCountEqual(all_diags, set(all_diags))

        for line in all_diags:
            x,y = zip(* [line.p0, line.p1] )
            p.plot(x, y, color='pink', linewidth=2)

        plot_path = self.cfg.output_dir / "diagonals_square.png" ;

        p.savefig(plot_path, dpi=300, bbox_inches='tight')


    #
    #
    #
    #
    def test_road_width_degenerecy( self ):

        point1 = (*self.eps_mesh.interpolate_xy_coords((0.5, 0.75))[:2], 0)
        point2 = (*self.eps_mesh.interpolate_xy_coords((0.5, 1.25))[:2], 0)

        points = [point1,point2]
        indx = [0,1]

        info = np.finfo(np.float32)

        start= self.eps_mesh.interpolate_xy_coords((1, 1))

        EPS = np.spacing( np.float32(start[1]) )

        points_CW, index_CW = self.eps_mesh.create_quad_from_points(points, indx, EPS)

        #
        # not degenerate clipping poly 
        # let's examine the intersections that get produced
        self.assertFalse(points_CW[0][0] == points_CW[1][0])

        self.eps_mesh.perform_clipping( [ points_CW ] )

        self.assertEqual( len(self.eps_mesh.all_mesh_edges[2]), 2 )

        inters = self.eps_mesh.all_mesh_edges[2]

        self.assertNotEqual(inters[1].point[0] , inters[0].point[0] )
        self.assertEqual( inters[1].clip_feature, FeatureType.EDGE )
        self.assertEqual( inters[0].clip_feature, FeatureType.EDGE )

        EPS = np.spacing( np.float32(start[1]) ) / 50

        #
        # DEGENERATE clipping poly 
        # 

        points_CW_2, index_CW_2 = self.eps_mesh.create_quad_from_points(points, indx, EPS)
        #self.assertTrue(points_CW[0][0] == points_CW[1][0])
        # 
        # what happens when we pass by eps of mesh vertex?
        #
        point1 = (*self.eps_mesh.interpolate_xy_coords((1, 0.9))[:2], 0)
        point2 = (*self.eps_mesh.interpolate_xy_coords((1, 1.5))[:2], 0)

        points = [point1, point2]
        
        vert_test_points_CW, index_CW = self.eps_mesh.create_quad_from_points(points, indx, 0.5)

        self.assertEqual( vert_test_points_CW[1][0]-0.5 , self.eps_mesh.buffer_array[1,0,0].item() )
        shft_verts = []
        EPS = np.spacing( np.float32(start[1]) )
        for vert in vert_test_points_CW:
            new_vert = (vert[0]-0.5-EPS.item(), vert[1], vert[2])
            shft_verts.append(new_vert)

        self.assertNotEqual( vert_test_points_CW[0][0], shft_verts[0][0] )

        clip_p0 = shft_verts[1]
        clip_p1 = shft_verts[2]

        poly_edge_inters = self.eps_mesh.clip_mesh_edges( clip_p0, clip_p1, (0,1), 4 )

        #
        # 
        self.assertEqual( poly_edge_inters[1].mesh_feature, FeatureType.EDGE )
        self.assertEqual( poly_edge_inters[0].mesh_feature, FeatureType.EDGE )


        '''
        [Intersection(point=(543263.9900000002, 4943493.0), mesh_feature=<FeatureType.EDGE: 2>, clip_feature=<FeatureType.VERTEX: 1>, mesh_edge_index=6, mesh_edge_order=None, mesh_vertex_index=None, letter=None, graticule_type=<Graticule.LON: 1>, grat_pix0=(1, 0), grat_pix1=(1, 2), clip_poly_id=(0, 1), clip_poly_order=0, mesh_t=0.5, clip_t=0.0),
        Intersection(point=(543264.0, 4943488.0), mesh_feature=<FeatureType.VERTEX: 1>, clip_feature=<FeatureType.EDGE: 2>, mesh_edge_index=3, mesh_edge_order=None, mesh_vertex_index=4, letter=None, graticule_type=<Graticule.LAT: 2>, grat_pix0=(0, 1), grat_pix1=(2, 1), clip_poly_id=(0, 1), clip_poly_order=1, mesh_t=0.0, clip_t=0.5),
        Intersection(point=(543264.0, 4943488.0), mesh_feature=<FeatureType.VERTEX: 1>, clip_feature=<FeatureType.EDGE: 2>, mesh_edge_index=15, mesh_edge_order=None, mesh_vertex_index=4, letter=None, graticule_type=<Graticule.DIAG: 0>, grat_pix0=(0, 0), grat_pix1=(2, 2), clip_poly_id=(0, 1), clip_poly_order=2, mesh_t=0.0, clip_t=0.5)]
        '''
        shft_verts = []
        EPS = ( np.spacing( np.float32(start[0]) ) / 50 ).item()

        pre_adj_clip_p0 = self.eps_mesh.interpolate_xy_coords((1, 0.5))[:2]
        pre_adj_clip_p1 = self.eps_mesh.interpolate_xy_coords((1, 1.5))[:2]

        clip_p0 = (pre_adj_clip_p0[0]-EPS , pre_adj_clip_p0[1])
        clip_p1 = (pre_adj_clip_p1[0]-EPS , pre_adj_clip_p1[1])

        #
        # DEGENERATE! 
        #
        poly_edge_inters = self.eps_mesh.clip_mesh_edges( clip_p0, clip_p1, (0,1), 4 )

        #
        #
        # further interogate
        #

        e0 = self.eps_mesh.get_point(1,0)[:2]
        e1 = self.eps_mesh.get_point(1,1)[:2]

        self.assertNotEqual( clip_p0[0] , e0[0] )

        p0_glm = glm.vec2(clip_p0)
        p1_glm = glm.vec2(clip_p1)
        e0_glm = glm.vec2(e0)
        e1_glm = glm.vec2(e1)

        #
        # want to get the smallest difference we can for 
        # a GLM vec2 which uses 32 bit float
        #
        self.assertNotEqual(glm.vec2(start[0] - np.spacing( np.float32(start[0]) ) ,1) ,\
                             glm.vec2(start[0],1) )

        self.assertEqual(glm.vec2(start[0] - np.spacing( np.float32(start[0]) )/2 ,1) ,\
                             glm.vec2(start[0],1) )

        clip_p0 = ((pre_adj_clip_p0[0] - np.spacing( np.float32(start[0]) )).item() \
                    , pre_adj_clip_p0[1])
        clip_p1 = ((pre_adj_clip_p1[0] - np.spacing( np.float32(start[0]) ) ).item() \
                    , pre_adj_clip_p1[1])

        shft_verts = [\
            (clip_p0[0] - 1, clip_p0[1]),
            clip_p0,
            clip_p1,
             (clip_p1[0] - 1, clip_p1[1])
        ]

        poly_edge_inters = self.eps_mesh.clip_mesh_edges( clip_p0, clip_p1, (0,1), 4 )

        #
        # DEGENERECY! 
        # two different mesh edge intersections share the same intersection point!
        #
        #self.assertNotEqual(poly_edge_inters[0].point, poly_edge_inters[1].point)


        p0_glm = glm.vec2(clip_p0)
        p1_glm = glm.vec2(clip_p1)

        e0 = self.eps_mesh.get_point(0,0)[:2]
        e1 = self.eps_mesh.get_point(1,1)[:2]

        e0_glm = glm.vec2(e0)
        e1_glm = glm.vec2(e1)

        mesh_dir = e1_glm - e0_glm;
        clip_dir = p1_glm - p0_glm;
        rel_dir = e0_glm - p0_glm
        eps = 1e-9

        diag_denominator = self.mesh.cross_2D(clip_dir, mesh_dir)
        diag_clip_t = self.mesh.cross_2D(rel_dir, mesh_dir) / diag_denominator

        e0 = self.eps_mesh.get_point(0,1)[:2]
        e1 = self.eps_mesh.get_point(1,1)[:2]

        e0_glm = glm.vec2(e0)
        e1_glm = glm.vec2(e1)

        mesh_dir = e1_glm - e0_glm;
        clip_dir = p1_glm - p0_glm;
        rel_dir = e0_glm - p0_glm

        diag_denominator = self.mesh.cross_2D(clip_dir, mesh_dir)
        diag_clip_t = self.mesh.cross_2D(rel_dir, mesh_dir) / diag_denominator

        #self.assertGreater( abs(self.mesh.cross_2D(clip_dir, mesh_dir )), eps   )
        self.assertGreater( abs( self.mesh.cross_2D(mesh_dir, rel_dir)), eps )

        add_letter( poly_edge_inters )
        #
        # [xmin, xmax, ymin, ymax]
        # want to avoid ymin == ymax
        window = 1.0
        micro_extent = [ self.eps_mesh.buffer_array[1,1,0] - np.spacing( np.float32(start[0]) ) * 2 ,\
                         self.eps_mesh.buffer_array[1,1,0] + np.spacing( np.float32(start[0]) ) * 2 ,\
                         self.eps_mesh.buffer_array[1,1,1] - window, \
                         self.eps_mesh.buffer_array[1,1,1] + window ]

        file_path =  self.cfg.output_dir / f"micro_diff_degenerecy_new.png"

        
        plt.close()
        plot_points( self.triang 
                    , self.mesh
                    , [clip_p0[:2], clip_p1[:2]]
                    , poly_edge_inters
                    , [ v[:2] for v in shft_verts ]
                    , file_path 
                    , [12,7,2]
                    , False
                    , micro_extent ) 

    def test_degenerate_intersection(self):

        start = self.mesh.interpolate_xy_coords((1, 1))

        pre_adj_clip_p0 = self.mesh.interpolate_xy_coords((1, 0.5))[:2]
        pre_adj_clip_p1 = self.mesh.interpolate_xy_coords((1, 1.5))[:2]
        e0 = self.mesh.get_point(0,0)[:2]
        e1 = self.mesh.get_point(1,1)[:2]

        clip_p0 = ((pre_adj_clip_p0[0] - np.spacing( np.float32(start[0]) )).item() \
                    , pre_adj_clip_p0[1])
        clip_p1 = ((pre_adj_clip_p1[0] - np.spacing( np.float32(start[0]) ) ).item() \
                    , pre_adj_clip_p1[1])
                
        e0_glm = glm.vec2(e0)
        e1_glm = glm.vec2(e1)

        p0_glm = glm.vec2(clip_p0)
        p1_glm = glm.vec2(clip_p1)

        mesh_dir = e1_glm - e0_glm;
        clip_dir = p1_glm - p0_glm;
        rel_dir = e0_glm - p0_glm

        diag_denominator = self.mesh.cross_2D(clip_dir, mesh_dir)
        diag_clip_t = self.mesh.cross_2D(rel_dir, mesh_dir) / diag_denominator
        diag_intersection = tuple(p0_glm + diag_clip_t * clip_dir)

        edge_index = self.mesh.get_mesh_edge_inter(
                                Graticule.DIAG, 
                                 diag_intersection )

        #
        # notice this result is "wrong" because this degenerate intersection
        # falls on the latitude line it gets routed to (1,0)
        #
        self.assertEqual(edge_index, 14) # 12 is first diag, 14 is second




    def test_degeneracy_resolution(self):

        shft_verts = []

        #
        # when a vertex is passed on the inside
        #
        start = self.eps_mesh.interpolate_xy_coords((1, 1))
        pre_adj_clip_p0 = self.eps_mesh.interpolate_xy_coords((1, 0.5))[:2]
        pre_adj_clip_p1 = self.eps_mesh.interpolate_xy_coords((1, 1.5))[:2]

        clip_p0 = ((pre_adj_clip_p0[0] - np.spacing( np.float32(start[0]) )).item() \
                    , pre_adj_clip_p0[1])
        clip_p1 = ((pre_adj_clip_p1[0] - np.spacing( np.float32(start[0]) ) ).item() \
                    , pre_adj_clip_p1[1])

        shft_verts = [\
            (clip_p0[0] - 1, clip_p0[1]),
            clip_p0,
            clip_p1,
             (clip_p1[0] - 1, clip_p1[1])
        ]

        poly_edge_inters = self.eps_mesh.clip_mesh_edges( clip_p0, clip_p1, (0,1), 4 )

        self.assertEqual(len(poly_edge_inters) , 1 )

        self.assertEqual(poly_edge_inters[0].mesh_feature, FeatureType.VERTEX)

        # plot_points_anim( self.triang 
        #             , self.mesh
        #             , [clip_p0[:2], clip_p1[:2]]
        #             , poly_edge_inters
        #             , [ v[:2] for v in shft_verts ]
        #             , file_path 
        #             , [12,7,2]
        #             , False
        #             , micro_extent ) 





if __name__ == "__main__":
    sys.stdout = sys.__stdout__
    unittest.main(argv=["first-arg-is-ignored"], exit=False)