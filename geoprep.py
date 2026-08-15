import numpy as np
import rasterio
from affine import Affine
from shapely.geometry import LineString


def create_roof_raster(
    file_name,
    n,
    peak_z=1000.0,
    slope=100.0,
    crs="EPSG:26910"
):
    """
    Create an n x n GeoTIFF where the highest elevations occur
    along a vertical ridge through the center of the raster.

    peak_z : elevation along ridge
    slope  : elevation loss per horizontal meter
    """

    transform = Affine(
        10.0, 0.0, 543254.0,
        0.0, -10.0, 4943498.0
    )

    width = n
    height = n

    # Geographic x coordinate of the raster's central vertical line.
    #
    # Raster extends from:
    # x = transform.c
    # to
    # x = transform.c + width * 10
    ridge_x = transform.c + (width * transform.a) / 2.0

    # Generate coordinates of PIXEL CENTERS
    cols = np.arange(width)

    x_coords = (
        transform.c
        + (cols + 0.5) * transform.a
    )

    # Distance each column is from ridge
    distance_from_ridge = np.abs(x_coords - ridge_x)

    # Roof profile:
    #
    #            ridge
    #              /\
    #             /  \
    #            /    \
    #
    row_z = peak_z - slope * distance_from_ridge

    # Same profile for every row
    data = np.tile(row_z, (height, 1)).astype(np.float32)

    # Make the test LineString run along the ridge,
    # stopping 10 m short of top and bottom.
    top = transform.f
    bottom = transform.f + height * transform.e

    line = LineString([
        (ridge_x, bottom + 10.0),
        (ridge_x, top - 10.0),
    ])

    with rasterio.open(
        file_name,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(data, 1)

    return line