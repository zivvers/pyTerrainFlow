import requests
import math
#from config import load_config
import py3dep
from osgeo import gdal, ogr, osr
import tempfile
import os

# globales
IMAGE_SERVER = "https://imagery.oregonexplorer.info/arcgis/rest/services/OSIP_2024/OSIP_2024_SL/ImageServer/exportImage"


def transform_point(x, y) :
    src = osr.SpatialReference()
    src.ImportFromEPSG(26910)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4326)   # or 4269

    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    transform = osr.CoordinateTransformation(src, dst)

    lon, lat, _ = transform.TransformPoint(float(x), float(y), 0.0)

    return lon, lat


def download_osip_tile(bbox: str, full_path: str) -> bool:

    #bbox = f"{x_min},{y_min},{x_min+tile_size},{y_min+tile_size}"
    # print(f"BBOX: {bbox}")

    #
    # Oregon spans two UTM zones:
    #W: Zone 10 , EPSG:26910
    #E: Zone 11 , EPSG:26911
    #

    params = {
        "bbox": bbox,
        "size": "1024,1024",
        "bboxSR": 26910,
        "imageSR": 26910,
        "format": "tiff",
        "pixelType": "U8",
        "f": "Image",
        "noData":-999,
        "noDataInterpretation" : "esriNoDataMatchAny",
        "interpolation" : "RSP_BilinearInterpolation",
        "adjustAspectRatio" : True,
        "validateExtent" : False
    }
    
    response = requests.get(IMAGE_SERVER, params=params)

    if response.status_code == 200:
        
        with open(full_path, "wb") as f:
            f.write(response.content)
        # print(f"Saved {full_path}")
        return True
    else:
        print(f"Failed tile download! w/ {response} for bbox: {bbox}")
        return False



if __name__ == "__main__":

    #cfg = load_config()

    TILE_SIZE = 1024

    # EPSG:26910 w Oregon approx box
    # [360000, 4650000, 620000, 5140000]

    western_oregon_bbox = [360000, 4650000, 620000, 5140000]

    global_xmin, global_ymin, global_xmax, global_ymax = western_oregon_bbox

    ncols = math.ceil((global_xmax - global_xmin) / TILE_SIZE)
    nrows = math.ceil((global_ymax - global_ymin) / TILE_SIZE)

    print(f"# cols: {ncols}, # rows: {nrows}")

    # arbitrary
    col_i = 12
    row_i = 12

    #input_dir = cfg.input_dir

    #assert(cfg.input_dir.exists())

    # tile_y_max = global_ymax - row_i * TILE_SIZE
    # tile_y_min = tile_y_max - TILE_SIZE
    # tile_x_min = global_xmin + col_i * TILE_SIZE
    # tile_x_max = tile_x_min + TILE_SIZE
    #txmax = min(xmax, txmin + step)

    #
    # !!
    #

    tile_x_min, tile_y_max = 540354+2000, 4942208+2000
    tile_x_max = tile_x_min + TILE_SIZE
    tile_y_min = tile_y_max - TILE_SIZE

    bounds = (tile_x_min, tile_y_min, tile_x_max, tile_y_max)



    #lon, lat = transform_point(tile_x_min, tile_y_min)

    #print(f"LON: {lon}, LAT: {lat}")

    bbox = ",".join([str(item) for item in bounds])

    #imagery_tif = f"osip_EPSG_26910_{col_i}_{row_i}_{int(tile_x_min)}_{int(tile_y_min)}.tif"
    imagery_tif = f"osip_EPSG_26910_{int(tile_x_min)}_{int(tile_y_max)}.tif"
    dem_tif = f"dem_EPSG_26910_{int(tile_x_min)}_{int(tile_y_max)}.tif"

    imagery_full_path =  imagery_tif

    dem_full_path =  dem_tif
    
    saved = download_osip_tile(bbox, imagery_full_path)

    if saved:
        print(f"successfully saved: {imagery_full_path}")
    else:
        print("could not fetch!!")
    

    print("get DEM")
    dem_bounds = (tile_x_min, tile_y_min - 6, tile_x_max + 6, tile_y_max)

    assert(dem_bounds[2] - dem_bounds[0] == 1030)
    assert(dem_bounds[3] - dem_bounds[1] == 1030)

    #
    # geo_crs for get_dem function
    # is for newer py3dep
    dem = py3dep.get_dem(
        dem_bounds,
        resolution=10,  # 10 meters
        crs=26910       
        )

    tmp_file = tempfile.NamedTemporaryFile(suffix=".tif", delete=False).name
    dem.rio.to_raster(tmp_file)

    print(f"raster width: {dem.rio.width}, height: {dem.rio.height}")

    raster_ds = gdal.Open(tmp_file)

    new_proj = r"C:\dev\repo\Accessibility\build\vcpkg_installed\x64-windows\share\proj"
    os.environ["PROJ_LIB"] = new_proj
    os.environ["PROJ_DATA"] = os.environ.get("PROJ_LIB")

    out = gdal.Warp( str(dem_full_path),
        raster_ds,
        dstSRS="EPSG:26910",
        outputBounds=dem_bounds,
        xRes=10, # meters per pixel
        yRes=10,
        resampleAlg=gdal.GRA_Bilinear ,
        format="GTiff")

    test_ds = gdal.Open( str( dem_full_path ) )
    print("size:", test_ds.RasterXSize, test_ds.RasterYSize)
    print("geotransform:", test_ds.GetGeoTransform())






