'''
BChoat 2025/04/28

functions to reproject rasters using xarray, rioxarray, and rasterio.
'''

# %% import libraries
#################################################

import os
from pathlib import Path
import rioxarray
from rasterio.enums import Resampling
import geopandas as gpd
from shapely.geometry import box
import numpy as np

# %% function to replace values with NAs
##################################################
def apply_na_filter(raster, na_val):
    """
    Apply NA value filtering to a raster based on various condition types (improved).

    Parameters:
    - raster: xarray.DataArray - The raster to filter
    - na_val: Various types:
        - Single value: replaces exact matches with NA
        - List/tuple/set: replaces any matching values with NA
        - Dict with single key-value pair: applies operator comparison
          e.g., {">": 253} replaces values > 253 with NA

    Returns:
    - Filtered raster with NA values applied
    """
    if na_val is None:
        return raster

    mask = None

    if isinstance(na_val, (list, tuple, set)):
        mask = ~raster.isin(list(na_val))
    elif isinstance(na_val, dict):
        if len(na_val) != 1:
            raise ValueError("Condition dictionary must have a single key-value pair.")
        operator, threshold = list(na_val.items())[0]
        if operator == ">":
            mask = raster <= threshold
        elif operator == ">=":
            mask = raster < threshold
        elif operator == "<":
            mask = raster >= threshold
        elif operator == "<=":
            mask = raster > threshold
        elif operator == "==":
            mask = raster != threshold
        elif operator == "!=":
            mask = raster == threshold
        else:
            raise ValueError(f"Unsupported operator: {operator}")
    else:  # Single value
        mask = raster != na_val

    if mask is not None:
        raster = raster.where(mask)

    # Preserve existing nodata value if it exists, otherwise set to NaN
    original_nodata = raster.rio.nodata
    if original_nodata is None or np.isnan(original_nodata):
        raster = raster.rio.write_nodata(np.nan)
    else:
        raster = raster.rio.write_nodata(original_nodata)

    return raster

# %% define function
#################################################
# I had Claude write this up quickly and made only small modifications
def reproject_tif(input_tif, output_tif=None, target_crs="EPSG:4326",
                  resampling_method=Resampling.nearest, resolution=None,
                  template_raster=None, na_val=None):
    """
    Read a GeoTIFF file and reproject it to a new coordinate system.

    Parameters
    ----------
    input_tif : str, Path, or xarray.DataArray
        Path to the input GeoTIFF file
    output_tif : str, optional
        Path to save the reprojected GeoTIFF. If None, won't save to file.
    target_crs : str, optional
        Target coordinate reference system as EPSG code, WKT or Proj4 string.
        Default is "EPSG:4326" (WGS84 lat/lon). 
        Ignored if template_raster is provided.
    resampling_method : rasterio.enums.Resampling, optional
        Resampling method to use, default is nearest neighbor
    resolution : tuple or float, optional
        Target resolution in target CRS units. Can be (x_res, y_res) tuple
        or single value for both dimensions. 
        Ignored if template_raster is provided.
    template_raster : str, optional
        Path to a template raster. If provided, the input will
        be reprojected to match the CRS, resolution,
        and extent of this raster.
    na_val : numeric, list of numeric, or dict with operator, optional
        Values in the raster to be replaced with NA. Can be:
        - A single value (e.g., 255)
        - A list of values (e.g., [255, -9999])
        - A dict with operator and threshold (e.g., {">": 253}, {"==": 0})
          Supported operators: >, >=, <, <=, ==, !=
        If None, no values are replaced.

    Returns
    -------
    xarray.DataArray
        The reprojected raster as a DataArray
    """
    # Open the GeoTIFF with rioxarray
    if isinstance(input_tif, str) or isinstance(input_tif, Path):
        print(f"Reading {input_tif}...")
        raster = rioxarray.open_rasterio(input_tif)
    else:
        raster = input_tif

    if na_val:
        print("applying NA filter")
        # Apply NA value filtering
        raster = apply_na_filter(raster, na_val)

    # Handle template raster if provided
    if template_raster:
        if not os.path.exists(template_raster):
            raise FileNotFoundError(
                f"Template raster not found: {template_raster}"
                )

        print(f"Using template raster: {template_raster}")
        template = rioxarray.open_rasterio(template_raster).load()

        # Reproject to match template
        print(f"Reprojecting to match template (CRS: {template.rio.crs})...")
        reprojected = raster.rio.reproject_match(
            template,
            resampling=resampling_method
        )

        # make sure file is closed
        template.close()
    else:
        # Standard reprojection with target_crs and resolution
        # Prepare resolution parameter
        kwargs = {}
        if resolution is not None:
            if isinstance(resolution, (int, float)):
                resolution = (resolution, resolution)
            kwargs['resolution'] = resolution

        # Reproject to the target CRS
        print(f"Reprojecting to {target_crs}...")
        reprojected = raster.rio.reproject(
            target_crs,
            resampling=resampling_method,
            **kwargs
        )

    # Save to a new file if output path is provided
    if output_tif:
        print(f"Saving reprojected raster to {output_tif}...")
        reprojected.rio.to_raster(output_tif)

    # make sure input raster is closed
    raster.close()

    return reprojected


# %% function to clip raster (xarray.DataArray w/rioxarray) object to bbox
####################################

# this function was quickly generated by chatgpt and I modified a bit
def clip_raster_to_latlon_bbox(raster, min_lon, min_lat,
                               max_lon, max_lat, output_tif=None,
                               na_val=None):
    """
    Clip a raster to a bounding box defined in EPSG:4326 (WGS84).

    Parameters:
    - raster: xarray.DataArray or rioxarray raster with CRS defined
        or str pointing to raster to read in (e.g., .tif, .asc...)
    - min_lon, min_lat, max_lon, max_lat: float -
        bounding box coordinates in WGS84
    - output_tif : str, optional
        Path to save the reprojected GeoTIFF. If None, won't save to file.
    - na_val : numeric, list of numeric, or dict with operator, optional
        Values in the raster to be replaced with NA. Can be:
        - A single value (e.g., 255)
        - A list of values (e.g., [255, -9999])
        - A dict with operator and threshold (e.g., {">": 253}, {"==": 0})
          Supported operators: >, >=, <, <=, ==, !=
        If None, no values are replaced.

    Returns:
    - Clipped raster as an xarray.DataArray
    """
    if isinstance(raster, str) or isinstance(raster, Path):
        print(f"Reading {raster}...")
        raster = rioxarray.open_rasterio(raster)

    if na_val:
        print("applying NA filter")
        # Apply NA value filtering
        raster = apply_na_filter(raster, na_val)

    # Create bounding box geometry in EPSG:4326
    bbox_wgs84 = box(min_lon, min_lat, max_lon, max_lat)
    gdf_wgs84 = gpd.GeoDataFrame(geometry=[bbox_wgs84], crs="EPSG:4326")

    # Ensure raster has CRS and get it
    if not raster.rio.crs:
        raise ValueError("Raster must have a defined CRS via rioxarray.")

    # Reproject bounding box to match raster CRS
    print(f"Projecting to rasterios crs of {raster.rio.crs}")
    gdf_projected = gdf_wgs84.to_crs(raster.rio.crs)

    assert raster.rio.crs == gdf_projected.crs
    print("Clipping raster")
    # Clip raster using reprojected bounding box
    raster_clipped = raster.rio.clip_box(*gdf_projected.total_bounds)
        
    # raster_clipped = raster.rio.clip(
    #     gdf_projected.geometry,
    #     gdf_projected.crs,
    #     drop=False
    # )
    print("Done clipping raster")
    raster.close()
    # Save to a new file if output path is provided
    if output_tif:
        print(f"Saving clipped raster to {output_tif}...")
        raster_clipped.rio.to_raster(output_tif)

    return raster_clipped
# %% ways in which this function can be used
###############################################

# # Basic usage - reproject to WGS84 (EPSG:4326)
# reprojected = reproject_tif("input.tif", "output_wgs84.tif")

# # Reproject to a different CRS
# reprojected = reproject_tif("input.tif", "output_utm.tif", target_crs="EPSG:32633")

# # With custom resolution (e.g., 0.01 degrees for lat/lon)
# reprojected = reproject_tif("input.tif", "output_wgs84.tif", resolution=0.01)

# # Using bilinear resampling instead of nearest
# from rasterio.enums import Resampling
# reprojected = reproject_tif("input.tif", "output_wgs84.tif", 
#                            resampling_method=Resampling.bilinear)

# # Just get the reprojected array without saving
# reprojected = reproject_tif("input.tif")


# %% execut if desired
#######################################

if __name__ == "__main__":
    from pathlib import Path

    # main path holding tiffs
    DIR_TIFFS = ('C:/Projects/OWRD/PrecipFreq/GIS/CONUS404/'
                 'CONUS404_tiffs/ForUseInPFA')
    # tif file name
    # FILE_TIFF = 'CONUS404_totalWarmSeasonPrec.tif'
    FILE_TIFF = 'NORM_6190_AHM.asc'
    # path to tiff to process/reproject
    INPUT_TIFF = Path(DIR_TIFFS, FILE_TIFF)

    # where to store output tiff
    DIR_OUT = DIR_TIFFS
    # filename to give to output tiff
    # OUTPUT_FILE = FILE_TIFF.replace(".tif", "_clipped.tif")
    OUTPUT_FILE = FILE_TIFF.replace(".asc", "_test2.tif")
    # path to saved output tiff (set to None if don't want to save)
    OUTPUT_TIFF = Path(DIR_OUT, OUTPUT_FILE)  #

    # target crs for output
    TARGET_CRS = "EPSG:4326"  # can also use wkt or proj4

    # resampling method to use
    RESAMPLING_METHOD = Resampling.bilinear  # bilinear  # nearest

    # resolution (tuple: (x_res, y_res) or float if same in both directions)
    RESOLUTION = None  # maintains res if set to None

    # if would rather match an existing raster, provide path to it here
    # None for no template raster
    TEMPLATE_RASTER = None

    # bbox to clip raster to (using to get rid of odd NA values)
    # defined as min_lon, min_lat, max_lon, max_lat
    bbox_in = (-125.302, 34.229, -113.93, 50.33)

    # test1 (clip then reproject)
    # clip
    # clipped_rast = clip_raster_to_latlon_bbox(
    #     raster=INPUT_TIFF,  # proj_rast,
    #     min_lon=bbox_in[0],
    #     min_lat=bbox_in[1],
    #     max_lon=bbox_in[2],
    #     max_lat=bbox_in[3],
    #     output_tif=None  # OUTPUT_TIFF
    # )
    # # reproject
    # proj_rast = reproject_tif(
    #     input_tif=clipped_rast,  # INPUT_TIFF,
    #     output_tif=OUTPUT_TIFF,  # None,  # OUTPUT_TIFF,
    #     target_crs=TARGET_CRS,
    #     resampling_method=RESAMPLING_METHOD,
    #     resolution=RESOLUTION,
    #     template_raster=TEMPLATE_RASTER
    # )

    # test2 (reproject then clip)
    # reproject
    proj_rast = reproject_tif(
        input_tif=INPUT_TIFF,
        output_tif=None,  # OUTPUT_TIFF,
        target_crs=TARGET_CRS,
        resampling_method=RESAMPLING_METHOD,
        resolution=RESOLUTION,
        template_raster=TEMPLATE_RASTER
    )

    # clip
    clipped_rast = clip_raster_to_latlon_bbox(
        raster=proj_rast,
        min_lon=bbox_in[0],
        min_lat=bbox_in[1],
        max_lon=bbox_in[2],
        max_lat=bbox_in[3],
        output_tif=OUTPUT_TIFF
    )

