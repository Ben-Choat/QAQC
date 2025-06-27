import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.mask import mask
from shapely.geometry import box
import geopandas as gpd
import numpy as np
from rasterio.io import MemoryFile
from pathlib import Path
from glob import glob

# # Define input and output file paths
# src_path = "input_lambert.asc"
# dst_path = "output_wgs84.tif"

# main working directory
MAIN_DIR = Path('//westfolsom/Projects/2023/OWRD/PrecipFrequency/Data/Task2')

# Specific folders in which to search for files
SEARCH_FOLDERS = [
    "NORM_6190_Bioclim_ASCII",
     'TreeCover'
]

names = ['NORM_6190_Tave_wt.asc', 'NORM_6190_PPT_wt.asc']

files_work = glob(f"{MAIN_DIR}/*/{names[0]}")

name = names[0]

src_path = Path(
    '//westfolsom/Projects/2023/OWRD/PrecipFrequency/Data/'
    f'Task2/NORM_6190_Bioclim_ASCII/{name}.asc'
    )
# dst_path = Path(
#     '//westfolsom/Projects/2023/OWRD/PrecipFrequency/Data/'
#     f'Task2/NORM_reprojected//{name}.tif'
# )
dst_path = Path(
    f'C:/Projects/OWRD/{name}.tif'
)

# Define the bounding box in WGS84 (min_lon, min_lat, max_lon, max_lat)
wgs84_bbox = [-125.302, 34.229, -113.93, 50.33]

# Open the source dataset
with rasterio.open(src_path) as src:
    # Calculate the transformation parameters for reprojection
    transform, width, height = calculate_default_transform(
        src.crs,                  # Source CRS (Lambert Conformal Conical)
        {'init': 'EPSG:4326'},    # Destination CRS (WGS84)
        src.width,                # Source width
        src.height,               # Source height
        *src.bounds               # Source bounds
    )
    
    # Set up the parameters for the reprojected raster
    kwargs = src.meta.copy()
    kwargs.update({
        'crs': {'init': 'EPSG:4326'},
        'transform': transform,
        'width': width,
        'height': height,
        'driver': 'GTiff'  # Needed for MemoryFile
    })
    
    # Create a memory file for the reprojected data
    with MemoryFile() as memfile:
        # Open as a dataset and write the reprojected data
        with memfile.open(**kwargs) as mem_dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(mem_dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs={'init': 'EPSG:4326'},
                    resampling=Resampling.bilinear
                )
        
        # Read from the memory file for clipping
        with memfile.open() as mem_src:
            # Create a geometry from the WGS84 bbox
            bbox_polygon = box(*wgs84_bbox)
            
            # Create a GeoDataFrame with the bbox geometry
            bbox_gdf = gpd.GeoDataFrame({'geometry': [bbox_polygon]}, crs="EPSG:4326")
            
            # Clip the raster using the bbox
            clipped_data, clipped_transform = mask(
                mem_src,
                bbox_gdf.geometry,
                crop=True,
                all_touched=True
            )
            
            # Update the metadata for the clipped data
            clipped_meta = mem_src.meta.copy()
            clipped_meta.update({
                'height': clipped_data.shape[1],
                'width': clipped_data.shape[2],
                'transform': clipped_transform
            })
            
            # Save the clipped data directly to the final output file
            with rasterio.open(dst_path, 'w', **clipped_meta) as dst:
                dst.write(clipped_data)
