'''
BChoat 2025/04/29

Script to loop through covariate rasters relevant to OWRD precip frequency/PMP
and reproject/clip/standardize them and save to a new location.
'''

# %% import libs/functions
##########################################

# assuming TransformRasters_wXarray_Functs is in same folder as this file
from pathlib import Path
from TransformRasters_wXarray_Functs import (
    reproject_tif,
    clip_raster_to_latlon_bbox
    )
from rasterio.enums import Resampling

# %% define dirs, vars, and such
##########################################

# base directory out of which analysis/processing is based
# WORK_DIR = Path('//westfolsom/Projects/2023/OWRD/PrecipFrequency/Data/Task2')
WORK_DIR = Path(
    '/mnt/West_Folsom_Projects/2023/OWRD/PrecipFrequency/Data/Task2'
    )

# directory holding input rasters
INPUT_DIR = WORK_DIR

# folders in INPUT_DIR that hold tiffs/ascs of interest and
# name of folder to place outputs in (matching L. Otto's previous convention)
INPUT_OUTPUT_FOLDERS = {
    # "NORM_6190_Bioclim_ASCII": "CovariateTiffs_20250429/NORM_tif",
    "NORM_6190_Monthly_ASCII": "CovariateTiffs_20250429/NORM_tif"  # ,
    # "TreeCover": "CovariateTiffs_20250429/TreeCover",
    # "rslopeaspect": "CovariateTiffs_20250429/ElevRSlopeAspect",
    # "CONUS404": "CovariateTiffs_20250429/CONUS404",
    # "HydroSheds": "CovariateTiffs_20250429/ElevRSlopeAspect"
}

# point to template raster to be used when repojrecting
TEMPLATE_RASTER = Path(WORK_DIR, "CONUS404/CONUS404_totalAnnualPrec.tif")

# bbox to clip raster to (using to get rid of odd NA values)
# defined as min_lon, min_lat, max_lon, max_lat
BBOX_IN = (-125.302, 34.229, -113.93, 50.33)

# list of file extensions to search for (include * for search pattern)
# SEARCH_PATTERNS = ['*.asc', '*.tif', '*.tiff']
# SEARCH_PATTERNS = ['*/*nlcd_tcc_conus_*-4.tif']  # just treecover
# SEARCH_PATTERNS = ['*aspect.tif', "*dem*.tif"]  # just aspsect and dem
SEARCH_PATTERNS = ['*.asc']  # just asc files (PRISM normals)

# list of files to apply mode resampling to, others will use bilinear
MODE_SAMPLE = ['aspect.tif']

# destination crs
TARGET_CRS = "EPSG:4326"  # wgs84

# define dict indicating which files include values for which NAs should be
# used in their place.
# - na_val : numeric, list of numeric, or dict with operator, optional
#     Values in the raster to be replaced with NA. Can be:
#     - A single value (e.g., 255)
#     - A list of values (e.g., [255, -9999])
#     - A dict with operator and threshold (e.g., {">": 253}, {"==": 0})
#       Supported operators: >, >=, <, <=, ==, !=
#     If None, no values are replaced.
NA_DICT = {
    'nlcd_tcc': {">": 100}  # replace 255 with NA in treecover rasters
}

# %% get list of files to be transformed
################################################
# first chck that INPUT_DIR exists
if not INPUT_DIR.exists():
    raise ValueError(f"{INPUT_DIR} does not exist.")
# list of input file names
INPUT_FOLDERS = list(INPUT_OUTPUT_FOLDERS.keys())
# list to hold all files
files_work = []

for FOLDER in INPUT_FOLDERS:
    temp_path = Path(INPUT_DIR, FOLDER)
    temp_files = [
        result for pat in SEARCH_PATTERNS for result in temp_path.glob(pat)
        ]
    files_work.extend(temp_files)


# %% execute transformations
################################################

for file in files_work:
    print(f"process: {file}")

    # get input folder
    input_folder = [f for f in INPUT_FOLDERS if f in str(file)][0]

    # get output folder to add to output dir
    OUTPUT_FOLDER = INPUT_OUTPUT_FOLDERS[input_folder]
    # define output dir
    OUTPUT_DIR = Path(WORK_DIR, OUTPUT_FOLDER)
    if not OUTPUT_DIR.exists():
        print(f'Making output dir: {OUTPUT_DIR}')
        OUTPUT_DIR.mkdir()

    # get output filename based in input file name
    OUTPUT_TIFF = Path(OUTPUT_DIR, file.name)
    # replace .asc with .tif
    if OUTPUT_TIFF.suffix == ".asc":
        OUTPUT_TIFF = Path(str(OUTPUT_TIFF).replace(".asc", ".tif"))


    if file.name in MODE_SAMPLE:
        resample_meth = Resampling.mode
    else:
        resample_meth = Resampling.bilinear

    # handle NAs if in NA_DICT
    na_list = [x for x in list(NA_DICT.keys()) if x in str(file)]
    if len(na_list) > 1:
        raise ValueError(
            f"More than one NA_DICT value matched to current file {file}."
        )
    if len(na_list) > 0:
        NA_VAL = NA_DICT[na_list[0]]
        print(f"Will use NA filter: {NA_VAL}")
    else:
        NA_VAL = None

    # clip with buffer so reprojection is faster/smoother
    buffer = 0.05  # buffer in degrees for WGS84
    # clip raster and save it
    print("Clipping with buffer w/No NA_VAL")
    clipped_rast_buffer = clip_raster_to_latlon_bbox(
        raster=file,
        min_lon=BBOX_IN[0] - buffer,
        min_lat=BBOX_IN[1] - buffer,
        max_lon=BBOX_IN[2] + buffer,
        max_lat=BBOX_IN[3] + buffer,
        output_tif=None,
        na_val=None
    )

    print("reprojecting to match template raster w/NA_VAL")
    # reproject raster
    proj_rast = reproject_tif(
        input_tif=clipped_rast_buffer,
        output_tif=None,
        target_crs=TARGET_CRS,
        resampling_method=resample_meth,
        resolution=None,
        template_raster=TEMPLATE_RASTER,
        na_val=NA_VAL
    )

    print("clipping to final domain w/NA_VAL")
    # clip raster and save it
    clipped_rast = clip_raster_to_latlon_bbox(
        raster=proj_rast,
        min_lon=BBOX_IN[0],
        min_lat=BBOX_IN[1],
        max_lon=BBOX_IN[2],
        max_lat=BBOX_IN[3],
        output_tif=OUTPUT_TIFF,
        na_val=NA_VAL
    )

# %% close any open arrays
#########################################

    clipped_rast_buffer.close()
    proj_rast.close()
    clipped_rast.close()