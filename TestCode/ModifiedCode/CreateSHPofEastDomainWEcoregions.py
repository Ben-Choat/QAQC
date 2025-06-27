'''
BChoat 2025/04/18

This script combines our east OR domain, with chosen ecoregion level-iiis
based on transposable regions identified in the PMP work in cooperation
with OWRD (Keith)

Using genenv311 conda environment

NOTE:
- from PMP work:
- North cascades and Cascades regions were merged
- Blue Mountains region was split into two regions (a more central Oregon
    region and northeastern region)
- Columbia Plateau
- Blue Mountains North
- Northern Basin and Range

'''

# %% import libs
#################################################

from pathlib import Path
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon


# %% define dirs, vars, and such
#################################################

# name of output shape file
SHP_OUT = "EastOR_Extended_EcoRegions_bbxd.shp"
# "EastOR_Extended_EcoRegions_bbox_clip.shp"

# main dir holding input files
INPUT_DIR = Path("E:/Projects/OWRD_PrecipFreq/GIS")

# list of filenames to be read in
FILES_IN = [
    "CascadesEastOregon_OverlapWest_Outer.shp",
    "ca_eco_l3/ca_eco_l3.shp",
    "nv_eco_l3/nv_eco_l3.shp",
    "wa_eco_l3/wa_eco_l3.shp",
    "or_eco_l3/or_eco_l3.shp",
    "id_eco_l3/id_eco_l3.shp"
]

# list of ecoregions that we want to keep
# - North cascades and Cascades regions were merged
# - Blue Mountains region was split into two regions (a more central Oregon
#     region and northeastern region)
# - Columbia Plateau
# - Blue Mountains North
# - Northern Basin and Range
ECO_LIST = [
    "OregonEast_Defined",  # region we defined as east Oregon (not ecoregion)
    "Blue Mountains",
    "Northern Basin and Range",
    "Columbia Plateau",
    # "Northern Rockies",  # too far north I think
    "Idaho Batholith",
    "Snake River Plain",
    "Eastern Cascades Slopes and Foothills",
    # "Cascades",
    "Central Basin and Range"
]

# define outer bounaries as bounding box
# min_lon = -122.66
# max_lon = -114.22
# min_lat = 36.15
# max_lat = 49
min_lon = -122.66
max_lon = -114.6995
min_lat = 40.8191
max_lat = 46.8913
# %% start reading in and processing
#################################################

#########
# code from development process
# test1 = gpd.read_file(
#     Path(INPUT_DIR, FILES_IN[0])
# )
# test2 = gpd.read_file(
#     Path(INPUT_DIR, FILES_IN[1])
# )
# test1.head()
# test2.head()
#########

# define empty dict to hold spatial data
dict_all = {}

for i, file in enumerate(FILES_IN):
    print(i, file)
    temp_gdf = gpd.read_file(Path(INPUT_DIR, file))
    if i == 0:
        temp_gdf = temp_gdf[["OBJECTID", "geometry"]]
        temp_gdf['OBJECTID'] = "OregonEast_Defined"
        temp_gdf = temp_gdf.rename({"OBJECTID": "Name"}, axis=1)
    else:
        temp_gdf = temp_gdf[["US_L3NAME", "geometry"]]
        temp_gdf = temp_gdf.rename({"US_L3NAME": "Name"}, axis=1)
        temp_gdf = temp_gdf[temp_gdf['Name'].isin(ECO_LIST)].reset_index(drop=True)
    temp_gdf = temp_gdf.to_crs("epsg:4326")
    dict_all[Path(file).name.replace(".shp", "")] = temp_gdf

gdf_all = pd.concat(dict_all).reset_index().drop("level_1", axis=1)
gdf_all = gdf_all.rename({"level_0": "region"}, axis=1)

# define empty dict to hold spatial data
dict_all = {}

for region in ECO_LIST:
    temp_gdf = gdf_all[gdf_all['Name'] == region]
    merged_geom = temp_gdf['geometry'].unary_union
    merged_gdf = gpd.GeoDataFrame(
        pd.DataFrame({"Name": region}, index=[0]),
        geometry=[merged_geom],
        crs=gdf_all.crs
        )
    dict_all[region] = merged_gdf

gdf_all = pd.concat(dict_all).reset_index(drop=True)

gdf_all.plot(column="Name", legend=True)

# gdf_all.to_file(Path(INPUT_DIR, SHP_OUT))


# create a bounding box and clip gdf_all to it
bbox_polygon = Polygon([
    (min_lon, min_lat), (max_lon, min_lat),
    (max_lon, max_lat), (min_lon, max_lat),
    (min_lon, min_lat)
])

gdf_box = gpd.GeoDataFrame(geometry=[bbox_polygon], crs="epsg:4326")
gdf_all_clipped = gpd.overlay(gdf_all, gdf_box, how="intersection")
gdf_all_clipped = gdf_all_clipped[~gdf_all_clipped['Name'].isna()]
gdf_all_clipped.plot(column="Name", legend=True)

# gdf_all_clipped.to_file(Path(INPUT_DIR, SHP_OUT.replace(".shp", "_clip.shp")))
gdf_all_clipped.to_file(Path(INPUT_DIR,"CascadesEastOregon_Extended_smaller.shp"))

