'''
BChoat 2025/01/29

Donny K. delineated many watersheds based on OWRD dam locations.
There are 6 files of varying levels of how 'nested' the basins in them
were.

I'm using my environment genenv311.
This repo shold be close:
https://github.com/Ben-Choat/PythonEnvFiles/blob/main/genenv311.yml

There are many columns that should not be needed.

This script is being written to read those files in and clean them.
- remove uneeded columns
- maybe project to Oregon Lambert coord system?
- maybe combine into single shape file?
    Did combine into single file
- maybe add new column indicating significant overlap with other basins
    e.g., > 98% overlap and > 98% same basin area?
    used interception-over-union > 0.98
- recalculate area? (areas were practically identical after recalcing)
TODO:
- add column or another way ofindicating whether the basin should be simulated
    or if another basins results should be used for it
- ID duplicated geometries for points so don't replicate
- ID duplicatd nidids and hand appropriately
'''

# %% import libraries
#######################################

from pathlib import Path
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, Point
import numpy as np

# %% define working dirs, vars, etc.
#######################################

# dir holding original shape files
DIR_SHP = Path('E:/Projects/OWRD_PrecipFreq/TestData/Final_Output_20250127')

# dir where to place output files
DIR_OUT = Path('E:/Projects/OWRD_PrecipFreq/Priority_Basins')
if not DIR_OUT.is_dir():
    DIR_OUT.mkdir()

# list of shapefile names to read in and clean
FILE_NAMES = [
    'sub10sqmi_as_points.shp',  # unlike the rest, this holds points
    'unnest_basin_1.shp',
    'unnest_basin_2.shp',
    'unnest_basin_3.shp',
    'unnest_basin_4.shp',
    'unnest_basin_5.shp',
    'unnest_basin_6.shp'
]

# define list of column names to keep
COLS_KEEP = [
    'nidid', 'Area_sqmi', '_merge', 'geometry'
    ]
# other potential columns of interest
# 'dam_name_a', 'dam_height', 'dam_heig_1', 'nidHeightI',
# 'nidHeight',
# 'max_storag', 'nid_storag', 'nid_stor_1', 'max_stor_1',
# 'normal_s_1',
# 'volume', 'surfaceAre', 'drainageAr', 'maxDischar', 'hydraulicH',
# 'structural',  'damLength',
# 'location_a', 'latitude_a', 'longitude_', 'latitude_n',
# 'longitude_',
#  'location_n',
#  'hazard_all', 'river_all_',
#  'source_age', 'owner_name', 'owner_na_1', 'is_active_',
#  'statutory_',
#  'yearComple',  'dataUpdate',
#  'usaceDivis', 'usaceDistr', 'privateDam',
#  'huc2', 'huc4', 'huc6', 'huc8', 'femaRegion',


# %% function to take input poly and get overlap with all polys in
# another geodataframe
#######################################

def find_matching_polygon(target_geometry: Polygon | Point | MultiPolygon,
                          comparison_gdf: gpd.GeoDataFrame,
                          return_field: str,
                          threshold: float = 0.98) -> list:
    """
    Find index of matching polygon in comparison_gdf that matches
    target_geometry.

    Args:
        target_geometry: The polygon or point to match
        comparison_gdf: GeoDataFrame containing polygons or points
            to compare against
        return_field: name of column with id to return
        threshold: Similarity threshold (0-1), higher means more similar

    Returns:
        List of length 2;
            1. Index of matching polygon or 'None' if no match found
            2. the iou score
    """

    # make sure com
    # make sure comparison_gdf has reset index for referencing
    comparison_gdf = comparison_gdf.reset_index(drop=True)

    # if Polygon, calculate intersection over union for each polygon
    if isinstance(target_geometry, (Polygon, MultiPolygon)):
        if not any(comparison_gdf.geometry.intersects(target_geometry)):
            return ['None', 0]
        intersections = (comparison_gdf
                         .geometry
                         .intersection(target_geometry)
                         .area)
        unions = comparison_gdf.geometry.union(target_geometry).area

        # Calculate IoU scores
        iou_scores = intersections / unions

        # Find best match above threshold
        best_match_idx = np.argmax(iou_scores)
        # get best 'score' i.e., overlapping portion
        best_iou_score = iou_scores[best_match_idx]
        # get best iou 'ID'
        best_iou_id = comparison_gdf.loc[best_match_idx, return_field]

        if best_iou_score >= threshold:
            return [best_iou_id, best_iou_score]
        return ['None', 0]

    if isinstance(target_geometry, Point):
        gdf_points = comparison_gdf[
            comparison_gdf['geometry'] == target_geometry
            ].reset_index(drop=True)
        if gdf_points.shape[0] > 0:
            if gdf_points.shape[0] > 1:
                print(
                    'Multiple matches of the point geometry, using the first.'
                    )
            best_id = gdf_points.at[0, return_field]

            return [best_id, 1]
        else:
            return ['None', 0]

    raise ValueError('target_geometry must be shapely Point or Polygon')


# %% read in shape files and store in dict
#######################################

# define dict and keynames based on input files
shp_dict = {x.replace('.shp', ''): [] for x in FILE_NAMES}

# read in shape files into dict
for i, key in enumerate(shp_dict.keys()):
    print(i, key)
    # read in shapefile
    gdf_temp = gpd.read_file(Path(DIR_SHP, FILE_NAMES[i]))
    # subset to specific columns
    gdf_temp = gdf_temp.loc[:, COLS_KEEP]
    # convert to Oregon Lamber crs
    gdf_temp = gdf_temp.to_crs('epsg:2992')  # OREGON_LAMBERT
    # calculate area based on polygons in sqmi
    gdf_temp['Area_sqmi_poly'] = round(gdf_temp.geometry.area / (5280**2), 2)
    # add column indicating which 'nested' or shp file data is from
    gdf_temp['origin_file'] = FILE_NAMES[i].replace('.shp', '')
    # store in dict
    shp_dict[key] = gdf_temp

# combine into single geodataframe
gdf_combined = pd.concat(shp_dict, ignore_index=True)

# plot
# gdf_combined.plot(column='origin_file')
# gdf_combined.plot(column='source')
# gdf_combined.source.hist()

# %% correct owrd nidids for duplicates for points only
# keep larger area where dulpicates exist - smaller values were just
# artifacts from processing
##########################################

gdf_work = shp_dict[list(shp_dict.keys())[0]]
gdf_work['Area_sqmi'] = pd.to_numeric(gdf_work['Area_sqmi'])
nidid_duplicates = gdf_work.nidid.duplicated()
nidid_dupl_ids = gdf_work.loc[nidid_duplicates, 'nidid']
for id_in in nidid_dupl_ids:
    # get largest area of duplicated nidid
    gdf_temp = gdf_work.query("nidid == @id_in")
    gdf_temp = gdf_temp.sort_values(by='Area_sqmi', ascending=False)
    gdf_temp = gdf_temp.reset_index(drop=True)
    series_temp = gdf_temp.iloc[0]

    # series_temp = gdf_temp.loc[gdf_temp['Area_sqmi'].idxmax()]
    gdf_temp = pd.DataFrame(series_temp).transpose()

    # remove from gdf_work
    gdf_work = gdf_work.query("nidid != @id_in")
    # add largest area back
    gdf_work = pd.concat([gdf_work, gdf_temp])

gdf_work = gdf_work.reset_index(drop=True)

nidid_duplicates = gdf_work.nidid.duplicated()
nidid_dupl_ids = gdf_work.loc[nidid_duplicates, 'nidid']
if len(nidid_dupl_ids) != 0:
    raise ValueError('Still duplicates present!')
else:
    # update shp_dict0
    shp_dict[list(shp_dict.keys())[0]] = gdf_work
    # combine into single geodataframe again
    gdf_combined = pd.concat(shp_dict, ignore_index=True)

# %% investigate overlap using intersection over union scores
##########################################
# add 'intersected_area' column to geodataframe to hold results
gdf_combined['SimBasin'] = 'None'
gdf_combined['iou_score'] = 0

for i, row_in in gdf_combined.iterrows():
    nidid_in = row_in['nidid']
    # only look forward since function is symmetric
    # temp_df = gdf_combined.iloc[i:gdf_combined.shape[0], :]
    temp_df = gdf_combined[gdf_combined['nidid'] != nidid_in]
    input_geom = row_in['geometry']
    intersect_basins = find_matching_polygon(
                        target_geometry=input_geom,
                        comparison_gdf=temp_df,
                        return_field='nidid',
                        threshold=0.98
                        )

    gdf_combined.loc[
        gdf_combined['nidid'] == nidid_in, ['SimBasin', 'iou_score']
        ] = [intersect_basins[0], intersect_basins[1]]

print('complete')

# %% for areas where pracitally the same basins were identified
# handle by keeping the indicator for the smaller basin, since
# its results will be used for the larger basin
##########################################
# get rows w/overlapping basins id'd
gdf_temp = gdf_combined.query("SimBasin != 'None'")
original_Nrow = gdf_temp.shape[0]

# loopthrough nidid's in subset and modify smaller basin
if original_Nrow > 0:
    for nid in gdf_temp['nidid'].unique():
        print(nid)
        gdf_work = gdf_temp.query("nidid == @nid or SimBasin == @nid")
        gdf_work = gdf_work.sort_values(by='Area_sqmi', ascending=False)
        if gdf_work.shape[0] == 1:
            continue  # assume similar already removed
        # nid_keep = gdf_work.at[gdf_work.index[0], 'nidid']
        # defining to include 1+ but should just be 1
        nid_drop = gdf_work.at[gdf_work.index[0], 'nidid']

        # replace values in SimBasin as appropriate
        gdf_combined.loc[
            gdf_combined['nidid'] == nid_drop, 'SimBasin'
            ] = 'None'

# check to make sure results are as expected
gdf_temp = gdf_combined.query("SimBasin != 'None'")
new_Nrow = gdf_temp.shape[0]
if original_Nrow > 0 and new_Nrow != original_Nrow/2:
    raise ValueError('Expected half of original rows to be reset to None',
                     ' but this was not observed.')


# %% update gdf_combined
##########################################
# rename _merge and update values
gdf_combined = gdf_combined.rename({'_merge': 'source'}, axis=1)
gdf_combined['source'] = (gdf_combined['source']
                          .replace({
                              'both': 'NID_OWRD',
                              'all_state': 'OWRD',
                              'ntid': 'NID'
                              }))

# add a unique identifer column that will also be used when users provide
# new areas
gdf_combined['ID'] = gdf_combined['source'] + '_' + gdf_combined['nidid']

# reorder columns
gdf_combined = gdf_combined[[
    'ID', 'nidid', 'Area_sqmi', 'SimBasin',
    'iou_score', 'origin_file', 'source', 'geometry'
]]
# update column names to be friendly with 10 or fewer characters (max for .shp)
gdf_combined.columns = ['ID', 'nidid', 'Area_sqmi', 'SimBasin',
                        'iou_score', 'nest_lev', 'source', 'geometry']

# perform final check for duplicates
nidid_duplicates = gdf_combined.nidid.duplicated()
nidid_dupl_ids = gdf_combined.loc[nidid_duplicates, 'nidid']
if len(nidid_dupl_ids) != 0:
    print(gdf_combined[gdf_combined['nidid'].isin(nidid_dupl_ids)])
    raise ValueError('Still duplicates present!')

# %% get into final format and save
##########################################

# # save points
gdf_point = gdf_combined[gdf_combined['geometry'].type == 'Point']
gdf_point = gdf_point.sort_values(by='nidid').reset_index(drop=True)
# modify nidids for duplicated
gdf_point.to_file(
    Path(DIR_OUT, 'OWRD_priority_basins_points.shp')
)
# # save polys
gdf_poly = gdf_combined[gdf_combined['geometry'].type != 'Point']
gdf_poly = gdf_poly.sort_values(by='nidid').reset_index(drop=True)
gdf_poly.to_file(
    Path(DIR_OUT, 'OWRD_priority_basins_polys.shp')
)

# %% scratch code to inspect a bit
#############################

# gdf_combined[gdf_combined['SimBasin'] != 'None']
# gdf_combined[gdf_combined['nidid'] == 'OR03739']
# gdf_combined[gdf_combined['nidid'] == 'OR00053']
# df_test = shp_dict[list(shp_dict.keys())[1]]
# df_test = df_test.loc[:, COLS_KEEP]
# df_test.loc[:, ['nidHeight', 'nidHeightI', 'dam_height', 'dam_heig_1']]
# df_test.loc[df_test['dam_heig_1'] != df_test['dam_height'],
#                       ['dam_height', 'dam_heig_1']]


# df_test.loc[:, ['max_storag', 'nid_storag', 'nid_stor_1', 'max_stor_1',
#                           'normal_s_1']]

# df_test.loc[:, ['owner_name', 'owner_na_1', ]]

# %% check for duplicated  nidids
##########################################
# bool_duplicated = gdf_combined.nidid.duplicated()
# tot_duplicated = sum(bool_duplicated)
# id_duplicated = gdf_combined.loc[bool_duplicated, 'nidid']
# gdf_duplicated = gdf_combined[gdf_combined['nidid'].isin(id_duplicated)]
# gdf_duplicated = gdf_duplicated.sort_values(by='nidid')
# gdf_duplicated.shape

# geom_dupl = gdf_combined.geometry.duplicated()
# id_geom_dupl = gdf_combined.loc[geom_dupl, 'nidid']
# gdf_geom_dupl = gdf_combined[gdf_combined['nidid'].isin(id_geom_dupl)]
# gdf_geom_dupl.shape
# gdf_geom_dupl.sort_values(by='nidid')
# gdf_geom_dupl.groupby('nidid').count()
# # save to file
# # gdf_duplicated.to_file()


# points_inv = gdf_combined[gdf_combined['geometry'].type == 'Point']
# points_inv.query("SimBasin != 'None'")