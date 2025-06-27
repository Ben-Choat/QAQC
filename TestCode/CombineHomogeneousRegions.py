'''
BChoat 2024/05/07

Using the homogenous regions that Caileen Y. identified using
GEV distribution paramters, here we combine huc8's into single
files and save .shp files.
'''

# %% Import libraries
#################################

import geopandas as gpd
import pandas as pd

# %% define dirs, vars, etc.
#######################################

# directory holding huc8s and csv associating each huc8 with a region
# output shape files will be saved here as well
dir_in = 'C:/Projects/OWRD/PrecipFreq/HomogenousRegions_GEV'

# shape file name
shp_in = 'ClippedHUC8Shapefile.shp'

# csv name
csv_in = 'HUC8s_Iter11.csv'



# %% read in data and prep
#######################################

# shape file
gdf_work = gpd.read_file(f'{dir_in}/{shp_in}',
                        dtype = {'HUC8': str})

# csv
df_work = pd.read_csv(f'{dir_in}/{csv_in}',
                        dtype = {'HUC8': str})

# join into one gdf
gdf_work = pd.merge(gdf_work, df_work, on = 'HUC8')

# convert to nad83 to match other spatial data being used for max-stable work
gdf_work = gdf_work.to_crs('epsg:4269')


# %% loop through each unique region, merge polygons, and save a shape file
#######################################

for region in gdf_work.SortedRegion.unique():
    print(f'\n working with {region}')

    # subset to region of interest
    gdf_temp = gdf_work.query("SortedRegion == @region")
    region_name = region.replace(' ', '_')
    gdf_temp.to_file(f'{dir_in}/HomogeneousRegion_{region_name}.shp')


