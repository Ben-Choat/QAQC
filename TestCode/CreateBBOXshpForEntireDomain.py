'''
BChoat 2024/04/30

Script to create bounding box for entire domain (based on gauges)
'''


# %% import libraries
####################################################

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import box
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


# %% define directories, variables, and such
#####################################################

# directory/file holding data with lat/lons of sites
file_in = 'P:/2023/OWRD/PrecipFrequency/Data/QC_AMPSite_all_v12c_annual/1hour/1hourInvData.csv'

# directory/file to write subset of data to
file_out = 'P:/2023/OWRD/PrecipFrequency/MaxStable/GIS/BBOX_EntireDomain.shp'

# states shape file for plotting
states = 'C:/Projects/GIS_General/states_21basic/states.shp'
st_shp = gpd.read_file(states)

# %% create bbox based on gauge locations and add buffer
##########################################

# define buffer 
# a value in decimal degrees by which bbox will be extended in all directions
buffer_in = 0.5

df_in = pd.read_csv(file_in)
df_in = df_in.drop_duplicates(['SiteID'])

bbox_in = [df_in.Long.min(), df_in.Lat.min(), df_in.Long.max(), df_in.Lat.max()]

bbox_plt = box(bbox_in[0] - buffer_in, 
               bbox_in[1] - buffer_in, 
               bbox_in[2] + buffer_in,
               bbox_in[3] + buffer_in)
gdf_bbox = gpd.GeoDataFrame(geometry=[bbox_plt], crs="EPSG:4326")

fig, ax = plt.subplots()
gdf_bbox.boundary.plot(ax = ax, color = 'orange')
st_shp.boundary.plot(ax = ax)

gdf_bbox.to_file(file_out)


