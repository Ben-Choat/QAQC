'''
BChoat 2024/03/21

Script to identify good bounding box to use for developing a small
dataset with which to learn and test code.

Will be used to learn and develop workflow for OWRD precip freq. 
analysis with max-stable.
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
file_in = 'P:/2023/OWRD/PrecipFrequency/Data/QC_AMPSite_all_v17_annual/1hour/1hourInvData.csv'

# directory/file to write subset of data to
file_out = 'P:/2023/OWRD/PrecipFrequency/MaxStable/Learning/ToyDataSites.shp'

# states shape file for plotting
states = 'C:/Projects/GIS_General/states_21basic/states.shp'
st_shp = gpd.read_file(states)
# %% read in data and visualize it
##########################################

df_in = pd.read_csv(file_in)
df_in = df_in.drop_duplicates(['SiteID'])
df_in = df_in.query("maxYear - minYear >= 20")
df_in = df_in.query("maxYear > 2000")

# convert to geodataframe using lat lon as geometry
geometry = [Point(xy) for xy in zip(df_in['Long'], df_in['Lat'])]
gdf_in = gpd.GeoDataFrame(df_in, geometry=geometry, crs="EPSG:4326")
gdf_in = gdf_in[['Station ID', 'SiteID', 'Long', 'Lat', \
                 'Name', 'Elevation', 'maxYear', 'minYear', 'geometry']]

# gdf_in.plot(markersize=1)
map = gdf_in.explore(marker_kwds=dict(radius=0.1))
map

# define bbox based on map
bbox_in = [-119.0, 44.26, -117.9, 44.8]
bbox_plt = box(bbox_in[0], bbox_in[1], bbox_in[2], bbox_in[3])
gdf_bbox = gpd.GeoDataFrame(geometry=[bbox_plt], crs="EPSG:4326")
gdf_bbox.explore(m=map, style_kwds=dict(fill_alpha=0))
gdf_bbox.to_file(file_out.replace(".shp", "_BBOX.shp"))

# us coordinate indexer (cx) to subset points to bbox
gdf_sub = gdf_in.cx[bbox_in[0]:bbox_in[2], bbox_in[1]:bbox_in[3]]

gdf_sub.explore()

gdf_sub.shape

gdf_sub.to_file(file_out)
gdf_sub.to_csv(file_out.replace('.shp', '.csv'), index=False)


# %% make plot to save
################################
fig, ax = plt.subplots(1, 1, # figsize=(4, 9), 
                       subplot_kw={'projection': ccrs.PlateCarree()})
extent = gdf_in.total_bounds
extent = [extent[0], extent[2], extent[1], extent[3]]
ax.set_extent(extent, ccrs.PlateCarree())
st_shp.boundary.plot(ax=ax)
gdf_in.plot(ax=ax, markersize=1, color='blue')
gdf_sub.plot(ax=ax, markersize=1, color='red')

plt.savefig(file_out.replace(".shp", ".png"), dpi=300)