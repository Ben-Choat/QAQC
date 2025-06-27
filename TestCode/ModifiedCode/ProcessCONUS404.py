'''
BChoat 2025/04/18

Script to download CONUS404 data from OSN.

Info about OSN and other data sources here:
    https://hytest-org.github.io/hytest/dataset_catalog/README.html#storage-locations

This script is very much based on example workflow here:
    https://hytest-org.github.io/hytest/dataset_access/conus404_explore.html

'''

# %% import libs
######################################################

import os
from pathlib import Path
# import fsspec
import xarray as xr
import hvplot.xarray
import holoviews as hv
import intake
import metpy
import cartopy.crs as ccrs
# import zarr
import rioxarray
import pandas as pd
import hvplot.pandas  
import matplotlib.pyplot as plt

os.environ['USE_PYGEOS'] = '0'
# from dask.distributed import Client, LocalCluster

# %% define any input vars
######################################################

# directory where to save output ratsers
DIR_TIFFS = Path("C:/Projects/OWRD/PrecipFreq/GIS/CONUS404")

lats = (34.229, 50.33)  # for PMP - entire domain
longs = (-125.302, -113.93)  # for PMP - entire domain

# lists of months indicating season(s) they belong to as integers
cool_months = [10, 11, 12, 1, 2, 3, 4]
warm_months = [5, 6, 7, 8, 9]
annual_months = list(range(1, 13))


# %% begin exploring dataset
######################################################

# open the hytest data intake catalog
hytest_cat = intake.open_catalog(
    "https://raw.githubusercontent.com/hytest-org/hytest/main/dataset_catalog/hytest_intake_catalog.yml"
    )
list(hytest_cat)


# open the conus404 sub-catalog
cat = hytest_cat['conus404-catalog']
list(cat)

# osn: open storage network (OSN) Pod
# onprem-hw: USGS on-premises supercomputer
# ba: bias-adjusted
# pgw: psuedo global warming
#      (one storage system for the Tallgrass/Denali supercomputers and another
#      for the Hovenweep supercomputer)
#      only accessible to USGS employees or collaborators who have been granted
#      access to USGS supercomputers.


## Select the dataset you want to read into your notebook and preview its metadata
dataset = 'conus404-monthly-osn'
cat[dataset]

# read in the dataset and use metpy to parse the crs information on the dataset
print(f"Reading {dataset} metadata...", end='')
ds = cat[dataset].to_dask().metpy.parse_cf()

# print([i for i in ds.data_vars])

# ds.SNOW
# ds.PREC_ACC_NC

# crs = ds['SNOW'].metpy.cartopy_crs
# crs

# subset ds to area of interest
# use drop=True, to drop the values - not just mask them as NaN
# ds_subset = ds.where(
#     ((ds.lat >= lats[0])
#      & (ds.lat <= lats[1])
#      & (ds.lon >= longs[0])
#      & (ds.lon <= longs[1])).compute(),
#     drop=True
#     )
ds_subset = ds.sel(
    lat=slice(lats[0], lats[1]),
    lon=slice(longs[0], longs[1])
)
print(f"Subsetting to: lat=({lats[0]}, {lats[1]}), "
      f"lon=({longs[0]}, {longs[1]})")

# 4. Access the variable of interest as a dask array.
ds_prec = ds_subset["PREC_ACC_NC"]

# subset to warm, cool, all months
ds_prec_cool = ds_prec.sel(time=ds_prec.time.dt.month.isin(cool_months))
ds_prec_warm = ds_prec.sel(time=ds_prec.time.dt.month.isin(warm_months))

# calculate monthly precip for each season
# ds_prec_cool["PREC_ACC_NC"].sum()
# using .persist() will cause info to be stored locally - This is a
# good idea when the values will be referenced more than once - It
# allows you for example to use an_mean.values and have data instantly
# available.
an_mean = (
    ds_prec
    .groupby("time.year")
    .mean(dim=('time', 'x', 'y'))
    .persist()
)
cool_mean = (
    ds_prec_cool
    .groupby("time.year")
    .mean(dim=('time', 'x', 'y'))
    .persist()
)
warm_mean = (
    ds_prec_warm
    .groupby("time.year")
    .mean(dim=('time', 'x', 'y'))
    .persist()
)

# PLOT TO CHECK
# Convert each DataArray to a pandas-friendly format
an_mean_box = an_mean.hvplot.box(y='PREC_ACC_NC', title='an_mean')
cool_mean_box = cool_mean.hvplot.box(y='PREC_ACC_NC', title='cool_mean')
warm_mean_box = warm_mean.hvplot.box(y='PREC_ACC_NC', title='warm_mean')

# convert to dataframe for plotting
df_annual = an_mean.to_dataframe(name="precip")
df_annual['season'] = "annual"
df_cool = cool_mean.to_dataframe(name="precip")
df_cool['season'] = 'cool'
df_warm = warm_mean.to_dataframe(name="precip")
df_warm['season'] = 'warm'

df_all = pd.concat([df_annual, df_cool, df_warm])
# Plot with grouping by season
df_all.hvplot.box(
    y='precip',
    by='season',
    title='Total precip over given season - plotted by season'
    )

# looks reasonable to me

# %% get spatial values
################################################
# NOTE: persist will make this computationally/data-transfer
# "expensive". Values will be stored locally.
an_mean_spat = (
    ds_prec
    .mean(dim='time')
    .persist()
)
cool_mean_spat = (
    ds_prec_cool
    .mean(dim='time')
    .persist()
)
warm_mean_spat = (
    ds_prec_warm
    .mean(dim='time')
    .persist()
)


season_plot = "cool"  # cool, warm, or annual
das = {
    "cool": cool_mean_spat,
    "warm": warm_mean_spat,
    "annual": an_mean_spat
}
da = das[season_plot]
# da = ds_prec_cool.groupby("time.year").mean(dim="time").sel(year=1999)
# da = cool_mean_spat
# da = warm_mean_spat

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
p = ax.pcolormesh(da['lon'], da['lat'], da, cmap='YlGnBu')
ax.coastlines()
plt.colorbar(p, ax=ax, label='Precipitation (mm)')
plt.title(f"Precipitation in {season_plot}")
plt.show()

#### 
# looks good - save data

# %% Save the results as GeoTIFF rasters using rioxarray.
#    -  Ensure the output directory exists.
# os.makedirs(out_dir, exist_ok=True)

# import rioxarray
# from pyproj import CRS

def reproject_dataarray(daskArrayIn):
    '''
    reprojrect daskarray to another crs
    '''  
  
    # Get your pyproj CRS from MetPy
    crs = daskArrayIn.metpy.pyproj_crs  # da = your DataArray

    # Write the CRS using rioxarray
    daskArrayIn.rio.write_crs(crs, inplace=True)

    # Reproject to EPSG:4326 (WGS84)
    daskArrayIn_wgs84 = daskArrayIn.rio.reproject("EPSG:4326")

    return daskArrayIn_wgs84


an_mean_wgs84 = reproject_dataarray(an_mean_spat)
cool_mean_wgs84 = reproject_dataarray(cool_mean_spat)
warm_mean_wgs84 = reproject_dataarray(warm_mean_spat)

an_mean_wgs84.rio.to_raster(
    Path(DIR_TIFFS, "CONUS404_MeanAnnualPrec.tif"), driver="GTiff"
    )
cool_mean_spat.rio.to_raster(
    Path(DIR_TIFFS, "CONUS404_MeanCoolSeasonPrec.tif"), driver="GTiff"
    )
warm_mean_spat.rio.to_raster(
    Path(DIR_TIFFS, "CONUS404_MeanWarmSeasonPrec.tif"), driver="GTiff"
    )
print(f"Saved mean and max rasters to {DIR_TIFFS}")

# 8. Close the dataset.
# ds.close()