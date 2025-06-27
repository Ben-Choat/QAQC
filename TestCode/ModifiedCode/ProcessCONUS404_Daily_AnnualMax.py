'''
BChoat 2025/04/18

Script to download CONUS404 data from OSN.

This script specifically, is used to process mean, max daily precip by year

Info about OSN and other data sources here:
    https://hytest-org.github.io/hytest/dataset_catalog/README.html#storage-locations

This script is very much based on example workflow here:
    https://hytest-org.github.io/hytest/dataset_access/conus404_explore.html

'''

# %% import libs
######################################################

# NOTE: many of these libraries appear unused because they are not directly
# called, but may still need to be loaded to be accessed from other libraries
# that are directly called
import os
from pathlib import Path
import time
import fsspec
import xarray as xr
import hvplot.xarray
import holoviews as hv
import intake
import metpy
import cartopy.crs as ccrs
import rioxarray
import pandas as pd
import hvplot.pandas
import matplotlib.pyplot as plt
from CONUS404_Functions import reproject_dataarray, max_with_time

# These two should not need to be imported, but there is currently
# an issue with new zarr-3 (https://github.com/pydata/xarray/issues/10032) that
# still needs to be worked out. for now using zarr_format=2 when writing seems
# to work
# import zarr
# from numcodecs.zarr3 import Blosc

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
suffix = "_Annual_Maxima"

# %% begin exploring dataset
######################################################

# open the hytest data intake catalog
hytest_cat = intake.open_catalog(
    "https://raw.githubusercontent.com/hytest-org/hytest/main/dataset_catalog/"
    "hytest_intake_catalog.yml"
    )
list(hytest_cat)


# open the conus404 sub-catalog
cat = hytest_cat['conus404-catalog']
list(cat)

## Select the dataset you want to read into your notebook and preview its metadata
dataset = 'conus404-daily-osn'
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
print(f"Subsetting to: lat=({lats[0]}, {lats[1]}), "
      f"lon=({longs[0]}, {longs[1]})")

# use drop=True, to drop the values - not just mask them as NaN
ds_subset = ds.where(
    ((ds.lat >= lats[0])
     & (ds.lat <= lats[1])
     & (ds.lon >= longs[0])
     & (ds.lon <= longs[1])).compute(),
    drop=True
    )

# 4. Access the variable of interest as a dask array.
ds_prec = ds_subset["PREC_ACC_NC"]

# subset to warm, cool, all months
ds_prec_cool = ds_prec.sel(time=ds_prec.time.dt.month.isin(cool_months))
ds_prec_warm = ds_prec.sel(time=ds_prec.time.dt.month.isin(warm_months))


# %% Calculate mean monthly precipitation for each season for each grid cell
################################################
# NOTE: persist will make this computationally/data-transfer
# "expensive". Values will be stored locally.
# annual
start_time = time.time()

an_max_spat = (
    ds_prec
    .groupby('time.year')
    .apply(max_with_time)
    .persist()
    )

end_time = time.time()
print("Processing annual daily maxima for annual data "
      f"took {end_time-start_time:.2f} seconds")

# cool season
start_time = time.time()

cool_max_spat = (
    ds_prec_cool
    .groupby('time.year')
    .apply(max_with_time)
    .persist()
    )
end_time = time.time()
print("Processing annual daily maxima for cool season data "
      f"took {end_time-start_time:.2f} seconds")

start_time = time.time()
# warm_max_spat = (
#     ds_prec_warm
#     .groupby('time.year')
#     .max(dim='time')
#     .persist()
# )
warm_max_spat = (
    ds_prec_warm
    .groupby('time.year')
    .apply(max_with_time)
    .persist()
    )
end_time = time.time()
print("Processing annual daily maxima for warm season data "
      f"took {end_time-start_time:.2f} seconds")

da = an_max_spat.sel(year=1979)
# da = ds_prec_cool.groupby("time.year").mean(dim="time").sel(year=1999)
# da = cool_mean_spat
# da = warm_mean_spat

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
p = ax.pcolormesh(da['lon'], da['lat'], da['max_val'], cmap='YlGnBu')
ax.coastlines()
plt.colorbar(p, ax=ax, label='Precipitation (mm)')
# plt.title(f"Precipitation in {season_plot}")
plt.show()

####
# looks good - save data

# %% save to zar (uncomment to run)
# NOTE: INITIAL SAVE OF .zarr TAKES QUITE A WHILE (several minutes)
# Currently, there is an issue saving with zarr 3 (new), there is an active
# issue (https://github.com/zarr-developers/zarr-python/issues/2964)
# for now save as zarr_format=2
####################################

# annual
# an_max_spat_save = an_max_spat.drop_vars('metpy_crs')
# an_max_spat_save.attrs['crs'] = str(an_max_spat.metpy_crs)
# an_max_spat_save.to_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_annual.zarr"),
#     zarr_format=2,
#     mode='w'
# )

# # cool season
# cool_max_spat_save = cool_max_spat.drop_vars('metpy_crs')
# cool_max_spat_save.attrs['crs'] = str(cool_max_spat.metpy_crs)
# cool_max_spat_save.to_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_coolSeason.zarr"),
#     zarr_format=2,
#     mode='w'
# )

# # warm season
# warm_max_spat_save = warm_max_spat.drop_vars('metpy_crs')
# warm_max_spat_save.attrs['crs'] = str(warm_max_spat.metpy_crs)
# warm_max_spat_save.to_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_warmSeason.zarr"),
#     zarr_format=2,
#     mode='w'
# )


# %% read in files just saved as test
#######################################
from CONUS404_Functions import read_zarr_parse_coords
annual_test = xr.open_zarr(
    Path(DIR_TIFFS, "CONUS404_AnnualMaxima_annual.zarr")
    )
# ds_for_crs = Path(DIR_TIFFS, "SubdailyMaxima.zarr")
# ds_for_crs = xr.open_zarr(ds_for_crs)


################################
# read in and modify crs
# fix crs of these datasets
# annual_test = xr.open_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_annual.zarr")
#     )
# cool_test = xr.open_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_coolSeason.zarr"),
#     )
# warm_test = xr.open_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_warmSeason.zarr"),
#     )

# # update crs metadata
# annual_test.attrs['crs'] = ds_for_crs.attrs['crs_cf']
# cool_test.attrs['crs'] = ds_for_crs.attrs['crs_cf']
# warm_test.attrs['crs'] = ds_for_crs.attrs['crs_cf']

# # parse coordiantes to dataset
# annual_test = read_zarr_parse_coords(annual_test, crs_attrs='crs')
# cool_test = read_zarr_parse_coords(cool_test, crs_attrs='crs')
# warm_test = read_zarr_parse_coords(warm_test, crs_attrs='crs')


# # update attrs
# annual_test.attrs['crs'] = annual_test.rio.crs.to_wkt()
# cool_test.attrs['crs'] = cool_test.rio.crs.to_wkt()
# warm_test.attrs['crs'] = warm_test.rio.crs.to_wkt()

# # save back to crs
# annual_test.to_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_annual.zarr"),
#     zarr_format=2,
#     mode='w'
#     )
# cool_test.to_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_coolSeason.zarr"),
#     zarr_format=2,
#     mode='w'
#     )
# warm_test.to_zarr(
#     Path(DIR_TIFFS, "CONUS404_AnnualMaxima_warmSeason.zarr"),
#     zarr_format=2,
#     mode='w'
#     )
################################

# annual_test.max_time.dt.month.plot(bins=12)

# # plot
# da = annual_test.sel(year=2021)
# # da = ds_prec_cool.groupby("time.year").mean(dim="time").sel(year=1999)
# # da = cool_mean_spat
# # da = warm_mean_spat

# fig = plt.figure(figsize=(10, 8))
# ax = plt.axes(projection=ccrs.PlateCarree())
# # p = ax.pcolormesh(da['lon'], da['lat'], da['max_val'], cmap='YlGnBu')
# p = ax.pcolormesh(da['lon'], da['lat'], da['max_time'].dt.month, cmap='YlGnBu')
# ax.coastlines()
# # plt.colorbar(p, ax=ax, label='Precipitation (mm)')
# plt.colorbar(p, ax=ax, label='Month')
# # plt.title(f"Precipitation in {season_plot}")
# plt.show()

# annual_test.year.max(skipna=True)


# %% get mean maxes across all years and save as geotiff
##########################################
start_time = time.time()
an_mean_max = (
    # an_max_spat
    annual_test
    .mean(dim='year')
    .persist()
    )
end_time = time.time()
print("Processing annual mean maxima for annual data "
      f"took {end_time-start_time:.2f} seconds")

start_time = time.time()
cool_mean_max = (
    # cool_max_spat
    cool_test
    .mean(dim='year')
    .persist()
    )
end_time = time.time()
print("Processing annual mean maxima for cool season data "
      f"took {end_time-start_time:.2f} seconds")

start_time = time.time()
warm_mean_max = (
    # warm_max_spat
    warm_test
    .mean(dim='year')
    .persist()
    )
# annual
end_time = time.time()
print("Processing annual mean maxima for warm season data "
      f"took {end_time-start_time:.2f} seconds")


# %% plot for check
#######################################

da = warm_mean_max
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
p = ax.pcolormesh(da['x'], da['y'], da['max_val'], cmap='YlGnBu')
ax.coastlines()
plt.colorbar(p, ax=ax, label='Precipitation (mm)')
# plt.title(f"Precipitation in {season_plot}")
plt.show()


# %% save to zar (uncomment to run)
# NOTE: INITIAL SAVE OF .zarr TAKES QUITE A WHILE (several minutes)
# Currently, there is an issue saving with zarr 3 (new), there is an active
# issue (https://github.com/zarr-developers/zarr-python/issues/2964)
# for now save as zarr_format=2
###########################################

# annual
annual_save = reproject_dataarray(an_mean_max)
an_mean_max.rio.to_raster(
    raster_path=Path(DIR_TIFFS, "CONUS404_MeanAnnual_max_annual.tif"),
    driver="GTiff",
    tiled=True,
    compress="lzw"  # ? Recommended compression ?
)

# cool season
cool_save = reproject_dataarray(cool_mean_max)
cool_save.rio.to_raster(
    raster_path=Path(DIR_TIFFS, "CONUS404_MeanAnnual_max_coolSeason.tif"),
    driver="GTiff",
    tiled=True,
    compress="lzw"  # ? Recommended compression ?
)

# warm season
warm_save = reproject_dataarray(warm_mean_max)
warm_mean_max.rio.to_raster(
    raster_path=Path(DIR_TIFFS, "CONUS404_MeanAnnual_max_warmSeason.tif"),
    driver="GTiff",
    tiled=True,
    compress="lzw"  # ? Recommended compression ?
)

# reproject data to new crs
annual_val = reproject_dataarray(annual_test.max_val)
# cool_val = reproject_dataarray(cool_test.max_val)
# warm_val = reproject_dataarray(warm_test.max_val)

# write back to file
annual_val.rio.to_raster(
    Path(DIR_TIFFS, "CONUS404_AnnualMaxima_annual.tiff"),
    driver="GTiff",
    tiled=True,
    compress='lzw'
)
