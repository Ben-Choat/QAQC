'''
BChoat 2025/04/28

Some hopefully helpful short example applications and guidance using various 
xarray, dask, and custom functions to work with .zarr and CONUS404 data.

NOTE: I've created a module CONUS_Functions that includes some helper functions
    related to GIS operations on xarray.DataArrays/DataSets and general functions
    for computation on those same object types.


The code hear may not be immediately executable since I point to some dummy
directories and files.
'''


# %% calculate maximum of various durations
############################################

'''
Use ProcessCONUS404_MaxStorms.py defining input parameters near the bottom under
if __name__ == "__main__".

I've provided some commentary there as well for guidance.
'''

# pretty much all of this functionality relies on:
# not direclty used here, but all DataArrays/DataSets are xarray
import xarray as xr

# %% read data in from .zarr and assign coordinates to dataset and dataarrays
# in the dataset.
############################################
from pathlib import Path
from CONUS404_Functions import read_zarr_parse_coords

# NOTE: crs_attrs points to attribute name for which the crs is stored as
    # metadata. If it is not == "crs_cf" you need to adjust this input.

DIR_TIFFS = Path("C:/Projects/OWRD/PrecipFreq/GIS/CONUS404")
ds_test = read_zarr_parse_coords(
    Path(DIR_TIFFS, "SupraDailyMaxima.zarr"),
    crs_attrs='crs_cf',
    out_crs_format='wkt'
    # out_crs_format='proj4'
    )


# %% reproject raster to new crs
# NOTE: this can be done witih xarray and rasterio, but this function adds
# some additional help in case crs is not clearly defined - 
# I recommend just using this function when reprojecting
##################################################
from CONUS404_Functions import reproject_dataarray

# make sure to pass a dataarray and not a dataset
da = reproject_dataarray(ds_test.max_1_hour, crs_target="EPSG:4326")

# %% resample a dataarray to match resolution of a raster or different
# dataarray
##################################################
from CONUS404_Functions import resample_dataarray_to_match

an_mean_resampled = resample_dataarray_to_match(
    ds_test.max_1_hour,
    "C:/Projects/OWRD/PrecipFreq/ExplRasters/treecover_average.tif",
    resampling_method='bilinear'
)


# %% resample a dataarray that has datetimes as values in each 'cell'
# helpful for plotting but I don't think you can save to .tiff for example
##################################################
from CONUS404_Functions import reproject_dataarray_w_datetime

da = reproject_dataarray_w_datetime(
    ds_test.max_1_hour,
    crs_target='EPSG:4326'
)


# %% idenfity the maximum value from a time series of data array data, and return
# the max values and the times associated with the max values
# e.g., I used this to find max 1-day precip and the time they occurred
###################################################
from CONUS404_Functions import max_with_time

an_max_spat = (
    ds_test.max_1_hour
    .groupby('time.year')
    .apply(max_with_time)  # here is the max_with_time function being applied
    .persist()
)


# %% similar to max_with_time, but find the maximum precip using rolling sums
# through time. e.g., I used this to find 1, 2, 3, 4, and 5 day maximas and 
# associated datetimes.
# see ProcessCONUS404_MaxStorms.py for an example application
###################################################
from CONUS404_Functions import rolling_sum_max_and_time

da_out = rolling_sum_max_and_time(
                    data=ds_test.max_1_hour,
                    window_duration=3,
                    time_label='hr'
            )

# %% save a xarray.DataArray to .tiff
###########################################
import rioxarray

# subset DataSet to DataArray
da = ds_test.max_1_hour

# reproject to desired crs
da = reproject_dataarray(da, crs_target="EPSG:4326")  # could be wkt crs
da.rio.to_raster(
    raster_path=Path("path/to/target/location.tiff"),
    driver="GTiff",
    tiled=True,
    compress="lzw"  # ? Recommended compression ?
)


# %% read in .zarr and save to .nc
###############################################

# read in
ds = xr.open_zarr('path/to/data.zarr')

ds.to_netcdf('path/to/target/data.nc')