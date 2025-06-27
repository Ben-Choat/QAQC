import os
os.environ['USE_PYGEOS'] = '0'

import fsspec
import xarray as xr
import hvplot.xarray
import intake
import metpy
import cartopy.crs as ccrs

# open the hytest data intake catalog
hytest_cat = intake.open_catalog("https://raw.githubusercontent.com/hytest-org/hytest/main/dataset_catalog/hytest_intake_catalog.yml")
list(hytest_cat)

# open the conus404 sub-catalog
cat = hytest_cat['conus404-catalog']
list(cat)

## Select the dataset you want to read into your notebook and preview its metadata
dataset = 'conus404-hourly-osn' 
# dataset = 'conus404-daily-osn' 
# dataset = 'conus404-monthly-osn' 
cat[dataset]

print(f"Reading {dataset} metadata...", end='')
ds = cat[dataset].to_dask().metpy.parse_cf()
ds