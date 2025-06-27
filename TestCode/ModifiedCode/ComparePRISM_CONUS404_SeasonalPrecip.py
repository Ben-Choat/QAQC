'''
BChoat 2025/04/21

Compare seasonal precip between CONUS404 and PRISM
'''

# %% load libs
#############################################

import xarray as xr
import rioxarray
from rasterio.enums import Resampling
from pathlib import Path
import CONUS404_Functions.GIS_Functs as GF
import matplotlib.pyplot as plt

# %% define input vars, dirs, and such
####################################################
# Define your lat/lon bounding box
lat_min, lat_max = 34.229, 50.33
lon_min, lon_max = -125.302, -113.93

# directory where to save output ratsers
DIR_CONUS404 = Path("C:/Projects/OWRD/PrecipFreq/GIS/CONUS404")
DIR_PRISM = Path("C:/Projects/OWRD/PrecipFreq/ExplRasters")

# input rasters
raster1_file = Path(DIR_CONUS404, "CONUS404_totalCoolSeasonPrec_PRISM_Seasons.tif")
    # CONUS404_totalWarmSeasonPrec_PRISM_Seasons.tif
raster2_file = Path(DIR_PRISM, "NORM_6190_PPT_wt.tif")  # NORM_6190_PPT_sm.tif

# where to save raster data
raster_out1 = Path(DIR_CONUS404, "WinterPrecipDifference_CONUS_PRISM.tif")
raster_out2 = Path(DIR_CONUS404, "WinterPrecipPercentDiff_CONUS_PRISM.tif")
# plot title for plotting difffernce between two rasters:
plot_title1 = "Dec Through Feb. Precip\nCONUS404 minus PRISM"
plot_title2 = "Dec Through Feb. Precip\nPercent Difference (CONUS404 vs PRISM)"


# %% load data and preprocess
##################################################
# Open both rasters
raster1 = rioxarray.open_rasterio(raster1_file, masked=True).squeeze()
raster2 = rioxarray.open_rasterio(raster2_file, masked=True).squeeze()

# Reproject both to EPSG:4326 (if not already)
raster1 = raster1.rio.reproject("EPSG:4326", resampling=Resampling.bilinear)
raster2 = raster2.rio.reproject("EPSG:4326", resampling=Resampling.bilinear)

# Clip to bounding box
raster1_clipped = raster1.rio.clip_box(minx=lon_min, miny=lat_min,
                                       maxx=lon_max, maxy=lat_max)
raster2_clipped = raster2.rio.clip_box(minx=lon_min, miny=lat_min,
                                       maxx=lon_max, maxy=lat_max)

# Align the rasters if needed
raster2_aligned = raster2_clipped.rio.reproject_match(raster1_clipped)


# %% difference
###########################################################
# Subtract raster1 from raster2
difference = raster1_clipped - raster2_aligned

# plot difference
fig, ax = plt.subplots()
difference.plot(
    ax=ax,
    cbar_kwargs={"label": "Difference (mm)"} 
    )
ax.set_title(label=plot_title1)
# Save the result (optional)
difference.rio.to_raster(raster_out1)

# Show quick summary
# print(difference)

# %% percent difference
###########################################################

perc_difference = difference/((raster1_clipped+raster2_aligned)/2) * 100


# plot difference
fig, ax = plt.subplots()
perc_difference.plot(
    ax=ax,
    cbar_kwargs={"label": "% Difference calculated as: \ndifference/mean (%)"}
    )
ax.set_title(label=plot_title2)
# Save the result (optional)
difference.rio.to_raster(raster_out2)
