"""
Created on Tue Sep 24 09:20:53 2024

@author: cyu

Modified by BChoat on 2025/04/16
"""
# import os  # , sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from scipy.stats import genextreme
from rasterio.enums import Resampling
import numpy as np
import xarray as xr
# import rioxarray

# %% Plot the return period of the PMP values.

# strings to id scenario (used in all output files)
duration = "1day"
domain = "EntireDomain"
season = "all"
version = "v12"
regr_model = "2step"  # spatGEV or 2step
fit_metr = "1se"  # best or 1se
minimum_por = "POR60"  # "" for none

indir = Path("C:/Projects/OWRD/PrecipFreq/PredictedGEVParams_NewExpl"
             f"/MarginalGEVs/{duration}_{season}{version}_{minimum_por}"
             f"_{fit_metr}")
# loc_fn = Path(f"{indir}/gev_loc_2step_EntireDomain_all_3day_v12.tif")
loc_fn = Path(f"{indir}/gev_loc_{regr_model}_{domain}_{season}_"
              f"{duration}_{version}.tif")
scale_fn = Path(f'{indir}/gev_scale_{regr_model}_{domain}_{season}_'
                f'{duration}_{version}.tif')
shape_fn = Path(f'{indir}/gev_shape_{regr_model}_{domain}_{season}_'
                f'{duration}_{version}.tif')
# gen_outdir = "P:/2023/OWRD/PMP/Data/PFE_Point"
gen_outdir = Path("C:/Projects/OWRD/PrecipFreq/PredictedGEVParams_NewExpl"
                  "/MarginalGEVs")
if not Path.exists(gen_outdir): Path.mkdir(gen_outdir)
outdir = Path(gen_outdir, indir.name)
if not Path.exists(outdir): Path.mkdir(outdir)


pmp_fn = Path("P:/2023/OWRD/PMP/Task2_HMRMethod/"
              "GTF_to_PMP_HUC8Method_ProposedRegionsKeithEdit/"
              "UseWRF/304.8m_6deg_gtf1.5_RP1000/UseMinGTF0.1/"
              "10_3D_PMPEst_FocusArea.nc")
outfig_RP = Path(f'{gen_outdir}/'
                 f"PMP_RP_{duration}_{domain}_{season}_{version}_"
                 f"{regr_model}_{minimum_por}{fit_metr}_FocusArea.png"
                 )


# %%Open tif files
loc_ds = xr.open_dataset(loc_fn)
loc_ds = loc_ds.rename({'band_data': 'gev_loc'})

scale_ds = xr.open_dataset(scale_fn)
scale_ds = scale_ds.rename({'band_data': 'gev_scale'})

shape_ds = xr.open_dataset(shape_fn)
shape_ds = shape_ds.rename({'band_data': 'gev_shape'})

params = xr.Dataset()
params['gev_scale'] = scale_ds.gev_scale
params['gev_shape'] = shape_ds.gev_shape*-1  # use the -1* shape since scipy.stats.genextreme uses a different sign convention than the SpatialExtremes R Package used to generate these results
params['gev_loc'] = loc_ds.gev_loc
params = params.isel(band=0)
param_rio = params.copy()
param_rio = param_rio.rio.write_crs(scale_ds.spatial_ref.crs_wkt)


pmp = xr.open_dataset(pmp_fn)
pmp.rio.write_crs("EPSG:2992", inplace=True)
pmp_reproj = pmp.rio.reproject(
    scale_ds.spatial_ref.crs_wkt,
    Resampling=Resampling.nearest
    )

# param_rio = params.copy()
# param_rio = param_rio.rio.write_crs(scale_ds.spatial_ref.crs_wkt) 
params_reproj = params.rio.reproject_match(
    pmp_reproj, Resampling=Resampling.nearest
    )
# NOTE: xarray/pandas "where" keeps values where first condition is
# True and replaces other values with second argument (e.g., np.nan)
params_reproj['gev_shape'] = params_reproj.gev_shape.where(
    params_reproj.gev_shape != 3.4028234663852886e+38, np.nan
    )

pmp_vals = pmp_reproj.pmp
yy_cdf = genextreme.cdf(
    pmp_vals,
    params_reproj.gev_shape,
    params_reproj.gev_loc,
    params_reproj.gev_scale
    )
rp_yy = 1/(1-yy_cdf)

pmp_vals_rp = xr.DataArray(
            rp_yy, coords=[('y', pmp_reproj.y.values),
                           ('x', pmp_reproj.x.values)]
                    ).to_dataset(name='pmp_RP')

# save pmp to tiff
pmp_vals_rp.rio.write_crs("EPSG:4326", inplace=True)
# make sure coords are interpreted correctly
# pmp_vals_rp.pmp_RP["x"] = ("x", pmp_vals_rp.coords["x"].values)
# pmp_vals_rp.pmp_RP["y"] = ("y", pmp_vals_rp.coords["y"].values)
# pmp_vals_rp['pmp_RP'] = pmp_vals_rp.pmp_RP.where(~np.isinf(pmp_vals_rp.pmp_RP), np.nan)
pmp_vals_rp['pmp_RP'].values[np.isinf(pmp_vals_rp['pmp_RP'].values)] = np.nan
pmp_vals_rp.pmp_RP.rio.to_raster(
    str(outfig_RP).replace(".png", ".tif"),
    driver="GTiff"
    )

# plot
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
im = pmp_vals_rp.pmp_RP.plot(ax=ax, norm=mpl.colors.LogNorm())
# pmp_vals_rp.pmp_RP.plot(norm=mpl.colors.LogNorm(), vmin=650, vmax=1e12)
# ax.imshow(rp_yy, norm = mpl.colors.LogNorm()); #,vmin = 1, vmax = 10000,
# fig.colorbar(im,ax = ax, label = 'Return Period (yrs)')
ax.set_title(f"Return Period of 10 sq mi {duration} PMP Values")
ax.gridlines()

ax.coastlines()
ax.add_feature(cartopy.feature.STATES)

fig.savefig(
    outfig_RP,
    dpi=300,
    bbox_inches='tight',  # Crop the figure tightly to avoid unnecessary white space
)
