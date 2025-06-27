'''
BChoat 2025/04/18

Script to download CONUS404 data from OSN.

This script specifically, is used to process max storms of various durations
across the entire period of record

Info about OSN and other data sources here:
    https://hytest-org.github.io/hytest/dataset_catalog/README.html#storage-locations

This script is very much based on example workflow here:
    https://hytest-org.github.io/hytest/dataset_access/conus404_explore.html

    
NOTE: (some notes here - bchoat - delete later)
zarr does not exlpicitly account for crs. Metpy, can be used to parse
crs and will infer what the rpoper crs is depneding on the coordiantes in the
.zarr.

Once parsed, it still doesn't really seem to mean much, other than having 
information stored.

read in zarr -> parse_cf with metpy -> store in .attrs for any subsequently
created .zarrs. 

Create function to read in .zarrs and convert to .tiff using the crs in the
attrs.
'''

# %% import libs
######################################################

# NOTE: many of these libraries appear unused because they are not directly
# called, but may still need to be loaded to be accessed from other libraries
# that are directly called
import os
from pathlib import Path
import logging
import time
import fsspec
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
from dask.distributed import Client, LocalCluster
from CONUS404_Functions import rolling_sum_max_and_time

os.environ['USE_PYGEOS'] = '0'
# Configure logging to write to a file
logging.basicConfig(
        # filename='ProcessCONUS404_MaxStorms_SupraDaily.log',
        filename='ProcessCONUS404_MaxStorms_SubDaily.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
        )


# %% define main code to be executed below
#########################################

def main(dir_tiffs: str,
         out_filename: str,
         durs_in: list,
         time_label_in: str,
         n_workers: int,
         memory_limit_in: str,
         threads_per_worker: int):
    '''
    when using dask cluster, it expects code to be called from within
    if __name__ == "__main__": call, so here we define bulk of code as a
    function and call it at bottom of file

    inputs:
    - dir_tiffs (Path object or str): Main directory where to store outputs
    - out_filename (str): string of file name to use when saving output
        e.g., "SubdailyMaxima.zarr"
    - durs_in (list): list of integers specifying for which durations maxima
        should be identified
    - time_label_in (str): string indicating the unit of measurement for
        durations - must be "hr", "hour" or "day"
    - n_workers (int): number of processors for building dask cluster
    - memory_limit_in (str): specify limit for memory per worker, e.g., "4GB"
    - threads_per_worker (int): number of threads to use per worker

    outputs:
    - saves .zarr with maxima for each duration in dir_tiffs
    '''
    # create directory to store outputs if not already exists
    if not dir_tiffs.exists():
        dir_tiffs.mkdir(parents=True)

    lats = (34.229, 50.33)  # for PMP - entire domain
    longs = (-125.302, -113.93)  # for PMP - entire domain
    
    logging.info("will save outputs to %s", dir_tiffs)


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
    # list(cat)

    hytest_cat.close()

    # Select the dataset you want to read into your notebook and preview its metadata
    if time_label_in in ["hr", "hour"]:
        dataset = 'conus404-hourly-osn'
    elif time_label_in == "day":
        dataset = 'conus404-daily-osn'
    else:
        raise ValueError("time_label_in must be 'hr', 'hour', or 'day', "
                         f"but you provided {time_label_in}")
    # cat[dataset]

    # read in the dataset and use metpy to parse the crs information on the dataset
    # print(f"Reading {dataset} metadata...", end='')
    start_time = time.time()
    logging.info("Begin reading %s metadata", dataset)

    ds = cat[dataset].to_dask().metpy.parse_cf()

    end_time = time.time()
    logging.info(
        "Reading %s took %.2f seconds", dataset, end_time - start_time
    )

    logging.info(
            "Subsetting to: lat=(%.2f, %.2f), lon=(%.2f, %.2f)",
            lats[0], lats[1], longs[0], longs[1]
        )
    start_time = time.time()
    # subset ds to area of interest
    # use drop=True, to drop the values - not just mask them as NaN
    ds_subset = ds.where(
        ((ds.lat >= lats[0])
         & (ds.lat <= lats[1])
         & (ds.lon >= longs[0])
         & (ds.lon <= longs[1])).compute(),
        drop=True
    )
    end_time = time.time()
    logging.info(
        "Subsetting to domain took %.2f seconds", end_time - start_time
        )

    # time_chunk_size = round(size_time/max(durs_in)/100)
    # shooting for chunk sizes between 100MB and 1GB
    time_chunk_size = 366  # 732  # for hourly 366 results in 223 MiB chunks

    # Apply the chunking
    ##########
    # NOTE: Negative defines number of chunks (+) is size of chunks
    ds_prec = (
        ds_subset["PREC_ACC_NC"]
        .chunk({'time': time_chunk_size, 'x': -1, 'y': -1})
        )

    logging.info(
            "Begin looping through durations and calculating rolling sums"
            )

    results = {}
    with LocalCluster(
        n_workers=n_workers,  # 30, 4
        threads_per_worker=threads_per_worker,  # 2
        memory_limit=memory_limit_in,  # '4GB',  # or '4GB', etc.
        nanny=True  # auto restart closed workers
    ) as cluster, Client(cluster) as client:
        script_start_time = time.time()

        logging.info("Dask cluster started.")
        logging.info(client)

        start_time = time.time()
        logging.info("Begin persisting data within cluster")
        # persist data within cluster
        ds_prec = ds_prec.persist()
        end_time = time.time()
        logging.info(
            "Persisting data within cluster took %.2f seconds",
            end_time-start_time
        )

        start_time = time.time()

        for dur in durs_in:
            logging.info("begin processing %d %s", dur, time_label_in)
            results[f"{dur}_{time_label_in}"] = rolling_sum_max_and_time(
                    data=ds_prec,
                    window_duration=dur,
                    time_label=time_label_in
            ).compute()
        
        end_time = time.time()
        logging.info(
                "Calculating maxima of all durations took %.2f seconds",
                end_time-start_time
        )

    start_time = time.time()
    # merge into single dataset
    combined_dataset = xr.merge(results.values())
    end_time = time.time()
    logging.info(
            "Merging results from all durations took: "
            " %.2f seconds", end_time-start_time
            )

    # %% Save to a .zarr file
    ###############################################

    start_time = time.time()

    output_path = Path(dir_tiffs, out_filename)

    # store crs in atttrs
    crs_out = ds.metpy_crs.item().to_dict()
    combined_dataset.attrs['crs_cf'] = crs_out

    # clean up a bit
    combined_dataset_save = combined_dataset.drop_vars('metpy_crs')

    ###################
    # make sure have uniform chunk sizes (except last chunk)
    combined_dataset_save = combined_dataset_save.chunk({
        "y": 100, "x": 110  # based on Dimensions: y: 484 x: 330
    })
    combined_dataset_save.to_zarr(
        output_path,
        zarr_format=2,
        mode="w"
    )
    print(f"Combined dataset saved to {output_path}")

    end_time = time.time()
    logging.info(
            "Saving to .zarr took %.2f seconds",
            end_time-start_time
            )

    script_end_time = time.time()
    logging.info(
            "The entire script took "
            "%.2f seconds", script_end_time - script_start_time
            )


if __name__ == "__main__":
    # %% define any input vars
    ######################################################

    # NOTE:
    # These input parameters such as N_WORKERS and MEMORY_LIMIT_IN were at
    # least somewhat optimized for Shared3 where we have 36 total processors
    # and ~120 GB of RAM.

    # If applying on a different machine these will likely need to be adjusted.
    # time_chunk_size within the script is currently set to chunksize of 366
    # in time dimension, and x and y are set to have only 1 chunk (i.e., using -1).
    # These values may also need to be changed if working on a different machine.

    # On shared 3, processing hourly data to get 1, 3, 6, and 12 hour maxes took
    # about 2.5 hours. Processing daily data to get 1, 2, 3, 4, and 5 day
    # maxes took about 0.5 hours.

    #####
    # subdaily inputs
    # TIME_LABEL_IN = 'hr'  # used in var labels and some conditionals
    # DURS_IN = [1, 3, 6, 12]  # what durations to id maxes for
    # OUT_FILENAME = "SubdailyMaxima.zarr"  # name of file to save outputs to
    # N_WORKERS = 20  # number of processors to use
    # MEMORY_LIMIT_IN = '6GB'  # how much memory (RAM) to allow per worker
    # THREADS_PER_WORKER = 1  # how many threads per worker
    #####

    #####
    # >=1day intputs
    TIME_LABEL_IN = 'day'  # used in var labels and some conditionals
    DURS_IN = [1, 2, 3, 4, 5]  # what durations to id maxes for
    OUT_FILENAME = "SupraDailyMaxima.zarr"  # name of file to save outputs to
    N_WORKERS = 30  # number of processors to use
    MEMORY_LIMIT_IN = '4GB'  # how much memory (RAM) to allow per worker
    THREADS_PER_WORKER = 1  # how many threads per worker
    #####

    # directory where to save output ratsers
    DIR_TIFFS = Path(
        "/mnt/West_Folsom_Projects/2023/OWRD/PrecipFrequency/"
        "Data/Task2/CONUS404"
    )

    # execute script
    main(
        dir_tiffs=DIR_TIFFS,
        out_filename=OUT_FILENAME,
        durs_in=DURS_IN,
        time_label_in=TIME_LABEL_IN,
        n_workers=N_WORKERS,
        memory_limit_in=MEMORY_LIMIT_IN,
        threads_per_worker=THREADS_PER_WORKER
    )


# %% read in and plot as a test
#####################################

# DIR_TIFFS = Path("C:/Projects/OWRD/PrecipFreq/GIS/CONUS404")

from CONUS404_Functions import (
    read_zarr_parse_coords,
    reproject_dataarray_w_datetime,
    reproject_dataarray
    )
from rasterio.enums import Resampling

# read in make sure crs assigned
ds_test = read_zarr_parse_coords(
    Path(DIR_TIFFS, "SubdailyMaxima.zarr"),
    # Path(DIR_TIFFS, "SupraDailyMaxima.zarr"),
    crs_attrs='crs_cf',
    out_crs_format='wkt'
    # out_crs_format='proj4'
    )

dict_hold = ds_test.attrs['crs_cf']

ds_test.attrs = {x: dict_hold[x] for x in dict_hold.keys()}

# ds_test.to_netcdf(Path(DIR_TIFFS, "SupraDailyMaxima.nc"))
# ds_test2 = xr.open_dataset(Path(DIR_TIFFS, "SupraDailyMaxima.nc"))
ds_test.to_netcdf(Path(DIR_TIFFS, "SubDailyMaxima.nc"))
ds_test2 = xr.open_dataset(Path(DIR_TIFFS, "SubDailyMaxima.nc"))
ds_test2.max_6_hr.plot()

# # save to raster as test
# # test_save = reproject_dataarray(ds_test.max_1_day)
# # test_save.rio.to_raster(
# #     raster_path=Path(DIR_TIFFS, "TESTOUT.tif"),
# #     driver="GTiff",
# #     tiled=True,
# #     compress="lzw"  # ? Recommended compression ?
# # )

# # plot
# # da = ds_test.max_4_day
# # da = reproject_dataarray(da)
# # da = da.rio.reproject("EPSG:4326", resampling=Resampling.nearest)
# da = ds_test.time_max_4_day
# da = reproject_dataarray_w_datetime(da, crs_target='EPSG:4326')
# fig = plt.figure(figsize=(10, 8))
# # lambert or wgs84
# ax = plt.axes(projection=ccrs.PlateCarree())
# # ax = plt.axes(projection=ccrs.LambertConformal())
# # time or value
# # p = ax.pcolormesh(da.coords['x'], da.coords['y'], da, cmap='YlGnBu')
# p = ax.pcolormesh(da.coords['x'], da.coords['y'], da.dt.year, cmap='YlGnBu')
# ax.coastlines()
# # plt.colorbar(p, ax=ax, label='Precipitation (mm)')
# plt.colorbar(p, ax=ax, label='Month')
# # plt.title(f"Precipitation in {season_plot}")
# plt.show()
