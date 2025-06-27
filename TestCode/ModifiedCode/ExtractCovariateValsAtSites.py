'''
BChoat 2025/05/01

Script to extract covariate values at specific lat-longs
'''

from pathlib import Path
import rasterio as rio
from pyproj import Transformer
import pandas as pd
import numpy as np


def extract_raster_values(raster_path, points_df, input_crs="EPSG:4326",
                          lat_col="latitude", lon_col="longitude"):
    """
    Extract raster values at given coordinates.

    Parameters:
    - raster_path (str):
        Path to the raster file (.tif)
    - input_crs (str): 
        CRS of input coordinates (default is 'EPSG:4326' for WGS84)
    - points_df : pandas.DataFrame
        DataFrame containing point coordinates
    - lat_col : str, default "latitude"
        Name of column containing latitude values
    - lon_col : str, default "longitude"
        Name of column containing longitude values
    - crs : str or rasterio.crs.CRS, default "EPSG:4326"
        Coordinate reference system of the input points

    Returns:
    - List of raster values at the given points
    """
    # Create a copy of the input DataFrame to avoid modifying the original
    result_df = points_df.copy()

    # Convert the DataFrame to a GeoDataFrame
    # geometry = [
    #     Point(xy) for xy in zip(result_df[lon_col], result_df[lat_col])
    #     ]
    coords = [(x, y) for (x, y) in zip(result_df[lon_col], result_df[lat_col])]

    # gdf = gpd.GeoDataFrame(result_df, geometry=geometry, crs=input_crs)

    with rio.open(raster_path) as src:
        # Transform coordinates to raster CRS if needed
        if input_crs != src.crs.to_string():
            transformer = Transformer.from_crs(
                input_crs, src.crs, always_xy=True
                )
            coords_proj = [
                transformer.transform(lon, lat) for lon, lat in coords
                ]
        else:
            coords_proj = coords

        # Sample the raster
        values = np.array(list([float(x[0]) for x in src.sample(coords_proj)]))

        # Replace nodata values with ""
        if src.nodata is not None:
            values[values == src.nodata] = np.nan

    return values


# Example usage
if __name__ == "__main__":

    # path to csv with coords
    COORD_FILE = Path(
        # "//westfolsom/Projects/2023/OWRD/PrecipFrequency/Data/"
        "/mnt/West_Folsom_Projects/2023/OWRD/PrecipFrequency/Data/"
        "StationInventory_all_wo6_fix3/Site13WY_v1All.csv"
        )

    # name of columsn with latitude and longitude
    LAT_COL = "Lat"
    LONG_COL = "Long"

    TIF_DIR = Path(
        '/mnt/West_Folsom_Projects/2023/OWRD/PrecipFrequency'
        '/Data/Task2/CovariateTiffs_20250429'
    )
    # TIF_DIR = Path(
    #         '//westfolsom/Projects/2023/OWRD/PrecipFrequency/Data/'
    #         'Task2/CovariateTiffs_20250429'
    #     )

    # destination at which to save output
    CSV_OUT = Path(TIF_DIR, "CoVDataFrame_wCONUS404.csv")

    # read in coord data
    df_coords = pd.read_csv(COORD_FILE)

    # Path to your raster file
    INPUT_FOLDERS = [
        "NORM_tif", "TreeCover",
        "ElevRSlopeAspect", "CONUS404"
    ]
    if not TIF_DIR.exists():
        raise ValueError(f"{TIF_DIR} does not exist.")

    # list to hold all files
    files_work = []

    for FOLDER in INPUT_FOLDERS:
        temp_path = Path(TIF_DIR, FOLDER)
        temp_files = temp_path.glob("*.tif")
        files_work.extend(temp_files)

    # define dict to hold outputs
    # result_dict = dict.fromkeys([x.stem for x in files_work])
    result_dict = {
        "SiteID": list(range(df_coords.shape[0])),
        "Long": [val for val in df_coords['Long']],
        "Lat": [val for val in df_coords['Lat']]
    }

    for file in files_work:
        RASTER_FILE = file

        # Extract values
        result = extract_raster_values(
            raster_path=RASTER_FILE,
            points_df=df_coords,
            input_crs="EPSG:4326",
            lat_col=LAT_COL,
            lon_col=LONG_COL
        )

        # update output dict
        result_dict[file.stem] = result


    # convert results_dict to dataframe and save
    df_result = pd.DataFrame(result_dict)

    # replace nans with blanks
    df_result = df_result.replace(np.nan, "")

    # clean up column names
    df_result.columns = df_result.columns.str.replace("NORM_6190", "Site")
    df_result.columns = df_result.columns.str.replace(
        "nlcd_tcc_conus", "treecover"
        )
    df_result.columns = df_result.columns.str.replace("v2021-4", "Site")
    df_result.columns = df_result.columns.str.replace("aspect", "aspect_Site")
    df_result.columns = df_result.columns.str.replace(
        "na_dem_15s1", "elevation_Site"
        )
    df_result.columns = df_result.columns.str.replace(
        "profilecurvature", "profilecurvature_Site"
        )
    df_result.columns = df_result.columns.str.replace(
        "slope", "slope_Site"
        )
    df_result.columns = df_result.columns.str.replace(
        "tangentialcurvature", "tangentialcurvature_Site"
        )

    # update siteID to start at 1 instead of 0
    df_result['SiteID'] = df_result['SiteID'] + 1

    # save to csv
    df_result.to_csv(CSV_OUT, index=False)
