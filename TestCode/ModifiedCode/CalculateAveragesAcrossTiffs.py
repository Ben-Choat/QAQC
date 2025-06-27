'''
BChoat 2025/05/01

Script to calculate cell by cell averages across several tifs located
in a single directory.

Asked Claude 3.7 to write, and made small modifications
'''
import numpy as np
import rioxarray
import xarray as xr
from pathlib import Path


def average_tifs(input_dir, output_path=None, verbose=True):
    """
    Calculate the cell-by-cell average across multiple TIF files in a
    directory.

    Parameters:
    -----------
    input_dir : str
        Path to directory containing TIF files
    output_path : str, optional
        Path to save the output TIF file. If None, the result is returned
        but not saved.
    verbose : bool, default=True
        Whether to print information about the process

    Returns:
    --------
    xarray.DataArray
        The averaged raster as a DataArray
    """
    # Find all TIF files in the directory
    tif_files = list(Path(input_dir).glob("*.tif"))

    if len(tif_files) == 0:
        raise ValueError(f"No TIF files found in {input_dir}")

    if verbose:
        print(f"Found {len(tif_files)} TIF files to process")

    # Open the first file to use as a reference
    reference = rioxarray.open_rasterio(tif_files[0])

    # Initialize array to store all data
    all_data = []

    # Check each file and collect data
    for i, file_path in enumerate(tif_files):
        if verbose and i % 10 == 0:
            print(f"Processing file {i+1}/{len(tif_files)}: {file_path.name}")

        # Open the raster
        raster = rioxarray.open_rasterio(file_path)

        # Verify the format matches our reference
        if not _check_raster_compatibility(
                reference, raster, file_path.name, verbose
        ):
            continue

        # Add to our collection
        all_data.append(raster)

    if len(all_data) == 0:
        raise ValueError("No compatible TIF files were found")

    # Stack and average all rasters
    if verbose:
        print(f"Calculating average of {len(all_data)} rasters...")

    # Combine into a dataset along a new dimension and calculate mean
    stacked = xr.concat(all_data, dim="file")
    average = stacked.mean(dim="file")

    # Save the result if output path is provided
    if output_path is not None:
        if verbose:
            print(f"Saving result to {output_path}")
        average.rio.to_raster(output_path)

    if verbose:
        print("Processing complete!")

    return average


def _check_raster_compatibility(reference, raster, filename, verbose=True):
    """
    Check if a raster has the same format as the reference raster.

    Parameters:
    -----------
    reference : xarray.DataArray
        The reference raster
    raster : xarray.DataArray
        The raster to check
    filename : str
        Name of the file being checked (for error messages)
    verbose : bool
        Whether to print warnings

    Returns:
    --------
    bool
        True if compatible, False otherwise
    """
    # Check dimensions
    if raster.rio.shape != reference.rio.shape:
        if verbose:
            print(f"Warning: {filename} has different dimensions and will be skipped")
        return False

    # Check CRS
    if raster.rio.crs != reference.rio.crs:
        if verbose:
            print(f"Warning: {filename} has different CRS and will be skipped")
        return False

    # Check transform/resolution
    if not np.allclose(raster.rio.transform(), reference.rio.transform()):
        if verbose:
            print(
                f"Warning: {filename} has different transform/resolution and "
                "will be skipped"
                )
        return False

    # Check nodata value
    if raster.rio.nodata != reference.rio.nodata:
        if verbose:
            print(f"Warning: {filename} has different nodata value. "
                   "Will proceed but this may affect results.")

    return True


if __name__ == "__main__":

    #####################
    # this chunk enables execution from command line using e.g.,
    # python average_tifs.py path/to/tif/directory --output output_average.tif
    # import argparse

    # parser = argparse.ArgumentParser(
    #     description="Calculate average of multiple TIF files"
    #     )
    # parser.add_argument("input_dir", help="Directory containing TIF files")
    # parser.add_argument("--output", "-o", help="Output file path")
    # parser.add_argument(
    #     "--quiet", "-q", action="store_true", help="Suppress verbose output"
    #     )

    # args = parser.parse_args()

    # average_tifs(args.input_dir, args.output, verbose=not args.quiet)
    #####################

    #####################
    # this chunk was used to get average treecover for OWRD pfa work
    MAIN_DIR = '/mnt/West_Folsom_Projects/2023/OWRD/PrecipFrequency/Data/Task2/'
    INPUT_PATH = Path(MAIN_DIR, 'CovariateTiffs_20250429/TreeCover')
    OUTPUT_PATH = Path(INPUT_PATH, "treecover_average.tif")

    average_tc = average_tifs(
        input_dir=INPUT_PATH,
        output_path=OUTPUT_PATH
    )