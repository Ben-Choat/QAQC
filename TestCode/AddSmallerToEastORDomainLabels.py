'''
BChoat 2025/06/23

We modified the "extended" east OR domain to be smaller due to issues with
fitting the bootstrapped models.

As such, the scripts will look for the correct domain name in files and
filenames.

This script updates existing anisotropy filenames to included the _smaller
suffix so that the bs_fit_models.R file will find them.

It also creates a .csv file logging which files are being modified.
'''

# %% import libraries
##################################

from pathlib import Path

# %% define dirs, vars and such
#####################################

# main working directory
DIR_WORK = Path(
    "//westfolsom/Projects/2023/OWRD/PrecipFrequency/"
    "MaxStable/Application/Anisotropy"
    )

# txt filename that will hold record of files being modified
FILE_LOG = "ModifiedFileNames_East_smaller.txt"


# %% perform
#######################################

# list to hold info for FILE_LOG
log_list = ["BChoat 2025/06/23 ",
            "This file holds filenames for which anisotropy results were ",
            "produced using EastOR_extended, but filename was changed to ",
            "include _smaller suffix so subsequent scripts can find them.",
            "We are assuming here that anistropy results will not differe ",
            "significantly between the _smaller, and previous domains."]


# get folders in DIR_WORK
folders_work = [f for f in DIR_WORK.glob("*") if f.name.endswith("v12")]
print(f'folders found are: {folders_work}')

with open(Path(DIR_WORK, FILE_LOG), 'w', encoding='utf-8') as f:
    for item in log_list:
        f.write(f"{item}\n")

    for folder in folders_work:
        print(f'processing folder: {folder} \n')
        # get files in folder
        temp_path = folder / "ModelCoefficients"
        files_temp = list(temp_path.glob("*East*Extended.shp*.csv"))

        # loop through files and rename
        for file in files_temp:
            new_name = file.name.replace("_Extended.shp", "_Extended_smaller.shp")
            new_file = file.with_name(new_name)
            if not new_file.exists():
                log_list.append(str(file))
                file.rename(new_file)
                f.write(f"{str(file)}\n")

print("\nJOB COMPLETE\n")
