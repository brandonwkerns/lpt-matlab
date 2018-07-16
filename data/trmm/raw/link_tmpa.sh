#!/bin/sh

## Link TMPA data on orca.
## Edit data directory and years below.
## NOTE: It is OK to link a longer time period than you will use for actual tracking.
data_dir=/home/orca/data/satellite/trmm_global_rainfall
year1=2017
year2=2018

############### End of editing. ##########################

# Clean up
rm -f *.HDF

# Link Files
ln -sv $data_dir/$year1/*/*/*.HDF .
ln -sv $data_dir/$year2/*/*/*.HDF .

exit 0

