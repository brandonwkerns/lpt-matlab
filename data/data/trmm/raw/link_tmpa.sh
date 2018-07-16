#!/bin/sh

## Link TMPA data on orca.
## Edit data directory and years below.

data_dir=/home/orca/data/satellite/trmm_global_rainfall
year1=2017
year2=2018


############### End of editing. ##########################

# Clean up
rm -f *.HDF


# Link Files
for MM in 06 07 08 09 10 11 12
do
    ln -sv $data_dir/$year1/$MM/*/*.HDF .
done

for MM in 01 02 03 04 05 06
do
    ln -sv $data_dir/$year2/$MM/*/*.HDF .
done



exit 0

