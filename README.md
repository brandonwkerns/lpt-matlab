# lpt-matlab
# Brandon Kerns
# bkerns@uw.edu
##############################################################################################

Large Scale Precipitation Tracking (Kerns and Chen 2016, JGR) code in Matlab.

Large Scale Precipitation Tracking (LPT) is tracking of large-scale (e. g., thousands of km)
	contiguos precipitation features that persist several days or weeks. It is an extension of
	mesoscale cloud cluster tracking (Chen et al. 1996, JAS; Chen and Houze 1997, QJRMS).
	Rainfall is first accumulated in to 3-day overlapping time periods, then subjected to a
	Gaussian smoothing function (standard deviation 5 deg. extending out 3 standard deviations).
	The tracking criteria including the accumulation time, filter size, and overlap criteria
	can be changed in the script "options.m."


Input files are either .mat files or NetCDF (.nc).

Output files are .txt, .mat, and .nc.

The directories are organized as follows:
config/		    	        Configuration files.
data/
data/raw/			Where input .mat files are copied or linked.
data/interim/accumulated/	The intermediate rain accumulation data.
data/interim/filtered/		The intermediate filtered rain accumulation data.
data/processed/features/	The processed LPT "shapshot" features.
data/processed/tracks/		The processed LPT tracks.
plots/				Where plots get stored.
src/				Matlab source code.

#############################################################################################

Here are the steps to run LPT. Recommend doing it for TRMM 3B42 (TMPA) data first.

Hint for viewing multiple files at once in ncview:
You may need to set "ulimit -n 4096" or similar to see all the files at once.
By default this may be set to 1024. If the user file limit is exceeded,
you will get errors like this:
fi_initialize: can't properly open file ../data/interim/filtered/rain_filtered_2017110512.nc.


1) Clone this repository to your system.

   git clone /home/orca/bkerns/lib/lpt/lpt-matlab -b master

2) Copy a template config file to "options.m", and edit it as needed.
   --> The template "config/options.KC16.m" file has the options used for KC16.
   --> The template "config/options.trmm_template.m" file has the updated MJO LPT identification.

3) get the accumulated rain data files.

   3a) It is recommended to link the original files into a subdirectory of "data/raw" directory
       so you remember where they came from. E.g., "data/raw/tmpa" for TMPA.

   3b) Scripts are provided for dealing with TRMM 3B42 (TMPA) data and WRF model output.
       Otherwise, you can write your own Matlab script.

       Scripts:
       src/preprocess_tmpa.m           --> .mat files
       src/calc_rain_accumulation.m    --> .nc files

       (Run the scripts from within the src/ directory.)


4) Run the spatial filtering script.

   Script: src/calc_rain_filter.m (master) --> .nc files
   Dependencies: src/gaussSmooth.m, src/gaussSmoothKernel.m (function dependency--don't run these on command line!)


5) Run the feature identification script.

   Script: connect_features.m
   --> Outputs will be named "LPT_features*"

6) Run the time tracking script.

   Script: run_tracking.m
   Dependency: connect_time_clusters.m
   --> Outputs will be named "LPT_tracks*"
