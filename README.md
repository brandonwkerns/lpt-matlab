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
config/		    	       		Configuration files.
data/
data/CASE_LABEL/raw/			Where input .mat files are copied or linked.
data/CASE_LABEL/interim/accumulated/	The intermediate rain accumulation data.
data/CASE_LABEL/interim/filtered/	The intermediate filtered rain accumulation data.
data/CASE_LABEL/processed/objects/	The processed LPT objects, e.g., "shapshots".
data/CASE_LABEL/processed/tracks/	The processed LPT tracks.
plots/CASE_LABEL/			Where plots get stored.
src/					Matlab source code.

#############################################################################################

Here are the steps to run LPT. Recommend doing it for TRMM 3B42 (TMPA) data first.
Run the scripts while inside of the "src" directory.

Hint for viewing multiple files at once in ncview:
You may need to set "ulimit -n 4096" or similar to see all the files at once.
By default this may be set to 1024. If the user file limit is exceeded,
you will get errors like this:
fi_initialize: can't properly open file ../data/interim/filtered/rain_filtered_2017110512.nc.


In the paths below, the following get replaced by whatever you specify in options.m:
CASE_LABEL    as specified in options.m.
FILTER        for example, g20 got 20 point Gaussian smoothing.
ACCUM         for example, 72h.
THRESH        for example thresh12 for 12 mm/day threshold.


1) Clone this repository to your system.

   git clone https://github.com/brandonwkerns/lpt-matlab.git 
	-- OR --
   On UW "orca" machine:
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
       preprocess_tmpa.m           --> .mat files
       calc_rain_accumulation.m    --> .nc files

       (Run the scripts from within the src/ directory.)


4) Run the spatial filtering script.

   Script: calc_rain_filter.m (master) --> .nc files
   Dependencies: gaussSmooth.m, gaussSmoothKernel.m (function dependency--don't run these on command line!)


5) Run the feature identification script.

   Script: identify_objects.m
   --> Outputs will be text, mat, and nc files named "data/CASE_LABEL/processed/FILTER_ACCUM/THRESH/objects_*"


6) Run the time tracking script.

   Script: run_tracking.m
   Dependency: connect_time_clusters.m
   --> Outputs will be text, mat, and nc files in "data/CASE_LABEL/processed/FILTER_ACCUM/THRESH/"

7) Group in to families, and identify MJO candidates if you want to use this version of track data.

   Scripts:
   -- identify_clumps_of_worms.m (needed to proceed with steps below)
   -- identify_mjo_candidates.m (optional, not needed for steps below)

### The following additional steps apply for Option 4, which is now the preferred option.
### These steps identify "families" of LPT systems and rejoin many of the tracks with
### short splits and mergers.

8) Run the first "rejoin" script, which handles LPT systems that split then re-merge back together.

   Scripts:
   -- recombine_split_n_merge_lpts.m (Generates TIMECLUSTERS*.rejoin.* files)
   -- identify_clumps_of_worms_rejoin.m (needed -- groupint in to LPT system families)
   -- identify_mjo_candidates_rejoin.m (optional -- only if you want MJO candidates from this step).

9) Run the second "rejoin" script, which handles short duration (e.g., < 3 day) mergers and splits.

   Scripts:
   -- recombine_split_n_merge_lpts2.m (Generates TIMECLUSTERS*.rejoin.* files)
   -- identify_clumps_of_worms_rejoin2.m (needed -- grouping on to LPT system families)
   -- identify_mjo_candidates_rejoin2.m (get MJO candidates from this step).

10) Plotting.

    A set of sample plotting scripts are provided with the package. These are set up for TRMM TMPA,
    and would need to be adapted to use with other data sets.

    calc_rainrate_hov_MULTI_15deg_3day_full_year.m and calc_rainrate_hov_MULTI_15deg_3h_full_year.m
    calculate the time-lon rain which is used in plotting scripts.

    The plot*.py and plot*.m scripts do the actual plotting.

    Mostly, plots are written to the plots/CASE_LABEL/ directory.
    