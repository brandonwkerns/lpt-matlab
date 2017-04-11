# lpt-matlab
# Brandon Kerns
# bkerns@miami.edu
######################################

Large Scale Precipitation Tracking (Kerns and Chen 2016, JGR) code in Matlab.

Large Scale Precipitation Tracking (LPT) is tracking of large-scale (e. g., thousands of km)
	contiguos precipitation features that persist several days or weeks. It is an extension of
	mesoscale cloud cluster tracking (Chen et al. 1996, JAS; Chen and Houze 1997, QJRMS). 
	Rainfall is first accumulated in to 3-day overlapping time periods, then subjected to a 
	Gaussian smoothing function (standard deviation 5 deg. extending out 3 standard deviations). 
	The tracking criteria including the accumulation time, filter size, and overlap criteria
	can be changed in the script "options.m."


The easiest way to run LPT is to:

1) Clone this repository to your system.

   git clone brandonwkerns/lpt-matlab

2) Edit the "options.m" file if needed.

3) get the accumulated rain data into .mat files.

   Scripts are provided for dealing with TRMM 3B42 (TMPA) data and WRF model output.
   Otherwise, use your own script. The 

4) Run the spatial filtering script.

5) Run the feature identification script.

   Outputs will be named "LPT_feature*.mat"

6) Run the tracking script.

   Outputs will be named "LPT_track*.mat"


