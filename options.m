% This file contains the options and tracking criteria used for the LPT
% tracking. Load this file in each of the LPT scripts.

%% Data Criteria

% Rainfall data resolution (in degrees lat/lon and hours).
DATA_RES_X = 0.25;
DATA_RES_Y = 0.25;
DATA_RES_T = 3.0;


%% Accumulation and filtering criteria

% Accumulation period, in hours.
ACCUMULATION_PERIOD = 72.0;

% Spatial filter (Gaussian bell shape) settings.
% Standard deviation is in lat/lon.
% Filter width is in terms of how many standard deviations.
FILTER_STANDARD_DEVIATION = 5.0;
FILTER_WIDTH = 3.0;

% For smaller domains, it may be good to use "ghost points" so that the LPT
% features are not truncated near the edges. Otherwise, the filter will use
% zeros at the boundary. Ghost points use the data reflected at the
% boundary for calculating the filter. The ghost points are discarded after
% calculating the filter.
FILTER_USE_GHOST_POINTS = false;
FILTER_NUMBER_OF_GHOST_POINTS = 10;


%% Feature identification

% Threshold value to define a feature
% If the feature is > (<) the value, set FEATURE_IS_GT_VALUE to true (false).
FEATURE_THRESHOLD_VALUE = 12; 
FEATURE_IS_GT_VALUE = true; 

% Minimum feature size is based on the numbe of points.
% LPT features smaller than this are discarded.
FEATURE_MINIMUM_SIZE = 4; 


%% Tracking criteria

% Overlap criteria. Minimum Minimum feature size is based on the number of points.
% If either the fraction or points overlap is satisfied, the LPT features
% will be combined into a LPT track.
TRACKING_MINIMUM_OVERLAP_FRAC = 0.5;
TRACKING_MINIMUM_OVERLAP_POINTS = 10;

% How to handle splits and mergers.


% Tracks that do not last at least this many hours will be discarded.
TRACKING_MINIMUM_DURATION = 168;

