% This file contains the options and tracking criteria used for the LPT
% tracking. Load this file at the top of each of the LPT scripts.

CASE_LABEL = 's2s_ctl_init2013122500';   % Data goes in data/case/[raw|interim|processed]

%%
%% Data Criteria
%%

% Rainfall data resolution (in degrees lat/lon and hours).
% This is the resultion of the regular grid the data are interpolate to,
% if applicable.
% NOTE: I have only done LPT with 0.25 deg. data. But for time intervals
%       I have done 3 h for TRMM, 1 h for WRF ARW, and 6 h for ECMWF.
%       Also, LON and LAT below aren't used until the spatial filter step.
USE_NATIVE_GRID = true;   % Next 2 lines ignored if false. (Use native grid)
NATIVE_LON_VARNAME = 'longitude';
NATIVE_LAT_VARNAME = 'latitude';

DX  = 0.25; % Grid resolution. This and next two lines ignored if USE_NATIVE_GRID = true.
LON =  39.875 : DX : 178.875;  % Set the longitude grid.
LAT = -25.875 : DX : 24.875;   % Set the latitude grid.

% Time settings
DN1 = datenum(2013,12,25,0,0,0);      % Set the starting time as a datenum.
DN2 = datenum(2014,2,9,0,0,0);    % Set the ending time as a datenum.
DT  = 6.0;                           % Set the time interval (hours).

%%
%% Accumulation and filtering criteria
%%

% Accumulation period, in hours.
ACCUMULATION_PERIOD = 72.0; % hours

% If COLD_START_MODE is specified, assume there is no rain data before time zero.
%   Calculate the accumulation as follows:
%   For the first COLD_START_CONST_PERIOD, use the average rain rate during
%    the period, scaled to a 3 day accumulation, for the period.
%   After the COLD_START_CONST_PERIOD, use the accumulation up to that point,
%   scaled to a 3 day accumulation.
% If COLD_START_MODE is set to False, there should be gridded rain files for
%   the ACCUMULATION_PERIOD time prior to the initial time.
COLD_START_MODE = true;
COLD_START_CONST_PERIOD = 24.0; % hours

% Spatial filter (Gaussian bell shape) settings.
% Standard deviation is in lat/lon.
% Filter (half) width is in terms of how many standard deviations.
FILTER_STANDARD_DEVIATION = 4;   % In terms of grid points
FILTER_WIDTH = 3.0;

% For WRF regional domains, it may be good to use "ghost points" so that the LPT
% features are not truncated near the edges. Otherwise, the filter will use
% zeros at the boundary. Ghost points use the data reflected at the
% boundary for calculating the filter. The ghost points are discarded after
% calculating the filter. This only applies to the spatial filter step.
% The number of ghost points is the FILTER_WIDTH in the x direction only.
FILTER_USE_GHOST_POINTS = true;
FILTER_CYCLICAL_GHOST_POINTS = true;

%%
%% Feature identification
%%

% Threshold value to define a feature
% If the feature is > (<) the threshold value, set FEATURE_IS_GT_VALUE to true (false).
FEATURE_THRESHOLD_VALUE = 12;
FEATURE_IS_GT_VALUE = true;

% Minimum feature size is based on the numbe of points.
% LPT features smaller than this are discarded.
% (Kerns and Chen [2016] used 400 for 0.25 deg. pixels.)
FEATURE_MINIMUM_SIZE = 4;

% Max centroid latitude off the equator to identify a feature.
% Centroid abs(lat) > FEATURE_MAX_LAT is ignored.
% (Kerns and Chen [2016] used 15.0.)
FEATURE_MAX_LAT = 90.0;

%%
%% Tracking criteria
%%

TRACKING_VERBOSE = true;
CALC_MASK = false; % Calcluate mask arrays. WARNING: Takes ALOT of memory and MUCH longer to run!!!

% Overlap criteria. Minimum Minimum feature size is based on the number of points.
% If either the fraction or points overlap is satisfied, the LPT features
% will be combined into a LPT track.
TRACKING_MINIMUM_OVERLAP_FRAC = 0.5;
TRACKING_MINIMUM_OVERLAP_POINTS = 1600; %10;
TRACKING_MINIMUM_FRAMES = 2; % Minimum frames to keep a track.
TRACKING_MAX_TIME_TO_CONNECT = ACCUMULATION_PERIOD; % Max time between consecutive frames to connect CEs.
TRACKING_MAX_DIST_TO_CONNECT = 20.0; % Max distance between consecutive frames centroids.
% How to handle splits and mergers.


% Tracks that do not last at least this many hours will be discarded.
% NOTE: the pyhsical time scale is TRACKING_MINIMUM_DURATION + ACCUMULATION_PERIOD
%       due to the accumulation period.
TRACKING_MINIMUM_DURATION = 168;

%% Track splitting and merging method.
%%
%% 1 = KC16 method. For mergers, choose the previous track with longest areal accumulation history.
%%     And for splits, follow the system with the largest area at the time of the split.
%%     (Results in more spatially continuous tracks).
%% 2 = Separate tracks, keep together the track with largest accumulated area in time.
%%     (Results in longer tracked systems, but tends to be more jumpy in time).
%% 3 = DO NOT split or merge. Combine multiple CEs at the same time.
%%     Same as the Chen et al. (1996) "timeclusters" method for 208 K cloud clusters.
%% 4 = DO NOT split or merge. Keep duplicate copies of the tracks where they merge/split.
%%
SPLITTING_AND_MERGING_METHOD = 1;
