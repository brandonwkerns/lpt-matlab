% This file contains the options and tracking criteria used for the LPT
% tracking. Load this file at the top of each of the LPT scripts.

CASE_LABEL = 'trmm';   % Data goes in data/case/[raw|interim|processed]

%%
%% Data Criteria
%%

% Rainfall data resolution (in degrees lat/lon and hours).
% This is the resultion of the regular grid the data are interpolate to,
% if applicable.
% NOTE: I have only done LPT with 0.25 deg. data. But for time intervals
%       I have done 3 h for TRMM, 1 h for WRF ARW, and 6 h for ECMWF.
%       Also, LON and LAT below aren't used until the spatial filter step.
USE_NATIVE_GRID = false;   % Next 2 lines ignored if false. (Use native grid)
NATIVE_LON_VARNAME = 'XLONG';
NATIVE_LAT_VARNAME = 'XLAT';

DX  = 0.25; % Grid resolution. This and next two lines ignored if USE_NATIVE_GRID = true.
LON =  39.875 : DX : 178.875;  % Set the longitude grid.
LAT = -25.875 : DX : 24.875;   % Set the latitude grid.

% Time settings.
DN1 = datenum(2011,10,1,0,0,0);      % Set the starting time as a datenum.
DN2 = datenum(2012,3,31,21,0,0);    % Set the ending time as a datenum.
DT  = 3.0;                           % Set the time interval (hours).

%%
%% Accumulation and filtering criteria
%%

% Accumulation period, in hours.
ACCUMULATION_PERIOD = 72.0;

% Spatial filter (Gaussian bell shape) settings.
% Standard deviation is in lat/lon.
% Filter (half) width is in terms of how many standard deviations.
FILTER_STANDARD_DEVIATION = 20;  % In terms of grid points
FILTER_WIDTH = 3.0;

% For WRF regional domains, it may be good to use "ghost points" so that the LPT
% features are not truncated near the edges. Otherwise, the filter will use
% zeros at the boundary. Ghost points use the data reflected at the
% boundary for calculating the filter. The ghost points are discarded after
% calculating the filter. This only applies to the spatial filter step.
FILTER_USE_GHOST_POINTS = false;
FILTER_NUMBER_OF_GHOST_POINTS = 10;

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
FEATURE_MINIMUM_SIZE = 400;

% Max centroid latitude off the equator to identify a feature.
% Centroid abs(lat) > FEATURE_MAX_LAT is ignored.
% (Kerns and Chen [2016] used 15.0.)
FEATURE_MAX_LAT = 15.0;

%%
%% Tracking criteria
%%

% Overlap criteria. Minimum Minimum feature size is based on the number of points.
% If either the fraction or points overlap is satisfied, the LPT features
% will be combined into a LPT track.
TRACKING_MINIMUM_OVERLAP_FRAC = 0.5;
TRACKING_MINIMUM_OVERLAP_POINTS = 10;
TRACKING_MINIMUM_FRAMES = 2; % Minimum frames to keep a track.
TRACKING_MAX_TIME_TO_CONNECT = 3.0; % Max time between consecutive frames to connect CEs.
TRACKING_MAX_DIST_TO_CONNECT = 20.0; % Max distance between consecutive frames centroids.
% How to handle splits and mergers.


% Tracks that do not last at least this many hours will be discarded.
% NOTE: the pyhsical time scale is TRACKING_MINIMUM_DURATION + ACCUMULATION_PERIOD
%       due to the accumulation period.
TRACKING_MINIMUM_DURATION = 96;
