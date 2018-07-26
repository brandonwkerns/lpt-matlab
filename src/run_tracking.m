clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options
save('temp.mat');
OPT = load('temp.mat');
eval('!rm temp.mat')

% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
PROCESSED_DATA_DIR_IN = ['../data/',CASE_LABEL,'/processed/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), ...
                        'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/objects'];

PROCESSED_DATA_DIR_OUT = ['../data/',CASE_LABEL,'/processed/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), ...
                        'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

dn0 = DN1;%f0.time;
dn9 = DN2;%f9.time;

[year0,month0,day0,hour0] = datevec(dn0);
YYYY0 = sprintf('%d', year0);
MM0 = sprintf('%02d', month0);
DD0 = sprintf('%02d', day0);
HH0 = sprintf('%02d', hour0);

[year9,month9,day9,hour9] = datevec(dn9);
YYYY9 = sprintf('%d', year9);
MM9 = sprintf('%02d', month9);
DD9 = sprintf('%02d', day9);
HH9 = sprintf('%02d', hour9);

ymd0_ymd9 = [YYYY0,MM0,DD0,HH0,'_',YYYY9,MM9,DD9,HH9];
allPixelList=[PROCESSED_DATA_DIR_IN,'/objects_',ymd0_ymd9,'.mat'] ;

connect_time_clusters(allPixelList, OPT, TRACKING_VERBOSE) ;

% Move to output directory.
eval(['!mkdir -p ',PROCESSED_DATA_DIR_OUT])
eval(['!mv LONGSTATS_lpt_',ymd0_ymd9,' ',PROCESSED_DATA_DIR_OUT])
eval(['!mv TIMECLUSTERS_lpt_',ymd0_ymd9,'.mat ',PROCESSED_DATA_DIR_OUT])
eval(['!mv TIMECLUSTERS_lpt_*.nc ',PROCESSED_DATA_DIR_OUT])
