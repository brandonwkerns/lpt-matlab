clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
INTERIM_DATA_DIR_IN = '../data/interim/accumulate_rain';
INTERIM_DATA_DIR_OUT = '../data/interim/filter_rain';

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

for dn = DN1:datenum(0,0,0,DT,0,0):DN2

    [year,month,day,hour] = datevec(dn);

    YYYY = sprintf('%d', year);
    MM = sprintf('%02d', month);
    DD = sprintf('%02d', day);
    HH = sprintf('%02d', hour);

    this_interim_file_in = [INTERIM_DATA_DIR_IN,...
        '/rain_accumulated_',YYYY,MM,DD,HH,'.mat'];

    this_interim_file_out = [INTERIM_DATA_DIR_OUT,...
        '/rain_filtered_',YYYY,MM,DD,HH,'.mat'];

    disp([this_interim_file_in, ' --> ', this_interim_file_out])

    f = load(this_interim_file_in) ;

    %This gets used for the subsetting by LON and LAT from options.m
    lonKeepIndx=find(f.lon > LON(1)-0.01 & f.lon < LON(end)+0.01) ;
    latKeepIndx=find(f.lat > LAT(1)-0.01 & f.lat < LAT(end)+0.01) ;

    %Generate the filter criteria.
    filter_nx = floor(2.0 * FILTER_WIDTH * FILTER_STANDARD_DEVIATION/DX) + 1;
    filter_std_points = FILTER_STANDARD_DEVIATION/DX;

    RAINFILTER=gaussSmooth(f.rain,filter_nx,filter_nx,filter_std_points,filter_std_points) ;

    fout.lon=f.lon(lonKeepIndx);
    fout.lat=f.lat(latKeepIndx);
    fout.time=dn;
    fout.rain=RAINFILTER(latKeepIndx,lonKeepIndx);
    fout.info.filter=['Gaussian, ',num2str(FILTER_STANDARD_DEVIATION),...
        ' deg std dev., extending ',num2str(FILTER_WIDTH),'X std dev.'];
    fout.info.accumulation='3 days prior' ;
    fout.info.units='mm/day' ;

    eval(['!mkdir -p ',INTERIM_DATA_DIR_OUT]) ;
    eval(['save ',this_interim_file_out,' -struct fout'])

end
