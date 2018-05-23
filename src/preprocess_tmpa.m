clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

% Set variables that pertain to this step.
% Note: ALL gridded rain rate files should be linked or copied into the "RAW_DATA_DIR"
RAW_DATA_DIR = '../data/raw/tmpa';
INTERIM_DATA_DIR = '../data/interim/gridded_rain_rates'

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

% NOTE: Start 3 days before DN1 to make sure a 3-day accumulation can be calculated.
for dn = DN1-3.0:datenum(0,0,0,DT,0,0):DN2

    [year,month,day,hour] = datevec(dn);

    YYYY = sprintf('%d',year);
    MM = sprintf('%02d',month);
    DD = sprintf('%02d',day);
    HH = sprintf('%02d', hour);

    %this_dir = [RAW_DATA_DIR,'/',YYYY,'/',MM,'/',YYYY,MM,DD];
    this_dir = [RAW_DATA_DIR];

    this_raw_file = [this_dir,'/3B42.',YYYY,MM,DD,'.',HH,'.7.HDF'];
    this_interim_file = [INTERIM_DATA_DIR,...
        '/gridded_rain_rates_',YYYY,MM,DD,HH,'.mat'];

    % Read in raw file.
    F=[] ;

    if ( ~ exist(this_raw_file) )
        this_raw_file= [this_dir,'/3B42.',YYYY,MM,DD,'.',HH,'.7A.HDF'];
    end

    if (exist(this_raw_file))
        F = read3b42(this_raw_file);
        this_rain = F.precip ;
        this_rain(this_rain < -0.01) = NaN ;

        disp(this_raw_file) ;


		       % Prepare struct and write interim output file.
	fout.lon=F.lon ;
	fout.lat=F.lat ;
	fout.time=datenum(year,month,day,hour,0,0) ;
	fout.rain=this_rain ;
	fout.info.smoothing='none' ;
	fout.info.accumulation='none' ;
	fout.info.units='mm/day' ;
	
	eval(['!mkdir -p ',INTERIM_DATA_DIR])
	eval(['save ',this_interim_file,' -struct fout'])

    else
        disp(['WARNING! No TMPA data found for ',YYYY,MM,DD,HH,'. Skipping.'])
    end


end
