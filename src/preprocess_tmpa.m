clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

% Set variables that pertain to this step.
% Note: ALL gridded rain rate files should be linked or copied into the "RAW_DATA_DIR"
RAW_DATA_DIR = ['../data/',CASE_LABEL,'/raw'];
INTERIM_DATA_DIR = ['../data/',CASE_LABEL,'/interim/gridded_rain_rates'];

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

% NOTE: Start "ACCUMULATION_PERIOD_ hours before DN1 to make sure the full accumulation can be calculated.
for dn = DN1-ACCUMULATION_PERIOD/24.0:datenum(0,0,0,DT,0,0):DN2

  [year,month,day,hour] = datevec(dn);

  YYYY = sprintf('%d',year);
  MM = sprintf('%02d',month);
  DD = sprintf('%02d',day);
  HH = sprintf('%02d', hour);

  this_dir = [RAW_DATA_DIR];

  this_raw_file = [this_dir,'/3B42.',YYYY,MM,DD,'.',HH,'.7.HDF'];
  this_interim_file = [INTERIM_DATA_DIR,...
      '/gridded_rain_rates_',YYYY,MM,DD,HH,'.nc'];

  % Read in raw file.
  F=[] ;

  if ( ~ exist(this_raw_file) )
    this_raw_file= [this_dir,'/3B42.',YYYY,MM,DD,'.',HH,'.7A.HDF'];
  end

  if (exist(this_raw_file))
    F = read3b42(this_raw_file);

    %% Longitude 0 - 360
    F.lon = [F.lon(721:1440), 360.0 + F.lon(1:720)];
    F.precip = [F.precip(:,721:1440), F.precip(:,1:720)];


    this_rain = F.precip' ;
    this_rain(this_rain < -0.01) = NaN ;

    disp(this_raw_file) ;


    %% NetCDF output.
    eval(['!mkdir -p ',INTERIM_DATA_DIR])

    % Define mode.
    % Dims
    cmode = netcdf.getConstant('CLOBBER');
    ncid = netcdf.create(this_interim_file, cmode);
    dimid_lon  = netcdf.defDim(ncid, 'lon', numel(F.lon));
    dimid_lat  = netcdf.defDim(ncid, 'lat', numel(F.lat));
    dimid_time = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));

    % Vars
    varid_lon  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', dimid_lon);
    varid_lat  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', dimid_lat);
    varid_time = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', dimid_time);
    varid_rain = netcdf.defVar(ncid, 'rain', 'NC_DOUBLE', [dimid_lon dimid_lat dimid_time]);

    netcdf.endDef(ncid)

    % Data Mode
    netcdf.putVar(ncid, varid_lon, F.lon);
    netcdf.putVar(ncid, varid_lat, F.lat);
    netcdf.putVar(ncid, varid_time, 0,1, 86400.0 * (datenum(year,month,day,hour,0,0) - datenum(1970,1,1,0,0,0)));
    netcdf.putVar(ncid, varid_rain, this_rain);

    netcdf.close(ncid)

    ncwriteatt(this_interim_file,'time','units','seconds since 1970-1-1 0:0:0');
    ncwriteatt(this_interim_file,'rain','units','mm h-1');
    ncwriteatt(this_interim_file,'/','creation_date',datestr(now));
    ncwriteatt(this_interim_file,'/','smoothing',['None.']);
    ncwriteatt(this_interim_file,'/','accumulation', ['None.']);

  else
    disp(['WARNING! No TMPA data found for ',YYYY,MM,DD,HH,'. Skipping.'])
  end %end if file existance.

end %end dn loop

disp(['Interim gridded rain files written to: ',INTERIM_DATA_DIR])
