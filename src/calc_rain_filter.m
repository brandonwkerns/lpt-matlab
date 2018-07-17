clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
INTERIM_DATA_DIR_IN = ['../data/',CASE_LABEL,'/interim/accumulated/',...
                        sprintf('%d',ACCUMULATION_PERIOD), 'h'];
INTERIM_DATA_DIR_OUT = ['../data/',CASE_LABEL,'/interim/filtered/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), 'h'];

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
      '/rain_accumulated_',...
      sprintf('%d',ACCUMULATION_PERIOD), 'h_',YYYY,MM,DD,HH,'.nc'];

  this_interim_file_out = [INTERIM_DATA_DIR_OUT,...
      '/rain_filtered_',...
      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
      sprintf('%d',ACCUMULATION_PERIOD), 'h_',YYYY,MM,DD,HH,'.nc'];

  disp([this_interim_file_in, ' --> ', this_interim_file_out])

  %f = load(this_interim_file_in) ;
  f0.lon = ncread(this_interim_file_in, 'lon');
  f0.lat = ncread(this_interim_file_in, 'lat');
  f0.rain = ncread(this_interim_file_in, 'rain')';

  %This gets used for the subsetting by LON and LAT from options.m
  %lonKeepIndx=find(f.lon > LON(1)-0.01 & f.lon < LON(end)+0.01) ;
  %latKeepIndx=find(f.lat > LAT(1)-0.01 & f.lat < LAT(end)+0.01) ;
  if USE_NATIVE_GRID == false
    f.lon = LON;
    f.lat = LAT;
    f.rain = area_conserve_remap(f0.lon,f0.lat,f0.rain,f.lon,f.lat);
  else
    f = f0;
  end

  %Generate the filter criteria.
  filter_nx = floor(2.0 * FILTER_WIDTH * FILTER_STANDARD_DEVIATION) + 1;
  filter_std_points = FILTER_STANDARD_DEVIATION;

  RAINFILTER=gaussSmooth(f.rain,filter_nx,filter_nx,filter_std_points,filter_std_points,...
                          FILTER_USE_GHOST_POINTS) ;
  RAINFILTER = RAINFILTER';
  %% NetCDF output.
  eval(['!mkdir -p ',INTERIM_DATA_DIR_OUT])

  % Define mode.
  % Dims
  cmode = netcdf.getConstant('CLOBBER');
  %cmode = bitor(cmode, netcdf.getConstant('64BIT_OFFSET'));
  ncid = netcdf.create(this_interim_file_out, cmode);
  dimid_lon  = netcdf.defDim(ncid, 'lon', numel(f.lon));
  dimid_lat  = netcdf.defDim(ncid, 'lat', numel(f.lat));
  dimid_time = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));

  % Vars
  varid_lon  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', dimid_lon);
  varid_lat  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', dimid_lat);
  varid_time = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', dimid_time);
  varid_rain = netcdf.defVar(ncid, 'rain', 'NC_DOUBLE', [dimid_lon dimid_lat dimid_time]);

  netcdf.endDef(ncid)

  % Data Mode
  netcdf.putVar(ncid, varid_lon, f.lon);
  netcdf.putVar(ncid, varid_lat, f.lat);
  netcdf.putVar(ncid, varid_time, 0,1, 86400.0 * (dn - datenum(1970,1,1,0,0,0)));
  netcdf.putVar(ncid, varid_rain, RAINFILTER);

  netcdf.close(ncid)

  ncwriteatt(this_interim_file_out,'time','units','seconds since 1970-1-1 0:0:0');
  ncwriteatt(this_interim_file_out,'rain','units','mm day-1');
  ncwriteatt(this_interim_file_out,'/','creation_date',datestr(now));
  ncwriteatt(this_interim_file_out,'/','smoothing',['Gaussian, ',num2str(FILTER_STANDARD_DEVIATION),...
      ' deg std dev., extending ',num2str(FILTER_WIDTH),'X std dev.']);
  ncwriteatt(this_interim_file_out,'/','accumulation', [num2str(ACCUMULATION_PERIOD), ' hours prior']);

end %End for loop.

disp(['Interim accumulated rain files written to: ',INTERIM_DATA_DIR_OUT])

disp('Done.')
