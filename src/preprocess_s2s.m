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

this_raw_file = [RAW_DATA_DIR,'/total_precip_accum_ecmf_mdate2017-12-25_hdate2013-12-25.nc'];
F.lon = double(ncread(this_raw_file,'longitude'));
F.lat = double(ncread(this_raw_file,'latitude'));
F.time = double(ncread(this_raw_file,'time'))/24.0 + datenum(1900,1,1,0,0,0);
F.precip_accum = ncread(this_raw_file,'tp'); %order of data is: lon, lat, time.
F.precip = 0.0 * F.precip_accum;
F.precip(:,:,2:end) = F.precip_accum(:,:,2:end) - F.precip_accum(:,:,1:end-1);
conversion_factor = 6.0 ;%kg / m2 (6h) to mm h-1
F.precip = F.precip * conversion_factor;
disp(this_raw_file) ;


for dn = DN1:datenum(0,0,0,DT,0,0):DN2

  [year,month,day,hour] = datevec(dn);

  YYYY = sprintf('%d',year);
  MM = sprintf('%02d',month);
  DD = sprintf('%02d',day);
  HH = sprintf('%02d', hour);

  this_interim_file = [INTERIM_DATA_DIR,...
      '/gridded_rain_rates_',YYYY,MM,DD,HH,'.nc'];

  disp(this_interim_file)

  % Get "this rain"
  this_time_indx = find(F.time == dn);
  this_rain = squeeze(F.precip(:,:,this_time_indx));

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
  netcdf.putVar(ncid, varid_time, 0,1, 24.0 * (datenum(year,month,day,hour,0,0) - datenum(1970,1,1,0,0,0)));
  netcdf.putVar(ncid, varid_rain, this_rain);

  netcdf.close(ncid)

  ncwriteatt(this_interim_file,'time','units','hours since 1970-1-1 0:0:0');
  ncwriteatt(this_interim_file,'rain','units','mm h-1');
  ncwriteatt(this_interim_file,'/','creation_date',datestr(now));
  ncwriteatt(this_interim_file,'/','smoothing',['None.']);
  ncwriteatt(this_interim_file,'/','accumulation', ['None.']);

end %end dn loop

disp(['Interim gridded rain files written to: ',INTERIM_DATA_DIR])
