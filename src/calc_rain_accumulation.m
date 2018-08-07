clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
INTERIM_DATA_DIR_IN = ['../data/',CASE_LABEL,'/interim/gridded_rain_rates'];
INTERIM_DATA_DIR_OUT = ['../data/',CASE_LABEL,'/interim/accumulated/',...
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

    RAINCOLLECT = [] ;


    if (COLD_START_MODE == true) & (24.0*(dn-DN1) < COLD_START_CONST_PERIOD-0.01)
      accumulation_time_list = DN1+datenum(0,0,0,DT,0,0):datenum(0,0,0,DT,0,0):(DN1+COLD_START_CONST_PERIOD/24.0);
      % accumulation_scale_factor = ACCUMULATION_PERIOD / COLD_START_CONST_PERIOD;
    elseif (COLD_START_MODE == true) & (24.0*(dn-DN1) < ACCUMULATION_PERIOD-0.01)
      hour_rel1 = -1*(24.0*(dn-DN1));
      accumulation_time_list = dn + (hour_rel1:DT:0)/24;
      % accumulation_scale_factor = ACCUMULATION_PERIOD/(24.0*(dn-DN1));
    else
      hour_rel1 = -1*ACCUMULATION_PERIOD;
      accumulation_time_list = dn + (hour_rel1:DT:0)/24;
      % accumulation_scale_factor = 1.0;
    end

    kkk=0 ;
    for dn0 = [accumulation_time_list]
        kkk=kkk+1 ;

        %dn0 = dn + hour_rel/24;
        [year0,month0,day0,hour0] = datevec(dn0);

        YYYY0 = sprintf('%d', year0);
        MM0 = sprintf('%02d', month0);
        DD0 = sprintf('%02d', day0);
        HH0 = sprintf('%02d', hour0);

        this_interim_file_in = [INTERIM_DATA_DIR_IN,...
            '/gridded_rain_rates_',YYYY0,MM0,DD0,HH0,'.nc'];

        this_interim_file_out = [INTERIM_DATA_DIR_OUT,...
            '/rain_accumulated_',sprintf('%d',ACCUMULATION_PERIOD), 'h_',YYYY,MM,DD,HH,'.nc'];

        if ( ~ exist(this_interim_file_in) )
            disp(['WARNING: Did not find ',this_interim_file_in])
            disp(['WARNING: ',this_interim_file_out,' may be wrong!'])
            
        else

            disp(this_interim_file_in)
            f.lon = ncread(this_interim_file_in,'lon');
            f.lat = ncread(this_interim_file_in,'lat');
            f.rain = ncread(this_interim_file_in,'rain')';

            THISRAIN = f.rain ;
            THISRAIN(THISRAIN < -0.01) = NaN;
            RAINCOLLECT(:,:,kkk) = THISRAIN;
        end

    end

    RAINAVG=nanmean(RAINCOLLECT,3)*24.0; %*accumulation_scale_factor; % units: mm/h --> mm/day

    %%
    %% NetCDF output.
    %%

    eval(['!mkdir -p ',INTERIM_DATA_DIR_OUT])
    disp(['--> ',this_interim_file_out])

    % Define mode.
    % Dims
    ncid = netcdf.create(this_interim_file_out, 'CLOBBER');
    dimid_lon  = netcdf.defDim(ncid, 'lon', numel(f.lon));
    dimid_lat  = netcdf.defDim(ncid, 'lat', numel(f.lat));
    dimid_time = netcdf.defDim(ncid, 'time', 1);

    % Vars
    varid_lon  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', dimid_lon);
    varid_lat  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', dimid_lat);
    varid_time = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', dimid_time);
    varid_rain = netcdf.defVar(ncid, 'rain', 'NC_DOUBLE', [dimid_lon dimid_lat dimid_time]);

    netcdf.endDef(ncid)

    % Data Mode
    netcdf.putVar(ncid, varid_lon, f.lon);
    netcdf.putVar(ncid, varid_lat, f.lat);
    netcdf.putVar(ncid, varid_time, 24.0 * (dn - datenum(1970,1,1,0,0,0)));
    netcdf.putVar(ncid, varid_rain, RAINAVG');

    netcdf.close(ncid)

    ncwriteatt(this_interim_file_out,'time','units','hours since 1970-1-1 0:0:0');
    ncwriteatt(this_interim_file_out,'rain','units','mm day-1');
    ncwriteatt(this_interim_file_out,'/','creation_date',datestr(now));
    ncwriteatt(this_interim_file_out,'/','smoothing','none');
    ncwriteatt(this_interim_file_out,'/','accumulation', [num2str(ACCUMULATION_PERIOD), ' hours prior']);

end

disp(['Interim accumulated rain files written to: ',INTERIM_DATA_DIR_OUT])

disp('Done.')
