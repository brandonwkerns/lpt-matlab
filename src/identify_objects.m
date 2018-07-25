clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
INTERIM_DATA_DIR_IN = ['../data/',CASE_LABEL,'/interim/filtered/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), 'h'];
INTERIM_DATA_DIR_UNFILTERED = ['../data/',CASE_LABEL,'/interim/gridded_rain_rates'];
PROCESSED_DATA_DIR_OUT = ['../data/',CASE_LABEL,'/processed/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), ...
                        'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/objects'];

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

% Format string for "ceareas" ascii files.
fmt='%7.2f%8.2f%7.1f%7.1f%20.1f   %04d%02d%02d%02d      %7.2f%7.2f%7.1f  %7.2f\n' ;

%% Calculate area arrays in km.
%% These vary with latitude only.

interim_file_list = dir([INTERIM_DATA_DIR_IN,'/*.nc']);

% First get the "grid" from the first file.
f0.lon = ncread([INTERIM_DATA_DIR_IN,'/',interim_file_list(1).name],'lon');
f0.lat = ncread([INTERIM_DATA_DIR_IN,'/',interim_file_list(1).name],'lat');
f0.rain = ncread([INTERIM_DATA_DIR_IN,'/',interim_file_list(1).name],'rain')';

AREA=zeros(1,numel(f0.lat)) ;
for jj=1:numel(f0.lat)
    AREA(jj) = (DX*111.195) * (DX*111.195*cos(pi*f0.lat(jj)/180.0)) ;
end
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
allPixelList=[PROCESSED_DATA_DIR_OUT,'/objects_',ymd0_ymd9,'.mat'] ;
netcdf_output_fn = [PROCESSED_DATA_DIR_OUT,'/objects_',ymd0_ymd9,'.nc'] ;
fid=fopen([PROCESSED_DATA_DIR_OUT,'/objects_',ymd0_ymd9],'w') ;
eval(['!mkdir -p ',PROCESSED_DATA_DIR_OUT,'/',ymd0_ymd9]);

% Initialize the master output struct, which will contain the
% index for all of the cluster elements (CEs) in the same order
% as the ascii file opened as "fid" above. It will go in the
% matlab .mat file in "allPixelList."
fout_all.time_range=[dn0,dn9];
fout_all.time=[] ;
fout_all.lon=[] ;
fout_all.lat=[] ;
fout_all.area=[] ;
fout_all.volrain=[] ;
fout_all.pixels=[] ;
fout_all.max_filtered_rain=[] ;
fout_all.lon_max=[] ;
fout_all.lat_max=[] ;
fout_all.lon_rain_weighted=[] ;
fout_all.lat_rain_weighted=[] ;
fout_all.lon_median=[] ;
fout_all.lat_median=[] ;

INDX_ALL = 0;
begin_pixels_3d_size = 20000;
pixels_3d = zeros(begin_pixels_3d_size, numel(f0.lat), numel(f0.lon)) ;


% Now loop over the times updating the master output struct
% and individual time structs as needed.
for dn = DN1:datenum(0,0,0,DT,0,0):DN2

    [year,month,day,hour] = datevec(dn);
    yyyy = sprintf('%d', year);
    mm = sprintf('%02d', month);
    dd = sprintf('%02d', day);
    hh = sprintf('%02d', hour);

    this_interim_file_in = [INTERIM_DATA_DIR_IN,...
        '/rain_filtered_',...
        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
        sprintf('%d',ACCUMULATION_PERIOD), 'h_',yyyy,mm,dd,hh,'.nc'];

    F.lon = ncread(this_interim_file_in, 'lon') ;
    F.lat = ncread(this_interim_file_in, 'lat') ;
    F.rain = ncread(this_interim_file_in, 'rain')' ;

    thisPixelList=[PROCESSED_DATA_DIR_OUT,'/',ymd0_ymd9,'/objects_',yyyy,mm,dd,hh,'.mat'];
    thisCeareas=[PROCESSED_DATA_DIR_OUT,'/',ymd0_ymd9,'/objects_',yyyy,mm,dd,hh];
    thisFID=fopen(thisCeareas,'w') ;

    fileRain_Unfiltered=[INTERIM_DATA_DIR_UNFILTERED, ...
                        '/gridded_rain_rates_',yyyy,mm,dd,hh,'.nc'] ;
    Funfiltered.lon = ncread(fileRain_Unfiltered, 'lon');
    Funfiltered.lat = ncread(fileRain_Unfiltered, 'lat');
    Funfiltered.rain = ncread(fileRain_Unfiltered, 'rain');
    Funfiltered.rain = Funfiltered.rain';

    %TODO: This needs to do area_conserve_remap if custom grid is selected.

    RAINFILTER=F.rain ;
    RAINFILTER(~isfinite(RAINFILTER)) = 0.0;
    RAINFILTER_BW=logical(0.0*RAINFILTER) ;
    if (FEATURE_IS_GT_VALUE)
        RAINFILTER_BW(RAINFILTER > FEATURE_THRESHOLD_VALUE) = 1 ;
    else
        RAINFILTER_BW(RAINFILTER < FEATURE_THRESHOLD_VALUE) = 1 ;
    end
    RAINFILTER_BW(~isfinite(RAINFILTER)) = 0 ;

    RAINFILTER_keep=RAINFILTER;
    RAINFILTER_BW_keep=RAINFILTER_BW;

    RAINUNFILTERED=Funfiltered.rain ;
    RAINUNFILTERED_keep=RAINUNFILTERED;
    RAINUNFILTERED_keep(~isfinite(RAINUNFILTERED_keep))=0.0 ;

    % Initialize some stuff.
    CC1=[] ;
    stats1=[] ;

    % Use bwconncomp and regionprops to get CE features.
    % Use 4-point connectivity, i.e., don't connect diagonally.
    CC1=bwconncomp(RAINFILTER_BW_keep,4) ;
    stats1=regionprops(CC1,'all') ;

    fout.time=[] ;
    fout.lon=[] ;
    fout.lat=[] ;
    fout.area=[] ;
    fout.volrain=[] ;
    fout.pixels=[] ;
    fout.max_filtered_rain=[] ;
    fout.lon_max=[] ;
    fout.lat_max=[] ;
    fout.lon_rain_weighted=[] ;
    fout.lat_rain_weighted=[] ;
    fout.lon_median=[] ;
    fout.lat_median=[] ;

    pixels_2d = [];

    INDX=0 ;

    for iii=1:numel(stats1) ;

        centroidY=stats1(iii).Centroid(2) ;
        centroidX=stats1(iii).Centroid(1) ;
        thisLon=interp1(1:numel(f0.lon),f0.lon,centroidX) ;
        thisLat=interp1(1:numel(f0.lat),f0.lat,centroidY) ;

        % CE area
        thisArea=0.0 ;
        thisPixelY=stats1(iii).PixelList(:,2) ;
        thisPixelX=stats1(iii).PixelList(:,1) ;

        for jjj=1:numel(thisPixelY) %area only depends on Y.
            thisArea=thisArea + AREA(thisPixelY(jjj)) ;
        end

        %CE Volumetric Rain ( 10^9 mm^3/day -or- km^3/day -or- mm-(1000 km)^2 / day )
        thisVolRain=0.0 ;
        for jjj=1:numel(thisPixelY) %area only depends on Y.

            thisVolRain=thisVolRain + ...
                (RAINUNFILTERED_keep(thisPixelY(jjj),thisPixelX(jjj))*...
                 AREA(thisPixelY(jjj)) ) ;

        end
        thisVolRain=thisVolRain/1e6 ; % mm-km^2 ==> mm-(1000 km)^2

        % median centers.
        medianX=median(thisPixelX);
        medianY=median(thisPixelY);

        thisLonMedian=interp1(1:numel(f0.lon),f0.lon,medianX) ;
        thisLatMedian=interp1(1:numel(f0.lon),f0.lon,medianY) ;

        % Center of the MAX filtered rainfall.
        filtered_rain_list=[];
        for xxxx=1:numel(thisPixelY);
            filtered_rain_list(xxxx)=RAINFILTER_keep(...
                thisPixelY(xxxx), thisPixelX(xxxx)) ;
        end
        [max_filtered_rain,max_filtered_rain_index]=max(filtered_rain_list);

        max_filtered_rain_X=thisPixelX(max_filtered_rain_index);
        max_filtered_rain_Y=thisPixelY(max_filtered_rain_index);
        max_filtered_rain_lon=F.lon(max_filtered_rain_X);
        max_filtered_rain_lat=F.lat(max_filtered_rain_Y);

        % Rainfall weighted center of mass.
        rain_weighted_X=sum( filtered_rain_list*thisPixelX) / ...
            sum(filtered_rain_list) ;

        rain_weighted_Y=sum( filtered_rain_list*thisPixelY) / ...
            sum(filtered_rain_list) ;

        rain_weighted_lon=interp1(1:numel(F.lon),F.lon,rain_weighted_X) ;
        rain_weighted_lat=interp1(1:numel(F.lat),F.lat,rain_weighted_Y) ;

        %% Apply the "ce" area and centroid latitude criteria here

        if ( stats1(iii).Area < FEATURE_MINIMUM_SIZE | abs(thisLat) > FEATURE_MAX_LAT)
            RAINFILTER_BW_keep(CC1.PixelIdxList{iii})=0;
        else
            INDX=INDX+1 ;
            INDX_ALL=INDX_ALL+1 ;

            % Define the single day output struct.
            fout.time(INDX)=dn ;
            fout.lon(INDX)=thisLon ;
            fout.lat(INDX)=thisLat ;
            fout.area(INDX)=thisArea ;
            fout.volrain(INDX)=thisVolRain ;
            fout.max_filtered_rain(INDX)=max_filtered_rain ;
            fout.lon_max(INDX)=max_filtered_rain_lon ;
            fout.lat_max(INDX)=max_filtered_rain_lat ;
            fout.lon_rain_weighted(INDX)=rain_weighted_lon ;
            fout.lat_rain_weighted(INDX)=rain_weighted_lat ;
            fout.lon_median(INDX)=thisLonMedian ;
            fout.lat_median(INDX)=thisLatMedian ;

            fout.pixels(INDX).x=thisPixelX ;
            fout.pixels(INDX).y=thisPixelY ;

            %%
            %% Update 3D Pixel array.
            %%

            % If the pixels_3d is too small, add more space.
            if (mod(INDX_ALL, begin_pixels_3d_size) == 0)
              pixels_3d = cat(1, pixels_3d, zeros(begin_pixels_3d_size, numel(F.lat), numel(F.lon))) ;
            end

            pixels_2d = zeros(numel(F.lat), numel(F.lon));
            for kkk = 1:numel(thisPixelX)
              pixels_2d(thisPixelY(kkk), thisPixelX(kkk)) = 1;
            end
            pixels_3d(INDX_ALL,:,:) = pixels_2d;

            %%
            %% Pring out text output.
            %%

            fprintf(fid,fmt, thisLat, thisLon, ...
                    centroidY, centroidX, ...
                    thisArea, year, month, day, hour, ...
                    stats1(iii).Eccentricity, ...
                    stats1(iii).MajorAxisLength/stats1(iii).MinorAxisLength,...
                    stats1(iii).Orientation,thisVolRain ) ;

            fprintf(thisFID,fmt, thisLat, thisLon, ...
                    centroidY, centroidX, ...
                    thisArea, year, month, day, hour, ...
                    stats1(iii).Eccentricity, ...
                    stats1(iii).MajorAxisLength/stats1(iii).MinorAxisLength,...
                    stats1(iii).Orientation,thisVolRain ) ;

        end
    end

    fclose(thisFID) ;

    % Update the "grand index" struct.
    fout_all.time=[fout_all.time,fout.time] ;
    fout_all.lon=[fout_all.lon,fout.lon] ;
    fout_all.lat=[fout_all.lat,fout.lat] ;
    fout_all.area=[fout_all.area,fout.area] ;
    fout_all.volrain=[fout_all.volrain,fout.volrain] ;
    fout_all.max_filtered_rain=[fout_all.max_filtered_rain, ...
                        fout.max_filtered_rain] ;
    fout_all.lon_max=[fout_all.lon_max,fout.lon_max] ;
    fout_all.lat_max=[fout_all.lat_max,fout.lat_max] ;

    fout_all.pixels=[fout_all.pixels,fout.pixels] ;
    fout_all.lon_rain_weighted=[fout_all.lon_rain_weighted,...
                        fout.lon_rain_weighted] ;
    fout_all.lat_rain_weighted=[fout_all.lat_rain_weighted,...
                        fout.lat_rain_weighted] ;

    fout_all.lon_median=[fout_all.lon_median,fout.lon_median] ;
    fout_all.lat_median=[fout_all.lat_median,fout.lat_median] ;

    fout.grid.lon = F.lon;
    fout.grid.lat = F.lat;
    fout.grid.area=AREA ;

    disp(thisPixelList)
    eval(['save ',thisPixelList,' -struct fout'])

end

fclose(fid) ;
fout_all.grid.lon = F.lon;
fout_all.grid.lat = F.lat;

fout_all.grid.area=AREA ;

pixels_3d = pixels_3d(1:numel(fout_all.lon),:,:);

disp(allPixelList)
eval(['save ',allPixelList,' -struct fout_all'])


% NetCDF Output.
eval(['!mkdir -p ',PROCESSED_DATA_DIR_OUT])

% Define mode.
% Dims

disp(['Writing NetCDF: ', netcdf_output_fn])
cmode = netcdf.getConstant('CLOBBER');
ncid = netcdf.create(netcdf_output_fn, cmode);
dimid_obs  = netcdf.defDim(ncid, 'ceid', numel(fout_all.lon));

dimid_lon_grid  = netcdf.defDim(ncid, 'lon_grid', numel(fout_all.grid.lon));
dimid_lat_grid  = netcdf.defDim(ncid, 'lat_grid', numel(fout_all.grid.lat));

% Vars
varid_ceid = netcdf.defVar(ncid, 'ceid', 'NC_INT', dimid_obs);
varid_lon  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', dimid_obs);
varid_lat  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', dimid_obs);
varid_time = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', dimid_obs);
varid_endtime = netcdf.defVar(ncid, 'end_of_accumulation_time', 'NC_DOUBLE', dimid_obs);
varid_area = netcdf.defVar(ncid, 'area', 'NC_DOUBLE', dimid_obs);
varid_volrain = netcdf.defVar(ncid, 'volrain', 'NC_DOUBLE', dimid_obs);

varid_lon_grid  = netcdf.defVar(ncid, 'lon_grid', 'NC_DOUBLE', dimid_lon_grid);
varid_lat_grid  = netcdf.defVar(ncid, 'lat_grid', 'NC_DOUBLE', dimid_lat_grid);
varid_pixels_3d = netcdf.defVar(ncid, 'pixels', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_obs]);

netcdf.endDef(ncid)


% Data Mode
netcdf.putVar(ncid, varid_ceid, 1:numel(fout_all.lon));
netcdf.putVar(ncid, varid_lon, fout_all.lon);
netcdf.putVar(ncid, varid_lat, fout_all.lat);
netcdf.putVar(ncid, varid_time, 24.0 * (fout_all.time - datenum(1970,1,1,0,0,0)) - 0.5*ACCUMULATION_PERIOD);
netcdf.putVar(ncid, varid_endtime, 24.0 * (fout_all.time - datenum(1970,1,1,0,0,0)));
netcdf.putVar(ncid, varid_area, fout_all.area);
netcdf.putVar(ncid, varid_volrain, fout_all.volrain);

netcdf.putVar(ncid, varid_lon_grid, fout_all.grid.lon);
netcdf.putVar(ncid, varid_lat_grid, fout_all.grid.lat);
netcdf.putVar(ncid, varid_pixels_3d, permute(pixels_3d,[3,2,1]));


netcdf.close(ncid)

% Attributes
ncwriteatt(netcdf_output_fn,'lon','units','degrees_east');
ncwriteatt(netcdf_output_fn,'lon_grid','units','degrees_east');
ncwriteatt(netcdf_output_fn,'lat','units','degrees_north');
ncwriteatt(netcdf_output_fn,'lat_grid','units','degrees_north');
ncwriteatt(netcdf_output_fn,'time','units','hours since 1970-1-1 0:0:0');
ncwriteatt(netcdf_output_fn,'time','description','MIDDLE of accumulation period time. Most representative time.');
ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','units','hours since 1970-1-1 0:0:0');
ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','description','END of accumulation period time. For use with accumulated fields.');
ncwriteatt(netcdf_output_fn,'area','units','km2');
ncwriteatt(netcdf_output_fn,'area','coordinates','time lat lon');
ncwriteatt(netcdf_output_fn,'volrain','units','mm - km2');
ncwriteatt(netcdf_output_fn,'volrain','coordinates','time lat lon');

ncwriteatt(netcdf_output_fn,'/','description', 'LP objects (i.e., snapshots) file.');
ncwriteatt(netcdf_output_fn,'/','creation_date',datestr(now));
ncwriteatt(netcdf_output_fn,'/','smoothing',['Gaussian, ',num2str(FILTER_STANDARD_DEVIATION),...
    ' deg std dev., extending ',num2str(FILTER_WIDTH),'X std dev.']);
ncwriteatt(netcdf_output_fn,'/','accumulation', [num2str(ACCUMULATION_PERIOD), ' hours prior']);


disp('Done.')
