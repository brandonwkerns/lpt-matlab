function TIMECLUSTERS=connect_time_clusters(master_ce_file, OPT, verbose)

%%
%% TIMECLUSTERS=connect_time_clusters_g10_3day_t14_connect_branches(year, verbose)
%%
%% year = beginning year (e.g., 2011 for Oct. 2011 - March 2012)
%% verbose = 0 (default for off) or 1 for on.
%%
%% Approach here is to first group the cluster elements (CEs)
%%  in to their respective time clusters (TCs), allowing
%%  for mergers and splits.
%% Then, from the database of CEs, the tracking parameters are
%%  obtained, e.g., lat, lon of centroid and branch tracks.
%% (This function uses sub-functions each of which access and
%%  modify the TIMECLUSTERS struct array. So it cannot be run
%%  as a script--only a function.
%%

%% Set the tracking parameters here.
maxTimeToConnect = OPT.TRACKING_MAX_TIME_TO_CONNECT;% 3 ; %Maximum time between clusters (hours) to
                     %connect them. Especially useful if there are
                     %some empty frames vs. the case when some
                     %times have data but do not have clusters.

minFramesToKeep = OPT.TRACKING_MINIMUM_FRAMES ; %Need this many entries in track to keep it.
minDuration = OPT.TRACKING_MINIMUM_DURATION; % 96 ; % Minimum duration to keep it (hours).

maxLat = OPT.FEATURE_MAX_LAT ; %Positive. Max distance off the equator to keep it.

maxDistToConnect = OPT.TRACKING_MAX_DIST_TO_CONNECT ;%20.0 %; %Max. distance between centroids to match them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ( nargin < 3 )
    verbose=0 ;
end

disp(master_ce_file)
CE=load(master_ce_file) ;


%%
%% Intitialize time cluster stuff.
%%

TIMECLUSTERS=[] ;
nextClusterID=1 ;  %Start time cluster ID count at 1.
allDates=[] ;


DN=OPT.DN1:datenum(0,0,0,OPT.DT,0,0):OPT.DN2;


%%
%% Loop forward in time.
%% -- Mergers will be represented as duplicate, overlapping TCs.
%%    These need to be re-merged.
%% -- Splits will have the splitting branches as the same TC.
%%
for dn=[DN]

    ceINDXthisTime=find(CE.time == dn) ;
    if ( numel(ceINDXthisTime) < 1 )
        continue
    end

    if (numel(allDates) > 0.1)
        prevDate=max(allDates(allDates < dn-0.01 ) ) ;
    else
        prevDate=datenum(1900,1,1,0,0,0) ;
    end

    % If there is a gap longer than maxTimeToConnect,
    %  force all of this time's clusters to be new tracks.
    if ( 24*(dn - prevDate) > maxTimeToConnect+0.1 )
            prevDate=datenum(1900,1,1,0,0,0) ;
    end

    % Get a set of matching time clusters for each
    %  cluster feature present at this time.

    for thisCE=ceINDXthisTime  % Loop over the cluster features
                                 % present at this time.

        istr=num2str(thisCE) ;

        matchingClusterID=matchingTimeCluster(thisCE,prevDate);


        if ( numel( matchingClusterID ) <  1   )


            % This is a new cluster.


            if ( verbose > 0 )
                disp([istr,': No matching time cluster in the database.'])
            end


            % Starting new cloud cluster.
            startNewTimeCluster(nextClusterID) ;
            if ( verbose > 0 )
                disp(['        + Created new time cluster: ID = ',num2str(nextClusterID)])
            end


            addCluster(nextClusterID,thisCE) ;

            nextClusterID = nextClusterID + 1 ;

        else

            if ( verbose > 0 )
                disp([istr,': Matches with current existing cluster id = ',num2str(matchingClusterID)])
            end




            %
            % This will create duplicate tracks when there are more than one
            %  matchingClusterID's. To fix this, use below either
            %  mergeDuplicateTCs() or splitDuplicateTCs().
            %
            %{
            for iii=1:numel(matchingClusterID)
                addCluster(matchingClusterID(iii),thisCE) ;
            end
            %}

            % Use this version instead to avoid splitting tracks.
            %  Only have the timecluster track follow the largest CE.
            for iii=1:numel(matchingClusterID)
                % Is there already a CE assigned to the matchingClusterID?
                findMatchingTime=find(CE.time(TIMECLUSTERS(matchingClusterID(iii)).ceid) == dn) ;

                if ( numel(findMatchingTime)==0)
                    addCluster(matchingClusterID(iii),thisCE) ;
                else
                    %findMatchingTime
                    %TIMECLUSTERS(matchingClusterID(iii)).ceid(findMatchingTime)

                    matchingArea=CE.area(TIMECLUSTERS(matchingClusterID(iii)).ceid(findMatchingTime));

                    %[matchingArea,CE.area(thisCE)]

                    if ( CE.area(thisCE) > matchingArea  )

                        startNewTimeCluster(nextClusterID) ;
                        if ( verbose > 0 )
                            disp(['        + Created new time cluster: ID = ',num2str(nextClusterID)])
                        end

                        addCluster(nextClusterID,TIMECLUSTERS(matchingClusterID(iii)).ceid(findMatchingTime)) ;
                        nextClusterID = nextClusterID + 1 ;
                        TIMECLUSTERS(matchingClusterID(iii)).ceid(findMatchingTime)=thisCE;

                    else

                        startNewTimeCluster(nextClusterID) ;
                        if ( verbose > 0 )
                            disp(['        + Created new time cluster: ID = ',num2str(nextClusterID)])
                        end

                        addCluster(nextClusterID,thisCE) ;
                        nextClusterID = nextClusterID + 1 ;

                    end

                end
            end

            % End of the version of code to avoid split tracks.


        end

    end   % END Loop over the cluster features
          % present at this time.

    % Update the list of clusters.
    allDates=unique([allDates,dn]) ;


end % END Loop over Dates





%%
%% Now combining the merging tracks which were duplicates.
%%


%mergeDuplicateTCs() ;
splitDuplicateTCs() ;



%%
%% Take out time clusters with insufficient duration.
%%

removeShortLivedTCs(minDuration) ;
combineCloseProximityTCs(10.0,3.0) ;

%%
%% Get tracking parameters from the CE database
%%

disp('Calculating tracking parameters.')
calcTrackingParameters() ;

disp('Calculating mask arrays.')
maskArrays = calcMaskArrays(TIMECLUSTERS, CE, DN, OPT) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Output  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dn0 = OPT.DN1;%f0.time;
dn9 = OPT.DN2;%f9.time;

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


disp('Writing Output.')
FMT='        %4d%02d%02d%02d %7d %10.2f %10.2f %1d\n' ;


% Ascii output
fileout=['LONGSTATS_lpt_',ymd0_ymd9] ;
disp(fileout)
fid=fopen(fileout,'w') ;

for ii=1:numel(TIMECLUSTERS)

    fprintf(fid,'Cl%i:\n',ii) ;

    for jj = 1:numel(TIMECLUSTERS(ii).time)

        [y,m,d,h]=datevec(TIMECLUSTERS(ii).time(jj)) ;

        fprintf(fid,FMT,y,m,d,h,...
                round(TIMECLUSTERS(ii).area(jj)),...
                TIMECLUSTERS(ii).lat(jj),...
                TIMECLUSTERS(ii).lon(jj),...
                TIMECLUSTERS(ii).nclusters(jj)) ;

    end

end



fclose(fid) ;

% .mat file output
fout.TIMECLUSTERS = TIMECLUSTERS ;
fout.grid = CE.grid ;
fileout_mat=['TIMECLUSTERS_lpt_',ymd0_ymd9,'.mat'] ;
disp(fileout_mat)
eval(['save ', fileout_mat, ' -struct fout'])




% NetCDF output.
%%
%% Global Stuff
%%

% Get "stitch" track data.
MISSING = -9999.0
stitch_id = MISSING;
stitch_lon = MISSING;
stitch_lat = MISSING;
stitch_time = MISSING;
stitch_end_accumulation_time = MISSING;
stitch_area = MISSING;
stitch_size = MISSING;
stitch_effective_radius = MISSING;
stitch_volrain = MISSING;
for indx = 1:numel(TIMECLUSTERS)
  stitch_id = [stitch_id, MISSING, indx + 0.0 * TIMECLUSTERS(indx).time];
  stitch_time = [stitch_time, MISSING, TIMECLUSTERS(indx).time - 0.5 * OPT.ACCUMULATION_PERIOD/24.0];
  stitch_lon = [stitch_lon, MISSING, TIMECLUSTERS(indx).lon];
  stitch_lat = [stitch_lat, MISSING, TIMECLUSTERS(indx).lat];
  stitch_area = [stitch_area, MISSING, TIMECLUSTERS(indx).area];
  stitch_size = [stitch_size, MISSING, TIMECLUSTERS(indx).size];
  stitch_effective_radius = [stitch_effective_radius, MISSING, TIMECLUSTERS(indx).effective_radius];
  stitch_volrain = [stitch_volrain, MISSING, TIMECLUSTERS(indx).volrain];
end

% Define mode.
% Dims
netcdf_output_fn=['TIMECLUSTERS_lpt_',ymd0_ymd9,'.lptALL.nc'] ;
disp(['Writing NetCDF: ', netcdf_output_fn])
cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(netcdf_output_fn, cmode);

dimid_lpt_id  = netcdf.defDim(ncid, 'lpt_id', numel(TIMECLUSTERS));
dimid_time  = netcdf.defDim(ncid, 'time', numel(DN));
dimid_lon_grid  = netcdf.defDim(ncid, 'lon', numel(CE.grid.lon));
dimid_lat_grid  = netcdf.defDim(ncid, 'lat', numel(CE.grid.lat));
dimid_obs  = netcdf.defDim(ncid, 'obs', numel(stitch_lon));


% Global Vars
% Coordinates
varid_lpt_id = netcdf.defVar(ncid, 'lpt_id', 'NC_INT', dimid_lpt_id);
varid_obs = netcdf.defVar(ncid, 'obs', 'NC_INT', dimid_obs);
varid_time = netcdf.defVar(ncid, 'alltime', 'NC_DOUBLE', dimid_time);
varid_end_time = netcdf.defVar(ncid, 'end_of_accumulation_alltime', 'NC_DOUBLE', dimid_time);
varid_lon_grid  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', dimid_lon_grid);
varid_lat_grid  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', dimid_lat_grid);

% Summary stats
varid_lpt_duration = netcdf.defVar(ncid, 'duration', 'NC_DOUBLE', dimid_lpt_id);
varid_lpt_max_area = netcdf.defVar(ncid, 'max_area', 'NC_DOUBLE', dimid_lpt_id);
varid_lpt_max_volrain = netcdf.defVar(ncid, 'max_volrain', 'NC_DOUBLE', dimid_lpt_id);
varid_lpt_max_size = netcdf.defVar(ncid, 'max_size', 'NC_DOUBLE', dimid_lpt_id);
varid_lpt_max_effective_radius = netcdf.defVar(ncid, 'max_effective_radus', 'NC_DOUBLE', dimid_lpt_id);
varid_lpt_zonal_propagation_speed = netcdf.defVar(ncid, 'zonal_propagation_speed', 'NC_DOUBLE', dimid_lpt_id);
varid_lpt_meridional_propagation_speed = netcdf.defVar(ncid, 'meridional_propagation_speed', 'NC_DOUBLE', dimid_lpt_id);

% Stitched individual time clusters
varid_stitch_id = netcdf.defVar(ncid, 'id', 'NC_INT', dimid_obs);
varid_stitch_time = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', dimid_obs);
varid_stitch_end_time = netcdf.defVar(ncid, 'end_of_accumulation_time', 'NC_DOUBLE', dimid_obs);
varid_stitch_lon = netcdf.defVar(ncid, 'centroid_lon', 'NC_DOUBLE', dimid_obs);
varid_stitch_lat = netcdf.defVar(ncid, 'centroid_lat', 'NC_DOUBLE', dimid_obs);
varid_stitch_area = netcdf.defVar(ncid, 'area', 'NC_DOUBLE', dimid_obs);
varid_stitch_size = netcdf.defVar(ncid, 'size', 'NC_DOUBLE', dimid_obs);
varid_stitch_effective_radius = netcdf.defVar(ncid, 'effective_radius', 'NC_DOUBLE', dimid_obs);
varid_stitch_volrain = netcdf.defVar(ncid, 'volrain', 'NC_DOUBLE', dimid_obs);
netcdf.defVarFill(ncid,varid_stitch_id,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_time,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_end_time,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_lon,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_lat,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_area,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_size,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_effective_radius,false,MISSING);
netcdf.defVarFill(ncid,varid_stitch_volrain,false,MISSING);

% Mask
varid_lpt_mask_by_id  = netcdf.defVar(ncid, 'mask_by_id', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_time]);
varid_lpt_mask_by_id_with_accumulation  = netcdf.defVar(ncid, 'mask_with_accumulation', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_time]);
varid_lpt_mask_by_id_with_filter  = netcdf.defVar(ncid, 'mask_with_filter', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_time]);
varid_lpt_mask_by_id_with_filter_and_accumulation  = netcdf.defVar(ncid, 'mask_with_filter_and_accumulation', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_time]);
netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id,true,true,1);
netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_filter,true,true,1);
netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_accumulation,true,true,1);
netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_filter_and_accumulation,true,true,1);

netcdf.endDef(ncid)

% Get summary time series data.
max_areas = zeros(1,numel(TIMECLUSTERS));
max_volrains = zeros(1,numel(TIMECLUSTERS));
max_sizes = zeros(1,numel(TIMECLUSTERS));
max_effective_radii = zeros(1,numel(TIMECLUSTERS));
durations = zeros(1,numel(TIMECLUSTERS));
zonal_propagation_speeds = zeros(1,numel(TIMECLUSTERS));
meridional_propagation_speeds = zeros(1,numel(TIMECLUSTERS));
for indx = 1:numel(TIMECLUSTERS)
  durations(indx) = TIMECLUSTERS(indx).duration;
  max_areas(indx) = TIMECLUSTERS(indx).maxarea;
  max_volrains(indx) = TIMECLUSTERS(indx).maxvolrain;
  max_sizes(indx) = TIMECLUSTERS(indx).maxsize;
  max_effective_radii(indx) = TIMECLUSTERS(indx).max_effective_radius;
  zonal_propagation_speeds(indx) = TIMECLUSTERS(indx).zonal_propagation_speed;
  meridional_propagation_speeds(indx) = TIMECLUSTERS(indx).meridional_propagation_speed;
end

% Data mode
% Coordinates
netcdf.putVar(ncid, varid_lpt_id, 1:numel(TIMECLUSTERS));
netcdf.putVar(ncid, varid_time, 86400.0 * (DN - datenum(1970,1,1,0,0,0) - 0.5*OPT.ACCUMULATION_PERIOD/24.0));
netcdf.putVar(ncid, varid_end_time, 86400.0 * (DN - datenum(1970,1,1,0,0,0)));
netcdf.putVar(ncid, varid_lon_grid, CE.grid.lon);
netcdf.putVar(ncid, varid_lat_grid, CE.grid.lat);

% LPT overall summary
netcdf.putVar(ncid, varid_lpt_duration, durations);
netcdf.putVar(ncid, varid_lpt_max_area, max_areas);
netcdf.putVar(ncid, varid_lpt_max_volrain, max_volrains);
netcdf.putVar(ncid, varid_lpt_max_size, max_sizes);
netcdf.putVar(ncid, varid_lpt_max_effective_radius, max_effective_radii);

netcdf.putVar(ncid, varid_lpt_zonal_propagation_speed, zonal_propagation_speeds);
netcdf.putVar(ncid, varid_lpt_meridional_propagation_speed, meridional_propagation_speeds);

% stitched data
netcdf.putVar(ncid, varid_stitch_id, stitch_id);
netcdf.putVar(ncid, varid_stitch_time, 86400.0 * (stitch_time - datenum(1970,1,1,0,0,0)));
netcdf.putVar(ncid, varid_stitch_end_time, 86400.0 * (stitch_time - datenum(1970,1,1,0,0,0) - 0.5*OPT.ACCUMULATION_PERIOD/24.0));
netcdf.putVar(ncid, varid_stitch_lon, stitch_lon);
netcdf.putVar(ncid, varid_stitch_lat, stitch_lat);
netcdf.putVar(ncid, varid_stitch_area, stitch_area);
netcdf.putVar(ncid, varid_stitch_size, stitch_size);
netcdf.putVar(ncid, varid_stitch_effective_radius, stitch_effective_radius);
netcdf.putVar(ncid, varid_stitch_volrain, stitch_volrain);

% masks
netcdf.putVar(ncid, varid_lpt_mask_by_id, ...
  permute(maskArrays.all.mask_by_id, [3,2,1]));
netcdf.putVar(ncid, varid_lpt_mask_by_id_with_accumulation, ...
  permute(maskArrays.all.mask_by_id_with_accumulation, [3,2,1]));
netcdf.putVar(ncid, varid_lpt_mask_by_id_with_filter, ...
  permute(maskArrays.all.mask_by_id_with_filter, [3,2,1]));
netcdf.putVar(ncid, varid_lpt_mask_by_id_with_filter_and_accumulation, ...
  permute(maskArrays.all.mask_by_id_with_filter_and_accumulation, [3,2,1]));

netcdf.close(ncid)

% Attributes
ncwriteatt(netcdf_output_fn,'lpt_id','units','1');
ncwriteatt(netcdf_output_fn,'lon','units','degrees_east');
ncwriteatt(netcdf_output_fn,'lat','units','degrees_north');
ncwriteatt(netcdf_output_fn,'centroid_lon','units','degrees_east');
ncwriteatt(netcdf_output_fn,'centroid_lat','units','degrees_north');
ncwriteatt(netcdf_output_fn,'alltime','units','seconds since 1970-1-1 0:0:0');
ncwriteatt(netcdf_output_fn,'alltime','description','Entire tracking period. Middle of accumulation period time.');
ncwriteatt(netcdf_output_fn,'end_of_accumulation_alltime','units','seconds since 1970-1-1 0:0:0');
ncwriteatt(netcdf_output_fn,'end_of_accumulation_alltime','description','Entire tracking period. End of accumulation period time.');
ncwriteatt(netcdf_output_fn,'time','units','seconds since 1970-1-1 0:0:0');
ncwriteatt(netcdf_output_fn,'time','description','Stitched time of LPTs. Middle of accumulation period time.');
ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','units','seconds since 1970-1-1 0:0:0');
ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','description','Stitched time of LPTs. End of accumulation period time.');

ncwriteatt(netcdf_output_fn,'zonal_propagation_speed','units','m s-1');
ncwriteatt(netcdf_output_fn,'meridional_propagation_speed','units','m s-1');
ncwriteatt(netcdf_output_fn,'zonal_propagation_speed','description','Eastward propagation from linear best fit.');
ncwriteatt(netcdf_output_fn,'meridional_propagation_speed','description','Northward propagation from linear best fit.');

ncwriteatt(netcdf_output_fn,'area','units','km2');
ncwriteatt(netcdf_output_fn,'area','description','Stitched area of LPTs.');
ncwriteatt(netcdf_output_fn,'size','units','km');
ncwriteatt(netcdf_output_fn,'size','description','Stitched size of LPTs. Square root of area.');
ncwriteatt(netcdf_output_fn,'effective_radius','units','km');
ncwriteatt(netcdf_output_fn,'effective_radius','description','Stitched effective radius of LPTs.');
ncwriteatt(netcdf_output_fn,'volrain','units','mm km2');
ncwriteatt(netcdf_output_fn,'volrain','description','Stitched volumetric rain of LPTs.');

ncwriteatt(netcdf_output_fn,'max_area','units','km2');
ncwriteatt(netcdf_output_fn,'max_size','units','km');
ncwriteatt(netcdf_output_fn,'max_size','description','Square root of area.');
ncwriteatt(netcdf_output_fn,'max_effective_radus','units','km');
ncwriteatt(netcdf_output_fn,'max_volrain','units','mm km2');




%%
%% Individual LPT stuff.
%%

for indx = 1:numel(TIMECLUSTERS)

  % Define mode.

  netcdf_output_fn=['TIMECLUSTERS_lpt_',ymd0_ymd9,'.lpt',sprintf('%03d',indx),'.nc'] ;
  disp(['Writing NetCDF: ', netcdf_output_fn])
  cmode = netcdf.getConstant('CLOBBER');
  cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
  ncid = netcdf.create(netcdf_output_fn, cmode);

  dimid_lon_grid  = netcdf.defDim(ncid, 'lon', numel(CE.grid.lon));
  dimid_lat_grid  = netcdf.defDim(ncid, 'lat', numel(CE.grid.lat));
  dimid_obs  = netcdf.defDim(ncid, 'time', numel(TIMECLUSTERS(indx).time));
  dimid_1  = netcdf.defDim(ncid, 'summary', 1);

  % Individual LPT Vars
  varid_1_duration = netcdf.defVar(ncid, 'duration', 'NC_DOUBLE', dimid_1);

  varid_time = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', dimid_obs);
  varid_end_time = netcdf.defVar(ncid, 'end_of_accumulation_time', 'NC_DOUBLE', dimid_obs);
  varid_lon  = netcdf.defVar(ncid, 'centroid_lon', 'NC_DOUBLE', dimid_obs);
  varid_lat  = netcdf.defVar(ncid, 'centroid_lat', 'NC_DOUBLE', dimid_obs);
  varid_lat_grid  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', dimid_lat_grid);
  varid_lon_grid  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', dimid_lon_grid);

  varid_area = netcdf.defVar(ncid, 'area', 'NC_DOUBLE', dimid_obs);
  varid_volrain = netcdf.defVar(ncid, 'volrain', 'NC_DOUBLE', dimid_obs);

  varid_lpt_mask_by_id  = netcdf.defVar(ncid, 'mask_by_id', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_obs]);
  varid_lpt_mask_by_id_with_filter  = netcdf.defVar(ncid, 'mask_with_filter', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_obs]);
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id,true,true,1);
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_filter,true,true,1);


  % Data mode.
  netcdf.endDef(ncid)

  netcdf.putVar(ncid, varid_1_duration, TIMECLUSTERS(indx).duration);
  netcdf.putVar(ncid, varid_lon_grid, CE.grid.lon);
  netcdf.putVar(ncid, varid_lat_grid, CE.grid.lat);
  netcdf.putVar(ncid, varid_time, 86400.0 * (TIMECLUSTERS(indx).time - datenum(1970,1,1,0,0,0) - 0.5*OPT.ACCUMULATION_PERIOD/24.0));
  netcdf.putVar(ncid, varid_end_time, 86400.0 * (TIMECLUSTERS(indx).time - datenum(1970,1,1,0,0,0)));
  netcdf.putVar(ncid, varid_lon, TIMECLUSTERS(indx).lon);
  netcdf.putVar(ncid, varid_lat, TIMECLUSTERS(indx).lat);
  netcdf.putVar(ncid, varid_area, TIMECLUSTERS(indx).area);
  netcdf.putVar(ncid, varid_volrain, TIMECLUSTERS(indx).volrain);


  netcdf.putVar(ncid, varid_lpt_mask_by_id, ...
    permute(maskArrays.individual(indx).mask_by_id, [3,2,1]));
  netcdf.putVar(ncid, varid_lpt_mask_by_id_with_filter, ...
    permute(maskArrays.individual(indx).mask_by_id_with_filter, [3,2,1]));


  netcdf.close(ncid)


  % Attributes
  ncwriteatt(netcdf_output_fn,'time','units','seconds since 1970-1-1 0:0:0');
  ncwriteatt(netcdf_output_fn,'time','description','Middle of accumulation period time.');
  ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','units','seconds since 1970-1-1 0:0:0');
  ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','description','End of accumulation period time.');
  ncwriteatt(netcdf_output_fn,'lon','units','degrees_east');
  ncwriteatt(netcdf_output_fn,'lat','units','degrees_north');
  ncwriteatt(netcdf_output_fn,'centroid_lon','units','degrees_east');
  ncwriteatt(netcdf_output_fn,'centroid_lat','units','degrees_north');
  ncwriteatt(netcdf_output_fn,'area','units','km2');
  ncwriteatt(netcdf_output_fn,'volrain','units','mm km2');


end

disp('Done.')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Local functions below here  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function startNewTimeCluster(tcid)

    TIMECLUSTERS(tcid).ceid = [] ;

end


function addCluster(tcid,ceid)

%
% tcid = index of TIMECLUSTERS struct array to use.
% pixels = "pixels" struct for this time
% pixelindx = which "pixels" entry to use.
%
% Modifies global function variable TIMECLUSTERS
%

    TIMECLUSTERS(tcid).ceid=[TIMECLUSTERS(tcid).ceid,ceid];

end




function  matchingClusterID=matchingTimeCluster(theCE,timeToSearch);

    matchingClusterID=[] ;

    if ( numel(TIMECLUSTERS) < 1 )
        return
    else


        %% Loop through previous time clusters, looking for a match.

        for ii=1:numel(TIMECLUSTERS)


            times=[] ;

            for  jj=[TIMECLUSTERS(ii).ceid]

                %% Check distance here.
                dLON1=CE.lon(jj)-CE.lon(theCE) ;
                dLAT1=CE.lat(jj)-CE.lat(theCE) ;

                if ( sqrt(dLON1.^2 + dLAT1.^2 ) > maxDistToConnect )
                    continue
                end


                if (CE.time(jj) == timeToSearch)

                    X1=CE.pixels(jj).x ;
                    Y1=CE.pixels(jj).y ;

                    X2=CE.pixels(theCE).x ;
                    Y2=CE.pixels(theCE).y ;

                    hits=[] ;

                    A=[X1,Y1] ;
                    B=[X2,Y2] ;


                    HITS=intersect(A,B,'rows') ;
                    nHits=size(HITS,1);
                    %nHits=numel(HITS)


                    if ( nHits/numel(X1) > 0.75 | ...
                         nHits/numel(X2) > 0.75 | ...
                         nHits > 1600  )
                        %% 10 deg.^2 = 1600 pixels

                        matchingClusterID=[matchingClusterID,ii] ;

                    end

                end

            end

        end

    end

end



function mergeDuplicateTCs()

%% Accesses main function variables: TIMECLUSTERS
%% Modifies TIMECLUSTERS

    allTCindices=1:numel(TIMECLUSTERS) ;
    throwAwayTCs=[] ;

    for iiii=[allTCindices]

        otherTCindices=setxor(iiii,allTCindices) ;

        thisClusterCEList=TIMECLUSTERS(iiii).ceid ;


        for jjjj=[otherTCindices]

            if ( numel(intersect(TIMECLUSTERS(iiii).ceid,...
                                 TIMECLUSTERS(jjjj).ceid))>0)


                TIMECLUSTERS(iiii).ceid=sort(union(TIMECLUSTERS(iiii).ceid,...
                                              TIMECLUSTERS(jjjj).ceid));

                TIMECLUSTERS(jjjj).ceid=[] ;
                throwAwayTCs=[throwAwayTCs,jjjj] ;


            end

        end

    end

    % Throw away any empty TIMECLUSTERS

    TIMECLUSTERS(throwAwayTCs)=[] ;

end


function splitDuplicateTCs()

%% Accesses main function variables: TIMECLUSTERS
%% Modifies TIMECLUSTERS
%%
%% When two timecluster tracks merge, only keep the one with the
%% longer history in terms of accumulated added area.

    allTCindices=1:numel(TIMECLUSTERS) ;
    throwAwayTCs=[] ;

    for iiii=[allTCindices]



        otherTCindices=setxor(iiii,allTCindices) ;

        thisClusterCEList=TIMECLUSTERS(iiii).ceid ;


        for jjjj=[otherTCindices]

            intersections=intersect(TIMECLUSTERS(iiii).ceid,...
                                    TIMECLUSTERS(jjjj).ceid);



            if ( numel(intersections)>0)

                %Get earliest intersecting time.

                intersectionTimes=CE.time(intersections) ;
                firstIntersectionTime=min(intersectionTimes) ;

                %Get accumulated area of thisCluster
                thisClusterAccumArea=0.0 ;
                for thisClusterCEID=TIMECLUSTERS(iiii).ceid
                    %if ( CE.time(thisClusterCEID) < firstIntersectionTime-0.01 & ...
                    %     CE.time(thisClusterCEID) > firstIntersectionTime-1.01)

                        thisClusterAccumArea=thisClusterAccumArea+...
                            CE.area(thisClusterCEID);

                        %end
                end

                %Get accumulated area of otherCluster
                otherClusterAccumArea=0.0 ;
                for otherClusterCEID=TIMECLUSTERS(jjjj).ceid
                    %if ( CE.time(otherClusterCEID) < firstIntersectionTime-0.01 & ...
                    %     CE.time(otherClusterCEID) > firstIntersectionTime-1.01)

                        otherClusterAccumArea=otherClusterAccumArea+...
                            CE.area(otherClusterCEID);

                        %end
                end


                if ( thisClusterAccumArea > otherClusterAccumArea )

                    TIMECLUSTERS(jjjj).ceid=setxor(TIMECLUSTERS(jjjj).ceid,...
                                                   intersections);


                else

                    TIMECLUSTERS(iiii).ceid=setxor(TIMECLUSTERS(iiii).ceid,...
                                                   intersections);


                end



            end

        end

    end


end





function removeShortLivedTCs(minduration) ;

%% Accesses main function variables: TIMECLUSTERS, CE
%% Modifies TIMECLUSTERS

    allTCindices=1:numel(TIMECLUSTERS) ;
    throwAwayTCs=[] ;


    for iiii=[allTCindices]

        thisTCtimeList=[] ;

        if ( numel(TIMECLUSTERS(iiii).ceid) < 1 )
            throwAwayTCs=[throwAwayTCs,iiii] ;
        end

        for jjjj=1:numel(TIMECLUSTERS(iiii).ceid)

            thisTCtimeList=[thisTCtimeList,...
                            CE.time(TIMECLUSTERS(iiii).ceid(jjjj))] ;

        end

        duration=24*(max(thisTCtimeList)-min(thisTCtimeList)) ;

        if ( duration < minduration-0.01 )
            throwAwayTCs=[throwAwayTCs,iiii] ;
        end


    end

    % Throw unwanted ones away
    TIMECLUSTERS(throwAwayTCs)=[] ;

end


function combineCloseProximityTCs(maxCombineLonDiff,maxCombineTimeDiff) ;

    TC_eliminate_list=[];

    for ii=1:numel(TIMECLUSTERS)


        otherClusters=setxor(1:numel(TIMECLUSTERS),ii);

        %Test for GG_before
        for ii_before=[otherClusters]

            jumpTime=CE.time(TIMECLUSTERS(ii).ceid(1)) - ...
                     CE.time(TIMECLUSTERS(ii_before).ceid(end));

            jumpLon=CE.lon(TIMECLUSTERS(ii).ceid(1)) - ...
                    CE.lon(TIMECLUSTERS(ii_before).ceid(end));


            if ( abs(jumpTime) < maxCombineTimeDiff+0.01 & ...
                 abs(jumpLon) < maxCombineLonDiff+0.01 )

                TIMECLUSTERS(ii).ceid=sort([TIMECLUSTERS(ii).ceid,...
                                    TIMECLUSTERS(ii_before).ceid]);

                disp(['Combined LPTs: ',num2str(ii_before),' in to ',...
                      num2str(ii),'.'])


                %TIMECLUSTERS(ii_before).ceid=[] ;
                TC_eliminate_list=[TC_eliminate_list,ii_before];

            end

        end

    end


    TIMECLUSTERS(TC_eliminate_list)=[];

end

function calcTrackingParameters() ;

%% Accesses main function variables: TIMECLUSTERS, CE
%% Modifies TIMECLUSTERS

    %% Grid.
    gridLon=CE.grid.lon ;
    gridLat=CE.grid.lat ;

    gridDX=gridLon(2)-gridLon(1) ;
    gridDY=gridLat(2)-gridLat(1) ;

    AREA=CE.grid.area ;

    allTCindices=1:numel(TIMECLUSTERS) ;

    for iiii=[allTCindices]

        thisTCtimeList=[] ;


        for jjjj=1:numel(TIMECLUSTERS(iiii).ceid)

            thisTCtimeList=[thisTCtimeList,...
                            CE.time(TIMECLUSTERS(iiii).ceid(jjjj))] ;


        end

        thisTCtimeListUniq=unique(thisTCtimeList) ;

        % Time and duration for time cluster over all.
        TIMECLUSTERS(iiii).time=thisTCtimeListUniq ;
        TIMECLUSTERS(iiii).duration=24*(max(thisTCtimeList)-min(thisTCtimeList)) ;
        TIMECLUSTERS(iiii).nentries=numel(thisTCtimeListUniq) ;


        for tttt=1:numel(thisTCtimeListUniq) ;


            thisTimeMatch=find(thisTCtimeList == thisTCtimeListUniq(tttt));

            thisTimeCEID=TIMECLUSTERS(iiii).ceid(thisTimeMatch);


            % The CEs at each time. TC.ce(tttt).ceid, ect.x

            TIMECLUSTERS(iiii).ce(tttt).ceid=thisTimeCEID ;
            TIMECLUSTERS(iiii).ce(tttt).time=CE.time(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).lon=CE.lon(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).lat=CE.lat(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).area=CE.area(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).size=sqrt(CE.area(thisTimeCEID)) ;
            TIMECLUSTERS(iiii).ce(tttt).effective_radius=...
                sqrt(CE.area(thisTimeCEID)/pi) ;
            TIMECLUSTERS(iiii).ce(tttt).volrain=CE.volrain(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).pixels=CE.pixels(thisTimeCEID) ;


            % Combine the CEs for the net TC time dependent properties.
            % ( full feature centroid tracks)
            TIMECLUSTERS(iiii).lon(tttt)=...
                sum(CE.lon(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));

            TIMECLUSTERS(iiii).lat(tttt)=...
                sum(CE.lat(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));

            TIMECLUSTERS(iiii).lon_median(tttt)=...
                sum(CE.lon_median(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));
            TIMECLUSTERS(iiii).lat_median(tttt)=...
                sum(CE.lat_median(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));


            TIMECLUSTERS(iiii).area(tttt)=sum(CE.area(thisTimeCEID));
            TIMECLUSTERS(iiii).size(tttt)=sqrt(TIMECLUSTERS(iiii).area(tttt));
            TIMECLUSTERS(iiii).effective_radius(tttt)=...
                sqrt(TIMECLUSTERS(iiii).area(tttt)/pi);

            TIMECLUSTERS(iiii).volrain(tttt)=sum(CE.volrain(thisTimeCEID));

            TIMECLUSTERS(iiii).nclusters(tttt)=numel(thisTimeMatch);


            % The "Max Area" tracks.
            [MA,MaxAreaI]=max(TIMECLUSTERS(iiii).ce(tttt).area) ;

            TIMECLUSTERS(iiii).max_area_track.lon(tttt)=...
                CE.lon(thisTimeCEID(MaxAreaI));
            TIMECLUSTERS(iiii).max_area_track.lat(tttt)=...
                CE.lat(thisTimeCEID(MaxAreaI));

            TIMECLUSTERS(iiii).max_area_track.area(tttt)=...
                CE.area(thisTimeCEID(MaxAreaI));
            TIMECLUSTERS(iiii).max_are_tracka.size(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxAreaI)));
            TIMECLUSTERS(iiii).max_area_track.effective_radius(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxAreaI))/pi);
            TIMECLUSTERS(iiii).max_area_track.volrain(tttt)=...
                CE.volrain(thisTimeCEID(MaxAreaI));

            TIMECLUSTERS(iiii).max_area_track.lon_max(tttt)=...
                CE.lon_max(thisTimeCEID(MaxAreaI));
            TIMECLUSTERS(iiii).max_area_track.lat_max(tttt)=...
                CE.lat_max(thisTimeCEID(MaxAreaI));



            % The "Max Volrain" tracks.
            [MA,MaxVolrainI]=max(TIMECLUSTERS(iiii).ce(tttt).volrain) ;

            TIMECLUSTERS(iiii).max_volrain_track.lon(tttt)=...
                CE.lon(thisTimeCEID(MaxVolrainI));
            TIMECLUSTERS(iiii).max_volrain_track.lat(tttt)=...
                CE.lat(thisTimeCEID(MaxVolrainI));

            TIMECLUSTERS(iiii).max_volrain_track.area(tttt)=...
                CE.area(thisTimeCEID(MaxVolrainI));
            TIMECLUSTERS(iiii).max_are_tracka.size(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxVolrainI)));
            TIMECLUSTERS(iiii).max_volrain_track.effective_radius(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxVolrainI))/pi);
            TIMECLUSTERS(iiii).max_volrain_track.volrain(tttt)=...
                CE.volrain(thisTimeCEID(MaxVolrainI));

            TIMECLUSTERS(iiii).max_volrain_track.lon_max(tttt)=...
                CE.lon_max(thisTimeCEID(MaxVolrainI));
            TIMECLUSTERS(iiii).max_volrain_track.lat_max(tttt)=...
                CE.lat_max(thisTimeCEID(MaxVolrainI));




        end



        % The "bulk" cloud cluster properties.
        [TIMECLUSTERS(iiii).maxarea,TIMECLUSTERS(iiii).maxarea_indx]=...
            max(TIMECLUSTERS(iiii).area);

        TIMECLUSTERS(iiii).maxsize=...
            max(TIMECLUSTERS(iiii).size);

        TIMECLUSTERS(iiii).max_effective_radius=...
            max(TIMECLUSTERS(iiii).effective_radius);

        TIMECLUSTERS(iiii).maxvolrain=...
            max(TIMECLUSTERS(iiii).volrain);

        [TIMECLUSTERS(iiii).minlon,TIMECLUSTERS(iiii).minlon_indx]=...
            min(TIMECLUSTERS(iiii).lon);

        [TIMECLUSTERS(iiii).maxlon,TIMECLUSTERS(iiii).maxlon_indx]=...
            max(TIMECLUSTERS(iiii).lon);

        [TIMECLUSTERS(iiii).minlat,TIMECLUSTERS(iiii).minlat_indx]=...
            min(TIMECLUSTERS(iiii).lat);

        [TIMECLUSTERS(iiii).maxlat,TIMECLUSTERS(iiii).maxlat_indx]=...
            max(TIMECLUSTERS(iiii).lat);

        [TIMECLUSTERS(iiii).minabslat,TIMECLUSTERS(iiii).minabslat_indx]=...
            min(abs(TIMECLUSTERS(iiii).lat));

        [TIMECLUSTERS(iiii).maxabslat,TIMECLUSTERS(iiii).maxabslat_indx]=...
            max(abs(TIMECLUSTERS(iiii).lat));

        %% Zonal and meridional propagation speed.
        [FIT0,S,MU]=polyfit(TIMECLUSTERS(iiii).time,TIMECLUSTERS(iiii).lon,1) ;
        TIMECLUSTERS(iiii).zonal_propagation_speed = ...
            (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ;

        [FIT0,S,MU]=polyfit(TIMECLUSTERS(iiii).time,TIMECLUSTERS(iiii).lat,1) ;
        TIMECLUSTERS(iiii).meridional_propagation_speed = ...
            (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ;


        b=robustfit(TIMECLUSTERS(iiii).time,TIMECLUSTERS(iiii).lon) ;
        TIMECLUSTERS(iiii).zonal_propagation_speed_robust_fit = ...
            (b(2)) * 111000.0/(24.0*3600.0) ;


        if ( TIMECLUSTERS(iiii).zonal_propagation_speed > 0 )

            [min_lon,min_lon_indx]=min(TIMECLUSTERS(iiii).lon) ;
            [max_lon,max_lon_indx]=max(TIMECLUSTERS(iiii).lon(min_lon_indx+1:end)) ;

            [FIT0,S,MU]=polyfit(TIMECLUSTERS(iiii).time(min_lon_indx:min_lon_indx+max_lon_indx),...
                                TIMECLUSTERS(iiii).lon(min_lon_indx:min_lon_indx+max_lon_indx),1) ;

            TIMECLUSTERS(iiii).zonal_eastward_propagation_speed = ...
                (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ;


        else

            TIMECLUSTERS(iiii).zonal_eastward_propagation_speed = NaN ;

        end

        keepLat=find(TIMECLUSTERS(iiii).lat > -8.0 & ...
                     TIMECLUSTERS(iiii).lat < 8.0 ) ;

        [FIT0,S,MU]=polyfit(TIMECLUSTERS(iiii).time(keepLat),...
                            TIMECLUSTERS(iiii).lon(keepLat),1) ;
        TIMECLUSTERS(iiii).zonal_propagation_speed_8s8n = ...
            (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ;


    end

end

end %% End of parent function



function maskArrays = calcMaskArrays(TIMECLUSTERS, CE, DN, OPT) %FILTER_STANDARD_DEVIATION, ACCUMULATION_PERIOD)
% function: Short description
%
% Extended description
% Single point masks are all points within the LPT contour.
% Filter width masks include points within the filter radius (e.g., standard deviation) from LPT contour.
% Uses TIMECLUSTERS struct array
%   to calculate 3-D (lat/lon/time) mask arrays.


np = round(OPT.FILTER_STANDARD_DEVIATION / OPT.DX);
grid_nx = numel(OPT.LON);
grid_ny = numel(OPT.LAT);

maskArrays.all.mask_by_id = -1+zeros(numel(DN), numel(CE.grid.lat), numel(CE.grid.lon));
maskArrays.all.mask_by_id_with_filter = -1+zeros(numel(DN), numel(CE.grid.lat), numel(CE.grid.lon));
maskArrays.all.mask_by_id_with_accumulation = -1+zeros(numel(DN), numel(CE.grid.lat), numel(CE.grid.lon));
maskArrays.all.mask_by_id_with_filter_and_accumulation = -1+zeros(numel(DN), numel(CE.grid.lat), numel(CE.grid.lon));
maskArrays.individual = [];

for indx = 1:numel(TIMECLUSTERS)

  this_mask_array.mask_by_id = -1+zeros(numel(TIMECLUSTERS(indx).time), numel(CE.grid.lat), numel(CE.grid.lon));
  this_mask_array.mask_by_id_with_filter = -1+zeros(numel(TIMECLUSTERS(indx).time), numel(CE.grid.lat), numel(CE.grid.lon));

  for iii = 1:numel(TIMECLUSTERS(indx).time)

    tindx = find(DN > TIMECLUSTERS(indx).time(iii) - 0.0001 & ...
      DN < TIMECLUSTERS(indx).time(iii) + 0.0001);

    tindx2 = find(DN > TIMECLUSTERS(indx).time(iii) - 0.0001 - OPT.ACCUMULATION_PERIOD / 24.0 & ...
      DN < TIMECLUSTERS(indx).time(iii) + 0.0001);

    for ce = TIMECLUSTERS(indx).ce(iii);
      for iiii = 1:numel(ce.pixels.x)

        % Single point masks. These are easy.
        this_mask_array.mask_by_id(iii, ce.pixels.y(iiii), ce.pixels.x(iiii)) = TIMECLUSTERS(indx).ceid(iii);
        maskArrays.all.mask_by_id(tindx, ce.pixels.y(iiii), ce.pixels.x(iiii)) = indx;
        maskArrays.all.mask_by_id_with_accumulation(tindx2, ce.pixels.y(iiii), ce.pixels.x(iiii)) = indx;

      end
    end % for loop over CEs for this LPT

  end % for loop over timecluster entries


  this_mask_array.mask_by_id_with_filter = ...
    feature_spread(this_mask_array.mask_by_id, np);

  maskArrays.individual = [maskArrays.individual, this_mask_array];

end  % for

% Filter width masks: expand out from single point masks.

% Filter width masks.
maskArrays.all.mask_by_id_with_filter = ...
    feature_spread(maskArrays.all.mask_by_id, np);

maskArrays.all.mask_by_id_with_filter_and_accumulation = ...
    feature_spread(maskArrays.all.mask_by_id_with_accumulation, np);

end  % function




function array_out = feature_spread(array_in, np)



s = size(array_in);
array_last = array_in;
array_out = array_in;


if (np > 0)
  for nn = 1:np

    % shift left
    array_shift = circshift(array_last, -1, 3);
    array_shift(:,:,s(3)) = -1;
    array_out = max(array_out, array_shift);

    % shift right
    array_shift = circshift(array_last, 1, 3);
    array_shift(:,:,1) = -1;
    array_out = max(array_out, array_shift);


    % shift down
    array_shift = circshift(array_last, -1, 2);
    array_shift(:,s(2),:) = -1;
    array_out = max(array_out, array_shift);

    % shift up
    array_shift = circshift(array_last, 1, 2);
    array_shift(:,1,:) = -1;
    array_out = max(array_out, array_shift);

    array_last = array_out ;

  end
end

end  % function
