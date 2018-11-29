function lpt_systems_output_netcdf(systems, objects, ncfile, OPT)

  eval(['!rm -f ',ncfile])
  %% NetCDF output for "systems" to the specified file "ncfile."
  %% lp_systems_output_netcdf(systems, ncfile)
  
  TIMECLUSTERS = systems;
  CE = objects;
  netcdf_output_fn = ncfile; %['TIMECLUSTERS_lpt_',ymd0_ymd9,'.lptALL.nc'] ;
  DN=OPT.DN1:datenum(0,0,0,OPT.DT,0,0):OPT.DN2;

  %%
  %% Global Stuff
  %%

				% Get "stitch" track data.
  MISSING = -9999.0;
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
  disp(['Writing NetCDF: ', netcdf_output_fn])
  cmode = netcdf.getConstant('CLOBBER');
  cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
  ncid = netcdf.create(netcdf_output_fn, cmode);
  
  dimid_lpt_id  = netcdf.defDim(ncid, 'lpt_id', numel(TIMECLUSTERS));
  dimid_alltime  = netcdf.defDim(ncid, 'alltime', numel(DN));
  dimid_lon_grid  = netcdf.defDim(ncid, 'lon', numel(CE.grid.lon));
  dimid_lat_grid  = netcdf.defDim(ncid, 'lat', numel(CE.grid.lat));
  dimid_obs  = netcdf.defDim(ncid, 'obs', numel(stitch_lon));
  

				% Global Vars
				% Coordinates
  varid_lpt_id = netcdf.defVar(ncid, 'lpt_id', 'NC_INT', dimid_lpt_id);
  varid_obs = netcdf.defVar(ncid, 'obs', 'NC_INT', dimid_obs);
  varid_alltime = netcdf.defVar(ncid, 'alltime', 'NC_DOUBLE', dimid_alltime);
  varid_end_time = netcdf.defVar(ncid, 'end_of_accumulation_alltime', 'NC_DOUBLE', dimid_alltime);
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

  %% Data mode
  % Coordinates
  netcdf.putVar(ncid, varid_lpt_id, 1:numel(TIMECLUSTERS));
  netcdf.putVar(ncid, varid_alltime, 24.0 * (DN - datenum(1970,1,1,0,0,0) - 0.5*OPT.ACCUMULATION_PERIOD/24.0));
  netcdf.putVar(ncid, varid_end_time, 24.0 * (DN - datenum(1970,1,1,0,0,0)));
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
  netcdf.putVar(ncid, varid_stitch_time, 24.0 * (stitch_time - datenum(1970,1,1,0,0,0) - 0.5*OPT.ACCUMULATION_PERIOD/24.0));
  netcdf.putVar(ncid, varid_stitch_end_time, 24.0 * (stitch_time - datenum(1970,1,1,0,0,0)));
  netcdf.putVar(ncid, varid_stitch_lon, stitch_lon);
  netcdf.putVar(ncid, varid_stitch_lat, stitch_lat);
  netcdf.putVar(ncid, varid_stitch_area, stitch_area);
  netcdf.putVar(ncid, varid_stitch_size, stitch_size);
  netcdf.putVar(ncid, varid_stitch_effective_radius, stitch_effective_radius);
  netcdf.putVar(ncid, varid_stitch_volrain, stitch_volrain);
  
  netcdf.close(ncid)

				% Attributes
  ncwriteatt(netcdf_output_fn,'lpt_id','units','1');
  ncwriteatt(netcdf_output_fn,'lon','units','degrees_east');
  ncwriteatt(netcdf_output_fn,'lat','units','degrees_north');
  ncwriteatt(netcdf_output_fn,'centroid_lon','units','degrees_east');
  ncwriteatt(netcdf_output_fn,'centroid_lat','units','degrees_north');
  ncwriteatt(netcdf_output_fn,'alltime','units','hours since 1970-1-1 0:0:0');
  ncwriteatt(netcdf_output_fn,'alltime','description','Entire tracking period. Middle of accumulation period time.');
  ncwriteatt(netcdf_output_fn,'end_of_accumulation_alltime','units','hours since 1970-1-1 0:0:0');
  ncwriteatt(netcdf_output_fn,'end_of_accumulation_alltime','description','Entire tracking period. End of accumulation period time.');
  ncwriteatt(netcdf_output_fn,'time','units','hours since 1970-1-1 0:0:0');
  ncwriteatt(netcdf_output_fn,'time','description','Stitched time of LPTs. Middle of accumulation period time.');
  ncwriteatt(netcdf_output_fn,'end_of_accumulation_time','units','hours since 1970-1-1 0:0:0');
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


  
