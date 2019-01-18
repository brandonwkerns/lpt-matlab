function lpt_system_mask_output_netcdf(netcdf_output_fn, mask_arrays)

  %% Open file
  cmode = netcdf.getConstant('CLOBBER');
  cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
  ncid = netcdf.create(netcdf_output_fn, cmode);
  
  %% Dimensions
  dimid_alltime  = netcdf.defDim(ncid, 'time', numel(mask_arrays.time));
  dimid_lon_grid  = netcdf.defDim(ncid, 'lon', numel(mask_arrays.grid.lon));
  dimid_lat_grid  = netcdf.defDim(ncid, 'lat', numel(mask_arrays.grid.lat));

  %% Variables
  varid_alltime  = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', [dimid_alltime]);
  varid_lon_grid  = netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', [dimid_lon_grid]);
  varid_lat_grid  = netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', [dimid_lat_grid]);

  
  varid_lpt_mask_by_id  = netcdf.defVar(ncid, 'mask_by_lpt_system_id', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_alltime]);
  varid_lpt_mask_by_ceid  = netcdf.defVar(ncid, 'mask_by_lp_object_id', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_alltime]);

  varid_lpt_mask_by_id_with_accumulation  = netcdf.defVar(ncid, 'mask_with_accumulation', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_alltime]);
  varid_lpt_mask_by_id_with_filter  = netcdf.defVar(ncid, 'mask_with_filter', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_alltime]);
  varid_lpt_mask_by_id_with_filter_and_accumulation  = netcdf.defVar(ncid, 'mask_with_filter_and_accumulation', 'NC_INT', [dimid_lon_grid, dimid_lat_grid, dimid_alltime]);

  
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id,true,true,1);
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_ceid,true,true,1);
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_filter,true,true,1);
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_accumulation,true,true,1);
  netcdf.defVarDeflate(ncid,varid_lpt_mask_by_id_with_filter_and_accumulation,true,true,1);

  netcdf.endDef(ncid)


  %% Values
  netcdf.putVar(ncid, varid_alltime, 24.0 * (mask_arrays.time - datenum(1970,1,1,0,0,0)));
  netcdf.putVar(ncid, varid_lon_grid, mask_arrays.grid.lon);
  netcdf.putVar(ncid, varid_lat_grid, mask_arrays.grid.lat);


  netcdf.putVar(ncid, varid_lpt_mask_by_id, ...
		permute(mask_arrays.mask_by_lptid, [3,2,1]));
  netcdf.putVar(ncid, varid_lpt_mask_by_ceid, ...
		permute(mask_arrays.mask_by_objid, [3,2,1]));
  
  netcdf.putVar(ncid, varid_lpt_mask_by_id_with_accumulation, ...
		permute(mask_arrays.mask_by_lptid_with_accumulation, [3,2,1]));
  netcdf.putVar(ncid, varid_lpt_mask_by_id_with_filter, ...
		permute(mask_arrays.mask_by_lptid_with_filter, [3,2,1]));
  netcdf.putVar(ncid, varid_lpt_mask_by_id_with_filter_and_accumulation, ...
		permute(mask_arrays.mask_by_lptid_with_filter_and_accumulation, [3,2,1]));

  netcdf.close(ncid)

  %% Attributes
  ncwriteatt(netcdf_output_fn,'/','creation_date',datestr(now));

  DISC = ['This file contains space - time masks for a single LPT system. Time starts from the beginning of the first accumulation period to the end of the last accumulation period. "mask_by_lpt_system_id" is the LPT contour, with time as the END of the accumulation period (hence, first few times do not have anything). "mask_by_lp_object_id" is the same as the previous one, but the values are the LP Object ID, not the LPT ID. "with_filter" has the filter width added to the edges of the blobs. Finally, "accumulation" indicates that the entire accumulation period contributing to the blob is in the blob. This can be used for example to identify all grid cells contributing to the LPT contour.'];
	  
  ncwriteatt(netcdf_output_fn,'/','description',DISC);
  
  
  ncwriteatt(netcdf_output_fn,'time','units','hours since 1970-1-1 0:0:0');
  ncwriteatt(netcdf_output_fn,'lon','units','degrees_east');
  ncwriteatt(netcdf_output_fn,'lat','units','degrees_north');

  ncwriteatt(netcdf_output_fn,'time','description','For mask_by_lpt_system_id, mask_by_lp_object_id, and mask_with_filter, this is the END of accumulation period time. For mask_with_accumulation and mask_with_filter_and_accumulation, it is the instantaneous time as the mask is "spread out" over the accumulation time.');

  
end

