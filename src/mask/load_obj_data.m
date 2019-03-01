function OBJ = load_obj_data(data_dir, objid)
  %%
  %% Extract the LP Object struct for a specific objid.
  %% Usage:
  %%   OBJ = load_obj_data(data_dir, objid)  
  %%
  
  objid_str = num2str(objid);
  YYYY = objid_str(1:4);
  MM = objid_str(5:6);
  DD = objid_str(7:8);
  HH = objid_str(9:10);
  OBJ_all = load([data_dir, '/', YYYY, '/', MM, '/objects_', YYYY,MM,DD,HH,'.mat']);
  idx = find(OBJ_all.id == objid);

  fields = fieldnames(OBJ_all);
  for i = 1:numel(fields)
    if (~strcmp(fields{i}, 'grid'))
      OBJ.(fields{i}) = OBJ_all.(fields{i})(idx);
    end
  end
  OBJ.grid = OBJ_all.grid;
end

