function NEWTIMECLUSTERS = put_tracks_in_order(TIMECLUSTERS, verbose);

  % Put TIMECLUSTERS tracks in order by starting objid.

  starting_objid_list = [];
  for ii = 1:numel(TIMECLUSTERS)
    starting_objid_list(ii) = TIMECLUSTERS(ii).objid(1);
  end

  [sort_objid_list, sort_objid_list_indx] = sort(starting_objid_list);

  if (verbose)
    disp(['Sorting by order of starting objid: ', num2str(sort_objid_list_indx)])
  end

  NEWTIMECLUSTERS = TIMECLUSTERS(sort_objid_list_indx);
end %local function

