function NEWTIMECLUSTERS = eliminate_overlapping_tracks(TIMECLUSTERS, verbose)

  if (nargin < 2)
    verbose = false;
  end

  if verbose
    disp('Identifying any remaining overlapping tracks... This may take a while.')
  end
  
  TC_eliminate_list = [];
  NEWTIMECLUSTERS = TIMECLUSTERS;

  for ii = 1:numel(NEWTIMECLUSTERS)

    %% Skip it if already eliminated.
    if (sum(TC_eliminate_list == ii) > 0)
      continue
    end
    
    for jj = setxor(ii, 1:numel(NEWTIMECLUSTERS))

      %% Skip it if already eliminated.
      if (sum(TC_eliminate_list == jj) > 0)
	continue
      end
      
      if (numel(NEWTIMECLUSTERS(ii).objid) > 0 & ...
          numel(NEWTIMECLUSTERS(jj).objid) > 0 )

	if numel(intersect(NEWTIMECLUSTERS(ii).objid, NEWTIMECLUSTERS(jj).objid)) < 1
	  continue
	end
	
        if all(ismember(NEWTIMECLUSTERS(ii).objid, NEWTIMECLUSTERS(jj).objid)) | all(ismember(NEWTIMECLUSTERS(jj).objid, NEWTIMECLUSTERS(ii).objid))

	  NEWTIMECLUSTERS(ii).objid = unique([NEWTIMECLUSTERS(ii).objid, NEWTIMECLUSTERS(jj).objid]);
          TC_eliminate_list = [TC_eliminate_list, jj];
          NEWTIMECLUSTERS(jj).objid = [-999];

        end
      end

    end
  end

  if verbose
    disp(['Eliminating the following overlapping tracks: ', num2str(unique(TC_eliminate_list))])
  end
  NEWTIMECLUSTERS(TC_eliminate_list)=[];

end
