function NEWTIMECLUSTERS = eliminate_overlapping_tracks(TIMECLUSTERS, verbose)

  if (nargin < 2)
    verbose = false;
  end
  
  TC_eliminate_list = [];
  NEWTIMECLUSTERS = TIMECLUSTERS;

  for ii = 1:numel(NEWTIMECLUSTERS)
    for jj = setxor(ii, 1:numel(NEWTIMECLUSTERS))

      if (numel(NEWTIMECLUSTERS(ii).ceid) > 0 & ...
          numel(NEWTIMECLUSTERS(jj).ceid) > 0 )

        if all(ismember(NEWTIMECLUSTERS(ii).ceid, NEWTIMECLUSTERS(jj).ceid)) | all(ismember(NEWTIMECLUSTERS(jj).ceid, NEWTIMECLUSTERS(ii).ceid))

	  NEWTIMECLUSTERS(ii).ceid = unique([NEWTIMECLUSTERS(ii).ceid, NEWTIMECLUSTERS(jj).ceid]);
          TC_eliminate_list = [TC_eliminate_list, jj];
          NEWTIMECLUSTERS(jj).ceid = [-999];

        end
      end

    end
  end

  if verbose
    disp(['Eliminating the following overlapping tracks: ', num2str(unique(TC_eliminate_list))])
  end
  NEWTIMECLUSTERS(TC_eliminate_list)=[];

end
