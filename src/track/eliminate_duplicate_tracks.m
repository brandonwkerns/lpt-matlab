function NEWTIMECLUSTERS = eliminate_duplicate_tracks(TIMECLUSTERS, verbose)

  if (nargin < 2)
    verbose = false;
  end
  
  TC_eliminate_list = [];
  NEWTIMECLUSTERS = TIMECLUSTERS;

  for ii = 1:numel(NEWTIMECLUSTERS)
    for jj = setxor(ii, 1:numel(NEWTIMECLUSTERS))

      if (numel(NEWTIMECLUSTERS(ii).ceid) > 0 & ...
          numel(NEWTIMECLUSTERS(jj).ceid) > 0 & ...
          numel(NEWTIMECLUSTERS(ii).ceid) == numel(NEWTIMECLUSTERS(jj).ceid))

        if sum(abs(NEWTIMECLUSTERS(ii).ceid - NEWTIMECLUSTERS(jj).ceid)) == 0

          TC_eliminate_list = [TC_eliminate_list, jj];
          NEWTIMECLUSTERS(jj).ceid = [-999];

        end
      end

    end
  end

  if verbose
    disp(['Eliminating the following duplicate tracks: ', num2str(unique(TC_eliminate_list))])
  end
  NEWTIMECLUSTERS(TC_eliminate_list)=[];

end
