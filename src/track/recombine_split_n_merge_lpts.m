clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../../config')
options
save('temp.mat');
OPT = load('temp.mat');
eval('!rm temp.mat')

do_plotting=0;


%% Directories
PROCESSED_DATA_DIR = ['../../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

OBJECTS_DATA_DIR = ['../../data/',CASE_LABEL,'/processed/',...
                    'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                    sprintf('%d',ACCUMULATION_PERIOD), ...
                    'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/objects'];


%% Clumps of Worms
clumps_file = [PROCESSED_DATA_DIR, '/clumps_of_worms.txt'];
clumps = dlmread(clumps_file,'',1,0);



%%
%% Set the east propagation identification criteria here.
%%

max_n_splitting_times = 999 ; % Max splitting times to recombine.
n_splitting_times_collect = [];

for year1 = [2006] %1998:2017  ;

    year2=year1+1 ;

    yyyy1=num2str(year1) ;
    yyyy2=num2str(year2) ;

    y1_y2=[yyyy1,'_',yyyy2] ;

    disp(['########### ',y1_y2, ' ###########']) ;
    
    %% Read LPT systems
    dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.mat']);
    fn_in = [PROCESSED_DATA_DIR,'/', dir0(1).name];
    disp(fn_in)
    G = load(fn_in) ;

    for iiii = 2:30
      if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
	eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
      end
    end
    Gori = G;
    
    lpts_to_eliminate = [];





    %% Initialize an struct array for LP object (id and time) storage.
    OBJ_storage.objid = [];
    OBJ_storage.time = [];
    

    
    
    %% Get "clumps of worms" for this year.
    clump_idx_this_year = find(clumps(:,1) == year1);
    lptid_this_year = clumps(clump_idx_this_year, 2)';
    clump_num_this_year = clumps(clump_idx_this_year, 3)';
    
    %for this_clump_num = [3] %[unique(clump_num_this_year)]
    for this_clump_num = [unique(clump_num_this_year)]

      disp(['----------- Clump #', num2str(this_clump_num), ' -----------'])
      
      lptid_for_this_clump = lptid_this_year(clump_num_this_year == this_clump_num);

      if (numel(lptid_for_this_clump) < 1)
        continue
      end

      if do_plotting
	figure;
	subplot(1,2,1);	  
	for ii = [sort(lptid_for_this_clump)]
	  GG=G.TIMECLUSTERS(ii) ;	
	  plot(GG.lon, GG.time, 'b-o')
	  hold on
	end
      end

      Gori = G;

      count=0;
      more_to_do = 1;
      while more_to_do == 1
	more_to_do = 0;
	count=count+1;
	
	for ii = [sort(lptid_for_this_clump)] 

          disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
	
          GG=G.TIMECLUSTERS(ii) ;	

	  GGobjid=[GG.objid];
	  GGtimes = [];
	  for id = 1:numel(GGobjid)
	    %% Use the storage struct array if it has already been read.
	    %% Otherwise, read the OBJ time and add it to storage.
	    if (sum(OBJ_storage.objid == GGobjid(id)) > 0)
	      GGtimes = [GGtimes, OBJ_storage.time(OBJ_storage.objid == GGobjid(id))];
	    else
	      OBJ = load_obj_data(OBJECTS_DATA_DIR, GGobjid(id));
	      GGtimes = [GGtimes, OBJ.time];
	      
	      %% Append to storage so I don't need to read this one again.
	      OBJ_storage.time = [OBJ_storage.time, OBJ.time];
	      OBJ_storage.objid = [OBJ_storage.objid, GGobjid(id)];	      
	    end
	  end
	  
	  %%
	  %% Search for combos of tracks that have everything
	  %% in common *except* for a time period
	  %% that's not the beginning or end of the tracks.
	  %% A split somewhere in the middle
	  %% of the track that comes
	  %% back together into a single track.
	  %%
	  
	  other_lptids_in_this_clump = sort(setxor(lptid_for_this_clump, [ii]));
	  
	  for jj = [other_lptids_in_this_clump]
	    
	    HH=G.TIMECLUSTERS(jj) ;
	    
	    if (numel(HH.objid) < 1)
	      continue
	    end

	    HHobjid=[HH.objid];
	    HHtimes = [];
	    for id = 1:numel(HHobjid)
	      %% Use the storage struct array if it has already been read.
	      %% Otherwise, read the OBJ time and add it to storage.
	      if (sum(OBJ_storage.objid == HHobjid(id)) > 0)
		HHtimes = [HHtimes, OBJ_storage.time(OBJ_storage.objid == HHobjid(id))];
	      else
		OBJ = load_obj_data(OBJECTS_DATA_DIR, HHobjid(id));
		HHtimes = [HHtimes, OBJ.time];
		
		%% Append to storage so I don't need to read this one again.
		OBJ_storage.time = [OBJ_storage.time, OBJ.time];
		OBJ_storage.objid = [OBJ_storage.objid, HHobjid(id)];	      
	      end
	    end

	    %% To be a candidate for recombining, the following must apply:
	    %% 1. The have a time period of splitting.
	    %% 2. The splitting period is bounded by the intersecting time period.
	    
	    intersections = intersect([GG.objid], [HH.objid]);
	    splitting_obj_ids = setxor([GG.objid], [HH.objid]);

	    intersection_times = [];
	    for this_objid = [intersections]
	      intersection_times = [intersection_times, HHtimes(HHobjid == this_objid)];
	    end

	    if (numel(intersections) ~= numel(intersection_times))
	      disp([numel(intersections),numel(intersection_times)])
	      disp('WARNING: intersection times.')
	    end


	    splitting_times = [];
	    for this_objid = [splitting_obj_ids]
	      if sum(GGobjid == this_objid) > 0
		splitting_times = [splitting_times, GGtimes(GGobjid == this_objid)];
	      else
		splitting_times = [splitting_times, HHtimes(HHobjid == this_objid)];
	      end
	    end

	    if (numel(splitting_obj_ids) ~= numel(splitting_times))
	      disp([numel(splitting_obj_ids),numel(splitting_times)])
	      disp('WARNING: splitting times.')
	    end
	    
	    %% In some cases with multiple iterations, a time can be both a splitting
	    %% and an intersecting time. Still consider these times for LPT branch merger.
	    if (numel(splitting_times) < 1 | numel(intersection_times) < 1)
	      continue
	    end

	    for ttt = [intersection_times]
	      if (sum(splitting_times == ttt) > 0)
		intersection_times(intersection_times == ttt) = [];
	      end
	    end
	    

	    
	    %% If this pair has no splits or mergers, then skip it.
	    if (numel(splitting_times) < 1 | numel(intersection_times) < 1)
	      continue
	    end
	    
	    	    
	    split_and_merge_times = splitting_times(splitting_times > min(intersection_times) & splitting_times < max(intersection_times));
	    
	    split_and_merge_ids = [];
	    for tttt = split_and_merge_times
	      split_and_merge_ids = [split_and_merge_ids, GGobjid(GGtimes == tttt)];
	      split_and_merge_ids = [split_and_merge_ids, HHobjid(HHtimes == tttt)];
	    end
	    split_and_merge_ids = unique(split_and_merge_ids);
	    
	    if (numel(split_and_merge_ids) < 1)
	      continue
	    end

	    %% OK, a split-n-merge case is detected!
	    %% We won't know how many LPT system branches are affected,
	    %% So, first figure this out.
	    %% Then, we can assign the relevant object ids to each of those LPT branches.
	    
		
	    for jjjj = [other_lptids_in_this_clump]
	      if ( numel(intersect([G.TIMECLUSTERS(jjjj).objid], [split_and_merge_ids])) > 0)
		more_to_do = 1; % Assume there will be more cases, and trigger the loop again.
		G.TIMECLUSTERS(ii).objid = unique([G.TIMECLUSTERS(ii).objid, split_and_merge_ids]);
		G.TIMECLUSTERS(jjjj).objid = unique([G.TIMECLUSTERS(jjjj).objid, split_and_merge_ids]);
	      end
	    end
	    
	    if more_to_do == 1
	      disp('Starting a new iteration.')
	      break
	    end
	    
	  end
	  
	  %G.TIMECLUSTERS = calc_tracking_parameters(G.TIMECLUSTERS, OBJECTS);
	  	  
	  if more_to_do == 1
	    if (count == 10)
	      disp('Stuck in an infinite loop! Breaking out.')
	      more_to_do = 0;
	    end
	    %more_to_do = 0;  % Over-ride.
	    break
	  end
	  
	  
	end
      end

      %G.TIMECLUSTERS = calc_tracking_parameters(G.TIMECLUSTERS, OBJECTS_DATA_DIR);

      if do_plotting
	subplot(1,2,2);	  
	for iip = [sort(lptid_for_this_clump)]
	  GG=G.TIMECLUSTERS(iip) ;	
	  plot(GG.lon, GG.time, 'b-o')
	  hold on
	end
	drawnow;
      end
      
    end


    
    %% Take out duplicates.
    Gnew.TIMECLUSTERS = eliminate_overlapping_tracks(G.TIMECLUSTERS, 1) ;
    Gnew.TIMECLUSTERS = put_tracks_in_order(Gnew.TIMECLUSTERS, 1);
    
    %% Recalculate the tracking parameters.
    disp('Updating tracking parameters....')
    Gnew.TIMECLUSTERS = calc_tracking_parameters(Gnew.TIMECLUSTERS, OBJECTS_DATA_DIR);

    %% Output
    fn_out_base = [fn_in(1:end-4), '.rejoin'];    
    lpt_systems_output_netcdf(Gnew.TIMECLUSTERS, [fn_out_base,'.nc'], OPT);
    lpt_systems_output_ascii(Gnew.TIMECLUSTERS, [fn_out_base,'.txt']);
    lpt_systems_output_mat(Gnew.TIMECLUSTERS, [fn_out_base,'.mat']);

end
