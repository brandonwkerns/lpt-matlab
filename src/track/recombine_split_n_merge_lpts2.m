clear all
close all

%% Run this script *after* recombine_split_n_merge_lpts.m and identify_clumps_of_worms_rejoin.m.

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../../config')
options
save('temp.mat');
OPT = load('temp.mat');
eval('!rm temp.mat')

do_plotting=1;

%% Clumps of Worms
%clumps_file = ['../../data/',CASE_LABEL,'/processed/',...
%               'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
%               sprintf('%d',ACCUMULATION_PERIOD), ...
%               'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/identify_eastward_propagation/',...
%               'clumps_of_worms.rejoin.txt'];
clumps_file = './clumps_of_worms.rejoin.txt';

clumps = dlmread(clumps_file,'',1,0);


%% Directories
PROCESSED_DATA_DIR = './';

%{
PROCESSED_DATA_DIR = ['../../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];
%}
OBJECTS_DATA_DIR = ['../../data/',CASE_LABEL,'/processed/',...
                    'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                    sprintf('%d',ACCUMULATION_PERIOD), ...
                    'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/objects'];

%%
%% Set the east propagation identification criteria here.
%%

min_time_to_maintain_split = 3.0; % days
max_n_splitting_times = 999 ; % Max splitting times to recombine.
n_splitting_times_collect = [];

for year1 = [2011] %1998:2017  ;

    year2=year1+1 ;

    yyyy1=num2str(year1) ;
    yyyy2=num2str(year2) ;

    y1_y2=[yyyy1,'_',yyyy2] ;

    disp(['########### ',y1_y2, ' ###########']) ;


    %% Read LP Objects
    dir0 = dir([OBJECTS_DATA_DIR,'/objects_',num2str(year1),'*.mat']);
    disp([OBJECTS_DATA_DIR,'/', dir0(1).name])
    OBJECTS = load([OBJECTS_DATA_DIR,'/', dir0(1).name]) ;

    
    %% Read LPT systems
    dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.rejoin.mat']);
    %dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.mat']);
    fn_in = [PROCESSED_DATA_DIR,'/', dir0(1).name];
    disp(fn_in)
    G = load(fn_in) ;


    for iiii = 2:20
      if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
	eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
      end
    end
    Gori = G;

    lpts_to_eliminate = [];

    
    %% Get "clumps of worms" for this year.
    clump_idx_this_year = find(clumps(:,1) == year1);
    lptid_this_year = clumps(clump_idx_this_year, 2)';
    clump_num_this_year = clumps(clump_idx_this_year, 3)';
    
    for this_clump_num = [3] %[unique(clump_num_this_year)]
    %for this_clump_num = [unique(clump_num_this_year)]

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
      %{
      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp('MERGING')
      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
      count=0;
      more_to_do = 1;
      while more_to_do == 1
	more_to_do = 0;
	count=count+1;
	
	for ii = [sort(lptid_for_this_clump)] 

          disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
	
          GG=G.TIMECLUSTERS(ii) ;	

	  GGtimes=[GG.obj.time];
	  GGobjid=[GG.obj.objid];
	  
	  %%
	  %% Search for combos of tracks that have everything
	  %% in common *except* for a time period
	  %% that's not the beginning or end of the tracks.
	  %% A split somewhere in the middle
	  %% of the track that quickly comes
	  %% back together into a single track.
	  %%
	  
	  other_lptids_in_this_clump = sort(setxor(lptid_for_this_clump, [ii]));
	  
	  for jj = [other_lptids_in_this_clump]

	    HH=G.TIMECLUSTERS(jj) ;

	    if (numel(HH.objid) < 1)
	      continue
	    end

	    HHtimes=[HH.obj.time];
	    HHobjid=[HH.obj.objid];
	    
	    
	    %% To be a candidate for recombining, the following must apply:
	    %% 1. The have a time period of splitting and/or merging.
	    %% 2. The splitting/merging period is either at the beginning or end of the track.
	    %% 3. The splitting/merging period should be < 3 days.
	    
	    intersections = intersect([GG.objid], [HH.objid]);
	    union_ids = union([GG.objid], [HH.objid]);
	    %union_times = OBJECTS.time(union_ids);

	    %intersection_times = OBJECTS.time(intersections);
	    intersection_times = [];
	    for this_objid = [intersections]
	      %this_OBJ = load_obj_data(OBJECTS_DATA_DIR, this_objid);
	      intersection_times = [intersection_times, HHtimes(HHobjid == this_objid)];
	    end
	    
	    union_times = [];
	    for this_objid = [union_ids]
	      %this_OBJ = load_obj_data(OBJECTS_DATA_DIR, this_objid);
	      if sum(HHobjid == this_objid) > 0
		union_times = [union_times, HHtimes(HHobjid == this_objid)];
	      else
		union_times = [union_times, GGtimes(GGobjid == this_objid)];
	      end
	    end
	    %disp([numel(union_ids),numel(union_times)])
	    if (numel(union_ids) ~= numel(union_times))
	      disp('WARNING')
	    end

	    %union_times = unique(union_times);
	    
	    splitting_times = [];
	    for ttt = [union_times]
	      %idx_GG = intersect([GG.objid], find(OBJECTS.time == ttt));
	      %idx_HH = intersect([HH.objid], find(OBJECTS.time == ttt));
	      idx_GG = intersect([GG.objid], union_ids);
	      idx_HH = intersect([HH.objid], union_ids);

	      
	      if (numel(idx_GG) < 1 | numel(idx_HH) < 1)
		continue
	      end
	      if numel(setxor([idx_GG], [idx_HH])) > 0
		splitting_times = [splitting_times, ttt];
	      end
	    end

	    %% -- Left off here. --

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
	    
	    %% Mergers	    
	    split_and_merge_times = splitting_times(splitting_times < min(intersection_times));
	    if (numel(split_and_merge_times) < 1)
	      continue
	    end
	    
	    if (max(split_and_merge_times) - min(split_and_merge_times) > min_time_to_maintain_split)
	      continue
	    end
	    
	    
	    split_and_merge_ids = [];
	    for tttt = split_and_merge_times
	      %split_and_merge_ids = [split_and_merge_ids, intersect([GG.objid],[find(OBJECTS.time == tttt)])];
	      %split_and_merge_ids = [split_and_merge_ids, intersect([HH.objid],[find(OBJECTS.time == tttt)])];
	      split_and_merge_ids = [split_and_merge_ids, intersect([GG.objid],union_ids(union_times == tttt))];
	      split_and_merge_ids = [split_and_merge_ids, intersect([HH.objid],union_ids(union_times == tttt))];
	    end
	    split_and_merge_ids = unique(split_and_merge_ids);
	    
	    if (numel(split_and_merge_ids) < 1)
	      continue
	    end

	    

	    % Make sure splitting_obj_ids has entries from *both* tracks.
	    if numel(intersect([GG.objid], [split_and_merge_ids])) < 1
	      continue
	    end
	    if numel(intersect([HH.objid], [split_and_merge_ids])) < 1
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
	    if (count == 50)
	      disp('Stuck in an infinite loop! Breaking out.')
	      more_to_do = 0;
	    end
	    %more_to_do = 0;

	    break
	  end
	  
	  
	end
      end
      %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
      G.TIMECLUSTERS = calc_tracking_parameters(G.TIMECLUSTERS, OBJECTS_DATA_DIR);

      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp('SPLITTING')
      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      count=0;
      more_to_do = 1;
      while more_to_do == 1
	more_to_do = 0;
	count=count+1;
	
	for ii = [sort(lptid_for_this_clump)] 


          disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
	
          GG=G.TIMECLUSTERS(ii) ;	

	  GGtimes=[GG.obj.time];
	  GGobjid=[GG.obj.objid];
	  
	  %%
	  %% Search for combos of tracks that have everything
	  %% in common *except* for a time period
	  %% that's not the beginning or end of the tracks.
	  %% A split somewhere in the middle
	  %% of the track that quickly comes
	  %% back together into a single track.
	  %%
	  
	  other_lptids_in_this_clump = sort(setxor(lptid_for_this_clump, [ii]));
	  
	  for jj = [other_lptids_in_this_clump]

	    HH=G.TIMECLUSTERS(jj) ;
	    
	    if (numel(HH.objid) < 1)
	      continue
	    end

	    HHtimes=[HH.obj.time];
	    HHobjid=[HH.obj.objid];

	    
	    %% To be a candidate for recombining, the following must apply:
	    %% 1. The have a time period of splitting and/or merging.
	    %% 2. The splitting/merging period is either at the beginning or end of the track.
	    %% 3. The splitting/merging period should be < 3 days.
	    
	    intersections = intersect([GG.objid], [HH.objid]);
	    union_ids = union([GG.objid], [HH.objid]);
	    %union_times = OBJECTS.time(union_ids);

	    intersection_times = [];
	    for this_objid = [intersections]
	      %this_OBJ = load_obj_data(OBJECTS_DATA_DIR, this_objid);
	      intersection_times = [intersection_times, HHtimes(HHobjid == this_objid)];
	    end
	    
	    union_times = [];
	    for this_objid = [union_ids]
	      %this_OBJ = load_obj_data(OBJECTS_DATA_DIR, this_objid);
	      if sum(HHobjid == this_objid) > 0
		union_times = [union_times, HHtimes(HHobjid == this_objid)];
	      else
		union_times = [union_times, GGtimes(GGobjid == this_objid)];
	      end
	    end
	    %disp([numel(union_ids),numel(union_times)])
	    if (numel(union_ids) ~= numel(union_times))
	      disp('WARNING')
	    end
	    
	    %intersection_times = OBJECTS.time(intersections);


	    splitting_times = [];
	    for ttt = [union_times]
	      %idx_GG = intersect([GG.objid], find(OBJECTS.time == ttt));
	      %idx_HH = intersect([HH.objid], find(OBJECTS.time == ttt));

	      idx_GG = intersect([GG.objid], union_ids);
	      idx_HH = intersect([HH.objid], union_ids);

	      if (numel(idx_GG) < 1 | numel(idx_HH) < 1)
		continue
	      end
	      if numel(setxor([idx_GG], [idx_HH])) > 0
		splitting_times = [splitting_times, ttt];
	      end
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
	    %% Splits
	    split_and_merge_times = splitting_times(splitting_times > max(intersection_times));
	    if (numel(split_and_merge_times) < 1)
	      continue
	    end

	    if (max(split_and_merge_times) - min(split_and_merge_times) >  min_time_to_maintain_split)
	      continue
	    end
	    
	    
	    split_and_merge_ids = [];
	    for tttt = split_and_merge_times
	      %split_and_merge_ids = [split_and_merge_ids, intersect([GG.objid],[find(OBJECTS.time == tttt)])];
	      %split_and_merge_ids = [split_and_merge_ids, intersect([HH.objid],[find(OBJECTS.time == tttt)])];
	      split_and_merge_ids = [split_and_merge_ids, intersect([GG.objid],union_ids(union_times == tttt))];
	      split_and_merge_ids = [split_and_merge_ids, intersect([HH.objid],union_ids(union_times == tttt))];
	    end
	    split_and_merge_ids = unique(split_and_merge_ids);
	    
	    if (numel(split_and_merge_ids) < 1)
	      continue
	    end

	 % Make sure splitting_obj_ids has entries from *both* tracks.
	    if numel(intersect([GG.objid], [split_and_merge_ids])) < 1
	      continue
	    end
	    if numel(intersect([HH.objid], [split_and_merge_ids])) < 1
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

		
		disp('-----------')
		datevec(max(intersection_times))
		'with'
		datevec(max(splitting_times))

		max(split_and_merge_times) - min(split_and_merge_times)
		
	      end
	    end

	    
	    if more_to_do == 1
	      disp('Starting a new iteration.')
	    end

	    break
	    
	  end
	  	  	  
	  if more_to_do == 1
	    if (count == 50)
	      disp('Stuck in an infinite loop! Breaking out.')
	      more_to_do = 0;
	    end
	    %more_to_do = 0;
	    break
	  end

	  break
	  
	end
      end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %{
      
      G.TIMECLUSTERS = calc_tracking_parameters(G.TIMECLUSTERS, OBJECTS);
      disp('SPLIT-THEN-MERGE')


      

      count=0;
      more_to_do = 1;
      while more_to_do == 1
	more_to_do = 0;
	count=count+1;
	
	for ii = [sort(lptid_for_this_clump)] 


          disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
	
          GG=G.TIMECLUSTERS(ii) ;	
	  
	  %%
	  %% Search for combos of tracks that have everything
	  %% in common *except* for a time period
	  %% that's not the beginning or end of the tracks.
	  %% A split somewhere in the middle
	  %% of the track that quickly comes
	  %% back together into a single track.
	  %%
	  
	  other_lptids_in_this_clump = sort(setxor(lptid_for_this_clump, [ii]));
	  
	  for jj = [other_lptids_in_this_clump]

	    HH=G.TIMECLUSTERS(jj) ;
	    
	    if (numel(HH.objid) < 1)
	      continue
	    end
	    
	    %% To be a candidate for recombining, the following must apply:
	    %% 1. The have a time period of splitting and/or merging.
	    %% 2. The splitting/merging period is NOT at the beginning or end of the track.
	    
	    intersections = intersect([GG.objid], [HH.objid]);
	    union_ids = union([GG.objid], [HH.objid]);
	    union_times = OBJECTS.time(union_ids);

	    intersection_times = OBJECTS.time(intersections);
	    splitting_times = [];
	    for ttt = [union_times]
	      idx_GG = intersect([GG.objid], find(OBJECTS.time == ttt));
	      idx_HH = intersect([HH.objid], find(OBJECTS.time == ttt));
	      
	      %if (numel(idx_GG) < 1 | numel(idx_HH) < 1)
	      %  continue
	      %end
	      if (ttt < min(intersection_times) | ttt > max(intersection_times))
		continue
	      end
	      
	      if numel(setxor([idx_GG], [idx_HH])) > 0
		splitting_times = [splitting_times, ttt];
	      end
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
	    %% Splits
	    split_and_merge_times = splitting_times(splitting_times < max(intersection_times) & splitting_times > min(intersection_times));
	    if (numel(split_and_merge_times) < 1)
	      continue
	    end

	    %if (max(split_and_merge_times) - min(split_and_merge_times) >  min_time_to_maintain_split)
	    %  continue
	    %end
	    
	    
	    split_and_merge_ids = [];
	    for tttt = split_and_merge_times
	      split_and_merge_ids = [split_and_merge_ids, intersect([GG.objid],[find(OBJECTS.time == tttt)])];
	      split_and_merge_ids = [split_and_merge_ids, intersect([HH.objid],[find(OBJECTS.time == tttt)])];
	    end
	    split_and_merge_ids = unique(split_and_merge_ids);
	    
	    if (numel(split_and_merge_ids) < 1)
	      continue
	    end

	 % Make sure splitting_obj_ids has entries from *both* tracks.
	    
	    %% OK, a split-n-merge case is detected!
	    %% We won't know how many LPT system branches are affected,
	    %% So, first figure this out.
	    %% Then, we can assign the relevant object ids to each of those LPT branches.
	    
		
	    more_to_do = 1; % Assume there will be more cases, and trigger the loop again.
	    G.TIMECLUSTERS(ii).objid = unique([G.TIMECLUSTERS(ii).objid, split_and_merge_ids]);
	    G.TIMECLUSTERS(jj).objid = unique([G.TIMECLUSTERS(jj).objid, split_and_merge_ids]);
	    
	    
		%disp('-----------')
		%datevec(max(intersection_times))
		%'with'
		%datevec(max(splitting_times))

		%max(split_and_merge_times) - min(split_and_merge_times)
	    
	    
	    
	    if more_to_do == 1
	      disp('Starting a new iteration.')
	    end
	    
	  end
	  	  	  
	  if more_to_do == 1
	    if (count == 50)
	      disp('Stuck in an infinite loop! Breaking out.')
	      more_to_do = 0;
	    end
	    %more_to_do = 0;
	    break
	  end
	  
	end
      end


      %}

      
      G.TIMECLUSTERS = calc_tracking_parameters(G.TIMECLUSTERS, OBJECTS_DATA_DIR);

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
    %Gnew.TIMECLUSTERS = eliminate_duplicate_tracks(Gnew.TIMECLUSTERS, 1) ;
    Gnew.TIMECLUSTERS = put_tracks_in_order(Gnew.TIMECLUSTERS, 1);

    %% Recalculate the tracking parameters.
    disp('Updating tracking parameters....')
    Gnew.TIMECLUSTERS = calc_tracking_parameters(Gnew.TIMECLUSTERS, OBJECTS);

    %% Output
    fn_out_base = [fn_in(1:end-4), '2'];    
    lpt_systems_output_netcdf(Gnew.TIMECLUSTERS, OBJECTS, [fn_out_base,'.nc'], OPT);
    lpt_systems_output_ascii(Gnew.TIMECLUSTERS, [fn_out_base,'.txt']);
    lpt_systems_output_mat(Gnew.TIMECLUSTERS, OBJECTS, [fn_out_base,'.mat']);

end
