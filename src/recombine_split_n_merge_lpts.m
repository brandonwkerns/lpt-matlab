clear all
close all

% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options
save('temp.mat');
OPT = load('temp.mat');
eval('!rm temp.mat')

%do_plotting = true;
do_plotting = false;



% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options


%% Clumps of Worms
clumps_file = ['../data/',CASE_LABEL,'/processed/',...
               'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
               sprintf('%d',ACCUMULATION_PERIOD), ...
               'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/identify_eastward_propagation/',...
               'clumps_of_worms.txt'];

clumps = dlmread(clumps_file,'',1,0);


%% Directories
PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

OBJECTS_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                    'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                    sprintf('%d',ACCUMULATION_PERIOD), ...
                    'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/objects'];

%%
%% Set the east propagation identification criteria here.
%%

max_n_splitting_times = 999 ; % Max splitting times to recombine.
n_splitting_times_collect = [];

for year1 = 1998:2017  ;

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
    dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.mat']);
    fn_in = [PROCESSED_DATA_DIR,'/', dir0(1).name];
    disp(fn_in)
    G = load(fn_in) ;


    for iiii = 2:20

      if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
	eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
      end
      
    end


    lpts_to_eliminate = [];

    
    %% Get "clumps of worms" for this year.
    clump_idx_this_year = find(clumps(:,1) == year1);
    lptid_this_year = clumps(clump_idx_this_year, 2)';
    clump_num_this_year = clumps(clump_idx_this_year, 3)';

    for this_clump_num = [unique(clump_num_this_year)]

      disp(['----------- Clump #', num2str(this_clump_num), ' -----------'])
      
      lptid_for_this_clump = lptid_this_year(clump_num_this_year == this_clump_num);

      %% Remove "duplicate" tracks.
      i_remove = [];

      iii = 0;
      for ii = [lptid_for_this_clump] %   1:numel(G.TIMECLUSTERS)
        iii = iii + 1;
        GG=G.TIMECLUSTERS(ii) ;

        [syear, smon, sday, shour] = datevec(GG.time(1));

        if (syear == year1 & smon == 6 & sday == 1 & shour == 0)
          i_remove = unique([i_remove, iii]);
        end

        if (syear == year1+1 & smon == 6)
          i_remove = unique([i_remove, iii]);
        end

      end
      lptid_for_this_clump(i_remove) = [];
      if (numel(lptid_for_this_clump) < 1)
        continue
      end

				%lptid_for_this_clump=[21];
      for ii = [lptid_for_this_clump] %   1:numel(G.TIMECLUSTERS)

        disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
	
        GG=G.TIMECLUSTERS(ii) ;

	%% It may have already been eliminated. If so, skip it.
	if (numel(GG.ceid) < 1)
	  disp("Allready eliminated. Moving on.")
	  continue
	end
	
	
        [GG.year,GG.month,GG.day]=datevec(GG.time) ;
        [GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
        [GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;
        
        disp([num2str(GG.year0),sprintf('%02d',GG.month0),sprintf('%02d',GG.day0),sprintf('%02d',GG.hour0), ...
              ' to ', num2str(GG.year1),sprintf('%02d',GG.month1),sprintf('%02d',GG.day1),sprintf('%02d',GG.hour1)])
      

	%%
	%% Search for combos of tracks that have everything
	%% in common *except* for a time period
	%% that's not the beginning or end of the tracks.
	%% A split somewhere in the middle
	%% of the track that quickly comes
	%% back together into a single track.
	%%

	other_lptids_in_this_clump = setxor(lptid_for_this_clump, ii);
	
	for jj = [other_lptids_in_this_clump]

	  HH=G.TIMECLUSTERS(jj) ;

	  if (numel(HH.ceid) < 1)
	    continue
	  end

	  %% To be a candidate for recombining, the following must apply:
	  %% 1. The tracks begin and end at the same time.
	  %% 2. The tracks begin and and with the same CEID.

	  if (GG.time(1) == HH.time(1) & GG.time(end) == HH.time(end))
	    if (GG.ceid(1) == HH.ceid(1) & GG.ceid(end) == HH.ceid(end))
	  
	      intersections = intersect(GG.ceid, HH.ceid); 

	      n_splitting_times = numel(GG.ceid) - numel(intersections);
	      n_splitting_times_collect = [n_splitting_times_collect, ...
					   n_splitting_times];
	      
	      disp([num2str(n_splitting_times), ' splitting times.'])
	      
	      if (n_splitting_times < max_n_splitting_times)
		
		disp([num2str(ii), ' and ', num2str(jj) ' overlap.'])
		

		G.TIMECLUSTERS(ii).ceid = unique([GG.ceid, HH.ceid]);
		G.TIMECLUSTERS(jj).ceid = [];
		lpts_to_eliminate	= [lpts_to_eliminate, jj];
	    
	      end
	    end
	    
	  end

	end
		
      end
    end

    %% Take out the identified duplicates.
    disp(['Eliminating: ', num2str(numel(lpts_to_eliminate)),...
	  ' of ', num2str(numel(G.TIMECLUSTERS)), ' LPTs.'])

    Gnew = G ;
    Gnew.TIMECLUSTERS(lpts_to_eliminate) = [];

    %% Recalculate the tracking parameters.
    disp('Updating tracking parameters....')
    Gnew.TIMECLUSTERS = calc_tracking_parameters(Gnew.TIMECLUSTERS, OBJECTS);

    %% Output
    fn_out_base = [fn_in(1:end-4), '.rejoin'];    
    lpt_systems_output_netcdf(Gnew.TIMECLUSTERS, OBJECTS, [fn_out_base,'.nc'], OPT);
    lpt_systems_output_ascii(Gnew.TIMECLUSTERS, [fn_out_base,'.txt']);
    lpt_systems_output_mat(Gnew.TIMECLUSTERS, OBJECTS, [fn_out_base,'.mat']);

end
