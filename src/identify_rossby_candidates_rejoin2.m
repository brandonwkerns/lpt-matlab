clear all
close all


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
	       'clumps_of_worms.rejoin2.txt'];

clumps = dlmread(clumps_file,'',1,0);


%% Directories
PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

EASTWARD_PROP_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/identify_eastward_propagation'];

eval(['!mkdir -p ', EASTWARD_PROP_DATA_DIR])

%%
%% Set the west propagation identification criteria here.
%%

zonal_propagation_threshold = 999.0 ; % full LPT track net speed, in m/s.
min_duration    = 3.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_westward_prop_zonal_speed    = 0.0 ; % Westward propagation portion, in m/s.
min_westward_prop_duration    = 3.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_westward_prop_duration_in_lat_band = 3.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_net_lon_propagation   = -999.0 ;%20.0 ; % in deg. longitude.
min_total_lon_propagation = -999.0 ;%20.0 ; % in deg. longitude.
min_total_westward_lon_propagation = 10.0 ;%20.0 ; % in deg. longitude.
max_abs_latitude = 30.0 ;% in deg. latitude. Westward propagation period must get this close to the Equator at some point.



% For output table files.
FMT=['%10d%10d%10d%10.2f%16.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d %10.2f ', ...
     '%20d%15d%11.2f%11.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d%20.1f%15.1f\n'];

header='      year     index     clump  duration  mean_zonal_spd       begin         end    volrain      wprop_begin_idx  wprop_end_idx  wprop_spd  wprop_dur wprop_begin   wprop_end     wprop_lon_begin  wprop_lon_end';

fn_rossby_lpt_list = [EASTWARD_PROP_DATA_DIR,'/rossby_lpt_list.rejoin2.txt'];
fid_rossby_lpt_list = fopen(fn_rossby_lpt_list,'w');
fprintf(fid_rossby_lpt_list, '%s\n', header);


count=0 ; % Keep count of west propagating systems.

for year1 = 1998:2018  ;
%for year1 = [2011]  ;

  
  year2=year1+1 ;

  yyyy1=num2str(year1) ;
  yyyy2=num2str(year2) ;
  
  y1_y2=[yyyy1,'_',yyyy2] ;
  
  disp(['########### ',y1_y2, ' ###########']) ;
  
  
  dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.rejoin2.mat']);
  G=load([PROCESSED_DATA_DIR,'/', dir0(1).name]) ;
  for iiii = 2:20    
    if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
      eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
    end
  end
  
    
  %% Get "clumps of worms" for this year.
  clump_idx_this_year = find(clumps(:,1) == year1);
  lptid_this_year = clumps(clump_idx_this_year, 2)';
  clump_num_this_year = clumps(clump_idx_this_year, 3)';
  
  for this_clump_num = [unique(clump_num_this_year)]

    disp(['----------- Clump #', num2str(this_clump_num), ' -----------'])
      
    lptid_for_this_clump = lptid_this_year(clump_num_this_year == this_clump_num);
    
    westward_propagation_metric = 0.0 * lptid_for_this_clump;
    westward_propagation_metric2 = 0.0 * lptid_for_this_clump;
    longest_west_prop_zonal_speed_all = 0.0 * lptid_for_this_clump;
    longest_west_prop_duration_all = 0.0 * lptid_for_this_clump;
    longest_west_prop_duration_in_lat_band_all = 0.0 * lptid_for_this_clump;
    idx_begin_west_prop_all = 0.0 * lptid_for_this_clump;
    idx_end_west_prop_all = 0.0 * lptid_for_this_clump;
    
    meets_rossby_criteria = 0.0 * lptid_for_this_clump;
    
    n = 0;
    
    
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
    
    rossby_candidates_list = [];
    rossby_candidates_list_count = 0;
    
				%lptid_for_this_clump=[21];
    
    for ii = [lptid_for_this_clump] %   1:numel(G.TIMECLUSTERS)
      
      
      n = n + 1;
      disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
      close all;
      
      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=GG.date(end)-GG.date(1) ;
      
      [GG.year,GG.month,GG.day]=datevec(GG.date) ;
      [GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
      [GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;
      
      disp([num2str(GG.year0),sprintf('%02d',GG.month0),sprintf('%02d',GG.day0),sprintf('%02d',GG.hour0), ...
	    ' to ', num2str(GG.year1),sprintf('%02d',GG.month1),sprintf('%02d',GG.day1),sprintf('%02d',GG.hour1)])
      
      
      %%
      %% Now I'm using duration of periods of west propagation.
      %% and "divide and conquer" eliminating the shorter duration periods with westward motion.
      %%
      
      [mask_net_eastward_propagation, spd_raw] = west_east_divide_and_conquer(G, year1, ii, do_plotting);
      
      
      %% Get the westward velocity and duration of each of the periods
      %% of westward propagation.
      
      statsW=regionprops(~mask_net_eastward_propagation, 'all') ;
      if (numel(statsW) > 0)
      	longest_west_propagating_indx = 0;
      				%longest_period = 0;
	
	west_prop = [];
      	for jjj = 1:numel(statsW)
	  
      	  west_prop(jjj).idx_begin_west_prop = min(statsW(jjj).SubarrayIdx{2});
      	  west_prop(jjj).idx_end_west_prop = max(statsW(jjj).SubarrayIdx{2});
	  
      	  west_prop(jjj).lon_west_prop = GG.lon(west_prop(jjj).idx_begin_west_prop:west_prop(jjj).idx_end_west_prop);
      	  west_prop(jjj).lat_west_prop = GG.lat(west_prop(jjj).idx_begin_west_prop:west_prop(jjj).idx_end_west_prop);
      	  west_prop(jjj).time_west_prop = GG.time(west_prop(jjj).idx_begin_west_prop:west_prop(jjj).idx_end_west_prop);
	  west_prop(jjj).time_west_prop_in_lat_band = west_prop(jjj).time_west_prop;
	  west_prop(jjj).time_west_prop_in_lat_band(west_prop(jjj).lat_west_prop < -1*max_abs_latitude | west_prop(jjj).lat_west_prop > max_abs_latitude) = [];
	  
	  
	  
      				%Linear Regression
      	  [FIT0,S,MU] = polyfit(west_prop(jjj).time_west_prop, west_prop(jjj).lon_west_prop,1) ;
          west_prop_zonal_speed = (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ; 
      	  west_prop_duration = max(west_prop(jjj).time_west_prop) - min(west_prop(jjj).time_west_prop);
	  west_prop_duration_in_lat_band = 0.125 * (sum(abs(west_prop(jjj).lat_west_prop) < max_abs_latitude) - 1);
	  
	  net_westward_lon_propagation = -1*(GG.lon(west_prop(jjj).idx_end_west_prop) - GG.lon(west_prop(jjj).idx_begin_west_prop));
	  min_dist_to_eq_west = min(abs(GG.lat(west_prop(jjj).idx_begin_west_prop:west_prop(jjj).idx_end_west_prop)));
	  %% Check for "Rossby" criteria
	  
	  total_lon_propagation = max(GG.lon)-min(GG.lon);
	  net_lon_propagation = GG.lon(end)-GG.lon(1);
	  
	  if (GG.duration > min_duration-0.001)
            if (  GG.zonal_propagation_speed < zonal_propagation_threshold & ...
		  total_lon_propagation > min_total_lon_propagation & ...
		  net_westward_lon_propagation > min_total_westward_lon_propagation & ...
		  min_dist_to_eq_west < max_abs_latitude & ...
		  west_prop_zonal_speed <= min_westward_prop_zonal_speed & ...
		  west_prop_duration >= min_westward_prop_duration & ...
		  west_prop_duration_in_lat_band >= min_westward_prop_duration_in_lat_band)
	      
	      disp(['Period ',num2str(jjj),' Meets ROSSBY criteria!'])
	      east_prop(jjj).meets_rossby_criteria = 1;
	      
	      %% Update rossby_candidates_list
	      rossby_candidates_list_count = rossby_candidates_list_count + 1;
	      rossby_candidates_list(rossby_candidates_list_count).lptid = ii;
	      rossby_candidates_list(rossby_candidates_list_count).idx_begin_west_prop = west_prop(jjj).idx_begin_west_prop;
	      rossby_candidates_list(rossby_candidates_list_count).idx_end_west_prop = west_prop(jjj).idx_end_west_prop;

	      rossby_candidates_list(rossby_candidates_list_count).west_prop_zonal_speed = west_prop_zonal_speed;
	      rossby_candidates_list(rossby_candidates_list_count).west_prop_duration = west_prop_duration;
	      rossby_candidates_list(rossby_candidates_list_count).west_prop_duration_in_lat_band = west_prop_duration_in_lat_band;
	      
	      rossby_candidates_list(rossby_candidates_list_count).westward_propagation_metric = west_prop_duration ;
	      rossby_candidates_list(rossby_candidates_list_count).westward_propagation_metric2 = sum(GG.area .* spd_raw); 
	      
 	      
	    else
	      west_prop(jjj).meets_rossby_criteria = 0;
	    end
	    
	  end
	  
	end %for jjj = 1:numel(statsE)

      end %if (numel(statsE) > 0)

	 	
    end % End loop over the lpts in this clump.

      
    if (rossby_candidates_list_count > 0)
      %% -----------------------------------------------
      %% Now I have a list of periods within the LPT system that satisfy
      %% the ROSSBY candidate criteria!
      %% However, many of these are duplicated in two or more LPT system branches.
      %% SO: For each ROSSBY candidate, we need to choose the "dominant" LPT system branch
      %% that contains it.
      %% -----------------------------------------------

      ARRAY=[];
      for num = 1:rossby_candidates_list_count
	ARRAY(num,1) = rossby_candidates_list(num).westward_propagation_metric;
	ARRAY(num,2) = rossby_candidates_list(num).westward_propagation_metric2;
	ARRAY(num,3) = num;
      end
      ARRAY_sorted = sortrows(ARRAY);
      candidate_number_sorted = fliplr(ARRAY_sorted(:,3)');
	

      %% Write out to the file in order of eastward propagation metrics, ignoring duplicates.
      date_begin_list = [];
      date_end_list = [];
      for num = [candidate_number_sorted]

	this_lptid = rossby_candidates_list(num).lptid ;
	this_idx_begin_west_prop = rossby_candidates_list(num).idx_begin_west_prop ;
	this_idx_end_west_prop = rossby_candidates_list(num).idx_end_west_prop ;
        west_prop_zonal_speed = rossby_candidates_list(num).west_prop_zonal_speed;  
      	west_prop_duration = rossby_candidates_list(num).west_prop_duration;
	west_prop_duration_in_lat_band = rossby_candidates_list(num).west_prop_duration_in_lat_band;
	
	GG=G.TIMECLUSTERS(this_lptid) ;
	GG.date=GG.time ;
	GG.size=sqrt(GG.area) ;
	GG.area=GG.area/1e4 ;
	GG.nentries=numel(GG.date) ;
	GG.duration=GG.date(end)-GG.date(1) ;
	
	[GG.year,GG.month,GG.day]=datevec(GG.date) ;
	[GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
	[GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;
	
	
	[GG.year0e,GG.month0e,GG.day0e,GG.hour0e]=datevec(GG.time(this_idx_begin_west_prop)) ;
	[GG.year1e,GG.month1e,GG.day1e,GG.hour1e]=datevec(GG.time(this_idx_end_west_prop)) ;

	this_date_begin = 1000000*GG.year0e + 10000*GG.month0e + 100*GG.day0e + GG.hour0e;
	this_date_end = 1000000*GG.year1e + 10000*GG.month1e + 100*GG.day1e + GG.hour1e;

	disp([num2str(this_date_begin), ' - ',num2str(this_date_end)])
	
	if (sum(this_date_begin >= date_begin_list & this_date_end <= date_end_list) == 0)
	  
	  date_begin_list = [date_begin_list,this_date_begin];
	  date_end_list = [date_end_list,this_date_end];

	  
	  disp(['--> ', num2str(this_lptid), ' from ', num2str(this_idx_begin_west_prop),' to ', num2str(this_idx_end_west_prop),' is chosen for clump #', num2str(this_clump_num),'. <--'])
      	  


	  %% Total longitude propagation
	  total_lon_propagation = max(GG.lon)-min(GG.lon);
	  net_lon_propagation = GG.lon(end)-GG.lon(1);
	  
	  disp(['Total lon span (westmost to westmost): ', num2str(total_lon_propagation)])
	  disp(['Net lon propagation (start to end): ', num2str(net_lon_propagation)])



      
	  disp([num2str(GG.year0),sprintf('%02d',GG.month0),sprintf('%02d',GG.day0),sprintf('%02d',GG.hour0), ' to ', num2str(GG.year1),sprintf('%02d',GG.month1),sprintf('%02d',GG.day1),sprintf('%02d',GG.hour1)])

      
	  disp(['West propagation from ', num2str(GG.year0e),sprintf('%02d',GG.month0e),sprintf('%02d',GG.day0e),sprintf('%02d',GG.hour0e), ' to ', num2str(GG.year1e),sprintf('%02d',GG.month1e),sprintf('%02d',GG.day1e),sprintf('%02d',GG.hour1e)])
		
	  net_westward_lon_propagation = GG.lon(this_idx_end_west_prop) - GG.lon(this_idx_begin_west_prop);
	  disp(['Net westward lon propagation: ', num2str(net_westward_lon_propagation)])
	  disp(['Net westward duration (days): ', num2str(west_prop_duration)])
	  disp(['Net westward duration in lat range: ', num2str(west_prop_duration_in_lat_band)])
	  disp(['Net westward zonal speed (m/s): ', num2str(west_prop_zonal_speed)])
	  
	  
	  
	  fprintf(fid_rossby_lpt_list,FMT,...
                  year1,this_lptid,this_clump_num,GG.duration,GG.zonal_propagation_speed,...
                  GG.year0,GG.month0,GG.day0,GG.hour0,...
                  GG.year1,GG.month1,GG.day1,GG.hour1,...
                  mean(GG.volrain)*GG.duration,...
                  this_idx_begin_west_prop, this_idx_end_west_prop, ...
            	  west_prop_zonal_speed, west_prop_duration, ...
                  GG.year0e,GG.month0e,GG.day0e,GG.hour0e,...
                  GG.year1e,GG.month1e,GG.day1e,GG.hour1e,...
		  GG.lon(this_idx_begin_west_prop),GG.lon(this_idx_end_west_prop));


	  count = count + 1;
	  
	end
      end

    end
  end	
end


fclose(fid_rossby_lpt_list);

disp('########## Done! ###########')
disp([num2str(count), ' ROSSBY candidates were found.'])
disp(fn_rossby_lpt_list)
