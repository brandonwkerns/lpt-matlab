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
%% Set the east propagation identification criteria here.
%%

min_zonal_speed = -999.0 ; % full LPT track net speed, in m/s.
min_duration    = 7.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_eastward_prop_zonal_speed    = 0.0 ; % Eastward propagation portion, in m/s.
min_eastward_prop_duration    = 7.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_eastward_prop_duration_in_lat_band = 7.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_net_lon_propagation   = -999.0 ;%20.0 ; % in deg. longitude.
min_total_lon_propagation = -999.0 ;%20.0 ; % in deg. longitude.
min_total_eastward_lon_propagation = 10.0 ;%20.0 ; % in deg. longitude.
max_abs_latitude = 15.0 ;% in deg. latitude. Eastward propagation period must get this close to the Equator at some point.



% For output table files.
FMT=['%10d%10d%10d%10.2f%16.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d %10.2f ', ...
     '%20d%15d%11.2f%11.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d%20.1f%15.1f\n'];

header='      year     index     clump  duration  mean_zonal_spd       begin         end    volrain      eprop_begin_idx  eprop_end_idx  eprop_spd  eprop_dur eprop_begin   eprop_end     eprop_lon_begin  eprop_lon_end';

fn_mjo_lpt_list = [EASTWARD_PROP_DATA_DIR,'/mjo_lpt_list.rejoin2.txt'];
fid_mjo_lpt_list = fopen(fn_mjo_lpt_list,'w');
fprintf(fid_mjo_lpt_list, '%s\n', header);


count=0 ; % Keep count of east propagating systems.

for year1 = 1998:2018  ;
%for year1 = [2018]  ;

  
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
    
    eastward_propagation_metric = 0.0 * lptid_for_this_clump;
    eastward_propagation_metric2 = 0.0 * lptid_for_this_clump;
    longest_east_prop_zonal_speed_all = 0.0 * lptid_for_this_clump;
    longest_east_prop_duration_all = 0.0 * lptid_for_this_clump;
    longest_east_prop_duration_in_lat_band_all = 0.0 * lptid_for_this_clump;
    idx_begin_east_prop_all = 0.0 * lptid_for_this_clump;
    idx_end_east_prop_all = 0.0 * lptid_for_this_clump;
    
    meets_mjo_criteria = 0.0 * lptid_for_this_clump;
    
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
    
    mjo_candidates_list = [];
    mjo_candidates_list_count = 0;
    
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
      %% Now I'm using duration of periods of east propagation.
      %% and "divide and conquer" eliminating the shorter duration periods with westward motion.
      %%
      
      [mask_net_eastward_propagation, spd_raw] = west_east_divide_and_conquer(G, year1, ii, do_plotting);
      
      
      %% Get the eastward velocity and duration of each of the periods
      %% of eastward propagation.
      
      statsE=regionprops(mask_net_eastward_propagation, 'all') ;
      if (numel(statsE) > 0)
      	longest_east_propagating_indx = 0;
      				%longest_period = 0;
	
	east_prop = [];
      	for jjj = 1:numel(statsE)
	  
      	  east_prop(jjj).idx_begin_east_prop = min(statsE(jjj).SubarrayIdx{2});
      	  east_prop(jjj).idx_end_east_prop = max(statsE(jjj).SubarrayIdx{2});
	  
      	  east_prop(jjj).lon_east_prop = GG.lon(east_prop(jjj).idx_begin_east_prop:east_prop(jjj).idx_end_east_prop);
      	  east_prop(jjj).lat_east_prop = GG.lat(east_prop(jjj).idx_begin_east_prop:east_prop(jjj).idx_end_east_prop);
      	  east_prop(jjj).time_east_prop = GG.time(east_prop(jjj).idx_begin_east_prop:east_prop(jjj).idx_end_east_prop);
	  east_prop(jjj).time_east_prop_in_lat_band = east_prop(jjj).time_east_prop;
	  east_prop(jjj).time_east_prop_in_lat_band(east_prop(jjj).lat_east_prop < -1*max_abs_latitude | east_prop(jjj).lat_east_prop > max_abs_latitude) = [];
	  
	  
	  
      				%Linear Regression
      	  [FIT0,S,MU] = polyfit(east_prop(jjj).time_east_prop, east_prop(jjj).lon_east_prop,1) ;
          east_prop_zonal_speed = (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2));  
      	  east_prop_duration = max(east_prop(jjj).time_east_prop) - min(east_prop(jjj).time_east_prop);
	  east_prop_duration_in_lat_band = 0.125 * (sum(abs(east_prop(jjj).lat_east_prop) < max_abs_latitude) - 1);
	  
	  net_eastward_lon_propagation = GG.lon(east_prop(jjj).idx_end_east_prop) - GG.lon(east_prop(jjj).idx_begin_east_prop);
	  min_dist_to_eq_east = min(abs(GG.lat(east_prop(jjj).idx_begin_east_prop:east_prop(jjj).idx_end_east_prop)));
	  %% Check for MJO criteria
	  
	  total_lon_propagation = max(GG.lon)-min(GG.lon);
	  net_lon_propagation = GG.lon(end)-GG.lon(1);
	  
	  if (GG.duration > min_duration-0.001)
            if (  GG.zonal_propagation_speed > min_zonal_speed & ...
		  total_lon_propagation > min_total_lon_propagation & ...
		  net_eastward_lon_propagation > min_total_eastward_lon_propagation & ...
		  min_dist_to_eq_east < max_abs_latitude & ...
		  east_prop_zonal_speed >= min_eastward_prop_zonal_speed & ...
		  east_prop_duration >= min_eastward_prop_duration & ...
		  east_prop_duration_in_lat_band >= min_eastward_prop_duration_in_lat_band)
	      
	      disp(['Period ',num2str(jjj),' Meets MJO criteria!'])
	      east_prop(jjj).meets_mjo_criteria = 1;
	      
	      %% Update mjo_candidates_list
	      mjo_candidates_list_count = mjo_candidates_list_count + 1;
	      mjo_candidates_list(mjo_candidates_list_count).lptid = ii;
	      mjo_candidates_list(mjo_candidates_list_count).idx_begin_east_prop = east_prop(jjj).idx_begin_east_prop;
	      mjo_candidates_list(mjo_candidates_list_count).idx_end_east_prop = east_prop(jjj).idx_end_east_prop;

	      mjo_candidates_list(mjo_candidates_list_count).east_prop_zonal_speed = east_prop_zonal_speed;
	      mjo_candidates_list(mjo_candidates_list_count).east_prop_duration = east_prop_duration;
	      mjo_candidates_list(mjo_candidates_list_count).east_prop_duration_in_lat_band = east_prop_duration_in_lat_band;
	      
	      mjo_candidates_list(mjo_candidates_list_count).eastward_propagation_metric = east_prop_duration ;
	      mjo_candidates_list(mjo_candidates_list_count).eastward_propagation_metric2 = sum(GG.area .* spd_raw); 
	      
 	      
	    else
	      east_prop(jjj).meets_mjo_criteria = 0;
	    end
	    
	  end
	  
	end %for jjj = 1:numel(statsE)

      end %if (numel(statsE) > 0)

	 	
    end % End loop over the lpts in this clump.

      
    if (mjo_candidates_list_count > 0)
      %% -----------------------------------------------
      %% Now I have a list of periods within the LPT system that satisfy
      %% the MJO candidate criteria!
      %% However, many of these are duplicated in two or more LPT system branches.
      %% SO: For each MJO candidate, we need to choose the "dominant" LPT system branch
      %% that contains it.
      %% -----------------------------------------------

      ARRAY=[];
      for num = 1:mjo_candidates_list_count
	ARRAY(num,1) = mjo_candidates_list(num).eastward_propagation_metric;
	ARRAY(num,2) = mjo_candidates_list(num).eastward_propagation_metric2;
	ARRAY(num,3) = num;
      end
      ARRAY_sorted = sortrows(ARRAY);
      candidate_number_sorted = fliplr(ARRAY_sorted(:,3)');
	

      %% Write out to the file in order of eastward propagation metrics, ignoring duplicates.
      date_begin_list = [];
      date_end_list = [];
      for num = [candidate_number_sorted]

	this_lptid = mjo_candidates_list(num).lptid ;
	this_idx_begin_east_prop = mjo_candidates_list(num).idx_begin_east_prop ;
	this_idx_end_east_prop = mjo_candidates_list(num).idx_end_east_prop ;
        east_prop_zonal_speed = mjo_candidates_list(num).east_prop_zonal_speed;  
      	east_prop_duration = mjo_candidates_list(num).east_prop_duration;
	east_prop_duration_in_lat_band = mjo_candidates_list(num).east_prop_duration_in_lat_band;
	
	GG=G.TIMECLUSTERS(this_lptid) ;
	GG.date=GG.time ;
	GG.size=sqrt(GG.area) ;
	GG.area=GG.area/1e4 ;
	GG.nentries=numel(GG.date) ;
	GG.duration=GG.date(end)-GG.date(1) ;
	
	[GG.year,GG.month,GG.day]=datevec(GG.date) ;
	[GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
	[GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;
	
	
	[GG.year0e,GG.month0e,GG.day0e,GG.hour0e]=datevec(GG.time(this_idx_begin_east_prop)) ;
	[GG.year1e,GG.month1e,GG.day1e,GG.hour1e]=datevec(GG.time(this_idx_end_east_prop)) ;

	this_date_begin = 1000000*GG.year0e + 10000*GG.month0e + 100*GG.day0e + GG.hour0e;
	this_date_end = 1000000*GG.year1e + 10000*GG.month1e + 100*GG.day1e + GG.hour1e;

	disp([num2str(this_date_begin), ' - ',num2str(this_date_end)])
	
	if (sum(this_date_begin >= date_begin_list & this_date_end <= date_end_list) == 0)
	  
	  date_begin_list = [date_begin_list,this_date_begin];
	  date_end_list = [date_end_list,this_date_end];

	  
	  disp(['--> ', num2str(this_lptid), ' from ', num2str(this_idx_begin_east_prop),' to ', num2str(this_idx_end_east_prop),' is chosen for clump #', num2str(this_clump_num),'. <--'])
      	  


	  %% Total longitude propagation
	  total_lon_propagation = max(GG.lon)-min(GG.lon);
	  net_lon_propagation = GG.lon(end)-GG.lon(1);
	  
	  disp(['Total lon span (westmost to eastmost): ', num2str(total_lon_propagation)])
	  disp(['Net lon propagation (start to end): ', num2str(net_lon_propagation)])



      
	  disp([num2str(GG.year0),sprintf('%02d',GG.month0),sprintf('%02d',GG.day0),sprintf('%02d',GG.hour0), ' to ', num2str(GG.year1),sprintf('%02d',GG.month1),sprintf('%02d',GG.day1),sprintf('%02d',GG.hour1)])

      
	  disp(['East propagation from ', num2str(GG.year0e),sprintf('%02d',GG.month0e),sprintf('%02d',GG.day0e),sprintf('%02d',GG.hour0e), ' to ', num2str(GG.year1e),sprintf('%02d',GG.month1e),sprintf('%02d',GG.day1e),sprintf('%02d',GG.hour1e)])
		
	  net_eastward_lon_propagation = GG.lon(this_idx_end_east_prop) - GG.lon(this_idx_begin_east_prop);
	  disp(['Net eastward lon propagation: ', num2str(net_eastward_lon_propagation)])
	  disp(['Net eastward duration (days): ', num2str(east_prop_duration)])
	  disp(['Net eastward duration in lat range: ', num2str(east_prop_duration_in_lat_band)])
	  disp(['Net eastward zonal speed (m/s): ', num2str(east_prop_zonal_speed)])
	  
	  
	  
	  fprintf(fid_mjo_lpt_list,FMT,...
                  year1,this_lptid,this_clump_num,GG.duration,GG.zonal_propagation_speed,...
                  GG.year0,GG.month0,GG.day0,GG.hour0,...
                  GG.year1,GG.month1,GG.day1,GG.hour1,...
                  mean(GG.volrain)*GG.duration,...
                  this_idx_begin_east_prop, this_idx_end_east_prop, ...
            	  east_prop_zonal_speed, east_prop_duration, ...
                  GG.year0e,GG.month0e,GG.day0e,GG.hour0e,...
                  GG.year1e,GG.month1e,GG.day1e,GG.hour1e,...
		  GG.lon(this_idx_begin_east_prop),GG.lon(this_idx_end_east_prop));


	  count = count + 1;
	  
	end
      end

    end
  end	
end


fclose(fid_mjo_lpt_list);

disp('########## Done! ###########')
disp([num2str(count), ' MJO candidates were found.'])
disp(fn_mjo_lpt_list)
