clear all
close all


do_plotting = true;
%do_plotting = false;



% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options


%% Clumps of Worms
clumps_file = ['../data/',CASE_LABEL,'/processed/',...
               'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
               sprintf('%d',ACCUMULATION_PERIOD), ...
               'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/identify_eastward_propagation/',...
	       'clumps_of_worms.rejoin.txt'];

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


%%
%%
%%

% For output table files.
FMT=['%10d%10d%10d%10.2f%16.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d %10.2f ', ...
     '%20d%15d%11.2f%11.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d%20.1f%15.1f\n'];

header='      year     index     clump  duration  mean_zonal_spd       begin         end    volrain      eprop_begin_idx  eprop_end_idx  eprop_spd  eprop_dur eprop_begin   eprop_end     eprop_lon_begin  eprop_lon_end';

fn_mjo_lpt_list = [EASTWARD_PROP_DATA_DIR,'/mjo_lpt_list.rejoin.txt'];
fid_mjo_lpt_list = fopen(fn_mjo_lpt_list,'w');
fprintf(fid_mjo_lpt_list, '%s\n', header);


count=0 ; % Keep count of east propagating systems.

for year1 = [2004] %1998:2017  ;

  
    year2=year1+1 ;

    yyyy1=num2str(year1) ;
    yyyy2=num2str(year2) ;

    y1_y2=[yyyy1,'_',yyyy2] ;

    disp(['########### ',y1_y2, ' ###########']) ;

    
    dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.rejoin.mat']);
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

    for this_clump_num = [20] %[unique(clump_num_this_year)]

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

	

	%% Get the eastward velocity and duration of the longest continuous period
	%% of eastward propagation.

	statsE=regionprops(mask_net_eastward_propagation, 'all') ;
	if (numel(statsE) > 0)
      	  longest_east_propagating_indx = 0;
      	  longest_period = 0;
      	  for (jjj = 1:numel(statsE))
      	    if (statsE(jjj).Area > longest_period)
      	      longest_period = statsE(jjj).Area;
      	      longest_east_propagating_indx = jjj;
      	    end
      	  end
	  
      	  idx_begin_east_prop = min(statsE(longest_east_propagating_indx).SubarrayIdx{2});
      	  idx_end_east_prop = max(statsE(longest_east_propagating_indx).SubarrayIdx{2});
	  
      	  lon_east_prop = GG.lon(idx_begin_east_prop:idx_end_east_prop);
      	  lat_east_prop = GG.lat(idx_begin_east_prop:idx_end_east_prop);
      	  time_east_prop = GG.time(idx_begin_east_prop:idx_end_east_prop);
	  time_east_prop_in_lat_band = time_east_prop;
	  %time_east_prop_in_lat_band(lat_east_prop < -1*max_abs_latitude | lat_east_prop > max_abs_latitude) = [];

	  
	  
      				%Linear Regression
      	  [FIT0,S,MU] = polyfit(time_east_prop,lon_east_prop,1) ;
          longest_east_prop_zonal_speed = (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2));
	  
      	  longest_east_prop_duration = max(time_east_prop) - min(time_east_prop);
	  longest_east_prop_duration_in_lat_band = 0.125 * sum(abs(lat_east_prop) < max_abs_latitude);
	  
	  net_eastward_lon_propagation = GG.lon(idx_end_east_prop) - GG.lon(idx_begin_east_prop);
	  min_dist_to_eq_east = min(abs(GG.lat(idx_begin_east_prop:idx_end_east_prop)));
	  
	  eastward_propagation_metric(n) = longest_east_prop_duration ;
	  eastward_propagation_metric2(n) = sum(GG.area .* spd_raw); %sum(GG.area(idx_end_east_prop:end).*spd_raw(idx_end_east_prop:end)) + sum(GG.area(1:idx_begin_east_prop).*spd_raw(1:idx_begin_east_prop)) ;

	else
      	  longest_east_prop_zonal_speed = -999.0;
      	  longest_east_prop_duration = -999.0;
	  longest_east_prop_duration_in_lat_band = -999.0;
	  
      	  idx_begin_east_prop = -1;
      	  idx_end_east_prop = -1;
	  net_eastward_lon_propagation = -9999.0;
	  min_dist_to_eq_east = 9999.0;
	  
	  eastward_propagation_metric(n) = -999.0 ;
	  eastward_propagation_metric2(n) = -999.0 ;

	end
	
	%% Assign east prop metric here.
	
      	longest_east_prop_zonal_speed_all(n) = longest_east_prop_zonal_speed;
      	longest_east_prop_duration_all(n) = longest_east_prop_duration;
      	longest_east_prop_duration_in_lat_band_all(n) = longest_east_prop_duration_in_lat_band;
	idx_begin_east_prop_all(n) = idx_begin_east_prop;
	idx_end_east_prop_all(n) = idx_end_east_prop;



	%% Does it meet MJO criteria?

	total_lon_propagation = max(GG.lon)-min(GG.lon);
	net_lon_propagation = GG.lon(end)-GG.lon(1);
      
	if (GG.duration > min_duration-0.001)
          if (  GG.zonal_propagation_speed > min_zonal_speed & ...
             total_lon_propagation > min_total_lon_propagation & ...
             net_eastward_lon_propagation > min_total_eastward_lon_propagation & ...
             min_dist_to_eq_east < max_abs_latitude & ...
             longest_east_prop_zonal_speed > min_eastward_prop_zonal_speed & ...
             longest_east_prop_duration > min_eastward_prop_duration & ...
	     longest_east_prop_duration_in_lat_band > min_eastward_prop_duration_in_lat_band)

	    disp('YES! Meets MJO criteria!')
	    meets_mjo_criteria(n) = 1.0;
	  else
	    disp('NOPE! Does not meet MJO criteria.')
	  end
	  
	end
	
	
      end % End loop over the lpts in this clump.
      

      %% -----------------------------------------------
      %% Choose the LPT object to represent the "clump."
      %% Each clump can only contribute ONE MJO!
      %% -----------------------------------------------
      
      
      [dum, imax] = max(eastward_propagation_metric .* meets_mjo_criteria);
      if (sum(eastward_propagation_metric  .* meets_mjo_criteria == dum) > 1)
				% Tie-breaker, but only if really needed.
	[dum, imax] = max(eastward_propagation_metric  .* (eastward_propagation_metric == dum) + ...
			  eastward_propagation_metric2 .* (eastward_propagation_metric == dum) );

      end

      
      indx_of_max_east_propagation = lptid_for_this_clump(imax);
      ii = indx_of_max_east_propagation;
      
      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Below determine the overall system stats and whether it is an MJO LPT.
      %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      disp(['--> ', num2str(indx_of_max_east_propagation), ...
	    ' is chosen for clump #', num2str(this_clump_num),'. <--'])
      

      longest_east_prop_zonal_speed = longest_east_prop_zonal_speed_all(imax);
      longest_east_prop_duration = longest_east_prop_duration_all(imax);
      longest_east_prop_duration_in_lat_band = longest_east_prop_duration_in_lat_band_all(imax);

      idx_begin_east_prop = idx_begin_east_prop_all(imax);
      idx_end_east_prop = idx_end_east_prop_all(imax);

      
      GG=G.TIMECLUSTERS(indx_of_max_east_propagation) ;
      GG.date=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=GG.date(end)-GG.date(1) ;
      
      [GG.year,GG.month,GG.day]=datevec(GG.date) ;
      [GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
      [GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;


				% Total longitude propagation
      total_lon_propagation = max(GG.lon)-min(GG.lon);
      net_lon_propagation = GG.lon(end)-GG.lon(1);

      disp(['Total lon span (westmost to eastmost): ', num2str(total_lon_propagation)])
      disp(['Net lon propagation (start to end): ', num2str(net_lon_propagation)])

      
      if (idx_begin_east_prop > -1 & idx_end_east_prop > -1)
	[GG.year0e,GG.month0e,GG.day0e,GG.hour0e]=datevec(GG.time(idx_begin_east_prop)) ;
	[GG.year1e,GG.month1e,GG.day1e,GG.hour1e]=datevec(GG.time(idx_end_east_prop)) ;

      
	disp([num2str(GG.year0),sprintf('%02d',GG.month0),sprintf('%02d',GG.day0),sprintf('%02d',GG.hour0), ...
	      ' to ', num2str(GG.year1),sprintf('%02d',GG.month1),sprintf('%02d',GG.day1),sprintf('%02d',GG.hour1)])

      
	disp(['East propagation from ',...
	      num2str(GG.year0e),sprintf('%02d',GG.month0e),sprintf('%02d',GG.day0e),sprintf('%02d',GG.hour0e), ...
	      ' to ', num2str(GG.year1e),sprintf('%02d',GG.month1e),sprintf('%02d',GG.day1e),sprintf('%02d',GG.hour1e)])

	net_eastward_lon_propagation = GG.lon(idx_end_east_prop) - GG.lon(idx_begin_east_prop);
	disp(['Net eastward lon propagation: ', num2str(net_eastward_lon_propagation)])
	disp(['Net eastward duration (days): ', num2str(longest_east_prop_duration)])
	disp(['Net eastward duration in lat range: ', num2str(longest_east_prop_duration_in_lat_band)])
	disp(['Net eastward zonal speed (m/s): ', num2str(longest_east_prop_zonal_speed)])
      
				% What fraction of time moving east?
      else
	disp('No eastward propagation periods detected.')
      end
      

      count_east_moving=0;
      dist_east_moving=0;
      dist_west_moving=0;

      for jj=2:GG.nentries
        if ( GG.lon(jj) > GG.lon(jj-1) )
          count_east_moving=count_east_moving+1 ;
          dist_east_moving=dist_east_moving+(GG.lon(jj) - GG.lon(jj-1)) ;
        else
          dist_west_moving=dist_west_moving+(GG.lon(jj-1) - GG.lon(jj)) ;
        end
      end
      frac_east_moving=1.0*count_east_moving / (GG.nentries-1);
      east_west_moving_ratio=dist_east_moving / dist_west_moving ;
      
      %% In this version, instead of just using eastward prop time, use entire LPT time.
      if (idx_begin_east_prop > 0)
        [min_west_lon, min_west_lon_idx]=min(GG.lon(1:idx_end_east_prop)) ;
        [max_east_lon, max_east_lon_idx]=max(GG.lon(idx_begin_east_prop:end)) ;
      else
        min_west_lon = NaN;
        max_east_lon = NaN;
      end
      
      min_dist_to_eq = min(abs(GG.lat));
      disp(['Minimum dist to EQ: ', num2str(min_dist_to_eq)])
      if (idx_begin_east_prop > -1 & idx_end_east_prop > -1)
	min_dist_to_eq_east = min(abs(GG.lat(idx_begin_east_prop:idx_end_east_prop)));
	disp(['Minimum dist to EQ during east propagation: ', num2str(min_dist_to_eq_east)])
      end

      if (GG.duration > min_duration-0.001)

        if (  GG.zonal_propagation_speed > min_zonal_speed & ...
             total_lon_propagation > min_total_lon_propagation & ...
             net_eastward_lon_propagation > min_total_eastward_lon_propagation & ...
             min_dist_to_eq_east < max_abs_latitude & ...
             longest_east_prop_zonal_speed > min_eastward_prop_zonal_speed & ...
             longest_east_prop_duration > min_eastward_prop_duration)
	  
	  disp('YES! This is an MJO candidate!')

          fprintf(fid_mjo_lpt_list,FMT,...
                  year1,ii,this_clump_num,GG.duration,GG.zonal_propagation_speed,...
                  GG.year0,GG.month0,GG.day0,GG.hour0,...
                  GG.year1,GG.month1,GG.day1,GG.hour1,...
                  mean(GG.volrain)*GG.duration,...
                  idx_begin_east_prop, idx_end_east_prop, ...
            	  longest_east_prop_zonal_speed, longest_east_prop_duration, ...
                  GG.year0e,GG.month0e,GG.day0e,GG.hour0e,...
                  GG.year1e,GG.month1e,GG.day1e,GG.hour1e,...
		  GG.lon(idx_begin_east_prop),GG.lon(idx_end_east_prop));
	  
          count=count+1;
	  
	  
	else
	  disp('Nope, not an MJO candidate.')
	end
	
      else
	disp('Nope, not an MJO candidate.')
      end
      
    end

end


fclose(fid_mjo_lpt_list);

disp('########## Done! ###########')
disp([num2str(count), ' MJO candidates were found.'])
disp(fn_mjo_lpt_list)
