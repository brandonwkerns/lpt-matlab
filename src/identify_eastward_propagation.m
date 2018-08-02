clear all
close all

%% NINO 3.4 stuff.
nino34=load('nino34.mat');


% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

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

min_zonal_speed = -999.0 ; % in m/s.
min_duration    = 7.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_eastward_prop_zonal_speed    = 0.0 ; % in m/s.
min_eastward_prop_duration    = 5.0 ; % in Days. Doesn't include 3-Day accumulation period.

min_net_lon_propagation   = -999.0 ;%20.0 ; % in deg. longitude.
min_total_lon_propagation = -999.0 ;%20.0 ; % in deg. longitude.

mc_lon_1 = 100.0 ; % West end of MC for MC crossing
mc_lon_2 = 130.0 ; % East end of MC for MC crossing

% Search area for initial time cluster.
search_area = [50.0, 180.0, -15.0, 15.0];

%%
%%
%%

% For output table files.
FMT=['%10d%10d%10.2f%16.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d %10.2f%15.2f ' ...
     '%10.1f%20d%15d%11.2f%11.2f\n'];

header='      year     index  duration  mean_zonal_spd       begin         end    volrain  %_tot_volrain    nino3.4     eprop_begin_idx  eprop_end_idx  eprop_spd  eprop_dur';

fid_mc_crossing_lpts=fopen([EASTWARD_PROP_DATA_DIR,'/list_mc_crossing_io_lpts.txt'],'w');
fid_non_mc_crossing_lpts=fopen([EASTWARD_PROP_DATA_DIR,'/list_non_mc_crossing_io_lpts.txt'],'w');
fid_wpac_lpts=fopen([EASTWARD_PROP_DATA_DIR,'/list_wpac_lpts.txt'],'w');
fid_non_east_propagating_lpts=fopen([EASTWARD_PROP_DATA_DIR,'/list_non_east_propagating_lpts.txt'],'w');

fprintf(fid_mc_crossing_lpts, '%s\n', header);
fprintf(fid_non_mc_crossing_lpts, '%s\n', header);
fprintf(fid_wpac_lpts, '%s\n', header);
fprintf(fid_non_east_propagating_lpts, '%s\n', header);


for year1=2004:2017  ;
%for year1=[2014]  ;

    year2=year1+1 ;
    
    yyyy1=num2str(year1) ;
    yyyy2=num2str(year2) ;
    
    y1_y2=[yyyy1,'_',yyyy2] ;
    
    disp(y1_y2) ;

    
    
    
    %VOLRAIN=load(['../../../total_volrain_full_year/',...
    %              'monthly_volrain_',y1_y2,'.mat']);
    
    %TOTALVOLRAIN=sum(VOLRAIN.volrain) ;
    TOTALVOLRAIN=100000.0 ; %HACK! This needs updated.

    
%    F=load(['/home/disk/manta8/bkerns/pubs/lpt/timelon/rain_with_rmm/rain_hov_',y1_y2,'_15deg_3day.mat']) ;

    dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.mat']);
    G=load([PROCESSED_DATA_DIR,'/', dir0(1).name]) ;
    
    
    count=0 ; % Keep count of east propagating systems.
    
    
    for ii=1:numel(G.TIMECLUSTERS)

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

      if (GG.lon(1) < search_area(1) | GG.lon(1) > search_area(2) | ...
	  GG.lat(1) < search_area(3) | GG.lat(1) > search_area(4))

	disp('Out of search area. Skipping.')
	continue
      end

      
      %%
      %% Daily least squares propagation speed.
      %% I was hoping this would help get rid of some of the "wiggles"
      %% but it did not! Even using 3 day periods, and I don't want to use longer periods.
      %% So now I'm using duration of periods of east propagation.
      %% and "eating" the short duration periods with westward motion.
      %%
      
      %% Try median filter on the 3h centroid propagation speeds.
      avg_interval = 3; % Number of time points.
      avg_interval_half = (avg_interval - 1) / 2;

      spd_raw = [];
      spd_raw(2:GG.nentries) = (GG.lon(2:end) - GG.lon(1:end-1)) * 110000.0 / 10800.0;
      spd_raw(1) = spd_raw(2);							  
      
      spd_median_filter = spd_raw ;
      for jj = avg_interval_half+1:GG.nentries- avg_interval_half
	x_fit = GG.time(jj-avg_interval_half:jj+avg_interval_half);
	y_fit = spd_raw(jj-avg_interval_half:jj+avg_interval_half);
	spd_median_filter(jj) = median(y_fit);
      end
      
				% plot
      %{
      figure('visible','off')
      subplot(4,1,1)
      plot(GG.time, spd_raw, 'k','linewidth',0.5);
      title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Zonal Propagation Speed)'])      
      datetick('x');
      ylabel('m/s')
      hold on
      plot(GG.time, spd_median_filter, 'r','linewidth',1.0);
      plot(GG.time, 0.0*spd_median_filter, 'k--');
      plot(GG.time(spd_median_filter > 0), spd_median_filter(spd_median_filter > 0), 'ro','markersize',3,'markerfacecolor','r');
      
      saveas(gcf, ['plots/median_filter_speed_', num2str(year1), '_', sprintf('%02d',ii), '.png'])
      
      set(gca,'ylim',[-10.0, 10.0])
      
      saveas(gcf, ['plots/median_filter_speed_zoomed_', num2str(year1), '_', sprintf('%02d',ii), '.png'])
      %}
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Try the consecutive periods of > 0.
      %% This is an iterative process!
      %% For the first step, eat all cases that are one single
      %% 3h period of eastward (westward) propagation surrounded by
      %% westward (eastward) propagation. Use the median filter for this.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      where_eastward_propagation = (spd_raw > 0.0);
      where_eastward_propagation2 = where_eastward_propagation;
      for jj = 2:numel(where_eastward_propagation)-1
	where_eastward_propagation2(jj) = median(where_eastward_propagation2(jj-1:jj+1));
      end
      where_eastward_propagation2(1) = where_eastward_propagation2(2);
      where_eastward_propagation2(end) = where_eastward_propagation2(end-1);
      
      
      mask_net_eastward_propagation = where_eastward_propagation2;
      
				%figure('visible','off')
      %{
      subplot(4,1,2)
      plot(GG.time,where_eastward_propagation,'bo-');
      set(gca,'ylim',[0,1],'ytick',[0,1],'yticklabel',{'West','East'})
      title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Unprocessed)'])
      datetick('x')

      
      subplot(4,1,3)
      plot(GG.time,mask_net_eastward_propagation,'bo-');
      set(gca,'ylim',[0,1],'ytick',[0,1],'yticklabel',{'West','East'})
      title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Outlier filter)'])
      datetick('x')
      %}
      keep_going = 1;
      niter = 0;
      maxiter = 10;
      while(keep_going)
	niter = niter + 1;
	old_mask_net_eastward_propagation = mask_net_eastward_propagation;
	
		  % Let easterly propagation periods eat westerly jogs
        statsE=regionprops(mask_net_eastward_propagation, 'all') ;
        statsW=regionprops(~mask_net_eastward_propagation, 'all') ;
	
	for jjj = 1:numel(statsW)
	  
	  indxW = statsW(jjj).SubarrayIdx{2};
	  
			%Find the adjacent eastward propagating areas.
	  indx_before = min(indxW) - 1;
	  if (indx_before < 1)
	    continue
	  end
	  for jjj_before = 1:numel(statsE)
	    if (sum(statsE(jjj_before).SubarrayIdx{2} == indx_before) > 0)
	      break
	    end	 
	  end
	  jjj_after = jjj_before + 1;
	  if (jjj_after > numel(statsE))
	    continue
	  end
	  
	      % If "area" (e.g., duration) of statsW is less than both
	      % "areas" of statsE, then eat it.
	  if (statsW(jjj).Area <= statsE(jjj_before).Area && ...
	      statsW(jjj).Area <= statsE(jjj_after).Area)
	    mask_net_eastward_propagation(indxW) = 1;
	  end
	end
	
	
	      % Let westerly propagation periods eat easterly jogs	  
        statsE=regionprops(mask_net_eastward_propagation, 'all') ;
        statsW=regionprops(~mask_net_eastward_propagation, 'all') ;
	
	for jjj = 1:numel(statsE)
	  
	  indxE = statsE(jjj).SubarrayIdx{2};
	  
			%Find the adjacent eastward propagating areas.
	  indx_before = min(indxE) - 1;
	  if (indx_before < 1)
	    continue
	  end
	  for jjj_before = 1:numel(statsW)
	    if (sum(statsW(jjj_before).SubarrayIdx{2} == indx_before) > 0)
	      break
	    end	 
	  end
	  jjj_after = jjj_before + 1;
	  if (jjj_after > numel(statsW))
	    continue
	  end
	  
	      % If "area" (e.g., duration) of statsW is less than both
	      % "areas" of statsE, then eat it.
	  if (statsE(jjj).Area <= statsW(jjj_before).Area && ...
	      statsE(jjj).Area <= statsW(jjj_after).Area)
	    mask_net_eastward_propagation(indxE) = 0;
	  end
	end
	
	check = (mask_net_eastward_propagation ~= old_mask_net_eastward_propagation);
	if (sum(check) == 0)
	  keep_going = 0;
	  disp(['Finished in ',num2str(niter),' iterations.'])
	end
	
	if (niter > maxiter)
	  disp(['Warning: exceeded ',num2str(maxiter),' iterations! Stopping.'])
	  break
	end
	
	
      end
      
      %subplot(4,1,4)
      %plot(GG.time,mask_net_eastward_propagation,'bo-');
      %set(gca,'ylim',[0,1],'ytick',[0,1],'yticklabel',{'West','East'})
      %title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Final)'])
      %datetick('x')
      
      %saveas(gcf, ['plots/eastward_westward_separation_', num2str(year1), '_', sprintf('%02d',ii),'.png'])


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
	time_east_prop = GG.time(idx_begin_east_prop:idx_end_east_prop);

	%Linear Regression
	[FIT0,S,MU]=polyfit(time_east_prop,lon_east_prop,1) ;
        longest_east_prop_zonal_speed = ...
            (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2));
	longest_east_prop_duration = max(time_east_prop) - min(time_east_prop);
	
      else
	longest_east_prop_zonal_speed = -999.0;
	longest_east_prop_duration = -999.0;
	idx_begin_east_prop = -1;
	idx_end_east_prop = -1;
      end


      
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Below determine the overall system stats and whether it is an MJO LPT.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Total longitude propagation
        total_lon_propagation=max(GG.lon)-min(GG.lon);
        net_lon_propagation=GG.lon(end)-GG.lon(1);

        % What fraction of time moving east?
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
	%if (idx_begin_east_prop > 0)
        %  [max_east_lon, max_east_lon_idx]=max(GG.lon(idx_begin_east_prop:idx_end_east_prop)) ;
        %  [min_west_lon, min_west_lon_idx]=min(GG.lon(idx_begin_east_prop:idx_end_east_prop)) ;
	%else
          [max_east_lon, max_east_lon_idx]=max(GG.lon) ;
          [min_west_lon, min_west_lon_idx]=min(GG.lon) ;
	%end
	
        
        %% ENSO
        [y,m,d,h]=datevec(GG.time(1));

        findThisMonth = find( ...
            nino34.three_months.time == datenum(y,m,1,0,0,0) ) ;

        SSTA=nino34.three_months.ssta(findThisMonth) ;
	if (numel(SSTA) == 0)
	  SSTA = NaN;
	end
        
        if (GG.duration > min_duration-0.001)

          if ( (max_east_lon < mc_lon_2 ) & ...
                 GG.zonal_propagation_speed > min_zonal_speed & ...
                 net_lon_propagation > min_total_lon_propagation & ...
		 longest_east_prop_zonal_speed > min_eastward_prop_zonal_speed & ...
		 longest_east_prop_duration > min_eastward_prop_duration)
                
                fprintf(fid_non_mc_crossing_lpts,FMT,...
                        year1,ii,GG.duration,GG.zonal_propagation_speed,...
                        GG.year0,GG.month0,GG.day0,GG.hour0,...
                        GG.year1,GG.month1,GG.day1,GG.hour1,...
                        mean(GG.volrain)*GG.duration,...
                        100*mean(GG.volrain)*GG.duration/TOTALVOLRAIN,...
                        SSTA, ...
			idx_begin_east_prop, idx_end_east_prop, ...
		        longest_east_prop_zonal_speed, longest_east_prop_duration);

                count=count+1;

          elseif ( ( min_west_lon > mc_lon_1) & ...
                 GG.zonal_propagation_speed > min_zonal_speed & ...
                 net_lon_propagation > min_total_lon_propagation & ...
		 longest_east_prop_zonal_speed > min_eastward_prop_zonal_speed & ...
		 longest_east_prop_duration > min_eastward_prop_duration)
                
                fprintf(fid_wpac_lpts,FMT,...
                        year1,ii,GG.duration,GG.zonal_propagation_speed,...
                        GG.year0,GG.month0,GG.day0,GG.hour0,...
                        GG.year1,GG.month1,GG.day1,GG.hour1,...
                        mean(GG.volrain)*GG.duration,...
                        100*mean(GG.volrain)*GG.duration/TOTALVOLRAIN,...
                        SSTA, ...
			idx_begin_east_prop, idx_end_east_prop, ...
		        longest_east_prop_zonal_speed, longest_east_prop_duration);

                count=count+1;

		
          elseif ( min_west_lon < mc_lon_1 & max_east_lon > mc_lon_2 & ...
                     GG.zonal_propagation_speed > min_zonal_speed & ...
                     net_lon_propagation > min_total_lon_propagation & ...
		     longest_east_prop_zonal_speed > min_eastward_prop_zonal_speed & ...
		     longest_east_prop_duration > min_eastward_prop_duration)
            
                fprintf(fid_mc_crossing_lpts,FMT,...
                        year1,ii,GG.duration,GG.zonal_propagation_speed,...
                        GG.year0,GG.month0,GG.day0,GG.hour0,...
                        GG.year1,GG.month1,GG.day1,GG.hour1,...
                        mean(GG.volrain)*GG.duration,...
                        100*mean(GG.volrain)*GG.duration/TOTALVOLRAIN,...
                        SSTA, ...
			idx_begin_east_prop, idx_end_east_prop, ...
		        longest_east_prop_zonal_speed, longest_east_prop_duration);
                        
                 count=count+1 ;
                
            else
                
                fprintf(fid_non_east_propagating_lpts,FMT,...
                        year1,ii,GG.duration,GG.zonal_propagation_speed,...
                        GG.year0,GG.month0,GG.day0,GG.hour0,...
                        GG.year1,GG.month1,GG.day1,GG.hour1,...
                        mean(GG.volrain)*GG.duration,...
                        100*mean(GG.volrain)*GG.duration/TOTALVOLRAIN,...
                        SSTA, ...
			idx_begin_east_prop, idx_end_east_prop, ...
		        longest_east_prop_zonal_speed, longest_east_prop_duration);
                
            end
                
        
        end        
    
    end
    
end


fclose(fid_mc_crossing_lpts);
fclose(fid_non_mc_crossing_lpts);
fclose(fid_wpac_lpts);
fclose(fid_non_east_propagating_lpts);






    
    
