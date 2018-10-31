function [mask_net_eastward_propagation, spd_raw] = west_east_divide_and_conquer(G, year1, ii, do_plotting)


  divide_and_conquer_longitude_cutoff = 20.0; % If propagation in east or west propagating region is larger than this, don't allow it to get "eaten".

  divide_and_conquer_days_cutoff = 7.0; % If propagation in east or west propagating region is at least this duration (days), don't allow it to get "eaten".


  GG=G.TIMECLUSTERS(ii) ;
  GG.date=GG.time-1.5 ;
  GG.size=sqrt(GG.area) ;
  GG.area=GG.area/1e4 ;
  GG.nentries=numel(GG.date) ;
  GG.duration=GG.date(end)-GG.date(1) ;
  
  [GG.year,GG.month,GG.day]=datevec(GG.date) ;
  [GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
  [GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;
  
  
  
  
  
  
  %% Make speed be a centered difference except forward (backward) diff at beginning (end).
  spd_raw = [];
  spd_raw(2:GG.nentries-1) = (GG.lon(3:end) - GG.lon(1:end-2)) * 110000.0 / (2.0*10800.0);
  spd_raw(1) = (GG.lon(2) - GG.lon(1)) * 110000.0 / 10800.0;
  spd_raw(GG.nentries) = (GG.lon(end) - GG.lon(end-1)) * 110000.0 / 10800.0;
  
  
  %% plot
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (do_plotting)
    close all;
    
    eval('!mkdir -p plots');
    figure('visible','off');
    
    subplot(5,1,1)
    plot(GG.time, spd_raw, 'k','linewidth',0.5);
    title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Zonal Propagation Speed)'])
    datetick('x');
    ylabel('m/s')
    hold on
    
    
    
    subplot(5,1,2)
    plot(GG.time, spd_raw, 'k','linewidth',0.5);
    title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Zonal Propagation Speed)'])
    datetick('x');
    ylabel('m/s')
    hold on
    
    
    set(gca,'ylim',[-10.0, 10.0])
    
  end
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Try the consecutive periods of > 0.
  %% This is an iterative process!
  %% For the first step, eat all cases that are one single
  %% 3h period of eastward (westward) propagation surrounded by
  %% westward (eastward) propagation. Use the median filter for this.
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  where_eastward_propagation = (spd_raw > 0.0);
  mask_net_eastward_propagation = where_eastward_propagation;
  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (do_plotting)
  				%figure('visible','off')
    
    subplot(5,1,3)
    plot(GG.time,where_eastward_propagation,'bo-');
    set(gca,'ylim',[0,1],'ytick',[0,1],'yticklabel',{'West','East'})
    title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Unprocessed)'])
    datetick('x')
    
    
    subplot(5,1,4)
    plot(GG.time,mask_net_eastward_propagation,'bo-');
    set(gca,'ylim',[0,1],'ytick',[0,1],'yticklabel',{'West','East'})
    title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Outlier filter)'])
    datetick('x')
  end
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
	      % Note: from here on, use mask_net_eastward_propagation.
  
  keep_going = 1;
  niter = 0;
  maxiter = 10;
  while(keep_going)
    niter = niter + 1;
    old_mask_net_eastward_propagation = mask_net_eastward_propagation;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Sort the eastward and westward propagation periods.
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
    statsE=regionprops(mask_net_eastward_propagation, 'all') ;
    statsW=regionprops(~mask_net_eastward_propagation, 'all') ;
    
    if ( numel(statsE) < 1 | numel(statsW) < 1)
      break
    end
    
    if (statsE(1).SubarrayIdx{2}(1) == 1) %First section is eastward propagating.
      statsALL = [statsE(1)];
      for kkk = 1:numel(statsW)
	statsALL = [statsALL, statsW(kkk)];
	if kkk + 1 <= numel(statsE)
	  statsALL = [statsALL, statsE(kkk+1)];
	end
      end
    else
      statsALL = [statsW(1)];
      for kkk = 1:numel(statsE)
	statsALL = [statsALL, statsE(kkk)];
	if kkk + 1 <= numel(statsW)
	  statsALL = [statsALL, statsW(kkk+1)];
	end	      
      end	    
    end
    
    areasALL = [statsALL.Area];
    [areas_sort, areas_sort_order] = sort(areasALL);
    
    for jjj = [areas_sort_order] %1:numel(statsW)
      
      indx_this_segment = statsALL(jjj).SubarrayIdx{2};
      lon_this_segment = GG.lon(indx_this_segment);
      
      %% If it is a westward segment
      if mean(mask_net_eastward_propagation(indx_this_segment)) < 0.5
	
	lon_begin_west_progation = max(lon_this_segment);
	lon_end_west_progation = min(lon_this_segment);
	west_lon_propagation = -1 * (lon_end_west_progation - lon_begin_west_progation);
	west_duration = max(GG.time(indx_this_segment)) - min(GG.time(indx_this_segment));
	
      			%Find the adjacent eastward propagating areas.
	jjj_before = jjj - 1;
      	if (jjj_before < 1)
      	  continue
      	end
      	jjj_after = jjj + 1;
      	if (jjj_after > numel(statsALL))
      	  continue
      	end
	
	lon_of_east_propagation_before = GG.lon(statsALL(jjj_before).SubarrayIdx{2});
	lon_begin_east_propagation_before = min(lon_of_east_propagation_before);%(1);
	lon_end_east_propagation_before = max(lon_of_east_propagation_before);%(end);
	east_lon_propagation_before = lon_end_east_propagation_before - lon_begin_east_propagation_before;
	
	lon_of_east_propagation_after = GG.lon(statsALL(jjj_after).SubarrayIdx{2});
	lon_begin_east_propagation_after = min(lon_of_east_propagation_after);%(1);
	lon_end_east_propagation_after = max(lon_of_east_propagation_after);%(end);
	east_lon_propagation_after = lon_end_east_propagation_after - lon_begin_east_propagation_after;
	
	
	%% Check how much westward propagation there was in StatsW.
	
	if (west_lon_propagation < divide_and_conquer_longitude_cutoff & west_duration < divide_and_conquer_days_cutoff)
	  
      	      % If "area" (e.g., duration) of statsW is less than both
      	      % "areas" of statsE, then eat it.
      	  if (statsALL(jjj).Area <= statsALL(jjj_before).Area && ...
      	      statsALL(jjj).Area <= statsALL(jjj_after).Area)
	    
      	    mask_net_eastward_propagation(indx_this_segment) = 1;
      	  end
	  
	  %% Also check propagation. If eastward propagation on each side is greater, eat it.      
      	  if (west_lon_propagation <= east_lon_propagation_before && ...
      	      west_lon_propagation <= east_lon_propagation_after)
	    
      	    mask_net_eastward_propagation(indx_this_segment) = 1;
      	  end
	  
	end
	
      else %% If it is a eastward segment
	
	
	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Let westerly propagation periods eat easterly jogs
	%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	lon_begin_east_progation = min(lon_this_segment) ; %(1);
	lon_end_east_progation = max(lon_this_segment) ; %(end);
	east_lon_propagation = (lon_end_east_progation - lon_begin_east_progation);
	east_duration = max(GG.time(indx_this_segment)) - min(GG.time(indx_this_segment));
	
	
      			%Find the adjacent eastward propagating areas.
	jjj_before = jjj - 1;
      	if (jjj_before < 1)
      	  continue
      	end
      	jjj_after = jjj + 1;
      	if (jjj_after > numel(statsALL))
      	  continue
      	end
	
	
	lon_of_west_propagation_before = GG.lon(statsALL(jjj_before).SubarrayIdx{2});
	lon_begin_west_propagation_before = max(lon_of_west_propagation_before);%(1);
	lon_end_west_propagation_before = min(lon_of_west_propagation_before);%(end);
	west_lon_propagation_before = lon_end_west_propagation_before - lon_begin_west_propagation_before;
	
	lon_of_west_propagation_after = GG.lon(statsALL(jjj_after).SubarrayIdx{2});
	lon_begin_west_propagation_after = max(lon_of_west_propagation_after);%(1);
	lon_end_west_propagation_after = min(lon_of_west_propagation_after);%(end);
	west_lon_propagation_after = lon_end_west_propagation_after - lon_begin_west_propagation_after;
	
	
	%% Check how much eastward propagation there was in StatsE.
	
	if (east_lon_propagation < divide_and_conquer_longitude_cutoff & east_duration < divide_and_conquer_days_cutoff)
	  
      	      % If "area" (e.g., duration) of statsW is less than both
      	      % "areas" of statsE, then eat it.
      	  if (statsALL(jjj).Area <= statsALL(jjj_before).Area && ...
	      statsALL(jjj).Area <= statsALL(jjj_after).Area)
	    mask_net_eastward_propagation(indx_this_segment) = 0;
      	  end
	  
	  %% Also check propagation. If eastward propagation on each side is greater, eat it.      
      	  if (east_lon_propagation <= west_lon_propagation_before && ...
      	      east_lon_propagation <= west_lon_propagation_after)
	    
      	    mask_net_eastward_propagation(indx_this_segment) = 0;
      	  end
	  
	end
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
  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (do_plotting)
    subplot(5,1,5)
    plot(GG.time,mask_net_eastward_propagation,'bo-');
    set(gca,'ylim',[0,1],'ytick',[0,1],'yticklabel',{'West','East'})
    title([num2str(year1), ' - ', num2str(year1+1), ' LPT #', num2str(ii), ' (Final)'])
    datetick('x')
    
    disp(['plots/eastward_westward_separation_', num2str(year1), '_', sprintf('%03d',ii),'.png'])
    saveas(gcf, ['plots/eastward_westward_separation_', num2str(year1), '_', sprintf('%03d',ii),'.png'])
  end
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
