clear all
close all

%% This script calculates rainfall contribution from: LPTs, MJO, and Rossby systems.
%% It saves a .mat file.

addpath('../../config')
options

PROCESSED_DATA_DIR = ['../../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

MASK_DIR = ['../../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/masks'];

PLOT_DIR = ['../../plots/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems/individual_clumps'];


corner_label={'5 deg. Filter','Threshold=12 mm/day'};
clumps=dlmread(['../../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'],'',1,0);
CLUMPS=clumps;
colors=hsv(12);


MJO=dlmread(['../../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'],'',1,0);

ROSSBY=dlmread(['../../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/rossby_lpt_list.rejoin2.txt'],'',1,0);

lon1 = 0.0;
lon2 = 360.0;
lat1 = -45.0;
lat2 = 45.0;

for year1=[1998]
%for year1=[1998:2017]

  %%%%%%%%%%%%%%%%%%%%%%%%

  year2=year1+1 ;

  yyyy1=num2str(year1) ;
  yyyy2=num2str(year2) ;

  y1_y2=[yyyy1,'_',yyyy2] ;
  if year1 == 2018
    y11_y22=[yyyy1,'060100_',yyyy1,'112721'] ;
  else
    y11_y22=[yyyy1,'060100_',yyyy2,'063021'] ;
  end
  disp(y1_y2) ;

  %% Use the "grand master mask" from the annual files to calculate the seasonal rain.

  fileInBase=['rain_map2_full_year_all_lpts__',y1_y2];
  f=load([PLOT_DIR,'/',fileInBase,'.mat']);

  grand_master_times = f.mask_time;
  [Y, M] = datevec(grand_master_times);
  keep_time = (M == 6 | M == 7 | M == 8);
  grand_master_times_keep = grand_master_times(keep_time);
  eliminate_time = ~keep_time;
  grand_master_mask = f.mask_mask_lpt(keep_time,:,:);
  grand_master_mask_mjo = f.mask_mask_mjo(keep_time,:,:);
  grand_master_mask_mjo_eastward_portion = f.mask_mask_mjo_eastward_portion(keep_time,:,:);
  grand_master_mask_rossby = f.mask_mask_rossby(keep_time,:,:);
  grand_master_mask_rossby_westward_portion = f.mask_mask_rossby_westward_portion(keep_time,:,:);

  %% Set the times outside the season of interest to ZERO.
  %grand_master_mask(eliminate_time,:,:) = 0;
  %grand_master_mask_mjo(eliminate_time,:,:) = 0;
  %grand_master_mask_mjo_eastward_portion(eliminate_time,:,:) = 0;
  %grand_master_mask_rossby(eliminate_time,:,:) = 0;
  %grand_master_mask_rossby_westward_portion(eliminate_time,:,:) = 0;
  
  this_mask.lon = f.lon;
  this_mask.lat = f.lat;
  
  %% Get coordinates of the lon/lat extent of the "grand master mask."
    
  max_mask = squeeze(max(grand_master_mask));
  max_mask_x = squeeze(max(max_mask));
  max_mask_y = squeeze(max(max_mask'));
  
  max_mask_mjo = squeeze(max(grand_master_mask_mjo));
  max_mask_mjo_eastward_portion = squeeze(max(grand_master_mask_mjo_eastward_portion));
  
  max_mask_rossby = squeeze(max(grand_master_mask_rossby));
  max_mask_rossby_westward_portion = squeeze(max(grand_master_mask_rossby_westward_portion));
  
  %% Use the mask to get the rain from this LPT system family.
  %count = 0; % Will divide to get mean.
  idx = 1;
  
    for dn_rain = grand_master_times_keep
      
      [y_rain, m_rain, d_rain, h_rain] = datevec(dn_rain);
      fileRain = ['../../data/trmm/interim/gridded_rain_rates/gridded_rain_rates_',...
		  num2str(y_rain), sprintf('%02d', m_rain), sprintf('%02d', d_rain), sprintf('%02d',h_rain),'.nc'];
				%disp(fileRain);
      F.precip=ncread(fileRain,'rain')' ;
      F.precip(~isfinite(F.precip)) = 0.0;
      F.precip(F.precip < 0.0) = 0.0;
      F.precip_mjo = F.precip;
      F.precip_mjo_eastward_portion = F.precip;
      F.precip_rossby = F.precip;
      F.precip_rossby_westward_portion = F.precip;
      
      F.lon=ncread(fileRain,'lon');
      F.lat=ncread(fileRain,'lat');

      this_time_eliminate_mask = (squeeze(grand_master_mask(idx,:,:)) < 1);
      F.precip(this_time_eliminate_mask) = 0.0;
      
      this_time_eliminate_mask_mjo = (squeeze(grand_master_mask_mjo(idx,:,:)) < 1);
      F.precip_mjo(this_time_eliminate_mask_mjo) = 0.0;
      
      this_time_eliminate_mask_mjo_eastward_portion = (squeeze(grand_master_mask_mjo_eastward_portion(idx,:,:)) < 1);
      F.precip_mjo_eastward_portion(this_time_eliminate_mask_mjo_eastward_portion) = 0.0;
      
      this_time_eliminate_mask_rossby = (squeeze(grand_master_mask_rossby(idx,:,:)) < 1);
      F.precip_rossby(this_time_eliminate_mask_rossby) = 0.0;
      
      this_time_eliminate_mask_rossby_westward_portion = (squeeze(grand_master_mask_rossby_westward_portion(idx,:,:)) < 1);
      F.precip_rossby_westward_portion(this_time_eliminate_mask_rossby_westward_portion) = 0.0;
      
      
      
      if idx == 1
	rain_sum = 3.0 * F.precip;
	rain_sum_mjo = 3.0 * F.precip_mjo;
	rain_sum_mjo_eastward_portion = 3.0 * F.precip_mjo_eastward_portion;
	rain_sum_rossby = 3.0 * F.precip_rossby;
	rain_sum_rossby_westward_portion = 3.0 * F.precip_rossby_westward_portion;
      else
	rain_sum = rain_sum + 3.0 * F.precip;
	rain_sum_mjo = rain_sum_mjo + 3.0 * F.precip_mjo;
	rain_sum_mjo_eastward_portion = rain_sum_mjo_eastward_portion + 3.0 * F.precip_mjo_eastward_portion;
	rain_sum_rossby = rain_sum_rossby + 3.0 * F.precip_rossby;
	rain_sum_rossby_westward_portion = rain_sum_rossby_westward_portion + 3.0 * F.precip_rossby_westward_portion;
      end
      
      idx = idx + 1;
      
    end
    
    rain_sum(rain_sum < 0.01) = NaN ;
    rain_sum_mjo(rain_sum_mjo < 0.01) = NaN ;
    rain_sum_mjo_eastward_portion(rain_sum_mjo_eastward_portion < 0.01) = NaN ;
    rain_sum_rossby(rain_sum_rossby < 0.01) = NaN ;
    rain_sum_rossby_westward_portion(rain_sum_rossby_westward_portion < 0.01) = NaN ;
    
    %n_days = grand_master_end_time - grand_master_start_time + DT/24.0;
    n_days = 3.0 * sum(keep_time) / 24.0;
    DATA=rain_sum/n_days ; %log(rain_sum) ;
    DATA_mjo=rain_sum_mjo/n_days ; %log(rain_sum) ;
    DATA_mjo_eastward_portion=rain_sum_mjo_eastward_portion/n_days ; %log(rain_sum) ;
    DATA_rossby=rain_sum_rossby/n_days ; %log(rain_sum) ;
    DATA_rossby_westward_portion=rain_sum_rossby_westward_portion/n_days ; %log(rain_sum) ;
    
    
    %%
    %% Begin plot
    %%
    close all
    
    figure('visible','off')
    set(gcf, 'position',[100,100,2500,1600])

    %%
    %% All LPTs
    %%
    subplot(321)


    plot_map_background([lon1, lon2, lat1, lat2]);
    hold on

    
    pcolor(F.lon, F.lat, DATA) ;
    hold on

    shading flat
    caxis([0,10.0])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    %colormap(flipud(gray()))
    colormap(cmap)


    %% Outline of spatial influence of the LPT family.
    contour(this_mask.lon, this_mask.lat, max_mask, [0.5,0.5],'r-', 'linewidth',2);


    hcb=colorbar('NorthOutside') ;

    set(hcb,'position',[0.15,0.95,0.8,0.008]) ;
    
        
    plot_map_background([lon1, lon2, lat1, lat2]);
    daspect([1,1,1]);
    
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['LPT Rainfall: ',yyyy1, ' - ', yyyy2])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')

    %%
    %% MJO LPTs
    %%
    subplot(323)


    plot_map_background([lon1, lon2, lat1, lat2]);
    hold on

    
    pcolor(F.lon, F.lat, DATA_mjo) ;
    hold on

    shading flat
    caxis([0,10.0])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    %colormap(flipud(gray()))
    colormap(cmap)


    %% Outline of spatial influence of the LPT family.
    contour(this_mask.lon, this_mask.lat, max_mask_mjo, [0.5,0.5],'r-', 'linewidth',2);
        
    plot_map_background([lon1, lon2, lat1, lat2]);
    daspect([1,1,1]);
    
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['MJO LPT Rainfall: ',yyyy1, ' - ', yyyy2])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')

    %%
    %% MJO LPTs eastward portion
    %%
    subplot(325)


    plot_map_background([lon1, lon2, lat1, lat2]);
    hold on

    
    pcolor(F.lon, F.lat, DATA_mjo_eastward_portion) ;
    hold on

    shading flat
    caxis([0,10.0])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    %colormap(flipud(gray()))
    colormap(cmap)


    %% Outline of spatial influence of the LPT family.
    contour(this_mask.lon, this_mask.lat, max_mask_mjo_eastward_portion, [0.5,0.5],'r-', 'linewidth',2);
        
    plot_map_background([lon1, lon2, lat1, lat2]);
    daspect([1,1,1]);
    
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['MJO LPT Eastward Portion Rainfall: ',yyyy1, ' - ', yyyy2])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')






    %%
    %% W-Prop LPTs
    %%
    subplot(324)


    plot_map_background([lon1, lon2, lat1, lat2]);
    hold on

    
    pcolor(F.lon, F.lat, DATA_rossby) ;
    hold on

    shading flat
    caxis([0,10.0])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    %colormap(flipud(gray()))
    colormap(cmap)


    %% Outline of spatial influence of the LPT family.
    contour(this_mask.lon, this_mask.lat, max_mask_rossby, [0.5,0.5],'r-', 'linewidth',2);
        
    plot_map_background([lon1, lon2, lat1, lat2]);
    daspect([1,1,1]);
    
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['Westward LPTs Rainfall: ',yyyy1, ' - ', yyyy2])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')

    %%
    %% W-Prop LPTs westward portion
    %%
    subplot(326)


    plot_map_background([lon1, lon2, lat1, lat2]);
    hold on

    
    pcolor(F.lon, F.lat, DATA_rossby_westward_portion) ;
    hold on

    shading flat
    caxis([0,10.0])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    %colormap(flipud(gray()))
    colormap(cmap)


    %% Outline of spatial influence of the LPT family.
    contour(this_mask.lon, this_mask.lat, max_mask_rossby_westward_portion, [0.5,0.5],'r-', 'linewidth',2);
        
    plot_map_background([lon1, lon2, lat1, lat2]);
    daspect([1,1,1]);
    
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['Westward LPTs Westward Portion Rainfall: ',yyyy1, ' - ', yyyy2])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')









    
    fileOutBase=['rain_map2_full_year_all_lpts__',y1_y2,'__jja'];

    
    eval(['!mkdir -p ',PLOT_DIR])
    disp([PLOT_DIR,'/',fileOutBase,'.png'])
    saveas(gcf,[PLOT_DIR,'/',fileOutBase,'.png'])



    %%
    %% .mat file output.
    %%

    
    fout=[];
    fout.lon = F.lon;
    fout.lat = F.lat;
    fout.rain_lpt = DATA;
    fout.rain_mjo = DATA_mjo;
    fout.rain_mjo_eastward_portion = DATA_mjo_eastward_portion;
    fout.rain_rossby = DATA_rossby;
    fout.rain_rossby_westward_portion = DATA_rossby_westward_portion;
    fout.info = "Rainfall from 3B42 V7 in mm/day averaged over 1 year.";
    fout.n_days = n_days;
    fout.mask_time = grand_master_times;
    fout.mask_lon = this_mask.lon;
    fout.mask_lat = this_mask.lat;
    fout.mask_mask_lpt = logical(grand_master_mask > 0.5);
    fout.mask_mask_mjo = logical(grand_master_mask_mjo > 0.5);
    fout.mask_mask_mjo_eastward_portion = logical(grand_master_mask_mjo_eastward_portion > 0.5);
    fout.mask_mask_rossby = logical(grand_master_mask_rossby > 0.5);
    fout.mask_mask_rossby_westward_portion = logical(grand_master_mask_rossby_westward_portion > 0.5);
    
    disp([PLOT_DIR,'/',fileOutBase,'.mat'])
    eval(['save ',PLOT_DIR,'/',fileOutBase,'.mat -struct fout'])
    
end
