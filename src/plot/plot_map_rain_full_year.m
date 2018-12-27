clear all
close all

addpath('../config')
options

PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters']

MASK_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/masks']

PLOT_DIR = ['../plots/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems/individual_clumps']


corner_label={'5 deg. Filter','Threshold=12 mm/day'};
clumps=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'],'',1,0);
CLUMPS=clumps;
colors=hsv(12);


MJO=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'],'',1,0);

%for year1=[2017]
for year1=[1998:2017]

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

  G=load([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',y11_y22,'.rejoin2.mat']) ;

  for iiii = 2:20
    if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
      eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
    end
  end

  %% Start with "grand master mask"
  grand_master_start_time = datenum(year1,6,1,0,0,0);
  grand_master_end_time = datenum(year2,5,31,21,0,0);
  grand_master_times = grand_master_start_time:0.125:grand_master_end_time;
  grand_master_mask = zeros(numel(grand_master_times), 400, 1440);
  grand_master_mask_mjo = zeros(numel(grand_master_times), 400, 1440);
  grand_master_mask_mjo_eastward_portion = zeros(numel(grand_master_times), 400, 1440);
    
  %% Get "clumps of worms" for this year.
  clump_idx_this_year = find(clumps(:,1) == year1);
  lptid_this_year = clumps(clump_idx_this_year, 2)';
  clump_num_this_year = clumps(clump_idx_this_year, 3)';
  
  for this_clump_num = [unique(clump_num_this_year)]
    
    disp(['----------- Clump #', num2str(this_clump_num), ' -----------'])

    clf
           

    lptid_for_this_clump = lptid_this_year(clump_num_this_year == this_clump_num);

    lon1 = 0.0;
    lon2 = 360.0;
    lat1 = -45.0;
    lat2 = 45.0;
    dn1 = datenum(2100,1,1,0,0,0);
    dn2 = datenum(1900,1,1,0,0,0);



    %% Rain
    
    %% Using mask.
    
    for ii=[lptid_for_this_clump]    
      this_mask.file = [MASK_DIR,'/lpt_system_mask_',yyyy1,'_',sprintf('%03d',ii),'.nc'] ;
      this_mask.lon = ncread(this_mask.file, 'lon');
      this_mask.lat = ncread(this_mask.file, 'lat');

      this_mask.time = ncread(this_mask.file, 'time');
      this_mask.dn = datenum(1970,1,1,0,0,0) + this_mask.time / 24.0 ;
      this_mask.mask = ncread(this_mask.file, 'mask_with_filter_and_accumulation');
      this_mask.mask = permute(this_mask.mask, [3,2,1]);

      this_beginning_time_indx = find(grand_master_times == max(grand_master_start_time,this_mask.dn(1)));
      this_end_time_indx = find(grand_master_times == min(grand_master_end_time,this_mask.dn(end)));
      this_time_indx_range = this_beginning_time_indx:this_end_time_indx;

      this_beginning_time_indx0 = find(this_mask.dn == max(grand_master_start_time,this_mask.dn(1)));
      this_end_time_indx0 = find(this_mask.dn == min(grand_master_end_time,this_mask.dn(end)));
      this_time_indx_range0 = this_beginning_time_indx0:this_end_time_indx0;

      
      grand_master_mask(this_time_indx_range,:,:) = max(double(grand_master_mask(this_time_indx_range,:,:)),double(this_mask.mask(this_time_indx_range0,:,:)));


      if ( sum(MJO(:,1) == year1 & ...
               MJO(:,2) == ii) > 0 )

    
	grand_master_mask_mjo(this_time_indx_range,:,:) = max(double(grand_master_mask_mjo(this_time_indx_range,:,:)),double(this_mask.mask(this_time_indx_range0,:,:)));%this_mask.mask(this_time_indx_range0,:,:);

	%% For eastward propagating portion, need to use the indices.
	idx11 = MJO((MJO(:,1) == year1 & ...
                     MJO(:,2) == ii),9);
	idx22 = MJO((MJO(:,1) == year1 & ...
                     MJO(:,2) == ii),10);

	for kkkk = 1:numel(idx11)

	  idx1 = idx11(kkkk);
	  idx2 = idx22(kkkk);
	  
	  this_beginning_time_indx = find(grand_master_times == max(grand_master_start_time,this_mask.dn(idx1)));
	  this_end_time_indx = find(grand_master_times == min(grand_master_end_time,this_mask.dn(idx2)));
	  this_time_indx_range = this_beginning_time_indx:this_end_time_indx;

	  this_beginning_time_indx0 = find(this_mask.dn == max(grand_master_start_time,this_mask.dn(idx1)));
	  this_end_time_indx0 = find(this_mask.dn == min(grand_master_end_time,this_mask.dn(idx2)));
	  this_time_indx_range0 = this_beginning_time_indx0:this_end_time_indx0;

	
	  grand_master_mask_mjo_eastward_portion(this_time_indx_range,:,:) = max(double(grand_master_mask_mjo_eastward_portion(this_time_indx_range,:,:)),double(this_mask.mask(this_time_indx_range0,:,:)));%this_mask.mask(this_time_indx_range0,:,:);

	end
	
      end % End eastward propagation portion.
      
    end %End check whether it is an MJO 

    disp([num2str(sum(sum(sum(grand_master_mask > 0.5)))), ', ',...
	  num2str(sum(sum(sum(grand_master_mask_mjo > 0.5)))),', ',...
	  num2str(sum(sum(sum(grand_master_mask_mjo_eastward_portion > 0.5))))])
    if (sum(sum(sum(grand_master_mask_mjo > 0.5))) > sum(sum(sum(grand_master_mask > 0.5))))
      disp('Warning! More MJO points than there are LPT points!!!')
    end
    
  end % End loop over lptid.
    
   %% Get coordinates of the lon/lat extent of the "grand master mask."
    
    max_mask = squeeze(max(grand_master_mask));
    max_mask_x = squeeze(max(max_mask));
    max_mask_y = squeeze(max(max_mask'));


    max_mask_mjo = squeeze(max(grand_master_mask_mjo));
    max_mask_mjo_eastward_portion = squeeze(max(grand_master_mask_mjo_eastward_portion));

    
    %% Use the mask to get the rain from this LPT system family.
    %count = 0; % Will divide to get mean.
    idx = 1;
    
    for dn_rain = grand_master_times

      [y_rain, m_rain, d_rain, h_rain] = datevec(dn_rain);
      fileRain = ['../data/trmm/interim/gridded_rain_rates/gridded_rain_rates_',...
		  num2str(y_rain), sprintf('%02d', m_rain), sprintf('%02d', d_rain), sprintf('%02d',h_rain),'.nc'];
      %disp(fileRain);
      F.precip=ncread(fileRain,'rain')' ;
      F.precip(~isfinite(F.precip)) = 0.0;
      F.precip(F.precip < 0.0) = 0.0;
      F.precip_mjo = F.precip;
      F.precip_mjo_eastward_portion = F.precip;
      
      F.lon=ncread(fileRain,'lon');
      F.lat=ncread(fileRain,'lat');

      this_time_eliminate_mask = (squeeze(grand_master_mask(idx,:,:)) < 1);
      F.precip(this_time_eliminate_mask) = 0.0;

      this_time_eliminate_mask_mjo = (squeeze(grand_master_mask_mjo(idx,:,:)) < 1);
      F.precip_mjo(this_time_eliminate_mask_mjo) = 0.0;

      
      this_time_eliminate_mask_mjo_eastward_portion = (squeeze(grand_master_mask_mjo_eastward_portion(idx,:,:)) < 1);
      F.precip_mjo_eastward_portion(this_time_eliminate_mask_mjo_eastward_portion) = 0.0;
      
      if idx == 1
	rain_sum = 3.0 * F.precip;
	rain_sum_mjo = 3.0 * F.precip_mjo;
	rain_sum_mjo_eastward_portion = 3.0 * F.precip_mjo_eastward_portion;
      else
	rain_sum = rain_sum + 3.0 * F.precip;
	rain_sum_mjo = rain_sum_mjo + 3.0 * F.precip_mjo;
	rain_sum_mjo_eastward_portion = rain_sum_mjo_eastward_portion + 3.0 * F.precip_mjo_eastward_portion;
      end
      
      idx = idx + 1;
      
    end
    
    rain_sum(rain_sum < 0.01) = NaN ;
    rain_sum_mjo(rain_sum_mjo < 0.01) = NaN ;
    rain_sum_mjo_eastward_portion(rain_sum_mjo_eastward_portion < 0.01) = NaN ;

    n_days = grand_master_end_time - grand_master_start_time + DT/24.0;
    DATA=rain_sum/n_days ; %log(rain_sum) ;
    DATA_mjo=rain_sum_mjo/n_days ; %log(rain_sum) ;
    DATA_mjo_eastward_portion=rain_sum_mjo_eastward_portion/n_days ; %log(rain_sum) ;

    

    %%
    %% Begin plot
    %%
    close all
    
    figure('visible','off')
    set(gcf, 'position',[100,100,800,1000])

    %% All LPTs
    subplot(311)


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

    
    %% MJO LPTs
    subplot(312)


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

    
    %% MJO LPTs eastward portion
    subplot(313)


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





    
    fileOutBase=['rain_map2_full_year_all_lpts__',y1_y2];

    
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
    fout.info = "Rainfall from 3B42 V7 in mm/day averaged over 1 year.";
    fout.n_days = n_days;
    fout.mask_time = grand_master_times;
    fout.mask_lon = this_mask.lon;
    fout.mask_lat = this_mask.lat;
    fout.mask_mask_lpt = logical(grand_master_mask > 0.5);
    fout.mask_mask_mjo = logical(grand_master_mask_mjo > 0.5);
    fout.mask_mask_mjo_eastward_portion = logical(grand_master_mask_mjo_eastward_portion > 0.5);
    
    disp([PLOT_DIR,'/',fileOutBase,'.mat'])
    eval(['save ',PLOT_DIR,'/',fileOutBase,'.mat -struct fout'])
    
end
