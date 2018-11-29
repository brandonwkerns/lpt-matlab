clear all
close all

addpath('../config')
options

PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

PLOT_DIR = ['../plots/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems/individual_clumps'];




figure('visible','off')
set(gcf,'position',[100,100,800,800])
set(gcf,'color','w')

corner_label={'5 deg. Filter','Threshold=12 mm/day'};
clumps=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'],'',1,0);
CLUMPS=clumps;
colors=hsv(12);


MJO=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'],'',1,0);

%for year1=[2001]
for year1=[2007:2017]

  clf

  %%%%%%%%%%%%%%%%%%%%%%%%

  year2=year1+1 ;

  yyyy1=num2str(year1) ;
  yyyy2=num2str(year2) ;

  y1_y2=[yyyy1,'_',yyyy2] ;
  if year1 == 2017
    y11_y22=[yyyy1,'060100_',yyyy2,'053121'] ;
  else
    y11_y22=[yyyy1,'060100_',yyyy2,'063021'] ;
  end
  disp(y1_y2) ;

  G=load([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',y11_y22,'.rejoin2.mat']) ;

  if isfield(G, 'TIMECLUSTERS2')
    G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS2];
  end
  if isfield(G, 'TIMECLUSTERS3')
    G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS3];
  end
  if isfield(G, 'TIMECLUSTERS4')
    G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS4];
  end



    
  %% Get "clumps of worms" for this year.
  clump_idx_this_year = find(clumps(:,1) == year1);
  lptid_this_year = clumps(clump_idx_this_year, 2)';
  clump_num_this_year = clumps(clump_idx_this_year, 3)';
  
  for this_clump_num = [unique(clump_num_this_year)]
    
    disp(['----------- Clump #', num2str(this_clump_num), ' -----------'])

    clf
           

    lptid_for_this_clump = lptid_this_year(clump_num_this_year == this_clump_num);

    lon1 = 999.0;
    lon2 = -999.0;
    lat1 = 999.0;
    lat2 = -999.0;
    dn1 = datenum(2100,1,1,0,0,0);
    dn2 = datenum(1900,1,1,0,0,0);


    for ii=[lptid_for_this_clump]%1:numel(G.TIMECLUSTERS)


      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=3.0*numel(GG.date)/24 ;

      lon1 = min([lon1, GG.lon]);
      lon2 = max([lon2, GG.lon]);
      lat1 = min([lat1, GG.lat]);
      lat2 = max([lat2, GG.lat]);
      dn1 = min([dn1, GG.date]);
      dn2 = max([dn2, GG.date]);
      
    end

    plot_map_background([lon1-5.0, lon2+5.0, lat1-5.0, lat2+5.0]);
    hold on



    %% Rain

    count = 0;

    for dn_rain = dn1:0.125:dn2

      [y_rain, m_rain, d_rain, h_rain] = datevec(dn_rain);
      fileRain = ['../data/trmm/interim/gridded_rain_rates/gridded_rain_rates_',...
		  num2str(y_rain), sprintf('%02d', m_rain), sprintf('%02d', d_rain), sprintf('%02d',h_rain),'.nc'];
      disp(fileRain);
      F.precip=ncread(fileRain,'rain')' ;
      F.precip(~isfinite(F.precip)) = 0.0;
      F.precip(F.precip < 0.0) = 0.0;
      
      F.lon=ncread(fileRain,'lon');
      F.lat=ncread(fileRain,'lat');


      if count < 1
	rain_sum = 3.0 * F.precip;
      else
	rain_sum = rain_sum + 3.0 * F.precip;
      end
      
      count = count + 1;

    end
    

    rain_sum = rain_sum;
				%F.rain(F.rain < 0.25) = NaN ;
    rain_sum(rain_sum < 10) = NaN ;
    
    hold on
    
    DATA=rain_sum ; %log(rain_sum) ;


    pcolor(F.lon, F.lat, DATA) ;
    hold on

    shading flat
				%caxis([log(0.2),log(5)])
    caxis([0,1000.0])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    colormap(flipud(gray()))











    
    
    for ii=[lptid_for_this_clump]%1:numel(G.TIMECLUSTERS)


      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=3.0*numel(GG.date)/24 ;
      
      
      this_clump_row = find(CLUMPS(:,1) == year1 & CLUMPS(:,2) == ii);
      this_clump_id = CLUMPS(this_clump_row,3);
      thisCol = colors(1+mod(this_clump_id, 10),:);
      

      plot(GG.lon, GG.lat, 'color',thisCol, 'linewidth', 6);
      plot(GG.lon(1), GG.lat(1), 'ro','markerfacecolor','r', 'markersize', 15);
      plot(GG.lon(end), GG.lat(end), 'mx','markerfacecolor','m', 'markersize', 15, 'linewidth',3);
      
      hold on
      
      
				% Beginning/Ending Track Labels
      text(GG.lon(1), GG.lat(1), num2str(ii),'clipping','on');
      text(GG.lon(end), GG.lat(end), num2str(ii),'clipping','on');
    end


    for ii=[lptid_for_this_clump]%1:numel(G.TIMECLUSTERS)

      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=3.0*numel(GG.date)/24 ;
      
      idx1 = -999;
      idx2 = -999;
      
      if ( sum(MJO(:,1) == year1 & ...
               MJO(:,2) == ii) > 0 )
	
	idx1 = MJO((MJO(:,1) == year1 & ...
	            MJO(:,2) == ii),9);
	idx2 = MJO((MJO(:,1) == year1 & ...
	            MJO(:,2) == ii),10);
	
      end
      
      
      if idx1(1) > -1
        for idxx = 1:numel(idx1)
          plot(GG.lon(idx1(idxx):idx2(idxx)), GG.lat(idx1(idxx):idx2(idxx)), 'k-', 'linewidth', 3.0);
        end
      end

   
    end
    

    axis([lon1 - 5.0,lon2 + 5.0,lat1 - 5.0, lat2 + 5.0]);
    daspect([1,1,1]);



    hcb=colorbar('NorthOutside') ;

    set(hcb,'position',[0.15,0.91,0.8,0.01]) ;
    set(gca,'position',[0.15,0.05,0.8,0.8]);
    
        
    plot_map_background([lon1-5.0, lon2+5.0, lat1-5.0, lat2+5.0]);
	
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['Rainfall and LPTs (N=',num2str(numel(lptid_for_this_clump)),'): ',yyyy1,',  clump ',num2str(this_clump_num)])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')
    
    fileOutBase=['rain_filter_track_map2_',y1_y2,'_clump', sprintf('%03d', this_clump_num)];
    
    eval(['!mkdir -p ',PLOT_DIR])
    disp([PLOT_DIR,'/',fileOutBase,'.png'])
    saveas(gcf,[PLOT_DIR,'/',fileOutBase,'.png'])
    
  end

end
