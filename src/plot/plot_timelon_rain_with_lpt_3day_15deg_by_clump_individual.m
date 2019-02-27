clear all
close all

addpath('../config')
options

PROCESSED_DATA_DIR = '../track' ; %['../data/',CASE_LABEL,'/processed/',...
                      %'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                      % sprintf('%d',ACCUMULATION_PERIOD), ...
                      % 'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

PLOT_DIR = '.' ; %['../plots/',CASE_LABEL,'/processed/',...
                  %    'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                  %     sprintf('%d',ACCUMULATION_PERIOD), ...
                  %     'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems/individual_clumps'];




figure('visible','off')
set(gcf,'position',[100,100,500,800])
set(gcf,'color','w')

corner_label={'5 deg. Filter','Threshold=12 mm/day'};
clumps=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin.txt'],'',1,0);
CLUMPS=clumps;
colors=hsv(12);


MJO=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin.txt'],'',1,0);

for year1=[2018]
%for year1=[1998:2017]

  clf

  %%%%%%%%%%%%%%%%%%%%%%%%

  year2=year1+1 ;

  yyyy1=num2str(year1) ;
  yyyy2=num2str(year2) ;

  y1_y2=[yyyy1,'_',yyyy2] ;
  % y11_y22=[yyyy1,'010400_',yyyy2,'063021'] ;
  if year1 == 2017
    y11_y22=[yyyy1,'060100_',yyyy2,'053121'] ;
  elseif year1 == 2018
    y11_y22=[yyyy1,'060100_',yyyy1,'112721'] ;
  else
    y11_y22=[yyyy1,'060100_',yyyy2,'063021'] ;
  end
  disp(y1_y2) ;

  F=load(['../data/trmm/interim/timelon/rain_hov_',y1_y2,'_15deg_3h_full_year.mat']) ;
  G=load([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',y11_y22,'.mat']) ;

  for iiii = 2:20

    if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
      eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
    end

  end


  F.rain=F.rain ;
  %F.rain(F.rain < 0.25) = NaN ;
  F.rain(F.rain < 0.35) = NaN ;


  hold on

  DATA=log(F.rain) ;


    
  %% Get "clumps of worms" for this year.
  clump_idx_this_year = find(clumps(:,1) == year1);
  lptid_this_year = clumps(clump_idx_this_year, 2)';
  clump_num_this_year = clumps(clump_idx_this_year, 3)';
  
  for this_clump_num = [unique(clump_num_this_year)]
    
    disp(['----------- Clump #', num2str(this_clump_num), ' -----------'])

    clf
    
    pcolor(F.lon,F.time-0.125,DATA) ;
    hold on


    shading flat
				%caxis([log(0.2),log(5)])
    caxis([log(0.2),log(2)])
    
    cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
    colormap(flipud(gray()))

       

    lptid_for_this_clump = lptid_this_year(clump_num_this_year == this_clump_num);

    lon1 = 999.0;
    lon2 = -999.0;
    dn1 = datenum(2100,1,1,0,0,0);
    dn2 = datenum(1900,1,1,0,0,0);
    
    for ii=[lptid_for_this_clump]%1:numel(G.TIMECLUSTERS)


      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.time=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=3.0*numel(GG.date)/24 ;

      lon1 = min([lon1, GG.lon]);
      lon2 = max([lon2, GG.lon]);
      dn1 = min([dn1, GG.date]);
      dn2 = max([dn2, GG.date]);
      
      
      this_clump_row = find(CLUMPS(:,1) == year1 & CLUMPS(:,2) == ii);
      this_clump_id = CLUMPS(this_clump_row,3);
      thisCol = colors(1+mod(this_clump_id, 10),:);
      
      
      longstats_plot_ts_circles(GG,[datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)],[40,200],[],thisCol,0.5) ;

				% Beginning/Ending Track Labels
      text(GG.lon(1), GG.time(1), num2str(ii),'clipping','on');
      text(GG.lon(end), GG.time(end), num2str(ii),'clipping','on');
    end

    %% Plot any MJO periods assoiated with this LPT system.
    for ii=[lptid_for_this_clump]%1:numel(G.TIMECLUSTERS)

      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.time=GG.time-1.5 ;
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
	  plot(GG.lon(idx1(idxx):idx2(idxx)), GG.time(idx1(idxx):idx2(idxx)), 'k-', 'linewidth', 2.0);
	end
      end
      
    end
    

    axis([lon1 - 5.0,lon2 + 5.0,dn1 - 5.0, dn2 + 5.0])


    TICKS=[datenum(year1,6,1,0,0,0):10:datenum(year2,6,1,0,0,0)];
    set(gca,'YTick',TICKS) ;
    datetick('y','mm/DD','keeplimits','keepticks')


    hcb=colorbar('NorthOutside') ;
    %    set(hcb,'xtick',log([.25,.5,1]),'xticklabel',{'6','12','24'})
    set(hcb,'xtick',log([.5,1,2]),'xticklabel',{'12','24','48'})
    
    set(hcb,'position',[0.15,0.91,0.8,0.01]) ;
    set(gca,'position',[0.15,0.05,0.8,0.8]);
    
    
    %set(gca,'xtick', 40:20:180) ;
    
    
    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on
    
    title(['Rainfall and LPTs (N=',num2str(numel(lptid_for_this_clump)),'): ',yyyy1,',  clump ',num2str(this_clump_num)])
    
    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')
    
    fileOutBase=['rain_filter_track_hov_',y1_y2,'_clump', sprintf('%03d', this_clump_num)];
    
    eval(['!mkdir -p ',PLOT_DIR])
    disp([PLOT_DIR,'/',fileOutBase,'.png'])
    saveas(gcf,[PLOT_DIR,'/',fileOutBase,'.png'])
    
  end

end
