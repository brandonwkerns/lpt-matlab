clear all
close all

YEAR1=1998;
YEAR2=1999;
addpath('../config')
options

PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

PLOT_DIR = ['../plots/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems'];




figure('visible','off')
set(gcf,'position',[100,100,500,800])
set(gcf,'color','w')

corner_label={'5 deg. Filter','Threshold=12 mm/day'};
CLUMPS=dlmread('../data/trmm_keep_overlapping_tracks/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.txt','',1,0);
colors=hsv(12);

MC_CROSSING=dlmread('../data/trmm_keep_overlapping_tracks/processed/g20_72h/thresh12/identify_eastward_propagation/list_mc_crossing_io_lpts.txt','',1,0);
NON_MC_CROSSING=dlmread('../data/trmm_keep_overlapping_tracks/processed/g20_72h/thresh12/identify_eastward_propagation/list_non_mc_crossing_io_lpts.txt','',1,0);
WPAC=dlmread('../data/trmm_keep_overlapping_tracks/processed/g20_72h/thresh12/identify_eastward_propagation/list_wpac_lpts.txt','',1,0);


% for year1=[2011]
% for year1=[1999:2016]
for year1=[1998:2017]

  clf

  %%%%%%%%%%%%%%%%%%%%%%%%

  year2=year1+1 ;

  yyyy1=num2str(year1) ;
  yyyy2=num2str(year2) ;

  y1_y2=[yyyy1,'_',yyyy2] ;
  % y11_y22=[yyyy1,'010400_',yyyy2,'063021'] ;
  if year1 == 2017
    y11_y22=[yyyy1,'060100_',yyyy2,'053121'] ;
  else
    y11_y22=[yyyy1,'060100_',yyyy2,'063021'] ;
  end
  disp(y1_y2) ;

  F=load(['../data/trmm/interim/timelon/rain_hov_',y1_y2,'_15deg_3day_full_year.mat']) ;
  G=load([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',y11_y22,'.mat']) ;

  if isfield(G, 'TIMECLUSTERS2')
    G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS2];
  end
  if isfield(G, 'TIMECLUSTERS3')
    G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS3];
  end
  if isfield(G, 'TIMECLUSTERS4')
    G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS4];
  end


  F.rain=F.rain/24 ;
  %F.rain(F.rain < 0.25) = NaN ;
  F.rain(F.rain < 0.35) = NaN ;


  hold on

  DATA=log(F.rain) ;

  pcolor(F.lon,F.time-3.0,DATA) ;
  hold on


  shading flat
  %caxis([log(0.2),log(5)])
  caxis([log(0.2),log(2)])

  cmap=dlmread('cmap_rain.dat') ;
  % colormap(cmap(9:end,:))
  colormap(flipud(gray()))


  axis([40,200,datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)])


  alreadyPlottedList=[-1] ;

  for ii=1:numel(G.TIMECLUSTERS)

    % if ( numel(find(alreadyPlottedList == ii)) > 0 )
    %     continue
    % end

    GG=G.TIMECLUSTERS(ii) ;
    GG.date=GG.time-1.5 ;
    GG.size=sqrt(GG.area) ;
    GG.area=GG.area/1e4 ;
    GG.nentries=numel(GG.date) ;
    GG.duration=3.0*numel(GG.date)/24 ;

    this_clump_row = find(CLUMPS(:,1) == year1 & CLUMPS(:,2) == ii);
    this_clump_id = CLUMPS(this_clump_row,3);
    thisCol = colors(1+mod(this_clump_id, 10),:);

    alreadyPlottedList=[alreadyPlottedList,ii];

    longstats_plot_ts_circles(GG,[datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)],[40,200],[],thisCol,0.2) ;

    % Beginning/Ending Track Labels
    text(GG.lon(1), GG.time(1), num2str(this_clump_id),'clipping','on');
    text(GG.lon(end), GG.time(end), num2str(this_clump_id),'clipping','on');
  end


  for ii=1:numel(G.TIMECLUSTERS)

    % if ( numel(find(alreadyPlottedList == ii)) > 0 )
    %     continue
    % end

    GG=G.TIMECLUSTERS(ii) ;
    GG.date=GG.time-1.5 ;
    GG.size=sqrt(GG.area) ;
    GG.area=GG.area/1e4 ;
    GG.nentries=numel(GG.date) ;
    GG.duration=3.0*numel(GG.date)/24 ;

    idx1 = -999;
    idx2 = -999;

    if ( sum(MC_CROSSING(:,1) == year1 & ...
             MC_CROSSING(:,2) == ii) > 0 )

      idx1 = MC_CROSSING((MC_CROSSING(:,1) == year1 & ...
	        MC_CROSSING(:,2) == ii),10);
      idx2 = MC_CROSSING((MC_CROSSING(:,1) == year1 & ...
	        MC_CROSSING(:,2) == ii),11);

    elseif ( sum(NON_MC_CROSSING(:,1) == year1 &  ...
                 NON_MC_CROSSING(:,2) == ii) > 0 )

      idx1 = NON_MC_CROSSING((NON_MC_CROSSING(:,1) == year1 & ...
	    NON_MC_CROSSING(:,2) == ii),10);
      idx2 = NON_MC_CROSSING((NON_MC_CROSSING(:,1) == year1 & ...
	    NON_MC_CROSSING(:,2) == ii),11);

    elseif ( sum(WPAC(:,1) == year1 &  ...
                 WPAC(:,2) == ii) > 0 )

      idx1 = WPAC((WPAC(:,1) == year1 & ...
	    WPAC(:,2) == ii),10);
      idx2 = WPAC((WPAC(:,1) == year1 & ...
	    WPAC(:,2) == ii),11);

    end

    if idx1 > -1
      plot(GG.lon(idx1:idx2), GG.time(idx1:idx2)-1.5, 'k-', 'linewidth', 1.5);
    end

  end



  TICKS=[datenum(year1,6,1,0,0,0):10:datenum(year2,6,1,0,0,0)];
  set(gca,'YTick',TICKS) ;
  datetick('y','mm/DD','keeplimits','keepticks')


  hcb=colorbar('NorthOutside') ;
  %    set(hcb,'xtick',log([.25,.5,1]),'xticklabel',{'6','12','24'})
  set(hcb,'xtick',log([.5,1,2]),'xticklabel',{'12','24','48'})

  set(hcb,'position',[0.15,0.91,0.8,0.01]) ;
  set(gca,'position',[0.15,0.05,0.8,0.8]);


  set(gca,'xtick', 40:20:180) ;
  % axis([40,200,datenum(year1,10,1,0,0,0),datenum(year2,4,1,0,0,0)])


  set(gca,'layer','top')
  set(gca,'FontSize',12)
  box on

  title(['Rainfall and LPT: June ',yyyy1,' - May ',yyyy2])

  text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')

  fileOutBase=['rain_filter_track_hov_',y1_y2,'_wide_15deg_3day_by_clump'];

  eval(['!mkdir -p ',PLOT_DIR])
  disp([PLOT_DIR,'/',fileOutBase,'.png'])
  saveas(gcf,[PLOT_DIR,'/',fileOutBase,'.png'])

end
