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
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems'];




figure('visible','off')
set(gcf,'position',[100,100,1000,800])
set(gcf,'color','w')

corner_label={'5 deg. Filter','Threshold=12 mm/day'};
CLUMPS=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'],'',1,0);
colors=hsv(12);


MJO=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'],'',1,0);

lon_range = [0, 360];
lon_ticks = 0:20:360;

%for year1=[2017]
for year1=[1998:2018]

  clf

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

  F=load(['../data/trmm/interim/timelon/rain_hov_',y1_y2,'_30deg_3h_full_year.mat']) ;
  G=load([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',y11_y22,'.rejoin2.mat']) ;

  for iiii = 2:20

    if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
      eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
    end

  end



  %F.rain=F.rain/24 ;
  %F.rain(F.rain < 0.25) = NaN ;
  F.rain(F.rain < 0.10) = NaN ;


  hold on

  DATA=log(F.rain) ;

  pcolor(F.lon,F.time-0.125,DATA) ;
  hold on


  shading flat
  %caxis([log(0.2),log(5)])
  caxis([log(0.1),log(2)])

  cmap=dlmread('cmap_rain.dat') ;
  % colormap(cmap(9:end,:))
  colormap(flipud(gray()))


  axis([lon_range(1),lon_range(2),datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)])


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
        plot(GG.lon(idx1(idxx):idx2(idxx)), GG.time(idx1(idxx):idx2(idxx)), 'b-', 'linewidth', 3.0);
      end
    end


    
  end



  TICKS=[datenum(year1,6,1,0,0,0):10:datenum(year2,6,1,0,0,0)];
  set(gca,'YTick',TICKS) ;
  datetick('y','mm/DD','keeplimits','keepticks')


  hcb=colorbar('NorthOutside') ;
  set(hcb,'xtick',log([0.25,.5,1,2]),'xticklabel',{'6','24','48'})

  set(hcb,'position',[0.15,0.91,0.8,0.01]) ;
  set(gca,'position',[0.15,0.05,0.8,0.8]);


  set(gca,'xtick', lon_ticks) ;


  set(gca,'layer','top')
  set(gca,'FontSize',12)
  box on

  title(['Rainfall and LPT: June ',yyyy1,' - May ',yyyy2])

  text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')

  fileOutBase=['rain_filter_track_hov_',y1_y2,'_wide_30deg_3h_mjo_only2'];

  eval(['!mkdir -p ',PLOT_DIR])
  disp([PLOT_DIR,'/',fileOutBase,'.png'])
  saveas(gcf,[PLOT_DIR,'/',fileOutBase,'.png'])

end
