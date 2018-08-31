clear all
close all

figure('visible','off')
set(gcf,'position',[100,100,500,800])
set(gcf,'color','w')

corner_label={'5 deg. Filter','Threshold=12 mm/day'};

%% Read TCs
ibtracs_file = 'Allstorms.ibtracs_wmo.v03r10.nc';
ibtracs.lon = ncread(ibtracs_file,'lon_wmo');
ibtracs.lat = ncread(ibtracs_file,'lat_wmo');
ibtracs.raw_time = ncread(ibtracs_file,'time_wmo'); % days since 1858-11-17 00:00:00
ibtracs.dn = ibtracs.raw_time + datenum(1858,11,17,0,0,0);

for year1=[1998]
% for year1=[1999:2016]
% for year1=[2017]

    clf

    %%%%%%%%%%%%%%%%%%%%%%%%

    year2=year1+1 ;

    yyyy1=num2str(year1) ;
    yyyy2=num2str(year2) ;

    y1_y2=[yyyy1,'_',yyyy2] ;
    y11_y22=[yyyy1,'010400_',yyyy2,'063021'] ;
    % y11_y22=[yyyy1,'060100_',yyyy2,'063021'] ;
    % y11_y22=[yyyy1,'060100_',yyyy2,'053121'] ;

    disp(y1_y2) ;

    F=load(['../data/trmm/interim/timelon/rain_hov_',y1_y2,'_15deg_3day_full_year.mat']) ;
    G=load(['../data/trmm/processed/g20_72h/thresh12/timeclusters/TIMECLUSTERS_lpt_',y11_y22,'.mat']) ;

    F.rain=F.rain/24 ;
    %F.rain(F.rain < 0.25) = NaN ;
    F.rain(F.rain < 0.35) = NaN ;


    hold on

    DATA=log(F.rain) ;

    pcolor(F.lon,F.time-3.0,DATA) ;
    hold on


    shading flat
    %caxis([log(0.2),log(5)])
    caxis([log(0.35),log(2)])

    cmap=dlmread('cmap_rain.dat') ;
    colormap(cmap(9:end,:))


    axis([40,180,datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)])


    alreadyPlottedList=[-1] ;

    for ii=1:numel(G.TIMECLUSTERS)

        if ( numel(find(alreadyPlottedList == ii)) > 0 )
            continue
        end

        GG=G.TIMECLUSTERS(ii) ;
        GG.date=GG.time-1.5 ;
        GG.size=sqrt(GG.area) ;
        GG.area=GG.area/1e4 ;
        GG.nentries=numel(GG.date) ;
        GG.duration=3.0*numel(GG.date)/24 ;

        if (GG.duration < 7-0.001)

            alreadyPlottedList=[alreadyPlottedList,ii];
            thisCol=[0.6,0.6,0.6] ;

            longstats_plot_ts_circles(GG,[datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)],[40,180],[],thisCol,0.2) ;

        else

          thisCol='k' ;

          longstats_plot_ts_circles(GG,[datenum(year1,6,1,0,0,0),datenum(year2,6,1,0,0,0)],[40,180],[],thisCol,0.2) ;

	end

	text(GG.lon(1), GG.time(1), num2str(ii),'clipping','on');

    end




    %% Draw TC Tracks
    plot(ibtracs.lon, ibtracs.dn, 'm^','markersize',3);


    TICKS=[datenum(year1,6,1,0,0,0):10:datenum(year2,6,1,0,0,0)];
    set(gca,'YTick',TICKS) ;
    datetick('y','mm/DD','keeplimits','keepticks')




    hcb=colorbar('NorthOutside') ;
    %    set(hcb,'xtick',log([.25,.5,1]),'xticklabel',{'6','12','24'})
    set(hcb,'xtick',log([.5,1,2]),'xticklabel',{'12','24','48'})

    set(hcb,'position',[0.15,0.91,0.8,0.01]) ;
    set(gca,'position',[0.15,0.05,0.8,0.8]);


    set(gca,'xtick', 40:10:180) ;
    % axis([40,180,datenum(year1,6,1,0,0,0),datenum(year1,8,1,0,0,0)])


    set(gca,'layer','top')
    set(gca,'FontSize',12)
    box on

    title(['Rainfall and LPT: June ',yyyy1,' - May ',yyyy2])

    text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')
    % plot(0.02, 0.94, 'b^', 'markerfacecolor','b','units','normalized');
    text(0.02,0.92,'Best Track TC','color','m','units','normalized', 'fontweight','bold')

    fileOutBase=['rain_filter_track_hov_',y1_y2,'_wide_15deg_3day_black_with_TCs'];

    saveas(gcf,[fileOutBase,'.png'])

end
