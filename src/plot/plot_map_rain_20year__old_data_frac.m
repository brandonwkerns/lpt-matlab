clear all
close all

addpath('../config')
options

CASE_LABEL = 'trmm_keep_overlapping_tracks'

PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

MASK_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/masks'];

PLOT_DIR = ['../plots/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems/individual_clumps'];


corner_label={'5 deg. Filter','Threshold=12 mm/day'};
clumps=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.txt'],'',1,0);
CLUMPS=clumps;
colors=hsv(12);


MJO=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.txt'],'',1,0);

lon1 = 0.0;
lon2 = 360.0;
lat1 = -45.0;
lat2 = 45.0;

DATA = [];
DATA_mjo = [];
DATA_mjo_eastward_portion = [];

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

  fileInBase=['rain_map2_full_year_all_lpts__',y1_y2,'__old'];

  fn = [PLOT_DIR,'/',fileInBase,'.mat'];
  F = load(fn);

  F.rain_lpt(~isfinite(F.rain_lpt)) = 0.0;
  F.rain_mjo(~isfinite(F.rain_mjo)) = 0.0;
  F.rain_mjo_eastward_portion(~isfinite(F.rain_mjo_eastward_portion)) = 0.0;
  
  if (numel(DATA) < 1)
    DATA = F.rain_lpt;
    DATA_mjo = F.rain_mjo;
    DATA_mjo_eastward_portion = F.rain_mjo_eastward_portion;
  else
    DATA = DATA + F.rain_lpt;
    DATA_mjo = DATA_mjo + F.rain_mjo;
    DATA_mjo_eastward_portion = DATA_mjo_eastward_portion + F.rain_mjo_eastward_portion;
  end

end


DATA = DATA/20.0;
DATA_mjo = DATA_mjo/20.0;
DATA_mjo_eastward_portion = DATA_mjo_eastward_portion/20.0;

min_to_plot = 0.1; % mm/day 
DATA(DATA < min_to_plot) = NaN;
DATA_mjo(DATA_mjo < min_to_plot) = NaN;
DATA_mjo_eastward_portion(DATA_mjo_eastward_portion < min_to_plot) = NaN;

%%
%% Begin plot
%%
close all

figure('visible','off')
set(gcf, 'position',[100,100,800,1000])

%% All LPTs
subplot(311)

%{
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
%contour(this_mask.lon, this_mask.lat, max_mask, [0.5,0.5],'r-', 'linewidth',2);


hcb=colorbar('NorthOutside') ;

set(hcb,'position',[0.15,0.95,0.8,0.008]) ;


plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',12)
box on

title(['LPT Rainfall: 1998 - 2018'])

text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')
%}

%% MJO LPTs
subplot(312)


plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, 100.0 * DATA_mjo ./DATA) ;
hold on

shading flat
caxis([0,100.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)

hcb=colorbar('NorthOutside') ;
set(hcb,'position',[0.15,0.95,0.8,0.008]) ;


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask_mjo, [0.5,0.5],'r-', 'linewidth',2);

plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',12)
box on

title(['MJO LPT Rainfall Fraction of LPT: 1998 - 2018'])

text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')


%% MJO LPTs eastward portion
subplot(313)


plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, 100.0 * DATA_mjo_eastward_portion ./ DATA) ;
hold on

shading flat
caxis([0,100.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask_mjo_eastward_portion, [0.5,0.5],'r-', 'linewidth',2);

plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',12)
box on

title(['MJO LPT Eastward Portion Rainfall Fraction of LPT: 1998 - 2018'])

text(0.02,0.97,corner_label,'units','normalized', 'fontweight','bold')





    
fileOutBase=['rain_map2_20yr__old_frac'];


eval(['!mkdir -p ',PLOT_DIR])
disp([PLOT_DIR,'/',fileOutBase,'.png'])
saveas(gcf,[PLOT_DIR,'/',fileOutBase,'.png'])

