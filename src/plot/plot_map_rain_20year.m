clear all
close all


addpath('../../config')
options

%PROCESSED_DATA_DIR = ['../../data/',CASE_LABEL,'/processed/',...
%                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
%                       sprintf('%d',ACCUMULATION_PERIOD), ...
%                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

%MASK_DIR = ['../../data/',CASE_LABEL,'/processed/',...
%                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
%                       sprintf('%d',ACCUMULATION_PERIOD), ...
%                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/masks'];

PLOT_DIR = ['../../plots/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/systems/individual_clumps'];

colors=hsv(12);
colors_percent=hsv(12);


%MJO=dlmread(['../data/',CASE_LABEL,'/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'],'',1,0);

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

  fileInBase=['rain_map2_full_year_all_lpts__',y1_y2];

  fn = [PLOT_DIR,'/',fileInBase,'.mat'];
  F = load(fn);

  F.rain_lpt(~isfinite(F.rain_lpt)) = 0.0;
  F.rain_mjo(~isfinite(F.rain_mjo)) = 0.0;
  F.rain_mjo_eastward_portion(~isfinite(F.rain_mjo_eastward_portion)) = 0.0;
  F.rain_rossby(~isfinite(F.rain_rossby)) = 0.0;
  F.rain_rossby_westtward_portion(~isfinite(F.rain_rossby_westward_portion)) = 0.0;
  
  if (numel(DATA) < 1)
    DATA = F.rain_lpt;
    DATA_mjo = F.rain_mjo;
    DATA_mjo_eastward_portion = F.rain_mjo_eastward_portion;
    DATA_rossby = F.rain_rossby;
    DATA_rossby_westward_portion = F.rain_rossby_westward_portion;
  else
    DATA = DATA + F.rain_lpt;
    DATA_mjo = DATA_mjo + F.rain_mjo;
    DATA_mjo_eastward_portion = DATA_mjo_eastward_portion + F.rain_mjo_eastward_portion;
    DATA_rossby = DATA_rossby + F.rain_rossby;
    DATA_rossby_westward_portion = DATA_rossby_westward_portion + F.rain_rossby_westward_portion;
  end

end

%total_rain_map = read_total_rain_map();

DATA = DATA/20.0;
DATA_mjo = DATA_mjo/20.0;
DATA_mjo_eastward_portion = DATA_mjo_eastward_portion/20.0;
DATA_rossby = DATA_rossby/20.0;
DATA_rossby_westward_portion = DATA_rossby_westward_portion/20.0;

min_to_plot = 0.1; % mm/day 
DATA(DATA < min_to_plot) = NaN;
DATA_mjo(DATA_mjo < min_to_plot) = NaN;
DATA_mjo_eastward_portion(DATA_mjo_eastward_portion < min_to_plot) = NaN;
DATA_rossby(DATA_rossby < min_to_plot) = NaN;
DATA_rossby_westward_portion(DATA_rossby_westward_portion < min_to_plot) = NaN;

%%
%% Begin plot
%%
close all

%% All LPTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('visible','off')
plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, DATA) ;
hold on

shading flat
caxis([0,8.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask, [0.5,0.5],'r-', 'linewidth',2);


hcb=colorbar('NorthOutside') ;
set(hcb,'position',[0.05,0.95,0.4,0.008]) ;


plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',10)
box on

title(['LPT Rainfall'])


fileOutBase=['rain_map2_20yr__lpt'];
eval(['!mkdir -p ',PLOT_DIR])
disp(['./',fileOutBase,'.png'])
saveas(gcf,['./',fileOutBase,'.png'])




%% MJO LPTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('visible','off')

plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, DATA_mjo_eastward_portion) ;
hold on

shading flat
caxis([0,8.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask_mjo, [0.5,0.5],'r-', 'linewidth',2);

hcb=colorbar('NorthOutside') ;
set(hcb,'position',[0.05,0.95,0.4,0.008]) ;

plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',10)
box on

title(['MJO LPT Rainfall'])


fileOutBase=['rain_map2_20yr__lpt'];
eval(['!mkdir -p ',PLOT_DIR])
disp(['./',fileOutBase,'.png'])
saveas(gcf,['./',fileOutBase,'.png'])



%% MJO LPTs (FRACTION) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('visible','off')

plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, 100.0*DATA_mjo_eastward_portion./DATA) ;
hold on

shading flat
caxis([0,100.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask_mjo_eastward_portion, [0.5,0.5],'r-', 'linewidth',2);

hcb=colorbar('NorthOutside') ;
set(hcb,'position',[0.05,0.95,0.4,0.008]) ;

plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',10)
box on

title(['MJO LPT Rainfall % of LPT'])

fileOutBase=['rain_map2_20yr__mjo_frac'];
eval(['!mkdir -p ',PLOT_DIR])
disp(['./',fileOutBase,'.png'])
saveas(gcf,['./',fileOutBase,'.png'])


%% Westward LPTs westward portion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('visible','off')

plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, DATA_rossby_westward_portion) ;
hold on

shading flat
caxis([0,8.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask_mjo_eastward_portion, [0.5,0.5],'r-', 'linewidth',2);

hcb=colorbar('NorthOutside') ;
set(hcb,'position',[0.05,0.95,0.4,0.008]) ;


plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',10)
box on

title(['Westward LPT Rainfall'])


fileOutBase=['rain_map2_20yr__westward'];
eval(['!mkdir -p ',PLOT_DIR])
disp(['./',fileOutBase,'.png'])
saveas(gcf,['./',fileOutBase,'.png'])





%% Westward LPTs westward portion (FRACTION) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('visible','off')

plot_map_background([lon1, lon2, lat1, lat2]);
hold on


pcolor(F.lon, F.lat, 100.0*DATA_rossby_westward_portion./DATA) ;
hold on

shading flat
caxis([0,100.0])

cmap=dlmread('cmap_rain.dat') ;
				% colormap(cmap(9:end,:))
				%colormap(flipud(gray()))
colormap(cmap)


%% Outline of spatial influence of the LPT family.
%contour(this_mask.lon, this_mask.lat, max_mask_mjo_eastward_portion, [0.5,0.5],'r-', 'linewidth',2);

hcb=colorbar('NorthOutside') ;
set(hcb,'position',[0.05,0.95,0.4,0.008]) ;

plot_map_background([lon1, lon2, lat1, lat2]);
daspect([1,1,1]);

set(gca,'layer','top')
set(gca,'FontSize',10)
box on

title(['Westward LPT Rainfall % of LPT'])

fileOutBase=['rain_map2_20yr__westward_frac'];
eval(['!mkdir -p ',PLOT_DIR])
disp(['./',fileOutBase,'.png'])
saveas(gcf,['./',fileOutBase,'.png'])


