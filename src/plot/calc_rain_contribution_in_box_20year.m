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

%% Set the box here.
box_label = 'MC';
lon1 = 100.0;
lon2 = 140.0;
lat1 = -15.0;
lat2 = 15.0;

DATA_all = [];
DATA_lpo = [];
DATA_lpt = [];
DATA_mjo = [];
DATA_mjo_eastward_portion = [];
DATA_rossby = [];
DATA_rossby_westward_portion = [];

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

 
  lon_keep_indx = find(F.lon > lon1 - 0.01 & F.lon < lon2 + 0.01);
  lat_keep_indx = find(F.lat > lat1 - 0.01 & F.lat < lat2 + 0.01);
  
  if (numel(DATA_lpt) < 1)
    DATA_lpt = nanmean(nanmean(F.rain_lpt(lat_keep_indx, lon_keep_indx)));
    DATA_mjo = nanmean(nanmean(F.rain_mjo(lat_keep_indx, lon_keep_indx)));
    DATA_mjo_eastward_portion = nanmean(nanmean(F.rain_mjo_eastward_portion(lat_keep_indx, lon_keep_indx)));
    DATA_rossby = nanmean(nanmean(F.rain_rossby(lat_keep_indx, lon_keep_indx)));
    DATA_rossby_westward_portion = nanmean(nanmean(F.rain_rossby_westward_portion(lat_keep_indx, lon_keep_indx)));
  else
    DATA_lpt = DATA_lpt + nanmean(nanmean(F.rain_lpt(lat_keep_indx, lon_keep_indx)));
    DATA_mjo = DATA_mjo + nanmean(nanmean(F.rain_mjo(lat_keep_indx, lon_keep_indx)));
    DATA_mjo_eastward_portion = DATA_mjo_eastward_portion + nanmean(nanmean(F.rain_mjo_eastward_portion(lat_keep_indx, lon_keep_indx)));
    DATA_rossby = DATA_rossby + nanmean(nanmean(F.rain_rossby(lat_keep_indx, lon_keep_indx)));
    DATA_rossby_westward_portion = DATA_rossby_westward_portion + nanmean(nanmean(F.rain_rossby_westward_portion(lat_keep_indx, lon_keep_indx)));
  end

end

DATA_lpt = DATA_lpt/20.0;
DATA_mjo = DATA_mjo/20.0;
DATA_mjo_eastward_portion = DATA_mjo_eastward_portion/20.0;
DATA_rossby = DATA_rossby/20.0;
DATA_rossby_westward_portion = DATA_rossby_westward_portion/20.0;



%% Total rain
total_rain_map_dir = '/home/orca/bkerns/projects/writingup/pmm_lpt/code/lpt_rain_maps';
total_rain_map_file_list = dir([total_rain_map_dir,'/total_rain_map_*.nc']);

icount = 0;
for ii = 1:numel(total_rain_map_file_list)
  icount = icount + 1;
  fn = [total_rain_map_dir,'/',total_rain_map_file_list(ii).name];
  disp(fn)
  
  total_rain.lon = ncread(fn, 'lon');
  total_rain.lat = ncread(fn, 'lat');
  total_rain.rain_all = ncread(fn, 'avg_rain_rate')';
  
  lon_keep_indx = find(total_rain.lon > lon1 - 0.01 & total_rain.lon < lon2 + 0.01);
  lat_keep_indx = find(total_rain.lat > lat1 - 0.01 & total_rain.lat < lat2 + 0.01);

  if (numel(DATA_all) < 1)
    DATA_all = nanmean(nanmean(total_rain.rain_all(lat_keep_indx, lon_keep_indx)));
  else
    DATA_all = DATA_all + nanmean(nanmean(total_rain.rain_all(lat_keep_indx, lon_keep_indx)));
  end
  
end
  
DATA_all = 24.0 * DATA_all / icount;





%% LPO Rain
objects_rain_map_dir = '/home/orca/bkerns/projects/writingup/pmm_lpt/code/lpt_rain_maps';
objects_rain_map_file_list = dir([objects_rain_map_dir,'/objects_rain_map_*.nc']);

icount = 0;
for ii = 1:numel(objects_rain_map_file_list)
  icount = icount + 1;
  fn = [objects_rain_map_dir,'/',objects_rain_map_file_list(ii).name];
  disp(fn)
  
  objects_rain.lon = ncread(fn, 'lon');
  objects_rain.lat = ncread(fn, 'lat');
  objects_rain.rain_all = ncread(fn, 'avg_rain_rate_from_lp_objects')';
  
  lon_keep_indx = find(objects_rain.lon > lon1 - 0.01 & objects_rain.lon < lon2 + 0.01);
  lat_keep_indx = find(objects_rain.lat > lat1 - 0.01 & objects_rain.lat < lat2 + 0.01);

  if (numel(DATA_lpo) < 1)
    DATA_lpo = nanmean(nanmean(objects_rain.rain_all(lat_keep_indx, lon_keep_indx)));
  else
    DATA_lpo = DATA_lpo + nanmean(nanmean(objects_rain.rain_all(lat_keep_indx, lon_keep_indx)));
  end
  
end
  
DATA_lpo = 24.0 * DATA_lpo / icount;





%% Print out results.

disp(['Time Period: 1998 - 2018, Full Year (June to May)']
disp(['Box: ',num2str(lon1), ' - ', num2str(lon2), ', ', num2str(lat1), ', ', num2str(lat2)])
disp(['Total rain                          : ', num2str(DATA_all), ' mm/day'])
disp(['Objects rain                        : ', num2str(DATA_lpo), ' mm/day'])
disp(['LPT rain                            : ', num2str(DATA_lpt), ' mm/day'])
disp(['MJO LPT rain                        : ', num2str(DATA_mjo), ' mm/day'])
disp(['MJO LPT (eastward portion) rain     : ', num2str(DATA_mjo_eastward_portion), ' mm/day'])
disp(['Westward LPT rain                   : ', num2str(DATA_rossby), ' mm/day'])
disp(['Westward LPT (westward portion) rain: ', num2str(DATA_rossby_westward_portion), ' mm/day'])
