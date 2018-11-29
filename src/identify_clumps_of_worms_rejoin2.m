clear all
close all


% do_plotting = true;
do_plotting = false;

%% NINO 3.4 stuff.
nino34=load('nino34.mat');


% Read in options that pertain to the entire tracking package.
% These settings are all in ../config/options.m
addpath('../config')
options

PROCESSED_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

EASTWARD_PROP_DATA_DIR = ['../data/',CASE_LABEL,'/processed/',...
                      'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                       sprintf('%d',ACCUMULATION_PERIOD), ...
                       'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/identify_eastward_propagation'];

eval(['!mkdir -p ', EASTWARD_PROP_DATA_DIR])

%%
%% Set the east propagation identification criteria here.
%%

min_zonal_speed = -999.0 ; % in m/s.
min_duration    = 7.0 ; % in Days. Doesn't include 3-Day accumulation period.
min_eastward_prop_zonal_speed    = 0.0 ; % in m/s.
min_eastward_prop_duration    = 7.0 ; % in Days. Doesn't include 3-Day accumulation period.

min_net_lon_propagation   = -999.0 ;%20.0 ; % in deg. longitude.
min_total_lon_propagation = 10.0 ;%20.0 ; % in deg. longitude.
max_abs_latitude = 7.5 ;% in deg. latitude. Must get this close to the Equator at some point.

mc_lon_1 = 100.0 ; % West end of MC for MC crossing
mc_lon_2 = 130.0 ; % East end of MC for MC crossing

% Search area for initial time cluster.
search_area = [50.0, 180.0, -15.0, 15.0];

%%
%%
%%

% For output table files.
FMT=['%10d%10d%10d%10.2f  %4d%0.2d%0.2d%0.2d  %4d%0.2d%0.2d%0.2d\n'];

header='      year     index     clump  duration       begin         end    ';

fid_clumps_of_worms=fopen([EASTWARD_PROP_DATA_DIR,'/clumps_of_worms.rejoin2.txt'],'w');

fprintf(fid_clumps_of_worms, '%s\n', header);

for year1 = 1998:2018  ;

  year2=year1+1 ;

  yyyy1=num2str(year1) ;
  yyyy2=num2str(year2) ;

  y1_y2=[yyyy1,'_',yyyy2] ;

  disp(y1_y2) ;


  dir0 = dir([PROCESSED_DATA_DIR,'/TIMECLUSTERS_lpt_',num2str(year1),'*.rejoin2.mat']);
  G=load([PROCESSED_DATA_DIR,'/', dir0(1).name]) ;

  for iiii = 2:20

    if isfield(G, ['TIMECLUSTERS', num2str(iiii)])
      eval(['G.TIMECLUSTERS = [G.TIMECLUSTERS, G.TIMECLUSTERS', num2str(iiii),'];'])
    end

  end


  count=0 ; % Keep count of east propagating systems.

  clump_ids = NaN * ones(1, numel(G.TIMECLUSTERS));
  clump_ids(1) = 1; % Start with clump ID #1
  next_clump_id = 2;
  for ii=2:numel(G.TIMECLUSTERS)

    disp([num2str(ii), ' of ', num2str(numel(G.TIMECLUSTERS))])
    found_a_match = false;

    for jj = setxor(ii,1:numel(G.TIMECLUSTERS))    %1:ii-1

      if numel(intersect(G.TIMECLUSTERS(ii).ceid, G.TIMECLUSTERS(jj).ceid)) > 0

        if (isfinite(clump_ids(jj)))
          clump_ids(ii) = clump_ids(jj);
          found_a_match = true;
        end
        % break;
      end

    end

    if ~found_a_match
      % If I made it here, there were no overlaps.
      clump_ids(ii) = next_clump_id;
      next_clump_id = next_clump_id + 1;
    end

  end %for ii=2:numel(G.TIMECLUSTERS)

  % Search clump id's that may still overlap.
  for this_clump_id = 1:nanmax(clump_ids)

    ceid1 = [];
    for ii = [find(clump_ids == this_clump_id)]
        ceid1 = unique([ceid1, G.TIMECLUSTERS(ii).ceid]);
    end

    for other_clump_id = setxor(this_clump_id, 1:nanmax(clump_ids))

      ceid2 = [];
      for ii = [find(clump_ids == other_clump_id)]
          ceid2 = unique([ceid2, G.TIMECLUSTERS(ii).ceid]);
      end

      if (numel(intersect(ceid1, ceid2)) > 0)

        clump_ids(clump_ids == other_clump_id) = this_clump_id;

      end
    end
  end



  %% Write out results.


  for this_clump_id = 1:nanmax(clump_ids)

    for ii = [find(clump_ids == this_clump_id)]

      GG=G.TIMECLUSTERS(ii) ;
      GG.date=GG.time-1.5 ;
      GG.size=sqrt(GG.area) ;
      GG.area=GG.area/1e4 ;
      GG.nentries=numel(GG.date) ;
      GG.duration=GG.date(end)-GG.date(1) ;

      [GG.year,GG.month,GG.day]=datevec(GG.date) ;
      [GG.year0,GG.month0,GG.day0,GG.hour0]=datevec(GG.time(1)) ;
      [GG.year1,GG.month1,GG.day1,GG.hour1]=datevec(GG.time(end)) ;

      fprintf(fid_clumps_of_worms,FMT,...
              year1,ii,this_clump_id, GG.duration,...
              GG.year0,GG.month0,GG.day0,GG.hour0,...
              GG.year1,GG.month1,GG.day1,GG.hour1);

    end

  end

end


fclose(fid_clumps_of_worms);
