%function TIMECLUSTERS=connect_time_clusters(master_ce_file, OPT, verbose)
clear all
close all


addpath('../../config')
options
save('temp.mat');
OPT = load('temp.mat');
eval('!rm temp.mat')

%verbose=0;
verbose=1;


% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
PROCESSED_DATA_DIR_IN = ['../../data/',CASE_LABEL,'/processed/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), ...
                        'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/objects'];

PROCESSED_DATA_DIR_OUT = ['../../data/',CASE_LABEL,'/processed/',...
                        'g',sprintf('%d',FILTER_STANDARD_DEVIATION), '_',...
                        sprintf('%d',ACCUMULATION_PERIOD), ...
                        'h/thresh',num2str(FEATURE_THRESHOLD_VALUE),'/timeclusters'];

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------


%%
%% TIMECLUSTERS=connect_time_clusters_g10_3day_t14_connect_branches(year, verbose)
%%
%% year = beginning year (e.g., 2011 for Oct. 2011 - March 2012)
%% verbose = 0 (default for off) or 1 for on.
%%
%% Approach here is to first group the cluster elements (CEs)
%%  in to their respective time clusters (TCs), allowing
%%  for mergers and splits.
%% Then, from the database of CEs, the tracking parameters are
%%  obtained, e.g., lat, lon of centroid and branch tracks.
%% (This function uses sub-functions each of which access and
%%  modify the TIMECLUSTERS struct array. So it cannot be run
%%  as a script--only a function.
%%

%% Set the tracking parameters here.
maxTimeToConnect = OPT.TRACKING_MAX_TIME_TO_RECONNECT;% 3 ; %Maximum time between clusters (hours) to
                     %connect them. Especially useful if there are
                     %some empty frames vs. the case when some
                     %times have data but do not have clusters.

minFramesToKeep = OPT.TRACKING_MINIMUM_FRAMES ; %Need this many entries in track to keep it.
minDuration = OPT.TRACKING_MINIMUM_DURATION; % Minimum duration to keep it (hours).

maxLat = OPT.FEATURE_MAX_LAT ; %Positive. Max distance off the equator to keep it.

maxDistToConnect = 9999.0 ; %20.0; %Max. distance between centroids to match them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% Intitialize time cluster stuff.
%%

TIMECLUSTERS=[] ;
nextClusterID=1 ;  %Start time cluster ID count at 1.
allDates=[] ;

DN=OPT.DN1:datenum(0,0,0,OPT.DT,0,0):OPT.DN2;

%%
%% Loop forward in time.
%% -- Mergers will be represented as duplicate, overlapping TCs.
%%    These need to be re-merged.
%% -- Splits will have the splitting branches as the same TC.
%%
for dn=[DN]

  [dnY,dnM,dnD,dnH] = datevec(dn);
  YYYY = num2str(dnY);
  MM = sprintf('%02d',dnM);
  DD = sprintf('%02d',dnD);
  HH = sprintf('%02d',dnH);
  
  if verbose
    disp(['--- ',YYYY, MM, DD, HH,' ---'])
  end

  %% Load Objects for this time.
  OBJ = load([PROCESSED_DATA_DIR_IN, '/', YYYY, '/', MM, '/objects_', YYYY,MM,DD,HH,'.mat']);
  
  objINDXthisTime=find(OBJ.time > dn - datenum(0,0,0,0,1,0) & OBJ.time < dn + datenum(0,0,0,0,1,0)) ;
  if ( numel(objINDXthisTime) < 1 )
    continue
  end

  if (numel(allDates) > 0.1)

    prevDates = allDates(allDates > dn - datenum(0,0,0,0,1,0) - OPT.DT/24.0 & ...
                         allDates < dn - datenum(0,0,0,0,1,0) ) ;

  else
    prevDates = -999.0;
  end

  already_matched_tc_list = [-999];


  % Get a set of matching time clusters for each
  %  cluster feature present at this time.

  for this_objid = OBJ.id(objINDXthisTime)  % Loop over the cluster features
                               % present at this time.
    istr=num2str(this_objid); %thisObj.id) ;

    for prevDate = prevDates
      matchingClusterID = matchingTimeCluster(OBJ, this_objid, TIMECLUSTERS, prevDate, maxDistToConnect, PROCESSED_DATA_DIR_IN, OPT);

      if ( numel( matchingClusterID ) <  1   )

        %% This is a new cluster.
	
        if ( verbose > 0 )
          disp([istr,': No matching time cluster in the database.'])
        end

        %% So I am starting new cloud cluster.
        TIMECLUSTERS = startNewTimeCluster(TIMECLUSTERS, nextClusterID);
        if ( verbose > 0 )
          disp(['        + Created new time cluster: ID = ',num2str(nextClusterID)])
        end

        TIMECLUSTERS = addCluster(TIMECLUSTERS, nextClusterID, this_objid) ;
        nextClusterID = nextClusterID + 1 ;

      else

	%% This is either a splitting LPT, or the continuation of a single LPT.
	
        if ( verbose > 0 )
          disp([istr,': Matches with current existing cluster id = ',num2str(matchingClusterID)])
        end

        for iii=1:numel(matchingClusterID)

	  if (sum(matchingClusterID(iii) == already_matched_tc_list) > 0)

	    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    %% ------- !!! SPLIT TRACKS !!! -----
            %% If two or more objects match to an LPT, I have a split.
            %% At this stage, split tracks will retain the full history
	    %% of the previous track before the split.
	    %%
            %% This will create duplicate tracks when there are more than one
            %%  matchingClusterID's. To fix this, use below either
            %%  mergeDuplicateTCs() or splitDuplicateTCs().
            %%
	    %% Depending on the tracking options, this will be refined below,
	    %% or it can be refined externally (e.g., using "rejoin" scripts).
	    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    
	    %% 1) Start the new time cluster
            TIMECLUSTERS = startNewTimeCluster(TIMECLUSTERS, nextClusterID) ;

	    %% 2) Copy the objids leading up to this point.
            for checkThisObjid = [TIMECLUSTERS(matchingClusterID(iii)).objid]
	      OBJprev = load_obj_data(PROCESSED_DATA_DIR_IN, checkThisObjid);
	      if (OBJprev.time < dn-0.001) % Only add record *before* the current time (dn).
		TIMECLUSTERS = addCluster(TIMECLUSTERS, nextClusterID, checkThisObjid);
		%disp(num2str(TIMECLUSTERS(nextClusterID).objid))
              end
            end

	    TIMECLUSTERS(nextClusterID).objid = unique(TIMECLUSTERS(nextClusterID).objid);
	    
	    %% 3) Add the current objid to the new LPT.
            TIMECLUSTERS = addCluster(TIMECLUSTERS, nextClusterID, this_objid) ;
            disp(['        + Split! Copied to new overlapping time cluster: ID = ',num2str(nextClusterID)])
	    %% 4) Set the counter to the next available LPT ID.
            nextClusterID = nextClusterID + 1 ;
	    
          else

	    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	    
	    %% In this case there was NO SPLIT TRACK.
	    %% Just attach to end of the previous track.
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            TIMECLUSTERS = addCluster(TIMECLUSTERS, matchingClusterID(iii),this_objid) ;
            already_matched_tc_list = [already_matched_tc_list, matchingClusterID(iii)];
          end
	  
        end
      end
    end

  end   % END Loop over the cluster features
        % present at this time.

  % Update the list of clusters.
  allDates=unique([allDates,dn]) ;

end % END Loop over Dates

%%
%% Some TCs were close enough to be considered a single track. Allow "center jumps".
%%
TIMECLUSTERS = combineCloseProximityTCs_by_area(TIMECLUSTERS, PROCESSED_DATA_DIR_IN, OPT.ACCUMULATION_PERIOD/24.0, OPT, maxDistToConnect, verbose);

%%
%% Now combining the merging tracks which were duplicates.
%%

if OPT.SPLITTING_AND_MERGING_METHOD == 1
  splitDuplicateTCs_split_by_instantaneous_area() ;
elseif OPT.SPLITTING_AND_MERGING_METHOD == 2
  splitDuplicateTCs_maximize_accumulated_area() ;
elseif OPT.SPLITTING_AND_MERGING_METHOD == 3
  mergeDuplicateTCs() ;
end
%% Otherwise (e.g., OPT.SPLITTING_AND_MERGING_METHOD == 4), keep each of the individual overlapping tracks.

%%
%% Take out time clusters with insufficient duration.
%%

TIMECLUSTERS = removeShortLivedTCs(TIMECLUSTERS, PROCESSED_DATA_DIR_IN, minDuration, verbose) ;
TIMECLUSTERS = eliminate_overlapping_tracks(TIMECLUSTERS, verbose);
TIMECLUSTERS = put_tracks_in_order(TIMECLUSTERS, verbose);


%%
%% Get tracking parameters from the CE database
%%

disp('Calculating tracking parameters.')
TIMECLUSTERS = calc_tracking_parameters(TIMECLUSTERS, PROCESSED_DATA_DIR_IN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Output  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dn0 = OPT.DN1;%f0.time;
dn9 = OPT.DN2;%f9.time;

[year0,month0,day0,hour0] = datevec(dn0);
YYYY0 = sprintf('%d', year0);
MM0 = sprintf('%02d', month0);
DD0 = sprintf('%02d', day0);
HH0 = sprintf('%02d', hour0);

[year9,month9,day9,hour9] = datevec(dn9);
YYYY9 = sprintf('%d', year9);
MM9 = sprintf('%02d', month9);
DD9 = sprintf('%02d', day9);
HH9 = sprintf('%02d', hour9);

ymd0_ymd9 = [YYYY0,MM0,DD0,HH0,'_',YYYY9,MM9,DD9,HH9];

disp('Writing Output.')
%FMT='        %4d%02d%02d%02d %7d %10.2f %10.2f %1d\n' ;

%% Make directory if it doesn't exist.
mkdir(PROCESSED_DATA_DIR_OUT);

%% Ascii output
fileout=[PROCESSED_DATA_DIR_OUT,'/LONGSTATS_lpt_',ymd0_ymd9] ;
disp(fileout)
lpt_systems_output_ascii(TIMECLUSTERS, fileout);

%% NetCDF Output
netcdf_output_fn=[PROCESSED_DATA_DIR_OUT,'/TIMECLUSTERS_lpt_',ymd0_ymd9,'.lptALL.nc'] ;
disp(netcdf_output_fn);
lpt_systems_output_netcdf(TIMECLUSTERS, netcdf_output_fn, OPT);

%% Mat file output
fileout_mat=[PROCESSED_DATA_DIR_OUT,'/TIMECLUSTERS_lpt_',ymd0_ymd9,'.mat'] ;
disp(fileout_mat);
lpt_systems_output_mat(TIMECLUSTERS, fileout_mat)

disp('Done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Local functions below here  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function TIMECLUSTERS=startNewTimeCluster(TIMECLUSTERS, tcid)
  TIMECLUSTERS(tcid).objid = [] ;
end


function TIMECLUSTERS = addCluster(TIMECLUSTERS, tcid,objid)
%
% tcid = index of TIMECLUSTERS struct array to use.
% pixels = "pixels" struct for this time
% pixelindx = which "pixels" entry to use.
%
% Modifies global function variable TIMECLUSTERS
%
  TIMECLUSTERS(tcid).objid=[TIMECLUSTERS(tcid).objid,objid];
end


function  matchingClusterID=matchingTimeCluster(OBJ, objid, TIMECLUSTERS, timeToSearch, maxDistToConnect, objects_data_dir, OPT)

  %%
  %% -- Function to match an LP Object with a collection of TIMECLUSTERS. --
  %% Usage: matchingClusterID=matchingTimeCluster(OBJ, objid, TIMECLUSTERS, timeToSearch)
  %%
  %% Inputs:
  %%   OBJ is the struct of objects containing the time "timeToSearch"
  %%   objid is the Object ID we are interestd in.
  %%   TIMECLUSTERS is the struct of timeclusters to look for matches.
  %%   timeToSearch is the time to look for matches in the timeclusters struct.
  %%   maxDistToConnect is the maximum distance between centroids (degrees) to connect them.
  %%   objects_data_dir is the parent directory of the objects files.
  %% Output:
  %%   matchingClusterID is an array with the index of the matching timecluster(s)
  %% Notes:
  %%   In merger cases or with overlapping TIMECLUSTERS, matchingClusterID will have multiple entries.
  %%
  
  matchingClusterID=[] ;

  %% If still at the initial time of timeclusters, nothing to search.
  if ( timeToSearch < -900 )
    return
  end
  
  %% If there are no TIMECLUSTERS yet, return empty array. This will force starting new TIMECLUSTERS.
  if ( numel(TIMECLUSTERS) < 1 )
    return
  else

    %% Load Objects for the timeToSearch.
    [dnY2,dnM2,dnD2,dnH2] = datevec(timeToSearch);
    YYYY2 = num2str(dnY2);
    MM2 = sprintf('%02d',dnM2);
    DD2 = sprintf('%02d',dnD2);
    HH2 = sprintf('%02d',dnH2);
    OBJ2 = load([objects_data_dir, '/', YYYY2, '/', MM2, '/objects_', YYYY2,MM2,DD2,HH2,'.mat']);

    %%
    %% Loop through previous time clusters, looking for a match.
    %%
    for ii=1:numel(TIMECLUSTERS)

      times=[] ;

      for  jj=[TIMECLUSTERS(ii).objid]

	%% If it's not at the right time, then no need to consider it further.
        if ( sum(OBJ2.id == jj) > 0)

          %% Check distance here.
	  %% I can skip it if it is beyond maxDistToConnect.
          dLON1 = OBJ2.lon(OBJ2.id == jj) - OBJ.lon(OBJ.id == objid) ;
          dLAT1 = OBJ2.lat(OBJ2.id == jj) - OBJ.lat(OBJ.id == objid) ;
          if ( sqrt(dLON1.^2 + dLAT1.^2 ) > maxDistToConnect )
            continue
          end

	  %% OK, the centroid is within maxDistToConnect, so determine how much overlap.
          X1=OBJ2.pixels(OBJ2.id == jj).x ;
          Y1=OBJ2.pixels(OBJ2.id == jj).y ;

          X2=OBJ.pixels(OBJ.id == objid).x ;
          Y2=OBJ.pixels(OBJ.id == objid).y ;

          hits=[] ;

          A=[X1,Y1] ;
          B=[X2,Y2] ;

          HITS=intersect(A,B,'rows') ;
          nHits=size(HITS,1);

          if ( nHits/numel(X1) > OPT.TRACKING_MINIMUM_OVERLAP_FRAC | ...
               nHits/numel(X2) > OPT.TRACKING_MINIMUM_OVERLAP_FRAC | ...
               nHits > OPT.TRACKING_MINIMUM_OVERLAP_POINTS )

            matchingClusterID=[matchingClusterID,ii] ;

          end

        end

      end

    end

  end

  matchingClusterID = unique(matchingClusterID) ;

end



function mergeDuplicateTCs()

  %% Accesses main function variables: TIMECLUSTERS
  %% Modifies TIMECLUSTERS

  still_work_to_be_done = true;
  throwAwayTCs=[] ;
  starting_ce_index = 1;
  iteration = 0;
  while still_work_to_be_done

    iteration = iteration + 1;
    disp(['----  iteration ',num2str(iteration),'  ----'])
    allTCindices = starting_ce_index:numel(TIMECLUSTERS) ;
    still_work_to_be_done = false;

    for iiii=[allTCindices]

      otherTCindices=setdiff(allTCindices,iiii) ;
      thisClusterCEList=TIMECLUSTERS(iiii).objid ;

      found_a_match = false;
      for jjjj=[otherTCindices]

        if ( numel(intersect(TIMECLUSTERS(iiii).objid,...
                             TIMECLUSTERS(jjjj).objid))>0)

          found_a_match = true;
          if verbose
            disp([num2str(jjjj), ' absorbed into ', num2str(iiii)])
          end

          TIMECLUSTERS(iiii).objid=sort(union(TIMECLUSTERS(iiii).objid,...
                                        TIMECLUSTERS(jjjj).objid));

          TIMECLUSTERS(jjjj).objid=[] ;
          throwAwayTCs=[throwAwayTCs,jjjj] ;

          still_work_to_be_done = true;
          break

        end
      end

      % If I get to this point, it means there were no overlapping tracks.
      % so, no need to check this track over and over again.
      if ~found_a_match
        starting_ce_index = iiii+1;
      end

      if (still_work_to_be_done)
        break
      end


    end

  end
  %%% Throw the flagged tracks away.
  throwAwayTCs = unique(throwAwayTCs);
  if verbose & numel(throwAwayTCs) > 0
    disp(['Removed the following due to being completely absorbed in merging step: ',num2str(throwAwayTCs)])
    disp(['(Database also re-ordered.)'])
  end
  TIMECLUSTERS(throwAwayTCs)=[] ;

end


function splitDuplicateTCs_split_by_instantaneous_area()

  %% Accesses main function variables: TIMECLUSTERS
  %% Modifies TIMECLUSTERS
  %%
  %% When two timecluster tracks merge, only keep the one with the
  %% longer history in terms of accumulated added area.

  still_work_to_be_done = true;
  throwAwayTCs = [];

  iteration = 0;
  starting_ce_index = 1;
  while still_work_to_be_done

    iteration = iteration + 1;
    disp(['----  iteration ',num2str(iteration),'  ----'])
    % allTCindices = 1:numel(TIMECLUSTERS) ;
    allTCindices = starting_ce_index:numel(TIMECLUSTERS) ;

    for iiii = [allTCindices]

      still_work_to_be_done = false;
      if (numel(TIMECLUSTERS(iiii).objid) < 1)
        starting_ce_index = iiii + 1;
        continue
      end
      otherTCindices=setdiff(allTCindices,iiii) ;

      intersectingTCindices = [];
      intersectingTCbeginning_times = [];
      intersecting_times = [];

      for kk = TIMECLUSTERS(iiii).objid
        intersecting_times = [intersecting_times, CE.time(TIMECLUSTERS(iiii).objid)];
      end
      for jjjj = [otherTCindices]
        intersections=intersect(TIMECLUSTERS(iiii).objid,...
                                TIMECLUSTERS(jjjj).objid);
        if (numel(intersections) > 0)
          intersectingTCindices = [intersectingTCindices, jjjj];
          intersectingTCbeginning_times = [intersectingTCbeginning_times, CE.time(TIMECLUSTERS(jjjj).objid(1))];
          intersecting_times = [intersecting_times, CE.time(TIMECLUSTERS(jjjj).objid)];
        end
      end

      intersecting_times = sort(intersecting_times);
      [intersectingTCbeginning_times, intersectingTCbeginning_times_indx] = sort(intersectingTCbeginning_times);
      intersectingTCindices = intersectingTCindices(intersectingTCbeginning_times_indx);
      if numel(intersectingTCindices) > 0
        if verbose
          disp([num2str(iiii), ' overlaps: ', num2str(intersectingTCindices)])
        end

        for jjjj = intersectingTCindices
          intersections=intersect(TIMECLUSTERS(iiii).objid,...
                                  TIMECLUSTERS(jjjj).objid);

      	  if (numel(intersections) < 1)
      	    continue
      	  end

          % Check if I have an overlapping/duplicate track here.
          if (numel(setxor(intersections, TIMECLUSTERS(iiii).objid)) < 1)
            if verbose
              disp(['Duplicate! Will throw away ID=',num2str(iiii)])
            end
            throwAwayTCs = [throwAwayTCs, iiii];
            TIMECLUSTERS(iiii).objid = [];
            continue
          end
          if (numel(setxor(intersections, TIMECLUSTERS(jjjj).objid)) < 1)
            if verbose
              disp(['Duplicate! Will throw away ID=',num2str(jjjj)])
            end
            throwAwayTCs = [throwAwayTCs, jjjj];
            TIMECLUSTERS(jjjj).objid = [];
            continue
          end

          %% "split" track portion, if it exists.
          %% (If no split, nothing happens here.)
          split_point_objid = -999;
          for ll = [TIMECLUSTERS(iiii).objid(TIMECLUSTERS(iiii).objid >= min(intersections))]
            if sum(intersections == ll) > 0
              split_point_objid = ll;
            else
              break;
            end
          end
      	  if split_point_objid < 1
            for ll = [TIMECLUSTERS(jjjj).objid(TIMECLUSTERS(jjjj).objid >= min(intersections))]
              if sum(intersections == ll) > 0
                split_point_objid = ll;
              else
                break;
              end
            end
      	  end
          indices1 = TIMECLUSTERS(iiii).objid(TIMECLUSTERS(iiii).objid > split_point_objid);
          indices2 = TIMECLUSTERS(jjjj).objid(TIMECLUSTERS(jjjj).objid > split_point_objid);


          if (numel(indices1) < 1 | numel(indices2) < 1)
            if verbose
              disp('No split. Moving on.')
            end
          else
            if CE.area(indices1(1)) > CE.area(indices2(1))
              startNewTimeCluster(nextClusterID) ;
              if ( verbose > 0 )
                disp(['        A track split off (ID=',num2str(jjjj),').'])
                disp(['        + Split portion is a new time cluster: ID = ',num2str(nextClusterID)])
              end
              TIMECLUSTERS(nextClusterID).objid=...
                  TIMECLUSTERS(jjjj).objid(TIMECLUSTERS(jjjj).objid > split_point_objid);
              nextClusterID = nextClusterID + 1 ;
              TIMECLUSTERS(jjjj).objid(TIMECLUSTERS(jjjj).objid > split_point_objid)=[];
              %break
            else
              startNewTimeCluster(nextClusterID) ;
              if ( verbose > 0 )
                disp(['        A track split off (ID=',num2str(iiii),').'])
                disp(['        + Split portion is a new time cluster: ID = ',num2str(nextClusterID)])
              end
              TIMECLUSTERS(nextClusterID).objid=...
                  TIMECLUSTERS(iiii).objid(TIMECLUSTERS(iiii).objid > split_point_objid);
              nextClusterID = nextClusterID + 1 ;
              TIMECLUSTERS(iiii).objid(TIMECLUSTERS(iiii).objid > split_point_objid)=[];
              %break
            end
          end

          %% "merge" track portion, if it exists.
          %% (If no merge, nothing happens here.)

          indices11 = TIMECLUSTERS(iiii).objid(TIMECLUSTERS(iiii).objid < min(intersections));
          indices22 = TIMECLUSTERS(jjjj).objid(TIMECLUSTERS(jjjj).objid < min(intersections));

          if (numel(indices11) < 1 | numel(indices22) < 1)
            if verbose
              disp('No merger. Moving on.')
            end
          else
            accum_area1 = 0.0;
            for idx = [indices11]
              accum_area1 = accum_area1 + CE.area(idx);
            end
            accum_area2 = 0.0;
            for idx = [indices22]
              accum_area2 = accum_area2 + CE.area(idx);
            end

            if accum_area1 > accum_area2
              if verbose
                disp(['  Merger --> ', num2str(iiii), ' wins,', num2str(jjjj),' cut short.'])
              end
              TIMECLUSTERS(jjjj).objid=indices22;
              %break
            else
              if verbose
                disp(['  Merger --> ', num2str(jjjj), ' wins,', num2str(iiii),' cut short.'])
              end
              TIMECLUSTERS(iiii).objid=indices11;
              %break
            end
          end
        end

        %% If these tracks are duplicates, mark the second entry for removal.
        still_work_to_be_done = true;
        break
      else
        starting_ce_index = iiii + 1;
      end %if numel(intersectingTCindices) > 0

      if still_work_to_be_done
        break
      end
    end %for iiii=[allTCindices]
  end %while still_work_to_0be_done

  %%% Throw the flagged tracks away.
  throwAwayTCs = unique(throwAwayTCs);
  if verbose & numel(throwAwayTCs) > 0
    disp(['Removed the following due to being completely absorbed in merging step: ',num2str(throwAwayTCs)])
    disp(['(Database also re-ordered.)'])
  end
  TIMECLUSTERS(throwAwayTCs)=[] ;

end


function splitDuplicateTCs_maximize_accumulated_area()

  %% Accesses main function variables: TIMECLUSTERS
  %% Modifies TIMECLUSTERS
  %%
  %% When two timecluster tracks merge, only keep the one with the
  %% longer history in terms of accumulated added area.

  still_work_to_be_done = true;
  throwAwayTCs = [];

  iteration = 0;
  while still_work_to_be_done

    iteration = iteration + 1;
    disp(['----  iteration ',num2str(iteration),'  ----'])
    allTCindices = 1:numel(TIMECLUSTERS) ;

    for iiii = [allTCindices]

      still_work_to_be_done = false;
      otherTCindices=setdiff(allTCindices,iiii) ;

      intersectingTCindices = [];
      for jjjj = [otherTCindices]
        intersections=intersect(TIMECLUSTERS(iiii).objid,...
                                TIMECLUSTERS(jjjj).objid);

        if (numel(intersections) > 0)
          intersectingTCindices = [intersectingTCindices, jjjj];
        end

      end

      if numel(intersectingTCindices) > 0
        if verbose
          disp([num2str(iiii), ' overlaps: ', num2str(intersectingTCindices)])
        end

        %Get accumulated area of thisCluster
        thisClusterAccumArea=0.0 ;
        for thisClusterOBJID=TIMECLUSTERS(iiii).objid
          thisClusterAccumArea=thisClusterAccumArea+...
              CE.area(thisClusterOBJID);
        end

        %Get accumulated area of otherCluster's
        otherClusterAccumArea = zeros(1, numel(intersectingTCindices));
        for kk = 1:numel(intersectingTCindices)
          for otherClusterOBJID=[TIMECLUSTERS(intersectingTCindices(kk)).objid]
            otherClusterAccumArea(kk)=otherClusterAccumArea(kk) + ...
                CE.area(otherClusterOBJID);
          end
        end

        if ( thisClusterAccumArea > max(otherClusterAccumArea) )

          if verbose
            disp(['  --> ', num2str(iiii), ' wins.'])
          end

          winner_index = iiii;
          loser_indices = sort(intersectingTCindices);

        else

          [maxAccumArea, maxAccumAreaIndx] = max(otherClusterAccumArea);

          if verbose
            disp(['  --> ', num2str(intersectingTCindices(maxAccumAreaIndx)), ' wins.'])
          end

          winner_index = intersectingTCindices(maxAccumAreaIndx);
          loser_indices = intersectingTCindices;
          loser_indices(maxAccumAreaIndx) = [];
          loser_indices = sort([loser_indices, iiii]);

        end

        %% Clip the tracks of each of the "loser" indices.
        %% Mark it for elimination, or split it up if needed.
        for this_loser_index = [loser_indices]

          % NOTE!!! The "winner" index is not necessarily the starting "iiii" indx!
          % Therefore, we need to check for intersections again.
          intersections=intersect(TIMECLUSTERS(winner_index).objid,...
                                  TIMECLUSTERS(this_loser_index).objid);
          if (numel(intersections) < 1)
            if verbose
              disp(['Skipping ', num2str(this_loser_index), ' since it does not intersect the Winner.'])
            end
            continue
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %% The "loser" track can either merge on to, or split off of the "winner" track,
          %%  OR it can do both! Handle each of these scanarios below.
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          %% "split" track portion, if it exists.
          %% (This does NOT get done if there is no split)
          if (sum(TIMECLUSTERS(this_loser_index).objid > max(intersections)) > 0)
            startNewTimeCluster(nextClusterID) ;
            if ( verbose > 0 )
              disp(['        A track split off (ID=',num2str(this_loser_index),').'])
              disp(['        + Split portion is a new time cluster: ID = ',num2str(nextClusterID)])
            end
            TIMECLUSTERS(nextClusterID).objid=...
                TIMECLUSTERS(this_loser_index).objid(TIMECLUSTERS(this_loser_index).objid > max(intersections));
            nextClusterID = nextClusterID + 1 ;
          end

          %% "merge" track portion, if it exists.
          %% (If no merge, the track gets flagged for removal.)
          TIMECLUSTERS(this_loser_index).objid=...
              TIMECLUSTERS(this_loser_index).objid(TIMECLUSTERS(this_loser_index).objid < min(intersections));

          if (numel(TIMECLUSTERS(this_loser_index).objid) < 1)
            if verbose
              disp(['Completely absorbed! Will throw away ID=',num2str(this_loser_index)])
            end
            throwAwayTCs = [throwAwayTCs, this_loser_index];
          end

        end

        still_work_to_be_done = true;
        break
      end %if numel(intersectingTCindices) > 0

    end %for iiii=[allTCindices]
  end %while still_work_to_be_done

  % Throw the flagged tracks away.
  throwAwayTCs = unique(throwAwayTCs);
  if verbose & numel(throwAwayTCs) > 0
    disp(['Removed the following due to being completely absorbed in merging step: ',num2str(throwAwayTCs)])
    disp(['(Database also re-ordered.)'])
  end
  TIMECLUSTERS(throwAwayTCs)=[] ;

end






function TIMECLUSTERS = removeShortLivedTCs(TIMECLUSTERS, obj_data_dir, minduration, verbose) ;

  %% Accesses main function variables: TIMECLUSTERS, CE
  %% Modifies TIMECLUSTERS

  allTCindices=1:numel(TIMECLUSTERS) ;
  throwAwayTCs=[] ;

  for iiii=[allTCindices]
    thisTCtimeList=[] ;

    if ( numel(TIMECLUSTERS(iiii).objid) < 1 )
      throwAwayTCs=[throwAwayTCs,iiii] ;
    end

    for jjjj=1:numel(TIMECLUSTERS(iiii).objid)
      OBJ = load_obj_data(obj_data_dir,TIMECLUSTERS(iiii).objid(jjjj));

      thisTCtimeList=[thisTCtimeList, OBJ.time];  %...
                      %CE.time(TIMECLUSTERS(iiii).objid(jjjj))] ;
    end

    duration=24*(max(thisTCtimeList)-min(thisTCtimeList)) ;

    if ( duration < minduration-0.01 )
      throwAwayTCs=[throwAwayTCs,iiii] ;
    end
  end

  % Throw unwanted ones away
  throwAwayTCs = unique(throwAwayTCs);
  if verbose
    disp(['Removed the following due to short duration: ',num2str(throwAwayTCs)])
    disp(['(Database also re-ordered.)'])
  end
  TIMECLUSTERS(throwAwayTCs)=[] ;

end


function TIMECLUSTERS = combineCloseProximityTCs_by_area(TIMECLUSTERS, obj_data_dir, maxCombineTimeDiff, OPT, maxDistToConnect, verbose) ;
  
  disp(' --- ') 
  disp(['Allow center jumps if they overlap within: ', num2str(24*maxCombineTimeDiff), ' hours.']) 
  disp(' --- ') 

  nextClusterID = numel(TIMECLUSTERS) + 1;
     
  TC_eliminate_list=[];
  
  imax = numel(TIMECLUSTERS);
  
  for ii=1:imax

    if (verbose)
      disp(['-- Working on ', num2str(ii), ' of ', num2str(imax), '.'])
    end
    
    otherClusters=setdiff(1:numel(TIMECLUSTERS),ii);
    
    %%Test for GG_before
    for ii_before=[otherClusters]

      %% Figure out how much time between the beginning of ii and end of ii_before.
      OBJ = load_obj_data(obj_data_dir,TIMECLUSTERS(ii).objid(1));
      OBJ_before = load_obj_data(obj_data_dir,TIMECLUSTERS(ii_before).objid(end));
      
      jumpTime=OBJ.time - OBJ_before.time;

      %% If it's not within the maxCombineTimeDiff, no need to consider it any further.
      if (jumpTime > 0.0 & jumpTime < maxCombineTimeDiff+0.01)

        dLON1=OBJ.lon - OBJ_before.lon; %CE.lon(TIMECLUSTERS(ii).objid(1)) - ...
              %CE.lon(TIMECLUSTERS(ii_before).objid(end));
	
        dLAT1=OBJ.lat - OBJ_before.lat; %CE.lat(TIMECLUSTERS(ii).objid(1)) - ...
              %CE.lat(TIMECLUSTERS(ii_before).objid(end));
	
        if ( sqrt(dLON1.^2 + dLAT1.^2 ) > maxDistToConnect )
	  continue
        end

	%% If it is within the time window, I need to calculate the overlap.
	X1=OBJ.pixels.x; %CE.pixels(TIMECLUSTERS(ii).objid(1)).x ;
	Y1=OBJ.pixels.y; %CE.pixels(TIMECLUSTERS(ii).objid(1)).y ;
	
	X2=OBJ_before.pixels.x; %CE.pixels(TIMECLUSTERS(ii_before).objid(end)).x ;
	Y2=OBJ_before.pixels.y; %CE.pixels(TIMECLUSTERS(ii_before).objid(end)).y ;
	
	hits=[] ;

	A=[X1,Y1] ;
	B=[X2,Y2] ;
	
	HITS=intersect(A,B,'rows') ;
	nHits=size(HITS,1);
	
	if (nHits/numel(X1) > OPT.TRACKING_MINIMUM_OVERLAP_FRAC | ...
            nHits/numel(X2) > OPT.TRACKING_MINIMUM_OVERLAP_FRAC | ...
            nHits > OPT.TRACKING_MINIMUM_OVERLAP_POINTS)
							  
	  if verbose
	    disp(['Combined LPTs: ',num2str(ii_before),' and ',...
		  num2str(ii),' in to ',num2str(nextClusterID), ...
		  ' (',num2str(24*jumpTime),' h, ',...
		  num2str(sqrt(dLON1.^2 + dLAT1.^2 )),' deg. jump).'])
	  end
	  
	  TIMECLUSTERS = startNewTimeCluster(TIMECLUSTERS, nextClusterID) ;
          TIMECLUSTERS(nextClusterID).objid=unique(sort([TIMECLUSTERS(ii).objid,...
							TIMECLUSTERS(ii_before).objid]));
          nextClusterID = nextClusterID + 1 ;
	  
          TC_eliminate_list=[TC_eliminate_list, ii];
          TC_eliminate_list=[TC_eliminate_list, ii_before];
	  
	end
	
      end
    end
  end
  
  TC_eliminate_list = unique(TC_eliminate_list);
  
  if verbose
    disp(['Removing the following due to having center jumps with other tracks: ',num2str(TC_eliminate_list)])
  end
  
  TIMECLUSTERS(TC_eliminate_list)=[];
  
end




function combineCloseProximityTCs_by_centroid(maxCombineDist,maxCombineTimeDiff) ;

  TC_eliminate_list=[];

  for ii=1:numel(TIMECLUSTERS)

    otherClusters=setdiff(1:numel(TIMECLUSTERS),ii);

    %Test for GG_before
    for ii_before=[otherClusters]

      jumpTime=CE.time(TIMECLUSTERS(ii).objid(1)) - ...
               CE.time(TIMECLUSTERS(ii_before).objid(end));

      jumpLon=CE.lon(TIMECLUSTERS(ii).objid(1)) - ...
              CE.lon(TIMECLUSTERS(ii_before).objid(end));

      jumpLat=CE.lat(TIMECLUSTERS(ii).objid(1)) - ...
              CE.lat(TIMECLUSTERS(ii_before).objid(end));

      jumpDist = sqrt(jumpLon^2 + jumpLat^2);

      if ( abs(jumpTime) < maxCombineTimeDiff+0.01 & ...
           abs(jumpDist) < maxCombineDist+0.01 )

        TIMECLUSTERS(ii).objid=unique(sort([TIMECLUSTERS(ii).objid,...
                            TIMECLUSTERS(ii_before).objid]));

        disp(['Combined LPTs: ',num2str(ii_before),' in to ',...
              num2str(ii),'.'])

        TC_eliminate_list=[TC_eliminate_list,ii_before];

      end
    end
  end

  TIMECLUSTERS(TC_eliminate_list)=[];

end


function NEWTIMECLUSTERS = put_tracks_in_order(TIMECLUSTERS, verbose);

  % Put TIMECLUSTERS tracks in order by starting objid.

  starting_objid_list = [];
  for ii = 1:numel(TIMECLUSTERS)
    starting_objid_list(ii) = TIMECLUSTERS(ii).objid(1);
  end

  [sort_objid_list, sort_objid_list_indx] = sort(starting_objid_list);

  if (verbose)
    disp(['Sorting by order of starting objid: ', num2str(sort_objid_list_indx)])
  end

  NEWTIMECLUSTERS = TIMECLUSTERS(sort_objid_list_indx);
end %local function

