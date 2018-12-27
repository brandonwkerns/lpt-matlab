function systems_out = calc_tracking_parameters(systems_in, objects)



  %% Accesses main function variables: TIMECLUSTERS, CE
  %% Modifies TIMECLUSTERS


  CE = objects;
  TIMECLUSTERS = systems_in;

    %% Grid.
    gridLon=CE.grid.lon ;
    gridLat=CE.grid.lat ;

    gridDX=gridLon(2)-gridLon(1) ;
    gridDY=gridLat(2)-gridLat(1) ;

    AREA=CE.grid.area ;

    allTCindices=1:numel(TIMECLUSTERS) ;

    for iiii=[allTCindices]

        thisTCtimeList=[] ;


        for jjjj=1:numel(TIMECLUSTERS(iiii).ceid)

            thisTCtimeList=[thisTCtimeList,...
                            CE.time(TIMECLUSTERS(iiii).ceid(jjjj))] ;


        end

        thisTCtimeListUniq=unique(thisTCtimeList) ;

        % Time and duration for time cluster over all.
        TIMECLUSTERS(iiii).time=thisTCtimeListUniq ;
        TIMECLUSTERS(iiii).duration=24*(max(thisTCtimeList)-min(thisTCtimeList)) ;
        TIMECLUSTERS(iiii).nentries=numel(thisTCtimeListUniq) ;

        for tttt=1:numel(thisTCtimeListUniq) ;

            thisTimeMatch=find(thisTCtimeList == thisTCtimeListUniq(tttt));
            thisTimeCEID=TIMECLUSTERS(iiii).ceid(thisTimeMatch);

            % The CEs at each time. TC.ce(tttt).ceid, ect.x
            TIMECLUSTERS(iiii).ce(tttt).ceid=thisTimeCEID ;
            TIMECLUSTERS(iiii).ce(tttt).time=CE.time(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).lon=CE.lon(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).lat=CE.lat(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).area=CE.area(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).size=sqrt(CE.area(thisTimeCEID)) ;
            TIMECLUSTERS(iiii).ce(tttt).effective_radius=...
                sqrt(CE.area(thisTimeCEID)/pi) ;
            TIMECLUSTERS(iiii).ce(tttt).volrain=CE.volrain(thisTimeCEID) ;
            TIMECLUSTERS(iiii).ce(tttt).pixels=CE.pixels(thisTimeCEID) ;


            % Combine the CEs for the net TC time dependent properties.
            % ( full feature centroid tracks)
            TIMECLUSTERS(iiii).lon(tttt)=...
                sum(CE.lon(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));

            TIMECLUSTERS(iiii).lat(tttt)=...
                sum(CE.lat(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));

            TIMECLUSTERS(iiii).lon_median(tttt)=...
                sum(CE.lon_median(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));
            TIMECLUSTERS(iiii).lat_median(tttt)=...
                sum(CE.lat_median(thisTimeCEID).*CE.area(thisTimeCEID))/...
                sum(CE.area(thisTimeCEID));

            TIMECLUSTERS(iiii).area(tttt)=sum(CE.area(thisTimeCEID));
            TIMECLUSTERS(iiii).size(tttt)=sqrt(TIMECLUSTERS(iiii).area(tttt));
            TIMECLUSTERS(iiii).effective_radius(tttt)=...
                sqrt(TIMECLUSTERS(iiii).area(tttt)/pi);

            TIMECLUSTERS(iiii).volrain(tttt)=sum(CE.volrain(thisTimeCEID));
            TIMECLUSTERS(iiii).nclusters(tttt)=numel(thisTimeMatch);

            % The "Max Area" tracks.
            [MA,MaxAreaI]=max(TIMECLUSTERS(iiii).ce(tttt).area) ;

            TIMECLUSTERS(iiii).max_area_track.lon(tttt)=...
                CE.lon(thisTimeCEID(MaxAreaI));
            TIMECLUSTERS(iiii).max_area_track.lat(tttt)=...
                CE.lat(thisTimeCEID(MaxAreaI));

            TIMECLUSTERS(iiii).max_area_track.area(tttt)=...
                CE.area(thisTimeCEID(MaxAreaI));
            TIMECLUSTERS(iiii).max_are_tracka.size(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxAreaI)));
            TIMECLUSTERS(iiii).max_area_track.effective_radius(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxAreaI))/pi);
            TIMECLUSTERS(iiii).max_area_track.volrain(tttt)=...
                CE.volrain(thisTimeCEID(MaxAreaI));

            TIMECLUSTERS(iiii).max_area_track.lon_max(tttt)=...
                CE.lon_max(thisTimeCEID(MaxAreaI));
            TIMECLUSTERS(iiii).max_area_track.lat_max(tttt)=...
                CE.lat_max(thisTimeCEID(MaxAreaI));

            % The "Max Volrain" tracks.
            [MA,MaxVolrainI]=max(TIMECLUSTERS(iiii).ce(tttt).volrain) ;

            TIMECLUSTERS(iiii).max_volrain_track.lon(tttt)=...
                CE.lon(thisTimeCEID(MaxVolrainI));
            TIMECLUSTERS(iiii).max_volrain_track.lat(tttt)=...
                CE.lat(thisTimeCEID(MaxVolrainI));

            TIMECLUSTERS(iiii).max_volrain_track.area(tttt)=...
                CE.area(thisTimeCEID(MaxVolrainI));
            TIMECLUSTERS(iiii).max_are_tracka.size(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxVolrainI)));
            TIMECLUSTERS(iiii).max_volrain_track.effective_radius(tttt)=...
                sqrt(CE.area(thisTimeCEID(MaxVolrainI))/pi);
            TIMECLUSTERS(iiii).max_volrain_track.volrain(tttt)=...
                CE.volrain(thisTimeCEID(MaxVolrainI));

            TIMECLUSTERS(iiii).max_volrain_track.lon_max(tttt)=...
                CE.lon_max(thisTimeCEID(MaxVolrainI));
            TIMECLUSTERS(iiii).max_volrain_track.lat_max(tttt)=...
                CE.lat_max(thisTimeCEID(MaxVolrainI));




        end



        % The "bulk" cloud cluster properties.
        [TIMECLUSTERS(iiii).maxarea,TIMECLUSTERS(iiii).maxarea_indx]=...
            max(TIMECLUSTERS(iiii).area);

        TIMECLUSTERS(iiii).maxsize=...
            max(TIMECLUSTERS(iiii).size);

        TIMECLUSTERS(iiii).max_effective_radius=...
            max(TIMECLUSTERS(iiii).effective_radius);

        TIMECLUSTERS(iiii).maxvolrain=...
            max(TIMECLUSTERS(iiii).volrain);

        [TIMECLUSTERS(iiii).minlon,TIMECLUSTERS(iiii).minlon_indx]=...
            min(TIMECLUSTERS(iiii).lon);

        [TIMECLUSTERS(iiii).maxlon,TIMECLUSTERS(iiii).maxlon_indx]=...
            max(TIMECLUSTERS(iiii).lon);

        [TIMECLUSTERS(iiii).minlat,TIMECLUSTERS(iiii).minlat_indx]=...
            min(TIMECLUSTERS(iiii).lat);

        [TIMECLUSTERS(iiii).maxlat,TIMECLUSTERS(iiii).maxlat_indx]=...
            max(TIMECLUSTERS(iiii).lat);

        [TIMECLUSTERS(iiii).minabslat,TIMECLUSTERS(iiii).minabslat_indx]=...
            min(abs(TIMECLUSTERS(iiii).lat));

        [TIMECLUSTERS(iiii).maxabslat,TIMECLUSTERS(iiii).maxabslat_indx]=...
            max(abs(TIMECLUSTERS(iiii).lat));

        %% Zonal and meridional propagation speed.
        if numel(TIMECLUSTERS(iiii).lon) > 1
          [FIT0,S,MU]=polyfit(TIMECLUSTERS(iiii).time,TIMECLUSTERS(iiii).lon,1);
          TIMECLUSTERS(iiii).zonal_propagation_speed = ...
              (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ;

          [FIT0,S,MU]=polyfit(TIMECLUSTERS(iiii).time,TIMECLUSTERS(iiii).lat,1) ;
          TIMECLUSTERS(iiii).meridional_propagation_speed = ...
              (FIT0(1)) * 111000.0/(24.0*3600.0*MU(2)) ;
        else
          TIMECLUSTERS(iiii).zonal_propagation_speed = NaN;
          TIMECLUSTERS(iiii).meridional_propagation_speed = NaN;
        end
    end
  
    systems_out = TIMECLUSTERS;
