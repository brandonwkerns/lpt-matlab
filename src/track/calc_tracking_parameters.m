function systems_out = calc_tracking_parameters(systems_in, objects_dir)



  %% Accesses main function variables: TIMECLUSTERS, CE
  %% Modifies TIMECLUSTERS


  %CE = objects;
  TIMECLUSTERS = systems_in;


  allTCindices=1:numel(TIMECLUSTERS) ;

  for iiii=[allTCindices]

    thisTCtimeList=[] ;
    OBJ_collect = []; 
    
    for objid = [TIMECLUSTERS(iiii).objid]
      OBJ = load_obj_data(objects_dir, objid);
      OBJ_collect = [OBJ_collect, OBJ];
      thisTCtimeList=[thisTCtimeList, OBJ.time];      
    end


    %% Grid.
    gridLon=OBJ.grid.lon ;
    gridLat=OBJ.grid.lat ;
    
    gridDX=gridLon(2)-gridLon(1) ;
    gridDY=gridLat(2)-gridLat(1) ;
    
    AREA=OBJ.grid.area ;


    
    thisTCtimeListUniq=unique(thisTCtimeList) ;

    %% Time and duration for time cluster over all.
    TIMECLUSTERS(iiii).time=thisTCtimeListUniq ;
    TIMECLUSTERS(iiii).duration=24*(max(thisTCtimeList)-min(thisTCtimeList)) ;
    TIMECLUSTERS(iiii).nentries=numel(thisTCtimeListUniq) ;

    %% Now for each of the times, gather the relevant objects and get bulk properties.
    for tttt=1:numel(thisTCtimeListUniq) ;

      thisTimeMatch=find(thisTCtimeList == thisTCtimeListUniq(tttt));
      thisTimeOBJID=TIMECLUSTERS(iiii).objid(thisTimeMatch);

      %% The OBJs at each time. TC.obj(tttt).objid, ect.x
      TIMECLUSTERS(iiii).obj(tttt).objid=thisTimeOBJID ;
      TIMECLUSTERS(iiii).obj(tttt).time = [OBJ_collect(thisTimeMatch).time] ;
      TIMECLUSTERS(iiii).obj(tttt).lon  = [OBJ_collect(thisTimeMatch).lon] ;
      TIMECLUSTERS(iiii).obj(tttt).lat  = [OBJ_collect(thisTimeMatch).lat] ;
      TIMECLUSTERS(iiii).obj(tttt).area = [OBJ_collect(thisTimeMatch).area] ;
      TIMECLUSTERS(iiii).obj(tttt).size = sqrt([OBJ_collect(thisTimeMatch).area]) ;
      TIMECLUSTERS(iiii).obj(tttt).effective_radius = sqrt([OBJ_collect(thisTimeMatch).area]/pi) ;
      TIMECLUSTERS(iiii).obj(tttt).volrain = [OBJ_collect(thisTimeMatch).volrain] ;
      %TIMECLUSTERS(iiii).obj(tttt).pixels  = [OBJ_collect(thisTimeMatch).pixels] ;

      %% Combine the OBJs for the net TC time dependent properties.
      %% ( full feature centroid tracks, weighted by area for multiple OBJs at a single time.)
      TIMECLUSTERS(iiii).lon(tttt)=...
      sum(TIMECLUSTERS(iiii).obj(tttt).lon .* TIMECLUSTERS(iiii).obj(tttt).area)/...
      sum(TIMECLUSTERS(iiii).obj(tttt).area);

      TIMECLUSTERS(iiii).lat(tttt)=...
      sum(TIMECLUSTERS(iiii).obj(tttt).lat .* TIMECLUSTERS(iiii).obj(tttt).area)/...
      sum(TIMECLUSTERS(iiii).obj(tttt).area);
      
      %{
      TIMECLUSTERS(iiii).lon_median(tttt)=...
      sum(CE.lon_median(thisTimeCEID).*CE.area(thisTimeCEID))/...
      sum(CE.area(thisTimeCEID));

      TIMECLUSTERS(iiii).lat_median(tttt)=...
      sum(CE.lat_median(thisTimeCEID).*CE.area(thisTimeCEID))/...
      sum(CE.area(thisTimeCEID));
      %}
      
      TIMECLUSTERS(iiii).area(tttt)=sum(TIMECLUSTERS(iiii).obj(tttt).area);
      TIMECLUSTERS(iiii).size(tttt)=sqrt(TIMECLUSTERS(iiii).area(tttt));
      TIMECLUSTERS(iiii).effective_radius(tttt)=...
      sqrt(TIMECLUSTERS(iiii).area(tttt)/pi);

      TIMECLUSTERS(iiii).volrain(tttt)=sum(TIMECLUSTERS(iiii).obj(tttt).volrain);
      TIMECLUSTERS(iiii).n_obj(tttt)=numel(TIMECLUSTERS(iiii).obj(tttt).objid);

      %% The "Max Area" tracks.
      %{
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
      %}
      %% The "Max Volrain" tracks.
      %{
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
      %}
        
    end



    %% The "bulk" cloud cluster properties.
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
