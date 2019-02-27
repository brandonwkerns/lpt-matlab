function lpt_systems_output_ascii(systems, asciifile)

  eval(['!rm -f ',asciifile])


%  CE = objects;
  TIMECLUSTERS = systems;

  FMT='        %4d%02d%02d%02d %7d %10.2f %10.2f %1d\n' ;


				% Ascii output
  %fileout=['LONGSTATS_lpt_',ymd0_ymd9] ;
  disp(asciifile)
  fid=fopen(asciifile,'w') ;

  for ii=1:numel(TIMECLUSTERS)

    fprintf(fid,'Cl%i:\n',ii) ;

    for jj = 1:numel(TIMECLUSTERS(ii).time)
      
      [y,m,d,h]=datevec(TIMECLUSTERS(ii).time(jj)) ;
      
      fprintf(fid,FMT,y,m,d,h,...
              round(TIMECLUSTERS(ii).area(jj)),...
              TIMECLUSTERS(ii).lat(jj),...
              TIMECLUSTERS(ii).lon(jj),...
              TIMECLUSTERS(ii).n_obj(jj)) ;
      
    end
    
  end
  
  fclose(fid) ;

  
