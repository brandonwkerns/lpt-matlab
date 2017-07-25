clear all
close all

% This script reads in the interim files from ../data/interim/gridded_rain_rates
% and creates accumulated rain files in ../data/interim/accumulate_rain
INTERIM_DATA_DIR_IN = '../data/interim/gridded_rain_rates';
INTERIM_DATA_DIR_OUT = '../data/interim/accumulate_rain'

file_list = dir([INTERIM_DATA_DIR_IN,'/*.mat'])

for ff=1:length(file_list)

    %% Open 3B42 rain rate
    fileIn=[mainDir,'/trmm_global_rainfall/',yyyy,'/',mm,'/',...
              yyyy,mm,dd, ...
              '/3B42.',yyyy,mm,dd,'.',hh,'.7.HDF'] ;



    RAINCOLLECT=[] ;


    kkk=0 ;

    for hour_rel=-72:3:0

        kkk=kkk+1 ;

        dn0=dateNum30(ii) + hour_rel/24 ;
        [yyyy,mm,dd,hh]=datevecstr(dn0) ;

        disp([yyyy,mm,dd,hh])


        fileRain=[mainDir,'/trmm_global_rainfall/',yyyy,'/',mm,'/',...
                  yyyy,mm,dd, ...
                  '/3B42.',yyyy,mm,dd,'.',hh,'.7.HDF'] ;

        if ( ~ exist(fileRain) )

            fileRain=[mainDir,'/trmm_global_rainfall/',yyyy,'/',mm,'/',...
                      yyyy,mm,dd, ...
                      '/3B42.',yyyy,mm,dd,'.',hh,'.7A.HDF'] ;

        end

        F=[] ;
        if ( exist(fileRain) )
            F=read3b42(fileRain) ;
            disp(fileRain) ;
        else

            disp(['No TMI data found for ',yyyy,mm,dd])

        end


        if ( length(F) > 0 )


            THISRAIN=F.precip ;
            THISRAIN(THISRAIN < -0.01) = NaN ;

            RAINCOLLECT(:,:,kkk)=THISRAIN ;


        end


    end


    RAINAVG=nanmean(RAINCOLLECT,3)*24.0 ;


    fout.lon=F.lon ;
    fout.lat=F.lat ;
    fout.time=dateNum30(ii) ;
    fout.rain=RAINAVG ;
    fout.info.smoothing='NO SMOOTHING' ;
    fout.info.accumulation='3 days prior' ;
    fout.info.units='mm/day' ;

    eval(['!mkdir -p ',outDir]) ;
    eval(['save ',outDir,'/',outFile,' -struct fout'])

end



end
