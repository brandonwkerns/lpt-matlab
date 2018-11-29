clear all
close all


for year1=[2018]


year2=year1+1 ;

yyyy1=num2str(year1) ;
yyyy2=num2str(year2) ;

y1_y2=[yyyy1,'_',yyyy2] ;

outdir=['../data/trmm/interim/timelon'];
outfile=[outdir,'/rain_hov_',y1_y2,'_15deg_3day_full_year.mat'] ;
mainDir='../data/trmm/interim/accumulated/72h' ;

lat1=-14.875 ;
lat2=14.875 ;

dateNumNow=now() ;
%dateNum30=datenum(year1,6,1,0,0,0):3:datenum(year2,6,1,0,0,0) ;
dateNum30=datenum(year1,6,1,0,0,0):3:datenum(year1,11,27,21,0,0) ;

%% Set up arrays
lon=0.875:.25:359.875 ;
time=dateNum30 ;
rain=NaN*ones(length(time),length(lon)) ;


lon1=lon(1) ;
lon2=lon(end) ;

for ii=1:length(dateNum30)
    
    
    [year,month,day,hour]=datevec(dateNum30(ii)) ;
    yyyy=num2str(year) ;
    mm=sprintf('%02d',month) ;
    dd=sprintf('%02d',day) ;
    
    dateOfYear=date2doy(datenum(year,month,day)) ;
    doy=sprintf('%03d',fix(dateOfYear)) ;
    
        
    
    hh=sprintf('%02d',hour) ;
    disp([yyyy,mm,dd,hh]) 
    
    %% Open 3B42 rain rate
    
    fileRain=[mainDir,'/rain_accumulated_72h_',yyyy,mm,dd,hh,'.nc'] ;
    
    F=[] ;
    if ( exist(fileRain) )
        %F=load(fileRain) ;
        disp(fileRain) ;
        F.precip=ncread(fileRain,'rain')' ;
        
        F.lon=ncread(fileRain,'lon');
	F.lat=ncread(fileRain,'lat');
    else
        disp(['No data found for ',yyyy,mm,dd])
    end
    
    if ( length(F) > 0 ) 
        
        %% Subset the data
        if ( length(F)>0 )
            keepLon=( F.lon > lon1-0.1 & F.lon < lon2 + 0.1 ) ;
            keepLat=( F.lat > lat1-0.1 & F.lat < lat2 + 0.1 ) ;
            
            keepRain=squeeze(F.precip(keepLat,keepLon)) ;
          
        else
            
            keepRain=NaN*ones(length(lon),80) ;
            
        end
 
        for jj=1:length(lon)
            
            thisLine=squeeze(keepRain(:, jj)) ;
            
            keepThisLine=find( isfinite(thisLine) & thisLine > -0.01 ) ;
            
            thisRain(jj)= mean(thisLine(keepThisLine)) ;
            
        end
        
        thisRainsmth=thisRain ;
        rain(ii,:)=thisRainsmth ;
        
    end
    
        
end




fout.lon=lon ;
fout.time=time ;
fout.rain=rain ;

eval(['!mkdir -p ',outdir])
eval(['save ', outfile,' -STRUCT fout'])  


end






