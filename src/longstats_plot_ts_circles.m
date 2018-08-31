function longstats_plot_ts_circles(fin, daterange, lonrange,latrange,color,circlefactor)

%
%  longstats_plot_ts_bw(fin, daterange, lonrange)
%
%  fin = structure from longstats2struct or longstats)extract
%        fin = longstats2struct(longstatsfile) ;
%

dateNumMin = daterange(1) ;
dateNumMax = daterange(2) ;

if ( numel(fin) < 1 )
    fin.lon=NaN ;
    fin.lat=NaN ;
    fin.date=NaN ;
    fin.size=NaN ;
end
    
fDateNum = fin.date ;
fSize=fin.size ;

if nargin < 3
    minlon = min(min(fin.lon(fin.lon > -200.0))) ;
    maxlon = max(max(fin.lon)) ;
else
    minlon = lonrange(1) ;
    maxlon = lonrange(2) ;
end


if nargin < 4
    latrange=[-90,90] ;
end



if nargin < 5
    color='k' ;
end

if nargin < 6
    circlefactor=1.5 ;
end

if length(latrange) < 2  % Putting [] gives default.
    latrange=[-90,90] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis([ minlon,maxlon,dateNumMin,dateNumMax ])
hold on ;

if (  numel(fin) > 0 )

for ii = 1:length(fin.nentries)

    if (fin.duration(ii) >= 0 & fin.lat(ii,1) > latrange(1) & ...
        fin.lat(ii,1) < latrange(2) )

        %
        % Plot the circles.
        %

        for jj=1:fin.nentries(ii)

            MS=circlefactor*sqrt(fin.area(ii,jj)) ;
            cx=fin.lon(ii,jj) ;
            cy=fin.date(ii,jj) ;
            
            %if ( fin.lat(ii,jj) > -8 & fin.lat(ii,jj) < 8.0 )
            h1 = plot(cx,cy,'o','color',color,'markersize',MS) ;
                %else
                %h1 = plot(cx,cy,'o','color',[0.6,0.6,0.6],'markersize',MS) ; 
                %end
            
        end  
        
        %
        % Plot the lines.
        %
 
        
        plot(fin.lon(ii,1:fin.nentries(ii)), fin.date(ii,1:fin.nentries(ii)), 'linewidth',1,'color',color);
        
        
    end
    
    
end


end


darr=datevec(dateNumMin:1:dateNumMax) ;
months=darr(:,2)' ; days=darr(:,3)' ;

ylabs=[] ;
for iii=1:length(months)
	ylabs=strvcat(ylabs,[num2str(months(iii)),'/',num2str(days(iii))]) ;	
end


set(gca,'ytick', dateNumMin:1:dateNumMax,...
    'yticklabel',  ylabs) ;
grid on ;



