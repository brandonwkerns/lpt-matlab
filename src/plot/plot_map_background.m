function mapplot = plot_map_background(area, s, conus,lw)

%
%   mapplot = plot_map_background(area, s, 0, lw)
%   mapplot = plot_map_background(area, s, conus)
%   mapplot = plot_map_background(area, s)
%   mapplot = plot_map_background(area)
%   mapplot = plot_map_background
%
%Plot a map background based on the range array specified.
%
%
%INPUT:
%area: of form [xmin,xmax,ymin,ymax]. defaults to [-100,-10,0,40]
%s: optional, if specified, feeds the S character string into plot. 
%     defaults to black ('k').
%
%OUTPUT:
%mapplot: the output from "plot" command
%
%---------------------------------------------
%Brandon Kerns, RSMAS
%bkerns@orca.rsmas.miami.edu
%8.8.2008
%---------------------------------------------
%CHANGES:
%
%BK 9.4.2008: Chenged default range to [-99.0,-10,0,40] so the labeling on 
%   the lower left side of the plot is not cluttered.
%
%BK 4.9.2009: If 
%
%

if nargin < 1 
  area=[-99.9,-10,0,40] ;
end	

if size(area,2) < 4 
  area=[-99.0,-10,0,40] ;
end

if nargin < 2 
  s = 'k' ;
end	

if nargin < 3 
  conus=0 ;
end	

if nargin < 4 
  lw=0.5 ;
end	


if conus == 0

%%call up world file
world = load_world ;

mapplot = plot(world.lon,world.lat,s,'linewidth',lw) ;
hold on
mapplot = plot(world.lon + 360.0,world.lat,s,'linewidth',lw) ;

axis(area) ;
grid on;

else

coast=load_coast ;
borders=load_borders ;

mapplot = plot(coast.lon, coast.lat, s) ;
hold on
mapplot = plot(coast.lon + 360.0, coast.lat, s) ;

axis(area) ;
hold on
plot(borders.lon, borders.lat, s) ;
plot(borders.lon + 360.0, borders.lat, s) ;

hold off


end

