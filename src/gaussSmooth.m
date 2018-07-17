function dataSmooth = gaussSmooth(data, stdx, stdy, nx, ny, ghost_points, cyclical)

% dataSmooth = gaussSmooth(data,nx,ny,stdx,stdy)
%    -- OR --
% dataSmooth = gaussSmooth(data,stdev)
%
% Data is a 2-D array
% The Gaussian kernel is defined by: nx,ny,stdx,stdy
%  using the function gaussSmoothKernel
%    nx = full width in x direction
%    ny = full width in y direction
%    stdx = standard deviation in x direction (gridpoints)
%    stdy = standard deviation in y direction (gridpoints)
%
%  NOTE:  Function gaussSmooth requires that nx and ny be odd.
%         If you give even numbers, 1 is added to make them odd.
%
%  ALSO: dataSmooth has the same dimensions as data. The way this works is
%        Matlab's built-in 'conv2' is called with the SHAPE option set to 'same'
%
% UPDATE 2012.11.12 you can call this in the simplified way:
%   dataSmooth = gaussSmooth(data,stdev)
%  and it will use stdev in both x and y direction, extending the
%  kernel out to 3 stdev's.

%% Manage optional args with default values.
[J,I] = size(data);
if ( nargin < 5 )
  nx = round(6*stdx + 1);
  ny = round(6*stdy + 1);
end

nx = round(nx);
ny = round(ny);
stdx = round(stdx);
stdy = round(stdy);

if ( mod(nx,2) == 0 ) , nx=nx+1 ; end
if ( mod(ny,2) == 0 ) , ny=ny+1 ; end

kernel = gaussSmoothKernel(nx,ny,stdx,stdy) ;

if ( nargin < 6 )
  ghost_points = false;
end

if ( nargin < 7 )
  cyclical = false;
end

if ghost_points == true
  if cyclical == true
    data_with_ghost_points = [data(:,I-nx:I), data, data(:,1:nx)];
  else
    data_with_ghost_points = [fliplr(data(:,1:nx)), data, fliplr(data(:,I-nx:I))];
  end
  dataSmooth0=conv2(data_with_ghost_points, kernel, 'same') ;
  dataSmooth = dataSmooth0(:,nx+1:end-(nx+1));
else
  dataSmooth=conv2(data, kernel, 'same') ;
end
