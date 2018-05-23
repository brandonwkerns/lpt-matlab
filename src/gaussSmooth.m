function dataSmooth = gaussSmooth(data,nx,ny,stdx,stdy) 

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


if ( nargin < 3 ) 
    
    kernel = gaussSmoothKernel(nx) ;
    
else

    if ( mod(nx,2) == 0 ) , nx=nx+1 ; end
    if ( mod(ny,2) == 0 ) , ny=ny+1 ; end


    kernel = gaussSmoothKernel(nx,ny,stdx,stdy) ;

end


dataSmooth=conv2(data, kernel, 'same') ;





