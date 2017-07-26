function F = gaussSmoothKernel(nx,ny,stdx,stdy) 
    
%
% F = gaussSmoothKernel(nx,ny,stdx,stdy)
%    -- or --
% F = gaussSmoothKernel(stdev)  
%
% Returns 2-D Gaussian smoothing kernel
% With x and y lengths specified.
% along with standard deviation on in the x and y direction.
% (standard deviations are in terms of gridpoints.)
%
% !!!! nx and ny must be odd, or else 1 will be added to
% them. !!!!
%    
% New simplified option: only specify stdx (which is used in x
% and y directions).  It is also in terms of gridpoints.
%     
% ----------------------------------------------
% 2010.11.11  BK  bkerns@rsmas.miami.edu
% 2012.11.12  BK  added simplified option to only specify the stdx
    
    if (nargin == 1)
        
        stdx=nx ;
        stdy=nx ;
        
        nx = 6*stdx + 1 ;
        ny = 6*stdy + 1 ;
      
    end
      
    if ( mod(nx,2) == 0  )
  
        disp('Warning: nx should be odd! Adding one to it.')
        nx=nx+1
        
    end
  
    if ( mod(ny,2) == 0  )
        
        disp('Warning: ny should be odd! Adding one to it.')
        ny=ny+1
        
    end
    
  
  
  
    F = ones(nx,ny) ;

    halfX=(nx-1)/2 ;
    x = [ -1*halfX : halfX ] ;
  
    halfY=(ny-1)/2 ;
    y = [ -1*halfY : halfY ] ;

    [X,Y]=meshgrid(x,y) ;

    F = exp(    -1.0*(X.^2)/(2*stdx^2) +   -1.0*(Y.^2)/(2*stdy^2)      ) ;

    %% Noramlize
    F = F/sum(sum(F)) ;


