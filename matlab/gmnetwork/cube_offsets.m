function offsets = cube_offsets ( template, cubesize ) 
%
% 1-dimensional offsets for 3x3x3 cubes in a 91x109x91 volume
%
% The cubes can start from x, y, z positions with x/y/z being {1, 2, 3}
% which leads to different 1D indexings. For each x, y, z position, 
% offsets{x, y, z} contains the indices for all cubes subsequently:
%
%  reshape (offsets{x, y, z}, [ n*n*n length(offsets{x, y, z})/n/n/n]) 
%  contains the offsets of every cube in separate columns, e.g.
%
%  >> test = offsets { 1, 1, 1 };
%  >> test = reshape ( test, [ n*n*n length(test)/n/n/n ] );
%  >> test ( :, 1:10 )  
%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%


% standard MNI image with 2x2x2 mm^3 voxels
Timg = load_untouch_header_only ( template );
dimx = Timg.dime.dim(2);
dimy = Timg.dime.dim(3);
dimz = Timg.dime.dim(4);

% for 3x3x3 voxel: cubesize = 3
n = cubesize;

% store the offsets in an easy to retrieve cell array
offsets = cell ( n, n, n );

% the offsets of one cube are spaced a 
% fixed amount of places apart that is 
% irrespective of 'cube grid' position
ind_  = (1:n)';                  % cube indices in the x direction'
ind_2 = ind_ + (0:(n-1)) * dimx;       
ind_2 = ind_2(:);                % cube indices in the y direction
ind_3 = ind_2 + (0:(n-1)) * ( dimx * dimy );
ind_3 = ind_3(:);                % cube indices in the z direction



for startx = 0:(n-1)             
  for starty = 0:(n-1)
    for startz = 0:(n-1)     % copy indices to 'allowed' positions
      
      % fprintf ( ' offsets for %d%d%d ... \r', startx, starty, startz );
      
      % on row: cube indices should not go beyond dimx
      rowx = startx:n:(dimx-n);
      ind_rows = ind_3 + rowx;
      ind_rows = ind_rows(:);
      
      % on column: cube indices should not go beyond dimy
      coly = starty:n:(dimy-n);
      ind_cols = ind_rows + coly * dimx;
      ind_cols = ind_cols(:);
      
      % in slice: cube indices should not go beyond dimz
      sliz = startz:n:(dimz-n);
      ind_slis = ind_cols + sliz * (dimx * dimy);
      ind_slis = ind_slis(:);
      
      % store the offsets for this starting position
      offsets { startx+1, starty+1, startz+1 } = ind_slis;
      
      % show the size of the grid
      % fprintf ( '%d %d %d\n', length(rowx), length(coly), length(sliz) );
      
    end;
  end;
end;



return;
