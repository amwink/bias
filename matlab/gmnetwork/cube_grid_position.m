function [ max_nonzero_cubes, bestoffset, greymatter ] = cube_grid_position( template, gm_volume, cubesize )
%
% Find optimal position of the cube grid in terms of non-zeros
%
% The grid with cubes of cubesize^3 voxels each can be placed at
% cubesize^3 different positions (starting at each of (1,1,1),
% (1,1,2), etc). Each position leads to a number of cubes that have
% no intensity variation (all intensities the same). This function
% returns the position that has the fewest zero-variance cubes.
%
% max_nonzero_cubes = number of cubes with nonzero variance
% bestoffset        = list of 1D cube offsets at the optimal position
% greymatter        = grey matter intensities sampled at these offsets
%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%



% generate the 1D offsets of the cubes for different cube grid positions
offsets   = cube_offsets ( template, cubesize );
gridstats = [];

% this is currently the best we have
max_nonzero_cubes=0;

% offsets{x, y, z} contains (cind = offsets{x, y, z};   clen = cubesize^3)
% cind (1:clen): first cube, cind (clen+1:2*clen): second cube, etc.
cubelen = cubesize * cubesize * cubesize;



% size offsets: number of starting positions in x, y and z (3D)
for xpos = 1:size(offsets, 1)
  for ypos = 1:size(offsets, 2)
    for zpos = 1:size(offsets, 3)
      
      % get the indices, resize to [cubesize #cubes]
      cubes = gm_volume ( offsets { xpos, ypos, zpos } );
      cubes = reshape ( cubes, [ cubelen, length(cubes) / cubelen ] );
      
      % each column is now a cube, so we can get mean and std, then
      % compute var after discarding NaNs and check >0
      cubes  ( ~cubes ) = nan;
      cvar = ( nanvar ( cubes ) > 0 );
      % select only those cubes without 0s and with nonzero variance
      cvar ( isnan ( cvar ) )   = 0;
      cvar                      = find   ( cvar );
      numgoodcubes              = length ( cvar ); 
      fprintf ( '%d\r', numgoodcubes );
      
      % if this is better than what we have -> replace
      % reshape the offsets to the same sizes as cubes
      % (only retain offsets of cubes we actually use)
      if ( numgoodcubes > max_nonzero_cubes ) 
	max_nonzero_cubes      = numgoodcubes;
	greymatter             = cubes ( :, cvar );
	bestoffset             = reshape ( offsets { xpos, ypos, zpos }, ...
		                 [ cubelen, length( offsets { xpos, ypos, zpos } ) / cubelen ] );
	bestoffset             = bestoffset ( :, cvar );
      end
      
    end
  end
end



return;
