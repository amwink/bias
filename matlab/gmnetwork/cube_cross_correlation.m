function [ ori_corr, ran_corr ] = cube_cross_correlation ( gm_incubes, gm_random, cubesize, add_diag )
%
% compute the correlations between every pair of cubes
%
% Given the grey matter densities sampled in the cubes, this
% scripts rotates the cube at all possible angles that are
% multiples of 45 (or 90) around each axis, and correlates them
% with all other cubes.
%
% ori_corr: result of the above procedure on the observed GM densities
% ran_corr: the same procedure on grey matter at permuted voxel locations
%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%

% different orderings cofficients for cubes rotated at different angles
angles = cube_rotations ( cubesize, add_diag );

% gm_incubes is size [ #percube #ofcubes ] so we need to visit #ofcubes places
ori_corr = zeros ( size ( gm_incubes, 2 ) );
ran_corr = ori_corr;

% subtract ( columm ) mean from each cube
t_mean = mean ( gm_incubes );
n_sum = gm_incubes - t_mean;

% l2 norm of each column ( = cube ) intensity variation
d_sum = sqrt ( sum ( n_sum.^2 ) );

% do the same for randomised cubes
rt_mean = mean ( gm_random );
rn_sum = gm_random - rt_mean;
rd_sum = sqrt ( sum ( rn_sum.^2 ) );

% tried method 2, but that is not faster
method = 1;

% now loop over de-meaned cubes and random cubes
if ( method == 1 )
  
  % fprintf ('method 1:\n');
  for i = 1:size( gm_incubes, 2 )
    
    % observed and random cubes and std
    cube =  n_sum ( :, i );
    rube = rn_sum ( :, i );
    cstd =  d_sum (    i );  
    rstd = rd_sum (    i );
    
    % all rotations for both
    crot = cube ( angles );
    rrot = rube ( angles );
    
    fprintf ( '  correlating cube %04d\r', i );
    
    % complete upper triangular correlation matrix
    for j = (i+1):size( gm_incubes, 2 )
      
      % target and random target
      tcube =  n_sum ( :, j );
      trube = rn_sum ( :, j );
      tcstd =  d_sum (    j );
      trstd = rd_sum (    j );
      
      % test correlations for all rotations both observed and random
      tr  = ( tcube' * crot ) ./ repmat ( cstd * tcstd, [ 1 size( angles,2) ] );
      rtr = ( trube' * rrot ) ./ repmat ( rstd * trstd, [ 1 size( angles,2) ] );
      
      % find maximum for both
      [  m,  mi ] = max (  tr );
      [ rm, rmi ] = max ( rtr );
      
      % store those
      ori_corr ( i, j ) =  m;
      ran_corr ( i, j ) = rm;
      
    end % for j
    
  end % for i
  
else
  
  fprintf ('method 2:\n');
  for i = 1:size( gm_incubes, 2 )
    
    % observed and random cubes and std
    cube =  n_sum ( :, i );
    rube = rn_sum ( :, i );
    cstd =  d_sum (    i );  
    rstd = rd_sum (    i );
    
    % all rotations for both
    crot = cube ( angles );
    rrot = rube ( angles );
    
    fprintf ( '  correlating cube %04d\r', i );
    
    % complete upper triangular correlation matrix - use cubes higher than i
    afteri = (i+1):size(gm_incubes, 2);
    nrot   = size   ( angles, 2 );
    clen   = length ( cube      );
    nleft  = length ( afteri    );
    
    % first the obesrved cubes
    % size [ clen nleft ]
    remain_cubes = n_sum ( :, afteri );
    % put all rotations in 1 column, add copies for all remaining cubes
    % size [ clen*nrot nleft ]
    allrotat_cub = crot(:) * ones ( 1, nleft );
    % reshape to be -- the size of the cube * number of all rotated remaining cubes  
    %               -- column:cube, row: <all cubes> per rotation
    %
    % this requires some juggling:
    %
    % (example with nrot=2; clen=3; nleft=4)
    % >> a = [ repmat ( [ 1:4 ], [ 3 1 ] ); repmat ( [ 11:14 ], [ 3 1 ] ) ];
    %    a =      1     2     3     4
    %             1     2     3     4
    %             1     2     3     4
    %            11    12    13    14
    %            11    12    13    14
    %            11    12    13    14
    % >> a = reshape ( permute ( reshape ( a', [ 4 3 2 ] ), [ 2 1 3 ] ), [ 3 4*2 ] )
    %    a =      1     2     3     4    11    12    13    14
    %             1     2     3     4    11    12    13    14
    %             1     2     3     4    11    12    13    14
    allrotat_cub = reshape ( permute ( reshape ( allrotat_cub', [ nleft clen nrot ] ), [ 2 1 3 ] ), [ clen ( nrot * nleft ) ] );
    % times the rotations of the current cube by enough copies of the remaining cubes and sum columns
    allrotat_cub = sum ( allrotat_cub .* repmat ( remain_cubes, [ 1 nrot ] ) );
    % now divide all these sums by the standard deviations
    allrotat_cub = allrotat_cub ./ repmat ( d_sum (i) * d_sum ( afteri ), [ 1 size(crot,2) ] );
    % reshape and take the maximum, size should be [ #angles #othercubes ],
    % after taking the maximum [ 1 #othercubes ] 1 line of the upper triangular matrix
    allrotat_cub = reshape ( allrotat_cub', [ size( remain_cubes, 2 ) size(crot,2) ] )';
    allrotat_cub = max ( allrotat_cub );
    ori_corr ( i, afteri ) = allrotat_cub;
    
  end % for i
  toc
  
end % if method



fprintf ( '\n  done.\n' );



return
