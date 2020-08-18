function all_rot = cube_rotations ( cubesize, add_diagonals )
%
% Takes the dimensions of a cube as argument and returns the array
% a which contains the indices to rotate a cube of that dimension.
%
% An example for a square whose indices are 1..9:
%
% a = [ 1 2 3; ...
%       4 5 6; ...
%       7 8 9 ];
%
% For nearest neighbour interpolation, rotating the square by 45
% degrees corresponds to shuffling the indices:
%
%           3
%         2   6
% a45 ~ 1   5   9 
%         4   8 
%           7  
%           
% a45 = [ 2 3 6; ...
%         1 5 9; ...
%         4 7 8 ];
%
% Because this is algebraically correct; a45(a45) = a90, and so on,
% this corresponds to rotating around the Z axis. others are harder
% to visualise. 
%
% This function returns the new positions for all rotions that are
% multiples of 45 degrees around the X, Y and Z axes.
%



% This works for 3x3x3 voxel-cubes (aka Rubik's cube)
if ( cubesize ~= 3 )
  error ( 'Sorry, this function is currently only implemented for cubes of size 3.' );
end



% Rotate "Rubik's cube" over the X axis
X90=[];
for i=cubesize:-1:1
  for t=0:(cubesize-1)
    X90=[X90, i+t*cubesize^2];
  end
  y=i; % Use i as starting value for next loop
  for j=1:(cubesize-1)
    y=i+j*cubesize;
    for t=0:(cubesize-1)
      X90=[X90, y+t*cubesize^2];
    end
  end
end
X180 = X90  ( X90 ); % multiples of 90 degrees
X270 = X180 ( X90 );

% add X rotation 45 degrees + multiples of 90
X45  = [2, 3, 12, 5, 6, 15, 8, 9, 18, 1, 11, 21, 4, 14, 24, 7, 17, 27, 10, 19, 20, 13, 22, 23, 16, 25, 26];
X135 =  X45 ( X90 );
X225 = X135 ( X90 );
X315 = X225 ( X90 );



% Rotate "Rubik's cube" over the Y axis
Y90=[];
for i=0:(cubesize-1)
  y=cubesize^2*(cubesize-1)+1+i*cubesize;
  for t=0:cubesize-1 	%Loop through other cubesizeensions	
    Y90=[Y90 (y-t*cubesize^2):(y-t*cubesize^2+cubesize-1)];
  end
end	
Y180 = Y90  ( Y90 ); % multiples of 90 degrees
Y270 = Y180 ( Y90 );

% add Y rotation 45 degrees + multiples of 90
Y45  = [10, 11, 12, 1, 2, 3, 4, 5, 6, 19, 20, 21, 13, 14, 15, 7, 8, 9, 22, 23, 24, 25, 26, 27, 16, 17, 18];
Y135 = Y45  ( Y90 );
Y225 = Y135 ( Y90 );
Y315 = Y225 ( Y90 );




% Rotate "Rubik's cube" over the Z axis
Z90=[];
for i=cubesize:cubesize^2:cubesize^3
  % loop through other dimensions
  for t=0:cubesize-1
    Z90=[Z90, i+t*cubesize];
  end
  % Loop through other dimensions
  for j=i-1:-1:(i-cubesize+1)
    for t=0:cubesize-1
      Z90=[Z90, j+t*cubesize];
    end
  end
end
Z180=Z90(Z90); % multiples of 90 degrees
Z270=Z180(Z90);

% add Z rotation 45 degrees + multiples of 90
Z45= [2, 3, 6, 1, 5, 9, 4, 7, 8, 11, 12, 15, 10, 14, 18, 13, 16, 17, 20, 21, 24, 19, 23, 27, 22, 25, 26];
Z135= Z45(Z90);
Z225= Z135(Z90);
Z315= Z225(Z90);



%% Reflection over X axis
cube=1:cubesize^3;
cube=reshape(cube, cubesize, cubesize, cubesize);

Xref=[];
for i=1:cubesize
	t=flipud(cube(:, :, i));
	Xref=[Xref, reshape(t, 1, cubesize^2)]; 
end

%% Reflect over Y axis
Yref=[];
for i=1:cubesize
	t=fliplr(cube(:, :, i));
	Yref=[Yref, reshape(t, 1, cubesize^2)]; 
end

%% Reflect over Z axis
Zref=[];
for i=cubesize:-1:1
	Zref=[Zref, reshape(cube(:, :, i), 1, cubesize^2)]; 
end



% inclusion of 45-degree turns is optional
if ( add_diagonals )
  all_rot = [ [1:cubesize^3]', ...
	      Y45', Y90', Y135', Y180', Y225', Y270', Y315', ...
	      X45', X90', X135', X180', X225', X270', X315', ...
	      Z45', Z90', Z135', Z180', Z225', Z270', Z315', ...
	      Yref', Xref', Zref' ];
else
  all_rot = [ [1:cubesize^3]', ...
	      Y90', Y180', Y270', ...
	      X90', X180', X270', ...
	      Z90', Z180', Z270', ...
	      Yref', Xref', Zref' ];
end



return
