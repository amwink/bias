function y = dsamp(x,dsf)
%
% This file is part of the software described in 
%
% Alle Meije Wink and Jos B. T. M. Roerdink (2004) 
% ``Denoising functional MR images: 
%   a comparison of wavelet denoising and Gaussian smoothing''
% IEEE Trans. Med. Im. 23 (3), pp. 374-387
%
% please refer to that paper if you use this software
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dsamp.m
%
% FORMAT y = DSAMP(x,dsf) 
%
% down sampling of array x by factor dsf.
%
%     x    = input signal (row vector)
%     dsf  = down-sampling factor
%
% Alle Meije Wink 2003/7/28

nrows = size(x,1);
ncols = size(x,2);

xmin = min(nrows,ncols);
xmax = max(nrows,ncols);

if xmin==1,	% 1D signal
  lox = xmax;
  if floor(lox/dsf)~=lox/dsf,
    error('Evenpart: signal length must be multiple of downsampling factor.');
  end
  list = 1:dsf:lox;
  y    = x(list);

else		% 2D signal
  if floor(nrows/dsf)~=nrows/dsf & floor(ncols/dsf)~=ncols/dsf,
    error('Evenpart: signal sizes must be multiples of downsampling factor.');
  end
  list1 = 1:dsf:nrows;
  list2 = 1:dsf:ncols;
  y     = x(list1,list2);

end

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%
