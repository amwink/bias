function y=bdyad(x);
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
%  bdyad -- Indexes before j-th dyad of 1-d wavelet xform
%   Usage
%     ix = bdyad(j);
%   Inputs
%     j    integer
%   Outputs
%     ix    list of all indices of wavelet coeffts after j-th level

y=dyad(x)-2^x;

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%
