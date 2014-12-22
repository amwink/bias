function [mat_out,s] = NormNoise2(mat_in)
%
% This file is part of the software described in 
%
% Alle Meije Wink and Jos B. T. M. Roerdink (2004) 
% ``Denoising functional MR images: 
%   a comparison of wavelet denoising and Gaussian smoothing''
% IEEE Trans. Med. Im. 23 (3), pp. 374-387
%
% please refer to this paper if you use this software
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NormNoise -- Normalize signal to noise level 1 for 2D matrix.
%              2D version of NormNoise.  
% 
%  Call:  mat_out = NormNoise(mat_in,qmf,action)
%                   
%                   mat_in: Input matrix (2D)
%                   an    : Analysis filters (in the freq. domain)
%                   sy    : Synthesis filters (in the freq. domain)
%                   lev   : Current level of decomposition
%                           (normalising is done on level 1)
%                
%                   An average scaling for the entire matrix is
%                   calculated by taking the mean value for the scaling
%                   over all rows (not columns!) of the matrix.
%
%    This is required pre-processing to use any of the DeNoising
%    tools on naturally-occurring data.
%
%  See Also: NormNoise, WaveShrink, CPDeNoise, WPDeNoise
%

% fetch level 1 wavelet coefficients
lev1wavs = [ mat_in(1:end/2,end/2+1:end) ...
             mat_in(end/2+1:end,1:end/2) ...
             mat_in(end/2+1:end,end/2+1:end) ];
lev1wavs = lev1wavs(:);

s=median(lev1wavs);

if ((s ~= 1) & (s ~= 0))
  mat_out =  .6745 .* mat_in ./s;
else
  mat_out = mat_in;
end

s = s/.6745;


return







