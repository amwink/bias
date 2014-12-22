function s = NormNoise1(vec_in)
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
% NormNoise -- Normalize signal to noise level 1 for 2D matrix.
%              2D version of NormNoise.  
% 
%  Call:  vec_out = NormNoise(vec_in,qmf,action)
%                   
%                   vec_in: Input vector (1D)
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
lev1wavs = vec_in(floor(end/2)+1:end,:);

if (size(vec_in,1)==1)
  s = median(lev1wavs);
else
  s = median(lev1wavs,1);
end
  
if (s ~= 0)
  s = s/.6745;
else
  s = 1;
end

return








