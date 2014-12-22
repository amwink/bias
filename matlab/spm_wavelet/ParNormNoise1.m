function sigmas = ParNormNoise1(wc)
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
% Parallel version of NormNoise1
%
% [ws] = ParNormNoise1(wc)
%        where wc = [timepoints x voxels]

% vectorised form of NormNoise1
i=1:size(wc,2);

sigmas(i)=median(wc(end/2+1:end,:));
sigmas=sigmas/.6745; 

% don't rescale zeros
sigmas(find(~sigmas))=1;

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%





