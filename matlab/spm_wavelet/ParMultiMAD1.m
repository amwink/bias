function ws = ParMultiMAD1(wc,L)
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
% MultiMAD -- Apply Shrinkage with level-dependent Noise level estimation
%  Usage 
%    s = MultiMAD(wc,L)
%  Inputs 
%    wc    Wavelet Transform of noisy sequence
%    L     low-resolution cutoff for Wavelet Transform
%  Outputs 
%    ws    result of applying VisuThresh to each wavelet level,
%          after scaling so MAD of coefficienst at each level = .6745 
%

ws=wc;
[n,J] = dyadlength(1:size(ws,1));


for j=J-1:-1:L    
    LengthDyadj=2^j;
    thr=sqrt(2*log(LengthDyadj));

    scale=median(ws(dyad(j),:))/.6745;

    ScaleNz=find(scale);
    LengthNz=length(ScaleNz);
    ScaleZ=find(~scale);
    LengthZ=length(ScaleZ);

    keyboard
    
    if (LengthNz)
      res = (ws(dyad(j),ScaleNz)./(ones(LengthDyadj,1)*scale(ScaleNz)) - thr);
      res = (res + abs(res))/2;
      ws(dyad(j),ScaleNz) = (ones(LengthDyadj,1)*scale(ScaleNz)).*sign(ws(dyad(j),ScaleNz)).*res;
    end
    if (LengthZ)
      ws(dyad(j),ScaleZ) = ws(dyad(j),ScaleZ) .* (ws(dyad(j),ScaleZ) > thr);
    end
end;


%
% Copyright (c) 1993-5.  Jonathan Buckheit, David Donoho and Iain Johnstone
%
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
