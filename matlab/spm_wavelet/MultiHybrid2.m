function ws = MultiHybrid2(wc,L)
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
% MultiHybrid -- Apply Shrinkage to Wavelet Coefficients
%  Usage 
%    ws = MultiHybrid(wc,L)
%  Inputs 
%    wc    Wavelet Transform of noisy sequence with N(0,1) noise
%    L     low-frequency cutoff for Wavelet Transform
%  Outputs 
%    ws    result of applying HybridThresh to each dyadic block
%

ws = wc;
[n,J] = dyadlength(ws(1,:));

for j=J-1:-1:L;
    l=length(dyad(j));

    tws=ws(dyad(j),dyad(j));
    ws(dyad(j),dyad(j))=reshape(HybridThresh(tws(:)'),l,l);
    tws=ws(bdyad(j),dyad(j));
    ws(bdyad(j),dyad(j))=reshape(HybridThresh(tws(:)'),l,l);
    tws=ws(dyad(j),bdyad(j));
    ws(dyad(j),bdyad(j))=reshape(HybridThresh(tws(:)'),l,l);
end;

return
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
    
