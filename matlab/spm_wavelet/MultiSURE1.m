function ws = MultiSURE1(wc,type,L)
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
% MultiSURE -- Apply Shrinkage to Wavelet Coefficients
%  Usage 
%    ws = MultiSURE(wc,L)
%  Inputs 
%    wc    Wavelet Transform of noisy sequence with N(0,1) noise
%    L     low-frequency cutoff for Wavelet Transform
%  Outputs 
%    ws    result of applying SUREThresh to each dyadic block
%


ws=wc;
[n,J] = dyadlength(ws);

for j=J-1:-1:L
    thr=length(dyad(j));
    thr=sqrt(2*log(thr));            % VisuShrink threshold
   
    thr=min([ValSUREThresh(ws(dyad(j))') thr]);
    switch(upper(type));
     case 'HARD'
      ws(dyad(j)) = HardThreshAbs(ws(dyad(j)),thr) ;
     case 'SOFT'
      ws(dyad(j)) = SoftThreshAbs(ws(dyad(j)),thr) ;
    end;
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
    








