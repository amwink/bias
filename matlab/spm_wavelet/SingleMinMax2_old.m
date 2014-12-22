function ws = SingleMinMax2_old(wc,L)
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
% MinMaxThresh -- Minimax Thresholding
%  Usage
%    x = MinMaxThresh(y)
%  Inputs
%    y   signal upon which to perform thresholding
%  Outputs
%    x   result
%
%  References
%    ``Ideal Spatial Adaptation via Wavelet Shrinkage''
%    by D.L. Donoho and I.M. Johnstone.
%

lamlist = [0 0 0 0 0 1.27 1.47 1.67 1.86 2.05 2.23 2.41 2.6 2.77 2.95 3.13];

lastc=max(bdyad(L));
len=size(wc,1);    

ws=wc(lastc+1:end,:);
ws=[wc(1:lastc,lastc+1:end)' ws];

wss=size(ws);
thr=prod(wss);            
  
l=ceil(log(thr)/log(2));
lam = lamlist(max([l 1]));
ws = HardThreshAbs(ws(:),lam) ;

ws=reshape(ws,wss);
ws=[ [ wc(1:lastc,1:lastc) ws(:,1:lastc)' ]; ...
       ws(:,lastc+1:end)                  ];

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
    
