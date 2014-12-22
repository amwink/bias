function ws = SingleMAD2_old(wc,L)
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

lastc=max(bdyad(L));
len=size(wc,1);    

thr=sqrt(2*log(len-lastc));            % VisuThresh threshold

ws=wc(lastc+1:end,:);
ws=[wc(1:lastc,lastc+1:end)' ws];

wss=size(ws);

scale=median(ws(:))/.6745;
if (scale)
  ws=scale*HardThreshAbs(ws(:)/scale,thr);
else
  ws=HardThreshAbs(ws(:),thr) ;
end;

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
    
