function ws = SingleSURE2_old(wc,type,L)
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

% structure
%
%   C[d D D]----
%   d d D D     | matrix transverse
%   D D D D     |
%   D D D D     |
% ^-------------
%
% is changed into
%
% C     and
%
% d d d D D
% D d d D D
% D D D D D  
% 
% one threshold is applied to the detail block 

lastc=max(bdyad(L));
len=size(wc,1);    

thr=sqrt(2*log(len-lastc));            % VisuThresh threshold

ws=wc(lastc+1:end,:);
ws=[wc(1:lastc,lastc+1:end)' ws];

wss=size(ws);
thr=min([ValSUREThresh(ws(:)) thr]);

switch(upper(type));
 case 'HARD'
  ws = HardThreshAbs(ws(:),thr) ;
 case 'SOFT'
  ws = SoftThreshAbs(ws(:),thr) ;
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
    


