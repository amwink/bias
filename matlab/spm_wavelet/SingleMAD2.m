function ws = SingleMAD2(wc,L)
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
[n,J] = dyadlength(ws(1,:));

ws1=[];
ws2=[];
ws3=[];

% collect the channels' coefficients
for j=J-1:-1:L 
    tl=length(dyad(j));

    ws1=[ws1 reshape(ws(dyad(j),dyad(j)),1,tl^2)];
    ws2=[ws2 reshape(ws(bdyad(j),dyad(j)),1,tl^2)];
    ws3=[ws3 reshape(ws(dyad(j),bdyad(j)),1,tl^2)];
end;
    
% compute the channels' thresholds
tl=length(ws1);
thr=sqrt(2*log(tl));  

% apply the threshold
scale=median(ws1)/.6745;
if (scale)
  ws1=scale*HardThreshAbs(ws1/scale,thr);
else
  ws1=HardThreshAbs(ws1,thr) ;
end;

scale=median(ws2)/.6745;
if (scale)
  ws2=scale*HardThreshAbs(ws2/scale,thr);
else
  ws2=HardThreshAbs(ws2,thr) ;
end;

scale=median(ws3)/.6745;
if (scale)
  ws3=scale*HardThreshAbs(ws3/scale,thr);
else
  w3=HardThreshAbs(ws3,thr) ;
end;

%restore the 2DFWT layout
for j=J-1:-1:L 
    tl=length(dyad(j));

    ws(dyad(j),dyad(j))=reshape(ws1(1:tl^2),tl,tl);
    ws1(1:tl^2)=[];
    ws(bdyad(j),dyad(j))=reshape(ws2(1:tl^2),tl,tl);
    ws2(1:tl^2)=[];
    ws(dyad(j),bdyad(j))=reshape(ws3(1:tl^2),tl,tl);
    ws3(1:tl^2)=[];
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
    
