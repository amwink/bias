function ws = SingleInvShrink2(wc,L)
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
% InvShrink -- Shrinkage with Exponential Factor applied
%  Usage
%    s = InvShrink(wc,L)
%  Inputs
%    wc      Wavelet Transform of noisy sequence with N(0,1) noise
%    L       low-frequency cutoff for Wavelet Transform
%    sa      noise level at highest resolution level
%    alpha   decay rate of noise level with resolution level
%  Outputs
%    ws      result of applying FixShrink to each wavelet level,
%            with appropriate weighting 
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
thr=sqrt(2*log(tl));      % VisuThresh threshold

scale1=std(ws1);
if (~scale1)
  scale1=1;
end;
scale2=std(ws2);
if (~scale2)
  scale2=1;
end;
scale3=std(ws3);
if (~scale3)
  scale3=1;
end;

% apply the threshold
ws1=scale1*SoftThreshAbs(ws1/scale1,thr);
ws2=scale2*SoftThreshAbs(ws2/scale2,thr);
ws3=scale3*SoftThreshAbs(ws3/scale3,thr);

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
    
