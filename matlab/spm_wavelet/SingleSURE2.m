function ws = SingleSURE2(wc,type,L)
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

thr1=min([ValSUREThresh(ws1) thr]);
thr2=min([ValSUREThresh(ws2) thr]);
thr3=min([ValSUREThresh(ws3) thr]);

% apply the threshold
switch(upper(type));
 case 'HARD'
  ws1 = HardThreshAbs(ws1,thr1) ;
  ws2 = HardThreshAbs(ws2,thr2) ;
  ws3 = HardThreshAbs(ws3,thr3) ;
 case 'SOFT'
  ws1 = SoftThreshAbs(ws1,thr1) ;
  ws2 = SoftThreshAbs(ws2,thr2) ;
  ws3 = SoftThreshAbs(ws3,thr3) ;
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
    
    
