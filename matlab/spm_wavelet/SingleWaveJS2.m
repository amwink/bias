function [ws] = SingleWaveJS2(wc,L)
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
% WaveJS -- James-Stein Shrinkage Applied to  Wavelet dyads
%  Usage
%    [xh,xwh] = WaveJS(y,L,qmf)
%  Inputs
%    y      array of dyadic length 2^J
%           NORMALIZED TO NOISE LEVEL 1 (See NOISENORM)
%    L      Low-Frequency cutoff for shrinkage (e.g. L=4)
%           SHOULD BE L << J
%    qmf    Quadrature Mirror Filters for Wavelet Transform
%  Outputs
%    xh   estimate, obtained by James-Stein Shrinkage
%         on Wavelet Coefficients
%    xwh  Wavelet Transform of estimate
%  See Also
%    FWT_PO, IWT_PO, JamesStein
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
     
% apply the threshold
ws1 = JamesStein(ws1);
ws2 = JamesStein(ws2);
ws3 = JamesStein(ws3);

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
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    






