function [ws] = SingleWaveJS2_old(wc,L)
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

ws=wc(lastc+1:end,:);
ws=[wc(1:lastc,lastc+1:end)' ws];

wss=size(ws);
ws=JamesStein(ws(:));

ws=reshape(ws,wss);
ws=[ [ wc(1:lastc,1:lastc) ws(:,1:lastc)' ]; ...
       ws(:,lastc+1:end)                  ];

return
       
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    






