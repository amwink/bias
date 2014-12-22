function ws = MultiMAD2(wc,L)
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

for j=J-1:-1:L    
    thr=length(dyad(j));
    thr=sqrt(2*log(thr));            % VisuShrink threshold

    tws=ws(dyad(j),dyad(j));
    scale1 = median(tws(:))/.6745;
    tws=ws(bdyad(j),dyad(j));
    scale2 = median(tws(:))/.6745;
    tws=ws(dyad(j),bdyad(j));
    scale3 = median(tws(:))/.6745;

    if (scale1)
       ws(dyad(j),dyad(j))=scale1*HardThreshAbs(ws(dyad(j),dyad(j))/scale1,thr) ;
    else
       ws(dyad(j),dyad(j))=HardThreshAbs(ws(dyad(j),dyad(j)),thr) ;
    end
    if (scale2)
       ws(bdyad(j),dyad(j))=scale2*HardThreshAbs(ws(bdyad(j),dyad(j))/scale2,thr) ;
    else
       ws(bdyad(j),dyad(j))=HardThreshAbs(ws(bdyad(j),dyad(j)),thr) ;
    end
    if (scale3)
       ws(dyad(j),bdyad(j))=scale3*HardThreshAbs(ws(dyad(j),bdyad(j))/scale3,thr) ;
    else
       ws(dyad(j),bdyad(j))=HardThreshAbs(ws(dyad(j),bdyad(j)),thr) ;
    end
end

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
    
