function ws = ParMultiSURE1(wc,type,L)
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
% Parallel version of MultiSURE1
%
% [ws] = ParMultiSURE1(wc,L)
%        where wc = [timepoints x voxels]

ws=wc;

[n,J] = dyadlength(1:size(ws,1));
spatdim=size(ws,2);

for j=J-1:-1:L
    LengthDyadj=2^j;
    thr=sqrt(2*log(LengthDyadj));
   
    % thr=min([ValSUREThresh(ws(dyad(j))') thr]) vectorised
    a = sort(ws(dyad(j),:)).^2 ;
    b = cumsum(a);
    n = 2^j;
    c = linspace(n-1,0,n);
    s = b+diag(c)*a;
    risk = (s+n-2.*((1:n)'*ones(1,spatdim)))/n;
    [guess,ibest] = min(risk);
    thr = sqrt(a(ibest));
    
    switch(upper(type));
     case 'HARD'
      ws(dyad(j),:) = ws(dyad(j),:) .* (ws(dyad(j),:) > ones(n,1)*thr);
     case 'SOFT'
      res = (ws(dyad(j),:) - ones(LengthDyadj,1)*thr);
      res = 0.5*(res + abs(res));
      ws(dyad(j),:) = sign(ws(dyad(j),:)).*res;
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
    








