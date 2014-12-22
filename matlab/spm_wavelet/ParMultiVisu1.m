function ws = ParMultiVisu1(wc,type,L)
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
% Parallel version of MultiVisu1
%
% [ws] = ParMultiVisu1(wc,L)
%        where wc = [timepoints x voxels]

ws=wc;

[n,J] = dyadlength(1:size(ws,1));

for j=J-1:-1:L
  
    LengthDyadj=2^j;
    thr=sqrt(2*log(LengthDyadj));
    
    switch(upper(type));
    case 'HARD'
         ws(dyad(j),:) = ws(dyad(j),:) .* (ws(dyad(j),:) > thr);
    case 'SOFT'
         res = (ws(dyad(j),:) - thr);
         res = (res + abs(res))/2;
         ws(dyad(j),:) = sign(ws(dyad(j),:)).*res;
    end;
end










