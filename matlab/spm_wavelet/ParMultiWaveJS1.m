function [ws] = ParMultiWaveJS1(wc,L)
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
% Parallel version of MultiWaveJs1
%
% [ws] = ParMultiWaveJS1(wc,L)
%        where wc = [timepoints x voxels]

ws = wc;

[n,J] = dyadlength(1:size(ws,1));

for j=J-1:-1:L;
  % ws(dyad(j))=JamesStein(ws(dyad(j))) vectorised
  % see JamesStein.m
  % keyboard
  L2DyadSum=sum(ws(dyad(j),:).^2);
  LengthDyadj=2^j;
  shrnk=max((L2DyadSum - (LengthDyadj-2))./L2DyadSum,0);
  ws(dyad(j),:)=(ones(LengthDyadj,1)*shrnk).*ws(dyad(j),:);
end;    




