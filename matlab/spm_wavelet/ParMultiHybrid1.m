function ws = ParMultiHybrid1(wc,L)
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
% Parallel version of MultiHybrid1
%
% [ws] = ParMultiHybrid1(wc,L)
%        where wc = [timepoints x voxels]

ws=wc;

[n,J] = dyadlength(1:size(ws,1));

spatdim=size(ws,2);

for j=J-1:-1:L;
    % ws(dyad(j))=HybridThresh(ws(dyad(j))) vectorised
    LengthDyadj=2^j;
    thr=sqrt(2*log(LengthDyadj));
    
    
        magic = sqrt(2*j);
        eta = (sum(ws(dyad(j),:).^2) - n)/n; % norm(a)=sum(a.^2).^(1/2)
        crit = j^(1.5)/sqrt(LengthDyadj);
	
	EtaLtCrit=find(eta<crit);
	LtDim=length(EtaLtCrit);
	EtaGeCrit=find(eta>=crit);
	GeDim=length(EtaGeCrit);
	
        % if eta < crit
	if (LtDim)
	  res = (ws(dyad(j),EtaLtCrit) - magic);
	  res = (res + abs(res))/2;
	  ws(dyad(j),EtaLtCrit) = sign(ws(dyad(j),EtaLtCrit)).*res;
	end;
        % if eta >= crit
	% thr = ValSUREThresh
	if (GeDim)
	  a = sort(ws(dyad(j),EtaGeCrit)).^2 ;
	  b = cumsum(a);
	  n = 2^j;
	  c = linspace(n-1,0,n);
	  s = b+diag(c)*a;
	  risk = (s+n-2.*((1:n)'*ones(1,GeDim)))/n;
	  [guess,ibest] = min(risk);
	  thr = sqrt(a(ibest));
	  thr = min(thr,magic);
	  % soft thresholding
	  res = (ws(dyad(j),EtaGeCrit) - ones(LengthDyadj,1)*thr);
	  res = (res + abs(res))/2;
	  ws(dyad(j),EtaGeCrit) = sign(ws(dyad(j),EtaGeCrit)).*res;
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
    











