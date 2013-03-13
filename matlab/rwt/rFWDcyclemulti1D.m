function [yc,yd] = rFWDcyclemulti1D(x, h, L, winfreq)
%
% rFWDcyclemulti1D.m
%
% Same as rFWD(), 
% returns the results in a nicer form
% (each level in a separate cell).
% 
% SYNTAX: [a, b] = rFWDcyclemulti1D(x, h, L)
% yw - wavelet coefficients.  NxNxLx3.
% ys - scaling coefficients.  NxN.
%
% (C) Alle Meije Wink (a.wink@vumc.nl) 
% 
  
  len=size(x,1);
  X=fft(x);
  
  if((nargin<4) | (~winfreq))
    H=[h(:) qmf(h(:))];
    cH=fft(H,len);
    H=conj(cH);
  else
    H=h;
  end
  
  [yc yd] = rFWD(X, H, L);
  
  yc=real(ifft(yc));
  yd=real(ifft(yd));
  
  myd=size(yd,1);
  nyd=size(yd,2);
 
  yd=mat2cell(yd,myd,ones(1,L)*nyd/L);

return 



