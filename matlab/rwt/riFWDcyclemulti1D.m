function r = riFWDcyclemulti1D(yc, yd, h, L, winfreq)
%
% riFWDcyclemulti1D.m
%
% Same as multi1Drdwt(), 
% returns the results in a nicer form
% (each level in a separate cell).
% 
% SYNTAX: [a, b] = riFWDcyclemulti1D(x, h, L)
% yw - wavelet coefficients.  NxNxLx3.
% ys - scaling coefficients.  NxN.
%
% (C) Alle Meije Wink (a.wink@vumc.nl)
% 

len=size(yc,1);
yd=cell2mat(yd);

YC=fft(yc);
YD=fft(yd);

if((nargin<5) | (~winfreq))
  H=[h(:) qmf(h(:))];
  H=fft(H,len);
else
  H=h;
end

r = riFWD(YC, YD, H, L);

r=real(ifft(r));

return 



