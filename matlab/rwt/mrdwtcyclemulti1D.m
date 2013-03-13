% mrdwtcyclemulti1D.m
%
% Same as multi1Drdwt(), 
% returns the results in a nicer form
% (each level in a separate cell).
% 
% SYNTAX: [a, b] = mrdwtcyclemulti1D(x, h, L)
% yw - wavelet coefficients.  NxNxLx3.
% ys - scaling coefficients.  NxN.
%
% (C) Alle Meije Wink (a.wink@vumc.nl)
% 

function [yc,yd] = mrdwtcyclemulti1D(x, h, L)

[yc,yd] = multi1Drdwt(x, h, L);

myd=size(yd,1);
nyd=size(yd,2);

yd=mat2cell(yd,myd,ones(1,L)*nyd/L);

return 



