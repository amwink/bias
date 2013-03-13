% mirdwtcyclemulti1D.m
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

function r = mirdwtcyclemulti1D(yc, yd, h, L)

yd=cell2mat(yd);

r = multi1Dirdwt(yc, yd, h, L);

return 



