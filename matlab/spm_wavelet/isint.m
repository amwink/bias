function k = isint(x)
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
%k = ISINT(x) detect integers
%
% Description:
%   isint(x) returns ones where the elements of x are integers and zeros
%   where they are not.
%
% Examples:
%   A = [pi NaN 0 0.001];
%   isint(A)
%   ans = 
%        0   0   1   0
%   any(isint(A))
%   ans =
%        1
%
% See also:
%   isinf, isnan, isstr, any
%
% M.D. van der Laan -- RuG/CS -- 10-dec-93 
%

k = abs(x - fix(x)) < eps;

return

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%
