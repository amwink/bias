function f = iffwt(h,w,k,unser)
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
% iffwt.m
%
% compute the ifwt in the frequency domain
%
% FORMAT w = iffwt(H,W,k,unser)
%
% H = the synthesis filters in the frequency domain
% W = the wavelet coefficients matrix
% k = the reconstruction level
% unser: filters defined in freq domain (0=no, 1=yes)
%
% f = the reconstructed signal
%
% Alle Meije Wink 2003/7/28

w = w(:);
low=length(w);

f=w.'; % If nothing happens, output = input

if ((nargin<3) | (~isint(low/2^k)))
  fprintf('level %d reconstruction useless\n',k);
  k=0;
else
  N=low/2^k;
end;

if((nargin>=4) & unser)
  Hfull=conj(h(1,:));
  Gfull=conj(h(2,:));
else
  h=[h(:) qmf(h(:))]; % make corresponding QMF 
  h=fft(h,low)';      % goto freq domain
  Hfull=conj(h(1,:));
  Gfull=conj(h(2,:));
end;

Hfull=conj(Hfull);
Gfull=conj(Gfull);

for j = 1:k,
  % Undecimate
  dsl = low/(2*N);		% downsampling level
  
  H = dsamp(Hfull,dsl);
  G = dsamp(Gfull,dsl);
  
  C = f(1:N);
  D = f(N+1:2*N);
  
  iZ = H .* [C C] +  G .* [D D];
  
  f(1:2*N) = iZ;
  
  % Setup next iteration
  N = 2*N;
end

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%











