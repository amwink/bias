function f = iffwt2(H,W,k,unser)
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
% iffwt2.m
%
% compute the ifwt2 in the frequency domain
%
% FORMAT w = iffwt2(H,W,k,unser)
%
% H = the synthesis filters in the frequency domain
% W = the wavelet coefficients matrix
% k = the reconstruction level
% unser: filters defined in freq domain (0=no, 1=yes)
%
% f = the reconstructed signal
%
% Alle Meije Wink 2003/7/28

% check for square image
if (size(W,1)~=size(W,2))
  disp('images must be square');
  k=0;
else
  low=size(W,1);
end;

if ((nargin<3) | (~isint(low/2^k)))
  fprintf('level %d reconstruction useless\n',k);
  k=0;
else
  N=low/2^k;
end;

% assign the scaling and wavelet filters
if((nargin>=4) & unser)
  G=conj(H(2,:));
  H=conj(H(1,:));
else
  H=[H(:) qmf(H(:))]; % make corresponding QMF 
  H=fft(H,low)';      % goto freq domain
  G=conj(H(2,:));
  H=conj(H(1,:));
end;

HH = H.' * H;
GH = G.' * H;
HG = H.' * G;
GG = G.' * G;

% Do the reconstruction
f = W; % If nothing happens, output = input

for j = 1:k,
  dsl = low/(2*N); % downsampling level

  HH1 = dsamp(HH,dsl);
  HG1 = dsamp(HG,dsl);
  GH1 = dsamp(GH,dsl);
  GG1 = dsamp(GG,dsl);

  C  = f(1:N,1:N);	   
  D1 = f(1:N,N+1:2*N);		
  D2 = f(N+1:2*N,1:N);		
  D3 = f(N+1:2*N,N+1:2*N);	
   
  iZ =   HH1 .* [C C ; C C]     + HG1 .* [D1 D1 ; D1 D1] ...
       + GH1 .* [D2 D2 ; D2 D2] + GG1 .* [D3 D3 ; D3 D3];

  f(1:2*N,1:2*N) = iZ;

  % set up next iteration
  N = 2*N;
end

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%


