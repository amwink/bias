function W = ffwt2(H,X,k,unser)
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
% ffwt2.m
%
% compute the fwt2 in the frequency domain
%
% FORMAT w = ffwt2(H,X,k,unser)
%
% H = the analysis filters in the frequency domain
% X = the matrix in the frequency domain
% k = the decomposition level
% unser (optional): if 1 -> wavelets defined in freq. domain
%
% w = the wavelet coefficients matrix
%
% Alle Meije Wink 2003/7/28

% check for square image
if (size(X,1)~=size(X,2))
  disp('images must be square');
  k=0;
else
  N=size(X,1);
end;

if ((nargin<3) | (~isint(N/2^k)))
  fprintf('decomposition to level %d useless\n',k);
end;
% Do the k-level decomposition, where N == blocksize*2^k

% initialise output matrix
W = zeros(N,N);

% define approximation and horizontal, vertical and diagonal detail
% filters
% assign the scaling and wavelet filters
if((nargin>=4) & unser)
  G=H(2,:);
  H=H(1,:);
else
  H=[H(:) qmf(H(:))]; % make corresponding QMF 
  H=fft(H,N)';      % goto freq domain
  G=H(2,:);
  H=H(1,:);
end;

HH = H.' * H;
GH = G.' * H;
HG = H.' * G;
GG = G.' * G;

for j=1:k
  % filter and downsample approximation and details
  C  = HH .* X;
  C  = .25*(C(1:end/2,1:end/2)+C(1:end/2,end/2+1:end)+...
	    C(end/2+1:end,1:end/2)+C(end/2+1:end,end/2+1:end));
  D1 = HG .* X;
  D1 = .25*(D1(1:end/2,1:end/2)+D1(1:end/2,end/2+1:end)+...
	    D1(end/2+1:end,1:end/2)+D1(end/2+1:end,end/2+1:end));
  D2 = GH .* X;
  D2 = .25*(D2(1:end/2,1:end/2)+D2(1:end/2,end/2+1:end)+...
	    D2(end/2+1:end,1:end/2)+D2(end/2+1:end,end/2+1:end));
  D3 = GG .* X;
  D3 = .25*(D3(1:end/2,1:end/2)+D3(1:end/2,end/2+1:end)+...
	    D3(end/2+1:end,1:end/2)+D3(end/2+1:end,end/2+1:end));

  W(1:N/2,N/2+1:N)   = D1; % hor.detail
  W(N/2+1:N,1:N/2)   = D2; % ver detail
  W(N/2+1:N,N/2+1:N) = D3; % dia detail

  % set up the next iteration
  X  = C;
  HH = HH(1:2:end,1:2:end);		
  HG = HG(1:2:end,1:2:end);	
  GH = GH(1:2:end,1:2:end);	
  GG = GG(1:2:end,1:2:end);	
  N = N/2;
end

W(1:N,1:N) = C;

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%
