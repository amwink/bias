function w = ffwt(h,X,k,unser)
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
% ffwt.m
%
% compute the fwt in the frequency domain
%
% FORMAT w = ffwt2(H,X,k,unser)
%
% H = the analysis filters in the frequency domain
% X = the signal in the frequency domain
% k = the decomposition level
% unser (optional): if 1 -> wavelets defined in freq. domain
%
% w = the wavelet coefficients
%
% Alle Meije Wink 2003/7/28
% See Also: FFWT

if (length(X)~=prod(size(X)))
  error('ffwt operates on 1D vectors only')
end;

oldxsize=size(X);
X=X(:).';
N=length(X);

if ((nargin<3) | (~isint(N/2^k)))
  fprintf('decomposition to level %d useless\n',k);
  k=0;
end;

if((nargin>=4) & unser)
  H=h(1,:);
  G=h(2,:);
else
  h=[h(:) qmf(h(:))];
  h=fft(h,N)';
  H = h(1,:);		% Full size fft
  G = h(2,:);		% Full size fft
end;

w = zeros(N,1);

for j=1:k,

  % filter
  C = H.*X;
  D = G.*X;

  % downsample
  C = C(1:N/2)+C(N/2+1:N);
  C = .5*C;
  D = D(1:N/2)+D(N/2+1:N);
  D = .5*D;
  
  % store the detail
  w(N/2+1:N) = D;

  % set up next iteration
  X = C;
  H = H(1:2:N);				
  G = G(1:2:N);
  N = N/2;				

end

% store approximation
w(1:N) = C;
w=reshape(w,oldxsize);

return;

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%




