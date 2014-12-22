function W=fwt2ffwt(FW,level)
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
% fwt2ffwt.m
%
% Convert an FWT into an FWD
%
% FORMAT: W=fwt2ffwt(FW,level)
%
% FW    = Wavelet decomposition of a signal s (1D or 2D)
% level = level of decomposition
%
% W     = FWD (frequency wavelet decomposition) of s
%
% Alle Meije Wink 30-11-2001

ps=prod(size(size(FW)));
if (ps>2)
  error('works only on 1D or 2D data')
end;

er=max(size(FW));
for i=1:level
  hr=er/2;
  switch(min(size(FW)))
   case 1,    W(hr+1:er)=fft(FW(hr+1:er));
   otherwise, W(1:hr,hr+1:er)=fft2(FW(1:hr,hr+1:er));
              W(hr+1:er,1:hr)=fft2(FW(hr+1:er,1:hr));
	      W(hr+1:er,hr+1:er)=fft2(FW(hr+1:er,hr+1:er));
  end;
  er=er/2;
end;

switch(min(size(FW)))
 case 1,    W(1:hr)=fft(FW(1:hr));
 otherwise, W(1:hr,1:hr)=fft2(FW(1:hr,1:hr));
end;

%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%
