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

function y=sinc(x)
y=ones(size(x));
i=find(x);
y(i)=sin(pi*x(i))./(pi*x(i));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% courtesy K. Bell for making this file available - Alle Meije Wink
% 
% http://ite.gmu.edu/DetectionandEstimationTheory/OAP/Figures/Ch2/sinc.m
%
% sinc.m
% same as Matlab's sinc function
% sinc(x) =  sin(pi*x)/(pi*x)
% Last updated 8/31/00 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
