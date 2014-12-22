%------------------------------------------------------------------------------------
% Function: rem_noise_adapt.m
% 
% Remove noise from a given wavelet subband using the method of [1].
% 
% Input: 
%
%       D..........- input wavelet subband (detail image) 
%       Dcoars.....- the corresponding detail image at the next coarser scale
%       sigma......- noise standard deviation in the subband D
%       Thr_factor - the threshold multiplication factor for the signal of interest
%       W..........- parameter which describes the window size: w=(window size-1)/2
% 
% Output: 
%       D_cl........- denoised subband D
%       Mask........- the estimated edge mask for the subband D
%
%[1] A. Pizurica, W. Philips, I. Lemahieu and M. Acheroy, 
%   "A Versatile Wavelet Domain Noise Filtration Technique for Medical Imaging," 
%    IEEE, Trans. on Medical Imaging, vol. 22, no. 3, pp. 323--331, March 2003. %
%
%  
% Author: Aleksandra Pizurica <Aleksandra.Pizurica@telin.UGent.be> 
%         TELIN/Ghent University
%
% Last modified: 7. May 2004.
%------------------------------------------------------------------------------------


function [D_cl,Mask]=rem_noise_adapt(D,Dcoars,sigma,Thr_factor,W);

[M,N]=size(D);

Mask=init_mask(D.*Dcoars,(Thr_factor*sigma)^2);

O=sum(sum(Mask));

if(O<M*N)
   r=O/(M*N-O);
else
   r=100;
end


if (r<1e-4)
   D_cl=zeros(M,N);
elseif (max(Mag_sector_row(1-Mask,D))<1)  % almost flat background - noise-free image
   D_cl=D;
else
   [arg,likel_ratio,confid]=approx_likel_ratio(Mask,D);
   [arg,prior_ratio,confid]=approx_prior_ratio(Mask,D,Dcoars,W);
   if(confid==1)
     SHRINK=adapt_energy(D,Dcoars,likel_ratio,prior_ratio,W,r);
     D_cl=D.*SHRINK;
   else
     D_cl=D;  % the likelihood_ratio or the prior_ratio were not reliably estimated from the nput data
   end
end



   
