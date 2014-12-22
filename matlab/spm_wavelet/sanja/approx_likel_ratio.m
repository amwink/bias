%------------------------------------------------------------------------------------
% Function: approx_likel_ratio.m
% 
% Piece-wise linear fitting of the log-likelihood ratio for the 
% wavelet coefficient magnitudes using empirical conditional densities according to [1].
% 
% Input: 
%       M - the detected edge mask for the given detail image
%       D - the input detail image  
%       ROI - region of interest from which the signal and noise statistics is estimated%
% 
% Output: 
%       arg    - range of values for which the likelihood ratio is estimated
%       ratio  - the fitted likelihood ratio
%       confid - confidence with respect to the estimated ratio: takes values {0,1} 
%                confid=0 means that no reliable estimate was found
%                confid=1 a reliable estimate is found
%
%[1] A. Pizurica, W. Philips, I. Lemahieu and M. Acheroy, 
%   "A Versatile Wavelet Domain Noise Filtration Technique for Medical Imaging," 
%    IEEE, Trans. on Medical Imaging, vol. 22, no. 3, pp. 323--331, March 2003. %
%
%  
% Author: Aleksandra Pizurica <Aleksandra.Pizurica@telin.UGent.be> TELIN/Ghent University
% Last modified: 7. May 2004.
%------------------------------------------------------------------------------------

function [arg,ratio,confid] = approx_likel_ratio(M,D)

% Compute p(x|1)

Max=ceil(max(max(abs(D))));
h=Mag_sector_row(M,D);
%step=round(Max/40);
step=4;
x=0:step:Max;
len=length(x);
[pdf_signal,x]=hist(h,x);
pdf_signal=pdf_signal/(sum(pdf_signal)*step);

%----------------------------------------------------
% Compute p(x|0)

h=Mag_sector_row(1-M,D);
pdf_noise=hist(h,x);
%pdf_noise(1)=pdf_noise(2);
pdf_noise=pdf_noise/(sum(pdf_noise)*step);

%--------------------------------------------------------------
arg=0:Max;

def_range=find(pdf_signal~=0 & pdf_noise~=0);
Rlog=log(pdf_signal(def_range))-log(pdf_noise(def_range));

if (length(def_range)<5)
   confid=0;
   ratio=zeros(1,Max);
else
   confid=1;

ind=1;
while ((Rlog(ind)<0) & (ind<length(Rlog)))
   ind=ind+1;
end

T=ind;  %boarder between two different fitting parts

range1=1:T;

if (T>2)
   Tl=T-1;
else
   Tl=1; 
end

%discard values at the very end of histogram tails
if(round(2*length(def_range)/3)>Tl)
   Tr=round(2*length(def_range)/3);
else
   Tr=length(def_range);
end

range2=Tl:Tr;


if(length(range1)>=2 & length(range2)>=2)

    fit_range_1=x(def_range(range1));
    fit_range_2=x(def_range(range2));

    [p1 s1]=polyfit(fit_range_1,Rlog(range1),1);
    [p2 s2]=polyfit(fit_range_2,Rlog(range2),1);

    Tx=x(def_range(T));

    rfit1=polyval(p1,0:Tx);
    rfit2=polyval(p2,(Tx+1):Max);

    Rlog_fit=cat(2,rfit1,rfit2);

else
   [p s]=polyfit(x(def_range),Rlog,1);
   Rlog_fit=polyval(p,0:Max);
end

ratio=exp(Rlog_fit);

%figure('units','normalized','position',[0.1 0.1 0.3 0.3]);
%plot(x,pdf_signal,x,pdf_noise);
%title('Magnitude likelihood');


%Title='Likelihood ratio';
%figure('units','normalized','position',[0.1 0.1 0.3 0.3]);
%plot(x(def_range),Rlog,'k.',0:Tx,rfit1,'r',(Tx+1):Max,rfit2,'b');
%title(Title);

hist_ratio=pdf_signal(def_range)./pdf_noise(def_range);

%figure('units','normalized','position',[0.1 0.1 0.3 0.3]);
%plot(x(def_range),hist_ratio,'k.',0:Max,ratio,'k');
%plot(0:Max,ratio,'k');
%title(Title);


end