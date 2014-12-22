
% Function y=Thr_univ(X) estimates the noise standard deviation
% using the median estimator of Donoho 
% y=median({X[i][j]-median({X[i][j]})}/0.6745
%
% Author: Aleksandra Pizurica/TELIN/Ghent University
% E-mail: Aleksandra.Pizurica@telin.UGent.be

function y=Thr_univ(X)

med=median(median(X));
%y=median(median(abs(X-med)))/0.6745;
y=median(median(abs(X)))/0.6745;