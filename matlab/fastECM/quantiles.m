function  q = quantile(x,p,method)
%QUANTILE Empirical quantile (percentile).
%
%         q = quantile(x,p)
%	  
%	  The empirical quantile is computed by one of three ways
%	  determined by a third input argument (with default 1).
%
%	  1. Interpolation so that F(X_(k)) == (k-0.5)/n.
%	  2. Interpolation so that F(X_(k)) == k/(n+1).
%	  3. Based on the empirical distribution.
%
%	  If input  x  is a matrix then the quantile is computed for 
%	  every column. Input  p  may be vector also. It even 
%	  accepts  x  being a vector and  p  a matrix!

%       Copyright (c) Anders Holtsberg, 1995, 1998

if nargin<3, method=1; end
if min(size(x)) == 1
   x = x(:);
   q = zeros(size(p));
else
   if min(size(p)) > 1 
      error('Not both matrix x and matrix p input')
   end
   q = zeros(length(p),size(x,2));
end
if any(any((p>1|p<0)))
   error('Input p is not probability')
end

x = sort(x); 
p = p(:);
n = size(x,1);
if method == 3
   qq1 = x(ceil(max(1,p*n)),:); 
   qq2 = x(floor(min(p*n+1,n)),:);
   qq = (qq1+qq2)/2;
else                         
   x = [x(1,:); x; x(n,:)];
   if method == 2
      % This method is from Hjort's "Computer
      % intensive statistical methods" page 102
      i = p*(n+1)+1;
   else % Metod 1
      i = p*n+1.5;
   end
   iu = ceil(i);
   il = floor(i);
   d = (i-il)*ones(1,size(x,2));
   qq = x(il,:).*(1-d)+x(iu,:).*d;
end

q(:) = qq;

