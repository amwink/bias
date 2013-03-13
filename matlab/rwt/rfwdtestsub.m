function theprofile=rfwdtestsub(lensig,numsig,levels);
%
% rfwdtestsub.m
%
% subroutine used in the test script multitestrfwd.m
% 
% decomposes a number <numsig> of <lensig>-point signals
% to <levels> wavelet decomposition levels, and reconstructs
% plots of the results will appear in new figures
%
% (C) Alle Meije Wink (a.wink@vumc.nl) 
%

% build time signal
[ts ots]=maketimesig(lensig,'white',10,0,'dho',1,32,0);
ts=ots+.1*randn(size(ots));
ts=ts(:)*ones(1,numsig);

% some constants
Ts=fft(ts);

% build othogonal spline wavelets in time and freq. domain
TYPES={'+ortho' '-ortho' '*ortho'};
typ=TYPES{3};

[an sy]=FFTfractsplinefilters(lensig,3,typ);
H=an'; %'
cH=conj(H);
h=real(ifft(conj(H(:,1))));

profile('on');

[rwts rwtd]=multi1Drdwt(ts,h,levels);
rwtr=multi1Dirdwt(rwts,rwtd,h,levels);

[fs fd]=rFWD(Ts,H,levels);
fr=riFWD(fs,fd,cH,levels);

[ffs ffd]=fsidwt(Ts,H,levels);
ffr=fisidwt(ffs,ffd,cH,levels);

theprofile=profile('info');

ifr=real(ifft(fr));
iffr=real(ifft(ffr));

mserror=mean((ifr(:)-ts(:)).^2);
msferror=mean((iffr(:)-ts(:)).^2);

if (nargout)
  fprintf('mse FWD = %0.3f mse fsidwt = %0.3f\n', ...
	  mserror,msferror);
end;

return
