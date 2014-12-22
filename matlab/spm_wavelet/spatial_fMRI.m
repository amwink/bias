function spatial_fMRI(method,files,level,degree,smooth,output_prefix,output_dir);
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
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%

% check the files that must be processed
fprintf('\rmapping volumes...\r');
V=spm_vol(files);
lV=length(V);

% set the output directory
newdname=spm_str_manip(V(1).fname,'h');
if (~exist('output_dir','var'))
  newdname=[newdname '/' method '/'];
else
  if (~isempty(output_dir))
    newdname=output_dir;
  end;
  if (newdname(end)~='/')
    newdname=[newdname '/'];
  end;
end;
if (~exist(newdname,'dir'))
  eval(['!mkdir ' newdname]);
end;

% set the output prefix
if (~exist('output_prefix','var'))
  output_prefix='';
end;
newdname=[newdname output_prefix];
  
% check if all images are equal in size
for i=1:lV;
  x(i)=V(i).dim(1); 
  y(i)=V(i).dim(2);
  z(i)=V(i).dim(3);
end;
if (prod(x-x(1))|prod(y-y(1))|prod(z-z(1)))
  error('images are not equal in size')
end;
xdim=x(1);
ydim=y(1);
zdim=z(1);

% call alternative routine 1: smoothing
if(strcmp(upper(method(1:3)),'SMO'))
  
  voxsize=[V(i).mat(1,1) V(i).mat(2,2) V(i).mat(3,3)];
  
  if(~exist('smooth','var'))
    if (length(method)>6)
      smooth=str2num(method(7:end));
      smooth=smooth*[1 1 1];
    else
      smooth=voxsize;
    end;
  end;
  
  kernel=smooth;       
  %kernel=smooth./voxsize; % (??? maybe inherited from spm_conv)
  %fprintf('smoothing kernel = [%0.4f %0.4f %0.4f] voxels\n',kernel(1),kernel(2),1)

  for i=1:lV,
    fprintf('\rprocessing %s with [%0.1f %0.1f %0.1f] mm...\r', ...
	    spm_str_manip(V(i).fname,'t'),kernel(1),kernel(2),kernel(3));
    pause(.1);
    newfname=[newdname spm_str_manip(V(i).fname,'t')];
    spm_smooth(V(i),newfname,kernel);
  end;

  fprintf('\r%-80s\n',sprintf('finished denoising with method %s',method))
  return;
  
end;

% make mean image
fprintf('\rcalculate mean image...\r');
mim=zeros(xdim,ydim,zdim);
if (lV>1)
  for i=1:lV,
    tmp=spm_read_vols(V(i));
    mim=mim+tmp;
  end;
else
  tmp=spm_read_vols(V(1));
  tmp=mean(tmp(:));
end
mim=mim/lV;

% call alternative routine 2: sanja's routine 
% (generalized likelihood based wavelet denoising)
if(strcmp(upper(method(1:3)),'GEN'))
  for i=1:lV,
    
    fprintf('\rprocessing %s...\r',spm_str_manip(V(i).fname,'t'));
    newfname=[newdname spm_str_manip(V(i).fname,'t')];
    v=spm_read_vols(V(i));
    mv=max(v(:));
    msk=(v>mv/8);
    
    v=v-mim;
    
    for j=1:zdim
      tmp=squeeze(v(:,:,j));
      tmpmsk=squeeze(msk(:,:,j));
      
      maxx=find(sum(tmpmsk,2));
      maxy=find(sum(tmpmsk,1));

      minx=min(maxx);
      maxx=max(maxx);
      miny=min(maxy);
      maxy=max(maxy);           
      
      tmp=tmp(minx:maxx,miny:maxy);
      currWdbz=warning('query','MATLAB:divideByZero');
      warning('off','MATLAB:divideByZero');

      %genlik_MRI(0,tmp,2,5);
      %genlik_BOLD(0,tmp,2,5);
      
      if(strcmp(upper(method(end-2:end)),'MRI'))
	fprintf('\rusing the MRI method\r');
	v(minx:maxx,miny:maxy,j)=genlik_MRI(0,tmp,2,5);
      else
	fprintf('\rusing the BOLD method\r');
	v(minx:maxx,miny:maxy,j)=genlik_BOLD(0,tmp,2,5);
      end
      warning(currWdbz.state,'MATLAB:divideByZero');
      v(:,:,j)=v(:,:,j).*tmpmsk;
    end
    
    v=v+mim;
    
    V(i).fname=[newdname spm_str_manip(V(i).fname,'t')];
    spm_write_vol(V(i),v);
    
  end;  

  fprintf('\r%-80s\n',sprintf('finished denoising with method %s',method))
  return;

end

% make spatial wavelet filters
[an sy]=FFTfractsplinefilters(xdim,degree,'*ortho');

% low-frequency cutoff for wavelet shrinkage
lfc=log(xdim)/log(2)-level;

% read image volumes
for i=1:lV,
  fprintf('\rprocessing %s...\r',spm_str_manip(V(i).fname,'t'));
  tmp = spm_read_vols(V(i));
  tmp = tmp-mim;
  
  for j=1:zdim
    tmp(:,:,j)=fft2(squeeze(tmp(:,:,j)));
    tmp(:,:,j)=ffwt2(an,squeeze(tmp(:,:,j)),level,1);    
    tmp(:,:,j)=ffwt2fwt(squeeze(tmp(:,:,j)),level);
  end;  
  for j=1:zdim
    plane=tmp(:,:,j);
    aplane=abs(plane);
    plane=sign(plane);
    [aplane sigma]=NormNoise2(aplane);
    
    % use MAD to check for autocorrelation
    MAD1=aplane(end/2+1:end,end/2+1:end);
    MAD2=aplane(end/4+1:end/2,end/4+1:end/2);
    if ( (level>1) & (mean(MAD2(:))/mean(MAD1(:))>1.25) )
      usemulti=1;
    else
      usemulti=0;
    end;

    usemulti=1;
    % usemulti=0 does not work properly
    %fprintf(' %d ',usemulti);
    
    if (usemulti)
      switch(upper(method))
       case 'HYBRID', aplane=MultiHybrid2(aplane,lfc);
       case 'INV',    aplane=MultiInvShrink2(aplane,lfc);
       case 'JAMES',  aplane=MultiWaveJS2(aplane,lfc);
       case 'MAD',    aplane=MultiMAD2(aplane,lfc);
       case 'MINMAX', aplane=MultiMinMax2(aplane,lfc);
       case 'SURH',   aplane=MultiSURE2(aplane,'Hard',lfc);
       case 'SURS',   aplane=MultiSURE2(aplane,'Soft',lfc);
       case 'VISH',   aplane=MultiVisu2(aplane,'Hard',lfc);
       case 'VISS',   aplane=MultiVisu2(aplane,'Soft',lfc);
       otherwise,     disp(sprintf('no denoising... %s unknown',method));
      end;
    else    
      switch(upper(method))
       case 'HYBRID', aplane=SingleHybrid2(aplane,lfc);
       case 'INV',    aplane=SingleInvShrink2(aplane,lfc);
       case 'JAMES',  aplane=SingleWaveJS2(aplane,lfc);
       case 'MAD',    aplane=SingleMAD2(aplane,lfc);
       case 'MINMAX', aplane=SingleMinMax2(aplane,lfc);
       case 'SURH',   aplane=SingleSURE2(aplane,'Hard',lfc);
       case 'SURS',   aplane=SingleSURE2(aplane,'Soft',lfc);
       case 'VISH',   aplane=SingleVisu2(aplane,'Hard',lfc);
       case 'VISS',   aplane=SingleVisu2(aplane,'Soft',lfc);
       otherwise,     disp(sprintf('no denoising... %s unknown',method));
      end;  
    end
      
    plane=(aplane.*plane)*sigma;
    tmp(:,:,j)=plane;
      
  end;  
      
  for j=1:zdim
    tmp(:,:,j)=fwt2ffwt(squeeze(tmp(:,:,j)),level);
    tmp(:,:,j)=iffwt2(sy,squeeze(tmp(:,:,j)),level,1);    
    tmp(:,:,j)=ifft2(squeeze(tmp(:,:,j)));
  end;  
  tmp=real(tmp)+mim;
  
  V(i).fname=[newdname spm_str_manip(V(i).fname,'t')];
  Vout(i)=spm_write_vol(V(i),tmp);
end;
clear tmp;

fprintf('\r%-80s\n',sprintf('finished denoising with method %s',method))

return
    
    
    
    
    















