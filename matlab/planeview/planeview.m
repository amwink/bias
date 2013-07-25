function axisnrs=planeview(volume,cmap)
%
% planeview.m 
%
% displays a volume plane-by-plane
%
% SYNTAX planeview(volume,cmap)
%
% volume = 3D data set
% cmap = color map (64x3-matrix or string)
%
% (C) Alle Meije Wink 10/10/2012
%

volume=squeeze(real(volume));
if(prod(size(size(volume)))<2|prod(size(size(volume)))>3)
  error('2D/3D data required')
end;

if (nargin<2)
  cmap='gray';
  if (nargin<1)
    error('argument <3D data set> missing');
  end;
end;

zdim=size(volume,3);

cols=ceil(sqrt(zdim));
rows=ceil(zdim/cols);

minval=double(min(volume(volume~=0)));
maxval=double(max(volume(volume~=0)));

if (minval==maxval)
  minval=minval-1;
end;
  
for i=1:zdim
    subplot(cols,rows,i);
    imagesc(rot90(volume(:,:,i)),[minval maxval]);
    set(get(gca,'title'),'string',num2str(i));
    axisnrs(i)=gca;
    set(gca,'Visible','off');
    axis square;
end;

colormap(cmap);

return







