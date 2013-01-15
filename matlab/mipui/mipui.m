function pos=mipui(v,mask);
%
% shows a MIP of a 3D data set, writes PNGs of views and returns the current co-ordinates
%
% format
% 
% POS = mipui ( volume, mask)
%
% where
%
% POS    == [x y z] == voxel coordinates of the 3d pointer
% volume == input volume
% mask   == mask volume
% 
% 'volume' and 'mask' can both either be 3D arrays or file names
% included are two nifti files, 'dmn.nii' and 'standard.nii'
% 
% 'dmn.nii' is the Default Mode Network from http://dx.crossref.org/10.1098%2Frstb.2005.1634
% 'standard.nii' is a brain mask in MNI standard space at the same resolution as 'dmn.nii'
% 
% given that
%
% >> D = nifti ( 'dmn.nii'      );
% >> d = D.dat ( : , : , :      );
% >> S = nifti ( 'standard.nii' );
% >> s = S.dat ( : , : , :      );
%
% the calls
%
% >> a = mipui ( 'dmn.nii' , 'standard.nii' );
% >> a = mipui ( d , 'standard.nii' );
% >> a = mipui ( 'dmn.nii' , s );
% >> a = mipui ( d , s );
%
% are equivalent.
%
% The button that shows the [X Y Z] position 
% exits the GUI, returning the position in POS.
% The button 'write mips' generates PNG files
% of the XY, XZ, and YZ planes of the mip
% with names mipxyX_Y_Z.png, mipxzX_Y_Z.png and mipyzX_Y_Z.png.
%
% (c) Alle Meije Wink, 10/10/2012
%
% requires: SPM5/8 
%

% check arguments
if (~nargin) 
  v=spm_select(1,'image','stat image','',pwd,'^*\.nii');
end

if ((isstr(v))&(~strcmp(v,'movecross')))
  v=nifti(v);
  v=v.dat(:,:,:);
end

% move the cross in an already existing window
if (strcmp(v,'movecross'))
  movecross;
  return
end

% remove other mip figures
mipfig=findobj('tag','mipfigure');
delete(mipfig);

mipfig=figure('tag','mipfigure');
  
if(size(size(v))~=3)
  error('mip is made for 3D data sets');
end

for i=1:size(v,3)
  mipyz(i,:)=-max(v(:,:,i));
  mipxz(:,i)=-max(v(:,:,i)')';
end

for i=1:size(v,1)
  mipxy(i,:)=-max(squeeze(v(i,:,:))');
end

if (nargin==2)
  if (isstr(mask))
    mask=nifti(mask);
    mask=mask.dat(:,:,:);
  end
  maskxy=rot90(squeeze(~sum(mask,3)));
  maskxz=rot90(squeeze(~sum(mask,2)));
  maskyz=squeeze(~sum(mask,1));
  maskyz=flipdim(maskyz,1);
  maskxz=flipdim(maskxz,1);
else
  maskxy=0;
  maskxz=0;
  maskyz=0;
end

medmip=(max(mipxy(:))+min(mipxy(:)))/2;

mipxy=rot90(mipxy);
mipxz=rot90(mipxz);
mipyz=rot90(mipyz);
mipxz=flipdim(mipxz,1);

subplot(2,2,1);
imagesc(mipxy+medmip*maskxy);     
xy=tidyaxis('x','y',1,'xy');

subplot(2,2,2);
imagesc(mipyz+medmip*maskyz);
yz=tidyaxis('z','y',2,'yz');

subplot(2,2,3);
imagesc(mipxz+medmip*maskxz);
xz=tidyaxis('x','z',3,'xz');

posbutton=initcrosses(size(v,1),size(v,2),size(v,3));

% wait until button is pressed
drawnow;

waitfor(posbutton,'userdata');
pos=get(posbutton,'string');
pos=str2num(pos);

delete(mipfig);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ax=tidyaxis(xs,ys,relpos,tagstr)

  colormap('gray');
  ax=gca;
  
  set(ax,'xtick',[],'ytick',[]);
  
  set(get(ax,'xlabel'),'String',xs);
  set(get(ax,'ylabel'),'String',ys);

  axis tight
  
  switch(relpos)
   case 1, set(ax,'position',[.1 .6 .3 .3])
   case 2, set(ax,'position',[.6 .6 .3 .3])
   case 3, set(ax,'position',[.1 .1 .3 .3])
  end
    
  set(ax,'tag',tagstr);  
  set(get(ax,'children'),'ButtonDownFcn','mipui(''movecross'')');
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function posbutton=initcrosses(xsi,ysi,zsi);
  
  x=round(xsi/2);
  y=round(ysi/2);
  z=round(zsi/2);
  
  % draw xy cross    
  ca=findobj('Tag','xy');
  set(gcf,'CurrentAxes',ca);
  setcross(x,y,'xycross'); 

   % redraw yz cross
  ca=findobj('Tag','yz');
  set(gcf,'CurrentAxes',ca);
  setcross(y,z,'yzcross'); 
  
  % redraw yz cross
  ca=findobj('Tag','xz');
  set(gcf,'CurrentAxes',ca);
  setcross(x,z,'xzcross'); 
 
  % make a button with the current position  
  posbutton=uicontrol (gcf, ...
		       'style','pushbutton', ...
		       'string',sprintf('%d %d %d',x,y,z), ...
		       'visible','on', ...
		       'units','normalized', ...
		       'tag','currentposition', ...
		       'Position',[.6 .3 .2 .1], ...
		       'callback','set(gco,''userdata'',get(gco,''string''))');

  wributton=uicontrol (gcf, ...
		       'style','pushbutton', ...
		       'string','write mips', ...
		       'visible','on', ...
		       'units','normalized', ...
		       'Position',[.6 .15 .2 .1], ...
		       'callback','savemips');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function movecross;
  
  ax=get(gca); 
  pos=gcp(ax);

  tag=get(gca,'tag');
  
  x=get(findobj('Tag','xycrossy'),'xdata');x=x(1);
  y=get(findobj('Tag','yzcrossx'),'ydata');y=y(1);
  z=get(findobj('Tag','xzcrossx'),'ydata');z=z(1);
  
  zmax=findobj('Tag','xz');
  zmax=get(zmax,'children');
  zmax=max(get(zmax(end),'ydata'));
  
  x=round(x);y=round(y);z=round(z);
  
  switch(tag)
   case 'xy',
    x=pos(1);
    y=pos(2);
   case 'yz',
    y=pos(2);
    z=pos(1);
   case 'xz',
    x=pos(1);
    z=pos(2);
  end

  x=round(x);y=round(y);z=round(z);

  % rewrite button
  set(findobj('tag','currentposition'),'string', ...
             sprintf('%d %d %d',x,y,z));
  
  % redraw xy cross
  ca=findobj('Tag','xy');
  set(gcf,'CurrentAxes',ca);
  setcross(x,y,'xycross'); 

   % redraw yz cross
  ca=findobj('Tag','yz');
  set(gcf,'CurrentAxes',ca);
  setcross(z,y,'yzcross'); 

  % redraw xz cross
  ca=findobj('Tag','xz');
  set(gcf,'CurrentAxes',ca);
  setcross(x,z,'xzcross'); 
  
return
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
function pos=gcp(axhandle)
  
  pos=get(gca,'currentpoint');
  pos=pos(1,1:2);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
function setcross(pos1,pos2,tagline); 

  delete(findobj('tag',[tagline 'x']))
  delete(findobj('tag',[tagline 'y']))
  
  line([pos1-2 pos1+2],[pos2 pos2],'color',[1 0 0], 'linewidth',2,'tag',[tagline 'x']);
  line([pos1 pos1],[pos2-2 pos2+2],'color',[1 0 0], 'linewidth',2,'tag',[tagline 'y']);
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function savemips;
 
  x=get(findobj('Tag','xycrossy'),'xdata');x=x(1);
  y=get(findobj('Tag','yzcrossx'),'ydata');y=y(1);
  z=get(findobj('Tag','xzcrossx'),'ydata');z=z(1);
  
  xyz=sprintf('%d_%d_%d',x,y,z);
  
  tmp=findobj('tag','xy');
  tmp=get(tmp,'children');
  tmp=tmp(ismember(tmp,findobj('type','image')));
  tmp=get(tmp,'cdata');
  tmp(y,:)=min(tmp(:));
  tmp(:,x)=min(tmp(:));  
  tmp=tmp-min(tmp(:));
  tmp=tmp/max(tmp(:));
  imwrite(tmp,['mipxy' xyz '.png'],'png');
  
  tmp=findobj('tag','yz');
  tmp=get(tmp,'children');
  tmp=tmp(ismember(tmp,findobj('type','image')));
  tmp=get(tmp,'cdata');
  tmp(y,:)=min(tmp(:));
  tmp(:,z)=min(tmp(:));  
  tmp=tmp-min(tmp(:));
  tmp=tmp/max(tmp(:));
  imwrite(tmp,['mipyz' xyz '.png'],'png');
  
  tmp=findobj('tag','xz');
  tmp=get(tmp,'children');
  tmp=tmp(ismember(tmp,findobj('type','image')));
  tmp=get(tmp,'cdata');
  tmp(z,:)=min(tmp(:));
  tmp(:,x)=min(tmp(:));  
  tmp=tmp-min(tmp(:));
  tmp=tmp/max(tmp(:));
  imwrite(tmp,['mipxz' xyz '.png'],'png');

return