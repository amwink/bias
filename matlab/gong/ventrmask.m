% File: ventrmask.m
%
% making a standard hull of a (standard space) brain
% pre:  the file standard.nii is present in the current 
% directory
% post: the file standardhull.nii contans the hull of 
% the brain
% 
% this image is thresholded from the outside inwards
% until the first above-threshold voxel is found, i.e.,
% it does not remove voxels inside the hull.
%
% (c) Alle Meije Wink 2015
% a.m.wink@gmail.com

N=nifti('standard.nii');
n=N.dat(:,:,:);

% make mask from outer brain contour
sn=size(n);
m=zeros(sn);

thresh=4000;

% fill inf to sup
for x=1:sn(1)
  for y=1:sn(2)
    z=1;
    while((z<sn(3))&(m(x,y,z)>=0)&(n(x,y,z)<thresh))
      m(x,y,z)=-1;
      z=z+1;
    end
  end
end

% fill sup to inf
for x=1:sn(1)
  for y=1:sn(2)
    z=sn(3);
    while((z>0)&(m(x,y,z)>=0)&(n(x,y,z)<thresh))
      m(x,y,z)=-1;
      z=z-1;
    end
  end
end

% fill ant to post
for x=1:sn(1)
  for z=1:sn(3)
    y=1;
    while((y<sn(2))&(m(x,y,z)>=0)&(n(x,y,z)<thresh))
      m(x,y,z)=-1;
      y=y+1;
    end
  end
end

% fill post to ant
for x=1:sn(1)
  for z=1:sn(3)
    y=sn(2);
    while((y>0)&(n(x,y,z)<thresh))
      m(x,y,z)=-1;
      y=y-1;
    end
  end
end

% fill left to right
for y=1:sn(2)
  for z=1:sn(3)
    x=1;
    while((x<sn(1))&(n(x,y,z)<thresh))
      m(x,y,z)=-1;
      x=x+1;
    end
  end
end

% fill right to left
for y=1:sn(2)
  for z=1:sn(3)
    x=sn(1);
    while((x>0)&(n(x,y,z)<thresh))
      m(x,y,z)=-1;
      x=x-1;
    end
  end
end

N.dat.fname='standardhull.nii';
N.dat(:,:,:)=(m+1);
create(N);