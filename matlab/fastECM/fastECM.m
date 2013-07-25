function fastECM( inputfile, rankmap, normmap, maxiter )
%
% function fastECM(string <inputfile>, bool <rankmap>, bool <normmap>)
% where - inputfile is the name of a 4D fMRI 
%         image time series in the NifTI format
%       - rankmap ~=0 produces a file rankECM.nii
%         with uniformly distributed [0,1] centrality (tied) ranks
%       - normmap ~=0 produces a file normECM.nii
%         with normally distributed N(0,1) centralities generated from ranks
%       - maxiter limits the number of iterations of the algorithm     
%  
% returns a file fastECM.nii in the 
%    same directory as <inputfile>
% which contains a 'fast eigenvector centrality mapping'
% whole-brain voxelwise functional connectivity analysis
% 
% centrality distributions can be reshaped using <rankmap> and <normmap>
%
% example:
% >> fastECM('/tmp/fmri4d.nii.gz',1,1,16);
% produces files
%    /tmp/fmri4d_fastECM.nii.gz
%    /tmp/fmri4d_rankECM.nii.gz
%    /tmp/fmri4d_normECM.nii.gz
% and iterates 16 times at most
%
% (c) Alle Meije Wink -- 16/03/2012
%     a.m.winkATgmail.com
%
% If you use this method in your research, remember to cite this paper
% in the journal "Brain Connectivity":
% 
% Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
% "Fast eigenvector centrality mapping of voxel-wise connectivity
% in functional MRI: implementation, validation and
% interpretation."
% Brain Connectivity 2012, Vol. 2 No. 5, pages 265-274
% URL: http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087
%

% use a local file if no input filename is supplied

if (~nargin)                  
  inputfile=[fileparts(which('fastECM')) filesep 'fmri4d.nii.gz'];
end

if (nargin<4)
  maxiter=100;
end

% check if nifti support exists, otherwise add supplied version

if (exist('load_nii')~=2)
  fprintf('Software for reading NifTI images not found.        \n');
  fprintf('Using Tools for Nifti/Analyze for NifTI file I / O. \n');
  fprintf('  www.mathworks.com/matlabcentral/fileexchange/8797 \n');
  usenifti=1;                    % use included nifti libraries 
  npath=[fileparts(which('fastECM')) filesep 'tools4nifti'];
  addpath(npath);                % (by adding them to the path)
else
  usenifti=0;                    % if not already on system
end

% load NifTI input file

fprintf('reading %s ...\n',inputfile);
if (inputfile(1)~=filesep)       % if path does not start with a separator
  if (    (isunix) | ...         % and in windos not with a drive letter or '\\'
       ( (~isunix) & (inputfile(2)~='\') & (inputfile(2)~=':') ) ...
     )                           % make it an absolute path
    inputfile=[pwd filesep inputfile];    
  end
end
M=load_nii(inputfile);           % read the information from the NifTI file
m=M.img;                         % get the 4d voxel data
M.img=[];                        % empty img after reading
M=rmfield(M,'img');              % and remove from M
tln=size(m,4);                   % store the time series length
m(find(~isfinite(m)))=0;         % get rid of NaNs

% make mask based on the positions of nonzeros in the 1st volume

msk=squeeze(m(:,:,:,1));         % get the 1st volume and make it 3D
msz=size(msk);                   % store the size of one image volume
msk=(msk>0);                     % make it a binary mask

% make mask 1d (:) and time series 2D [tln msk]

msk=msk(:);                      % re-shape the mask as 1D 
m=reshape(m,[prod(msz) tln])';   % reshape the 4d voxels as 2d: [tln msz]

% only select m inside mask

m=double(m(:,find(msk)));        % continue with only the nonzero voxels
np=size(m,2);                    % store the number of nonzero voxels

% compute mean and var

mav=mean(m);                     % compute the time series mean
mvr=std(m);                      % compute the time series variance
mvr=mvr+eps;                     % prevent divisions by 0

% prepare matrix

m=(m-(ones(tln,1)*mav))./(ones(tln,1)*mvr); 
                                 % mean 0, std 1 to make  covariance matrix
m=m/sqrt(tln-1);                 % make correlations instead (diagonal 1)

% initialise eigenvector estimate v

vprev=0;                         % 'initialise' previous ECM estimate
vcurr=ones(np,1)/sqrt(np);       % initialise estimate with L2-norm == 1

iter=0;                          % reset iteration counter
dnorm=1;                         % initial value for difference L2-norm
cnorm=0;                         % initial value for estimate L2-norm
while ((iter<maxiter) & (dnorm>cnorm))
  vprev=vcurr;                   % start with previous estimate
  prevsum=sum(vprev);            % sum of estimate
  vcurr_1=m*vprev;               % part one of M*v
  vcurr_2=m'*vcurr_1;            % part two of M*v
  vcurr_3=vcurr_2+prevsum;       % adding sum -- same effect as [M+1]*v
  vcurr=vcurr_3/norm(vcurr_3,2); % normalise L2-norm
  iter=iter+1;                   % increase iteration counter
  dnorm=norm(vcurr-vprev,2);     % L2-norm of difference prev-curr estimate
  cnorm=norm(vcurr*eps,2);       % L2-norm of current estimate
  fprintf('iteration %02d, || v_i - v_(i-1) || / || v_i * epsilon || = %0.16f / %0.16f\n', ...
	  iter,dnorm,cnorm)          % some stats for the users
end

% write the eigenvactor centrality map to the file fastECM.nii

write_map(inputfile, M, msk, vcurr, 'fastECM');

% produce a map of EC ranks if requested and write to rankECM.nii

if ( ( (nargin>1) & (rankmap ~= 0) ) | ( (nargin>2) & (normmap~=0) ) )
    
  rvcurr=tiedrank(vcurr);        % tied ranks: equal values lead to equal ranks
  rvcurr=(rvcurr)/(length(rvcurr)+1);
                                 % go from uniform [1,N] to uniform ]0,1[
  if (rankmap ~= 0) 
    write_map(inputfile, M, msk, rvcurr, 'rankECM');
  end

  % produce a map of gaussianised EC ranks if requested and write to normECM

  if ( (nargin>2) & (normmap ~= 0) )
   
    mu=0;                        % mean
    sig=1;                       % std dev
      
    nvcurr=mu+sqrt(2)*sig*erfinv(2*rvcurr-1);
                                 % uniform ]0,1[ to N(0,1) via 
                                 % inverse transform sampling
    write_map(inputfile, M, msk, nvcurr, 'normECM');
    
  end
  
end

% if nifti support was added before, remove it again to not leave prints

if (usenifti)
  rmpath(npath);                 % was own nifti added to path, remove now
end

return

function write_map(inputfile, M, msk, vcurr, mapfile)
%
% writes a map based on a template file
% - inputfile = filename of the template
% - M         = nifti record of the template
% - msk       = brain mask of the template
% - vcurr     = map values inside brain mask
% - mapfile   = filename of the map file
%
% this routine uses the following external code
% 
% " Tools for Nifti / Analyze " by Jimmy Shen
%   http://www.mathworks.com/matlabcentral/fileexchange/8797
% " quantile.m " by Anders Holtsberg
%   www.spatial-econometrics.com/distrib/quantile.m
%

% create output array mout
mout=zeros(size(msk));           % initialise output with 0 (background)
mout(find(msk))=vcurr;           % fill the nonzero voxels with map values
qnts=quantiles(vcurr(vcurr>0),[.05 .95]);  
                                 % quantiles to set contrast                               
                                 
% write output nifti file based on mout
[fd,fn,fx2]=fileparts(inputfile);
[fd1,fn,fx]=fileparts(fn);       % base file name and 1st extension (2nd may be .gz)
fx=[fx fx2];                     % in which case, concatenate by adding fx2
outputfile=[fd filesep fn '_' mapfile fx];
fprintf('writing %s ...\n',outputfile);
M.hdr.hist.descrip='fastECM';    % leave our name tag in the nifti record
M.hdr.dime.dim(5)=1;             % map is 3D not 4D
M.hdr.dime.pixdim(5)=0;          % so are its voxels
M.hdr.dime.cal_min=qnts(1);      % min of range (for win/lev)
M.hdr.dime.cal_max=qnts(2);      % max
M.hdr.dime.glmin=qnts(1);        % min of range
M.hdr.dime.glmax=qnts(2);        % max
M.img=mout;                      % add voxel data to map
save_nii(M,outputfile);          % write the file

return
