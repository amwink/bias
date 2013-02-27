function fastECM(inputfile)
%
% function fastECM(string <inputfile>)
%  - where inputfile is the name of a 4D fMRI 
%    image time series in the NifTI format
%  
% returns a file fastECM.nii in the 
%    same directory as <inputfile>
%
% which contains a 'fast eigenvector centrality mapping'
% whole-brain voxelwise functional connectivity analysis
%
% example:
% >> fastECM('/tmp/fmri4d.nii');
% produces a file
%    /tmp/fastECM.nii
%
% (c) Alle Meije Wink -- 16/03/2012
%     a.m.winkATgmail.com
%
% If you use this method in your research, remember to cite this paper
% in the journal "Brain Connectivity":
% 
% Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
% "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI: implementation, validation and interpretation."
% URL: http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087
%

% use a local file if no input filename is supplied

if (~nargin)                  
  inputfile=[fileparts(which('fastECM')) filesep 'fmri4d.nii'];
end

% check if nifti support exists, otherwise add supplied version

if (exist('nifti')~=2)
  fprintf('Software for reading NifTI images not found.\n');
  fprintf('Using niftimatlib-1.2 for NifTI file I / O. \n');
  fprintf('      (see http://niftilib.sourceforge.net) \n');
  usenifti=1;                    % use included nifti libraries 
  npath=[fileparts(which('fastECM')) filesep 'niftimatlib-1.2' filesep 'matlab'];
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
M=nifti(inputfile);              % read the information from the NifTI file
m=M.dat(:,:,:,:);                % read the 4d voxel data
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

m=m(:,find(msk));                % continue with only the nonzero voxels
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
while ((iter<100) & (dnorm>cnorm))
  vprev=vcurr;                   % start with previous estimate
  prevsum=sum(vprev);            % sum of estimate
  vcurr_1=m*vprev;               % part one of M*v
  vcurr_2=m'*vcurr_1;            % part two of M*v
  vcurr_3=vcurr_2+prevsum;       % adding sum -- same effect as [M+1]*v
  vcurr=vcurr_3/norm(vcurr_3,2); % normalise L2-norm
  iter=iter+1;                   % increase iteration counter
  dnorm=norm(vcurr-vprev,2);     % L2-norm of difference prev-curr estimate
  cnorm=norm(vcurr*eps,2);       % L2-norm of current estimate
  fprintf('iteration %02d, || v_i - v_(i-1) || / || v_i * epsilon || = %0.16f / %0.16f\r', ...
	  iter,dnorm,cnorm)      % some stats for the users
end
fprintf('\n');

% create output array mout

mout=zeros(msz);                 % initialise output with 0 (background)
mout(find(msk))=vcurr;           % fill the nonzero voxels with ECM values
qnts=quantiles(vcurr(vcurr>0),[.05 .95]);  
                                 % quantiles to set contrast                               
                                 
% write output nifti file based on mout

outputfile=[fileparts(M.dat.fname) filesep 'fastECM.nii'];
fprintf('writing %s ...\n',outputfile);
Mout=nifti;                      % empty nifti file record
Mout.mat=M.mat;                  % copy coordinate transforms from input
Mout.mat0=Mout.mat;              % extra mat0
Mout.dat=file_array(outputfile,size(mout),'float32',352,1,0);
                                 % initialise file array                               
Mout.dat(:,:,:)=mout;            % put the output data in the file array
Mout.cal=[qnts];                 % set min and max
create(Mout);                    % write the file

% if nifti support was added before, remove it again to not leave prints

if (usenifti)
  rmpath(npath);                 % was own nifti added to path, remove now
end

return
