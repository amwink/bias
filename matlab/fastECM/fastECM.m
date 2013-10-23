function fastECM( inputfile, rankmap, normmap, degmap, maxiter, maskfile, atlasfile )
%
% function fastECM(string <inputfile>, bool <rankmap>, bool <normmap>)
% where - inputfile is the name of a 4D fMRI
%         image time series in the NifTI format
%       - rankmap ~=0 produces a file rankECM.nii
%         with uniformly distributed [0,1] centrality (tied) ranks
%       - normmap ~=0 produces a file normECM.nii
%         with normally distributed N(0,1) centralities generated from ranks
%       - degmap ~=0 produces a node network power map (similar to degree)
%         this is the column sum of connection strengths
%       - maxiter limits the number of iterations of the algorithm
%       - maskfile selects voxels inside a mask
%       - atlasfile groups signals by pre-defined regions
%
% returns a file fastECM.nii in the
%    same directory as <inputfile>
% which contains a 'fast eigenvector centrality mapping'
% whole-brain voxelwise functional connectivity analysis
%
% centrality distributions can be reshaped using <rankmap> and <normmap>
%
% example:
% >> fastECM('/tmp/fmri4d.nii.gz',1,1,1,16);
% produces files
%    /tmp/fmri4d_fastECM.nii.gz
%    /tmp/fmri4d_rankECM.nii.gz
%    /tmp/fmri4d_normECM.nii.gz
%    /tmp/fmri4d_degCM.nii.gz
% and iterates 16 times at most
%
% (c) Alle Meije Wink -- 16/03/2012
%     a.m.winkATgmail.com
%
% If you use this method in your research, remember to cite this paper
% in the journal "Brain Connectivity":
%
% Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
% "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI:
%  implementation, validation and interpretation."
% Brain Connectivity 2012, Vol. 2 No. 5, pages 265-274
% URL: http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087
%

if (~nargin)
    
    % If called w/o agruments,
    % use demo file fmri4d.nii.gz in fastECM's own directory
    
    inputfile = [fileparts(which('fastECM')) filesep 'fmri4d.nii.gz'              ];
    maskfile  = [fileparts(which('fastECM')) filesep 'mask_csf.nii.gz'            ];
    atlasfile = [fileparts(which('fastECM')) filesep 'aal_MNI_V4_4mm_gong.nii.gz' ];
    
    if (exist(inputfile)~=2)
        
        inputfile=inputfile(1:(end-3));
        if (exist(inputfile)~=2)
            fprintf('warning: demo file %s[.gz] not found, exiting\n',inputfile);
            return;
        end
        
    end
    
    fprintf(['\nsyntax:\nfastECM(\n\t<str  filename>,\n\t<bool rankmap>,  <bool normmap>, '...
        '<bool degmap>,\n\t<int  maxiter>,\n\t<str  maskname>, <str atlasname>\n\t\b)\n\n']);
    
    % demo call that uses all options
    
    fprintf('inputfile = ''%s'';\n',inputfile);
    fprintf('maskfile  = ''%s'';\n',maskfile);
    fprintf('atlasfile = ''%s'';\n',atlasfile);
    
    democall='fastECM(inputfile,1,1,1,15,maskfile,atlasfile);';
    fprintf('\n%% demo call:\n>> %s\n\n',democall);
    eval(democall);
    fprintf('\n');
    
    return;
    
end % if nargin

% check for options for masking, region selection and max. iterations
if (nargin<7)
    atlasfile=0;
    if (nargin<6)
        maskfile=0;
        if (nargin<5)
            maxiter=99;
        end % if
    end % if
end %if

% check if nifti support exists, otherwise add supplied version

if (exist('load_untouch_nii')~=2)
    
    fprintf('Software for reading NifTI images not found.        \n');
    fprintf('Using Tools for Nifti/Analyze for NifTI file I / O. \n');
    fprintf('  www.mathworks.com/matlabcentral/fileexchange/8797 \n');
    usenifti=1;                    % use included nifti libraries
    npath=[fileparts(which('fastECM')) filesep 'tools4nifti'];
    addpath(npath);                % (by adding them to the path)
    
else
    
    usenifti=0;                    % if not already on system
    
end % if exist

% load NifTI input file

fprintf('reading %s ...\n',inputfile);

if (inputfile(1)~=filesep)       % if path does not start with a separator
    
    if (    (isunix) | ...         % and in windos not with a drive letter or '\\'
            ( (~isunix) & (inputfile(2)~='\') & (inputfile(2)~=':') ) ...
            )                           % make it an absolute path
        inputfile=[pwd filesep inputfile];
    end % if isunix
    
end % if inputfile

M=load_untouch_nii(inputfile);           % read the information from the NifTI file
m=M.img;                         % get the 4d voxel data
M.img=[];                        % empty img after reading
M=rmfield(M,'img');              % and remove from M
msz=size(m);
tln=msz(4);                      % store the time series length
msz=msz(1:3);                    % store the size of one image volume
m(find(~isfinite(m)))=0;         % get rid of NaNs

% check whether a mask filename has been given and if dimensions are OK

msk=1;

if ( maskfile ~=0 )
    
    if (maskfile(1)~=filesep)      % if path does not start with a separator
        
        if (    (isunix) | ...       % and in windos not with a drive letter or '\\'
                ( (~isunix) & (maskfile(2)~='\') & (maskfile(2)~=':') ) ...
                )                    % make it an absolute path
            maskfile=[pwd filesep maskfile];
        end % if isunix
        
    end % if maskfile
    
    Msk=load_untouch_nii(maskfile);        % read the information from the NifTI file
    msk=double(Msk.img);           % get the 3d voxel data
    Msk.img=[];                    % empty img after reading
    
    if (size(msk) ~= msz)          % test if mask is the right size
        
        warning (sprintf(['dimensions of %s incompatible with %s\n' ...
            'continuing with nonzeroes of first volume'], ...
            maskfile, inputfile ));
        mskfile=0;
        msk=1;
        
    end % if size
    
end % if maskfile

% if not then make mask based on the positions of nonzeros in the 1st volume

msk=msk.*squeeze(min(m,[],4));   % take the minimum of each timeseries
msk=(msk>0);                     % make mask binary

% if atlas is used, generate regional time series

atl=1;

if ( atlasfile ~=0 )
    
    if (atlasfile(1)~=filesep)     % if path does not start with a separator
        
        if (    (isunix) | ...       % and in windos not with a drive letter or '\\'
                ( (~isunix) & (atlasfile(2)~='\') & (atlasfile(2)~=':') ) ...
                )                    % make it an absolute path
            atlasfile=[pwd filesep atlasfile];
        end % if isunix
        
    end % if atlasfile
    
    Atl=load_untouch_nii(atlasfile);       % read the information from the NifTI file
    atl=msk.*double(Atl.img);      % get the 3d voxel data of atlas regions inside the mask
    Atl.img=[];                    % empty img after reading
    
    if (size(atl) ~= msz)          % test if mask is the right size
        
        warning (sprintf(['dimensions of %s incompatible with %s\n' ...
            'continuing by using all voxel time series'], ...
            atlasfile, inputfile ));
        atlasfile=0;
        atl=1;
        
    else
        
        msk=(atl>0);
        
    end % if size
    
end % if atlasfile

% make time series 2D [tln msk]

m=reshape(m,[prod(msz) tln])';   % reshape the 4d voxels as 2d: [tln msz]
m=double(m(:,find(msk)));        % continue with only the nonzero voxels

% if atlasfile found -> make regional time series
if ( ~atlasfile )
    
    np=size(m,2);                  % store the number of nonzero voxels
    
else
    
    mreg=max(atl(:));              % highest region label
    atl=atl(find(msk));            % atl now the same size as msk
    reg=unique(atl);               % region labels: nonzero values found inside mask/atl
    rsiz=length(reg);              % number of region labels found
    
    mm=zeros(tln,rsiz);
    for r=1:rsiz                   % construct regional means
        mm(:,r) = mean( m(:,find(reg==reg(r))),2 );
    end
    m=mm;                          % continue with regional mean time series
    clear mm;
    
    np=r;                          % store the number of included regiosnfac
    
end % if ~atlasfile

% compute mean and var

mav=mean(m);                     % compute the time series mean
% mvr=var(m);
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

% efficient power iteration

while ((iter<maxiter) & (dnorm>cnorm))
    
    vprev=vcurr;                   % start with previous estimate
    prevsum=sum(vprev);            % sum of estimate
    vcurr_1=m*vprev;               % part one of M*v
    vcurr_2=m'*vcurr_1;            % part two of M*v
    vcurr_3=vcurr_2+prevsum;       % adding sum -- same effect as [M+1]*v
    vcurr=vcurr_3/norm(vcurr_3,2); % normalise L2-norm
    
    if (~iter)
        
        if ( (nargin>3) & (degmap ~= 0) )
            dvcurr=vcurr_2;
            write_map(inputfile, M, msk, atl, vcurr_2, 'degCM', 'weighted degree centrality');
        end % if nargin
        
    end % if iter
    
    iter=iter+1;                   % increase iteration counter
    dnorm=norm(vcurr-vprev,2);     % L2-norm of difference prev-curr estimate
    cnorm=norm(vcurr,2)*eps;       % L2-norm of current estimate
    fprintf('iteration %02d, || v_i - v_(i-1) || / || v_i * epsilon || = %0.16f / %0.16f\n', ...
        iter,dnorm,cnorm)          % some stats for the users
    
end % while

% write the eigenvactor centrality map to the file fastECM.nii

write_map(inputfile, M, msk, atl, vcurr, 'fastECM', sprintf('ECM [iterations: %d]',iter) );

% produce a map of EC ranks if requested and write to rankECM.nii

if ( ( (nargin>1) & (rankmap ~= 0) ) | ( (nargin>2) & (normmap~=0) ) )
    
    rvcurr=tiedrank(vcurr);        % tied ranks: equal values lead to equal ranks
    rvcurr=(rvcurr)/(length(rvcurr)+1);
    % go from uniform [1,N] to uniform ]0,1[
    if (rankmap ~= 0)
        write_map(inputfile, M, msk, atl, rvcurr, 'rankECM', 'ECM [converted to <0,1> ranks]');
    end
    
    % produce a map of gaussianised EC ranks if requested and write to normECM
    
    if ( (nargin>2) & (normmap ~= 0) )
        
        mu=0;                        % mean
        sig=1;                       % std dev
        
        nvcurr=mu+sqrt(2)*sig*erfinv(2*rvcurr-1);
        % uniform ]0,1[ to N(0,1) via
        % inverse transform sampling
        write_map(inputfile, M, msk, atl, nvcurr, 'normECM', 'ECM [converted to N(0,1) values]');
        
    end % if nargin2
    
end % if nargin1

% for atlas-based analyses, create a text file of regional values

if (prod(size(atl)) ~= 1)
    
    regs=1:mreg;
    regs(end+1,reg)=vcurr;
    
    if(exist('dvcurr'))
        regs(end+1,reg)=dvcurr;
    end % if
    
    if(exist('rvcurr'))
        regs(end+1,reg)=rvcurr;
    end % if
    
    if(exist('nvcurr'))
        regs(end+1,reg)=nvcurr;
    end % if
    
    [fd,fn,fx2]=fileparts(inputfile);
    [fd1,fn,fx]=fileparts(fn);     % base file name and 1st extension (2nd may be .gz)
    fx=[fx fx2];                   % in which case, concatenate by adding fx2
    txtfile=[fd filesep fn '_fastECM.txt'];
    delete(txtfile);
    dlmwrite(txtfile,regs,'precision',4,'delimiter','\t')
    txtfile=[fd filesep fn '_reg_avg.txt'];
    delete(txtfile);
    dlmwrite(txtfile,[reg';m],'precision',4,'delimiter','\t')
    
end % if prod

% if nifti support was added before, remove it again to not leave prints

if (usenifti)
    rmpath(npath);                 % was own nifti added to path -- if yes then remove now
end % if

return % fastECM

function write_map(inputfile, M, msk, atl, vcurr, mapfile, descrip)
%
% writes a map based on a template file
% - inputfile = filename of the template
% - M         = nifti record of the template
% - msk       = brain mask of the template
% - atl       = atlas volume (if used, otherwise 1)
% - vcurr     = map values inside brain mask
% - mapfile   = filename of the map file
% - descrip   = description of contents
%
% this routine uses the following external code
%
% " Tools for Nifti / Analyze " by Jimmy Shen
%   http://www.mathworks.com/matlabcentral/fileexchange/8797
% " quantile.m " by Anders Holtsberg
%   www.spatial-econometrics.com/distrib/quantile.m
%

% if atlas used, make a map of regional centralities

if (prod(size(atl)) ~= 1)
    
    vcurr_2=zeros(size(atl));
    reg=unique(atl);               % region labels: nonzero values found inside atl
    
    for r=1:length(reg)
        nvox=find(atl==reg(r));
        vcurr_2(nvox)=vcurr(r)/(length(nvox));
    end % for
    
    vcurr=vcurr_2;
    clear vcurr_2;
    
end % if prod

% create output array mout

mout=zeros(size(msk));           % initialise output with 0 (background)
mout(find(msk))=vcurr;           % fill the nonzero voxels with map values
qnts=quantiles(vcurr(vcurr>0),[.01 .99]);    % quantiles to set contrast

% write output nifti file based on mout

[fd,fn,fx2]=fileparts(inputfile);
[fd1,fn,fx]=fileparts(fn);       % base file name and 1st extension (2nd may be .gz)
fx=[fx fx2];                     % in which case, concatenate by adding fx2
outputfile=[fd filesep fn '_' mapfile fx];
fprintf('writing %s ...\n',outputfile);
M.hdr.hist.descrip=sprintf('generated by fastECM - %s',descrip);
% leave our name tag in the nifti record
M.hdr.dime.dim(5)=1;             % map is 3D not 4D
M.hdr.dime.pixdim(5)=0;          % so are its voxels
M.hdr.dime.cal_min=qnts(1);      % min of range (for win/lev)
M.hdr.dime.cal_max=qnts(2);      % max
M.hdr.dime.glmin=qnts(1);        % min of range
M.hdr.dime.glmax=qnts(2);        % max
M.img=mout;                      % add voxel data to map
save_untouch_nii(M,outputfile);          % write the file

return
