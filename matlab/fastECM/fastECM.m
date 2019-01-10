function err = fastECM ( inputfile, rankmap, normmap, degmap, maxiter, maskfile, atlasfile, wholemat, dynamics, correlate ) 
%
% function fastECM ( string <inputfile>, 
%          bool <rankmap>, bool <normmap>, bool <degmap>, 
%          int maxiter, 
%          string <maskfile>, string <atlasfile>, 
%          bool wholemat, 
%          int dynamics ) 
% where - inputfile is the name of a 4D fMRI
%         image time series in the NifTI format
%       - rankmap ~= 0 produces a file rankECM.nii
%         with uniformly distributed [0, 1] centrality (tied) ranks
%       - normmap ~= 0 produces a file normECM.nii
%         with normally distributed N (0, 1) centralities generated from ranks
%       - degmap ~= 0 produces a node power map (similar to degree) 
%         this is the column sum of connection strengths
%       - maxiter limits the number of iterations of the algorithm
%       - maskfile selects voxels inside a mask
%       - atlasfile groups signals by pre-defined regions
%       - wholemat: do not use the 'fast' recipe, write out
%         connectivity matrix and compute graph measures (costly!) 
%       - dynamics: number of ECM to extract from a 4D fMRI
%       - correlate: way to approximate the eigenvector of 
%         correlated time series: 
%           'fast' -- default, correlation + 1
%           'relu' -- ReLU, rectified linear unit
%
% returns a file fastECM.nii in the same directory as <inputfile>
%     which contains a 'fast eigenvector centrality mapping'
%     whole-brain voxelwise functional connectivity analysis;
%     error code 0 if all is well, 1 otherwise
%
%
%
% Example 1
%
%     If the following files are present:
%        fmri4d.nii.gz       -- containing an fMRI time series
%        mask_csf.nii.gz      -- mask with non-brain tissue and CSF removed
%        aal_MNI_V4_4mm_gong.nii.gz -- containing a volume with atlas labels
%     the call
%        >> fastECM;
%     runs a demo analysis highlighting all these options
%
% Example 2
%
%     For using some, but not all options, the easiest way to
%     call fastECM is to pass these options in a struct:
%        >> op.inputfile = 'fmri4d.nii.gz';
%        >> op.dynamics = 25;
%        >> fastECM (op);
%
%
%
% (c) Alle Meije Wink -- 16/03/2012
%   a.m.wink AT gmail.com
%
%
%
% If you use this method in your research, remember to cite this paper
% in the journal "Brain Connectivity":
%
% Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
%     "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI:
%      implementation, validation and interpretation."
%      Brain Connectivity 2012, Vol. 2 No. 5, pages 265-274
% URL: http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087 
% or:  http://dare.ubvu.vu.nl/handle/1871/48750
%
%
%



% If called w/o agruments, produce a demo call
if (~nargin)  
  
  inputfile = [fileparts (which ('fastECM.m') ) filesep 'fmri4d.nii.gz'    ];
  maskfile  = [fileparts (which ('fastECM.m') ) filesep 'mask_csf.nii.gz'   ];
  atlasfile = [fileparts (which ('fastECM.m') ) filesep 'aal_MNI_V4_4mm_gong.nii.gz' ];
  
  if ( (exist (inputfile) *exist (maskfile) *exist (atlasfile) ) ~= 8) 
    
    inputfile = inputfile (1: (end-3) );
    maskfile  = maskfile  (1: (end-3) );
    atlasfile = atlasfile (1: (end-3) );
    
    if (exist (inputfile) ~= 2) 
      
      fprintf ('warning: not all demo files found, exiting\n');

      err = 1;
      return;
      
    end % if ~input files (uncompressed) 
    
  end % if ~inputfile or ~maskfile or ~atlasfiles
  
  fprintf (['\nsyntax:\nfastECM (\n\t<str filename>, \n\t<bool rankmap>, <bool normmap>, ' ...
	    '<bool degmap>, \n\t<int maxiter>, \n\t<str maskname>, <str atlasname>\n\t' ...
	    '<bool wholemat>\n\t<int dynamics>\n\t<str correlate>\n\t\b) \n\n']);
  
  correlate = 'fast';
  
  fprintf ('inputfile = ''%s'';\n', inputfile);
  fprintf ('maskfile = ''%s'';\n', maskfile);
  fprintf ('atlasfile = ''%s'';\n', atlasfile);
  fprintf ('correlate = ''%s'';\n', correlate);
  
  democall = ['fastECM ( inputfile, 1, 1, 1, 20, maskfile, atlasfile, 0, 25, correlate );'];
  fprintf ('\n%% demo call:\n>> %s\n\n', democall);
  err = eval (democall);
  fprintf ('\n'); 
  
  return; % call fastECM inside this call
  
else
  
  err = 0; % continue without errors 
  
end % if ~nargin



% check for options for no. of dynamics, 
%      whole matrix computation, 
%      atlas region selection, 
%      mask selection
%      max. iterations
%      write degree (node power) map
%      write normally distributed map
%      write ranks (uniformly distributed) map
% and set defaults
if (nargin<10) 
  correlate = 'fast';
  if (nargin<9) 
    dynamics = 1;
    if (nargin<8) 
      wholemat = 0;
      if (nargin<7) 
	atlasfile = 0;
	if (nargin<6) 
	  maskfile = 0;
	  if (nargin<5) 
	    maxiter = 25;
	    if (nargin<4) 
	      degmap = 0;
	      if (nargin<3) 
		normmap = 0;
		if (nargin<2) 
		  rankmap = 0;		
		end  % if rankmap
	      end   % if normmap
	    end	   % if degmap
	  end    % if maxiter
	end    % if maskfile
      end     % if atlasfile
    end     % if wholemat
  end      % if dynamics
end      % if correlate 



% if options were given as a struct
% use the existing fields
if ( isstruct (inputfile) ) 
  
 input_parameters = inputfile;
 fieldnames = { 'rankmap' 'normmap' 'degmap' 'maxiter' 'maskfile' 'atlasfile' 'wholemat' 'dynamics' 'correlate' 'inputfile' };
 
 for f = 1:length (fieldnames) 
   
   if isfield ( input_parameters, fieldnames{f} ) 
     eval ([ fieldnames{f} ' = input_parameters.' fieldnames{f} ';' ]);
   end % if isfield
   
 end % for f
 
 if isstruct (inputfile) 
   fprintf ('error: ''inputfile'' not defined in struct, exiting\n');
   err = 1;
   return;
 end % if isstruct
 
end % if isstruct



% check if nifti support exists
% otherwise add fastECM-supplied version
if (exist ('load_untouch_nii') ~= 2) 
  
  fprintf ('Software for reading NifTI images not found.  \n');
  fprintf ('Using Tools for Nifti/Analyze for NifTI file I / O. \n');
  fprintf (' www.mathworks.com/matlabcentral/fileexchange/8797 \n');
  usenifti = 1;     
  npath = [fileparts(which('fastECM.m')) filesep 'tools4nifti'];
  addpath (npath);     
  
else 
  
  usenifti = 0;      
  
end % if exist



% load NifTI input file
fprintf ('reading %s ...\n', inputfile);

if (inputfile (1) ~= filesep)    % if path does not start with a separator
  
  if (  (isunix) | ...           % and -in windows- not with a drive letter or '\\'
	( (~isunix) & (inputfile (2) ~= '\') & (inputfile (2) ~= ':') ) ...
	)                        % make it an absolute path
    
    inputfile = [pwd filesep inputfile];
    
  end % if isunix
  
end % if inputfile

% get the directory, base name and extension (s) 
[fd, fn, fx2] = fileparts (inputfile);  % dir and base file name + extension (may be combined) 
[fd1, fn, fx] = fileparts (fn);         % base file name and 1st extension (2nd may be .gz) 

% get the data and clear data from file record
M = load_untouch_nii (inputfile);       % read the information from the NifTI file
m = double (M.img);                     % get the 4d voxel data
M.img = [];                             % empty img after reading
M = rmfield (M, 'img');                 % and remove from M
msz = size (m);
tln = msz (4);                          % store the time series length
msz = msz (1:3);                        % store the size of one image volume
m (find (~isfinite (m) ) ) = 0;         % get rid of NaNs

% check whether a mask filename has been given and if dimensions are OK
msk = 1;                                % default value for no user-supplied mask
if ( maskfile ~= 0 & ~isempty (maskfile) ) 
  
  if (maskfile (1) ~= filesep)          % if path does not start with a separator
    
    if (  (isunix) | ...                % and -in windows- not with a drive letter or '\\'
	  ( (~isunix) & (maskfile (2) ~= '\') & (maskfile (2) ~= ':') ) ...
	  )                             % make it an absolute path
      maskfile = [pwd filesep maskfile];
    end % if isunix
    
 end % if maskfile
 
 Msk = load_untouch_nii (maskfile);     % read the information from the NifTI file
 msk = double (Msk.img);                % get the 3d voxel data
 Msk.img = [];                          % empty img after reading
 
 if (size (msk) ~= msz)                 % test if mask volume is the right size
   
   warning (sprintf (['dimensions of %s incompatible with %s\n' ...
		      'continuing with nonzero time series.'], ...
		     maskfile, inputfile ) );
   mskfile = 0;
   msk = 1;
   
 end % if size
 
else
  
  mskfile = 0;
  
end % if maskfile

% if the mask provided is not OK then the positions of nonzeros in the 1st volume
msk = msk.*squeeze (abs (min (m, [], 4) ) );  % take the minimum of each timeseries
msk = (msk>0);                          % make mask binary

% if atlas is used, generate regional time series
atl = 1;                                % default value for no user-supplied atlas
if ( atlasfile ~= 0 & ~isempty (atlasfile) ) 
  
  if (atlasfile (1) ~= filesep)         % if path does not start with a separator
    
    if (  (isunix) | ...                % and in windos not with a drive letter or '\\'
	  ( (~isunix) & (atlasfile (2) ~= '\') & (atlasfile (2) ~= ':') ) ...
	  )                             % make it an absolute path
      atlasfile = [pwd filesep atlasfile];
    end % if isunix
    
  end % if atlasfile
  
  Atl = load_untouch_nii (atlasfile);   % read the information from the NifTI file
  atl = double (Atl.img);               % get the 3d voxel data of atlas regions
  Atl.img = [];                         % empty img after reading
  
  if (size (atl) ~= msz)                % test if atlas volume is the right size
    
    warning (sprintf (['dimensions of %s incompatible with %s\n' ...
		       'continuing by using all voxel time series'], ...
		      atlasfile, inputfile ) );
    atlasfile = 0;
    atl = 1;
    
  else
    
    msk = msk.* (atl>0);
    
  end % if size
  
else
  
  atlasfile = 0;
  
end % if atlasfile

fprintf ('\nloaded and masked\n');



% make time series 2D [tln msk]
m = reshape(m, [prod(msz) tln])' ;      % reshape the 4d voxels as 2d: [tln msz]
m = double(m (:, find (msk) ) );        % continue with only the nonzero voxels



% if atlasfile found -> make regional time series
if ( ~atlasfile ) 
  
  np = size (m, 2);                     % store the number of nonzero voxels
  
else
  
  fprintf ('atlasing\n');
  
  mreg = max (atl (:) );                % highest region label
  atl = atl (find (msk) );              % atl now the same size as msk
  reg = unique (atl);                   % region labels: nonzero values found inside mask/atl
  
  mm = zeros (tln, length (reg) );
  for r = reg (:)'                      % construct regional means  
    mm (:, r) = mean ( m (:, find (atl == r) ), 2 );
  end % for r
  m = mm;                               % continue with regional mean time series
  clear mm;
  
  txtfile = [fd filesep fn '_fastECMtseries.txt'];
  dlmwrite (txtfile, m, ' ');           % write regional mean time series to text file
  
  matfile = [fd filesep fn '_fastECMtseries.mat']; 
  save (matfile, 'm');
  
  np = r;                               % store the number of included regios
  
  wholemat = 1;                         % regional matrices are small -> provide all info
  
end % if ~atlasfile



% compute mean and var
mav = mean (m);                         % compute the time series mean
mvr = std (m);                          % compute the time series standard deviation
mvr = mvr+eps;                          % prevent divisions by 0



% prepare matrix
m = (m - (ones (tln, 1) *mav) ) ...
    ./ (ones (tln, 1) *mvr);            % mean 0, std 1 to make covariance matrix
m = m/sqrt (tln-1);                     % make correlations instead (diagonal 1) 



% compute ECM on intervals given by the number of dynamics
ddiff = tln-dynamics;
for d = 1:dynamics
  
  % initialise eigenvector estimate v 
  vprev = 0;                            % 'initialise' previous ECM estimate
  vcurr = ones (np, 1) /sqrt (np);      % initialise estimate with L2-norm == 1
  
  iter = 0;                             % reset iteration counter
  dnorm = 1;                            % initial value for difference L2-norm
  cnorm = 0;                            % initial value for estimate L2-norm
  
  m0 = m (d: (ddiff+d), :);             % use the interval for the current dynamic
  
  if (strcmp (correlate, 'relu') )      % use ReLU equation to guarantee positive M
    
    m0 = [m0; abs(m0) ];                % concatenate absolute value of time series
    
  end % if strcmp
  
  % efficient power iteration for correlated time series 
  while ( (iter<maxiter) & (dnorm>cnorm) ) 
    
    vprev = vcurr;                      % start with previous estimate
    prevsum = sum (vprev);              % sum of estimate
    
    if (~wholemat)                      % A. 'fast recipe' -> cunning re-ordering of computations
			
      vcurr_1 = m0*vprev;               % 1. part one of M*v
      vcurr_2 = m0'*vcurr_1;            % 2. part two of M*v
      
      % if ReLU correlation is used, then values are positive (modulo 
      % some precision errors). Otherwise, this needs to be enforced
      if (strcmp (correlate, 'relu') ) 
	
	vcurr_3 = vcurr_2 .* (vcurr_2>0); % 3. remove remaining small negative values
	
      else
	
	vcurr_3 = vcurr_2+prevsum;      % 3. adding sum -- same effect as [M+1]*v
	
      end % if strcmp	
      
    else                                % B. original recipe -> provide full connectivity info
		
      vcurr_1 = (m0'*m0);               % 1. M = correlations
      vcurr_2 = vcurr_1*vprev;          % 2. M * v

      % if ReLU correlation is used, then values are positive (modulo 
      % some precision errors). Otherwise, this needs to be enforced
      if (strcmp (correlate, 'relu') ) 		     
	
	vcurr_3 = vcurr_2+1;            % 3. remove remaining small negative values	

      else
	
	vcurr_3 = vcurr_2+1;            % 3. [M+1]*v
	
      end
      
      if (~iter)                        % compute & write graph measures only in iteration 1 
	
	vcurr_1 = vcurr_1-diag (diag (vcurr_1) );
	connmat_out (:, :, d) = vcurr_1;
	dist = 1./vcurr_1;              % distance matrix <-> 1/connectivity
	
	for c = 1:np                    % shortest distances www.ee.columbia.edu/~marios/matlab/tips.pdf
	  dist = min (dist, repmat (dist (:, c), [1 np]) +repmat (dist (c, :), [np 1]) );
	end % for c
	
	mpl = mean (dist (dist ~= inf) );  % average shortest path to other nodes -> global path length
	
	% threshold & binarise for other graph measures 
	lv = length (vcurr_1);
	indi = triu (reshape (1: (lv*lv), [lv lv]), 1);
	
	if (exist ('backbone_wu') ~= 2);  % get 'backbone' i.e. MST + some extra
	  
	  eval (['!wget -P ' fileparts (which ('fastECM.m') ) ...
		 ' -O backbone_wu.m ' ...
		 'sites.google.com/site/bctnet/Home/functions/backbone_wu.m']); 	 
	  
	end % if exist
	
	no0corr = vcurr_1;                             % make correlation matrix where 0-pairs have very low, non-zero correlation
	no0corr (~no0corr) = min (no0corr (:) ) -eps;  % otherwise the BCT script will not add these nodes to the MST
	no0corr = no0corr+1;
	no0corr (find (diag (diag (no0corr) ) ) ) = 0;
	
	% compute MST and add nodes up to degree sqrt (#nodes) to the 'backbone'
	[mst clus] = backbone_wu (no0corr, fix (sqrt (lv) ) );	
	vcurr_bin (:, :, d) = sign (clus);
	vcurr_mst (:, :, d) = sign (mst);
	
	% apply bct measures to vcurr_bin: communities (index per node) 
	if (exist ('community_louvain') ~= 2);
	  
	  eval (['!wget -P ' fileparts (which ('fastECM.m') ) ...
		 ' -O community_louvain.m ' ...    
		 'sites.google.com/site/bctnet/Home/functions/community_louvain.m ']); 	 
	  
	end % if exist
	communities (:, d) = community_louvain (vcurr_bin (:, :, d) );
	
	% apply bct measures to vcurr_bin: betweenness (per node) 
	if (exist ('betweenness_bin') ~= 2);
	  
	  eval (['!wget -P ' fileparts (which ('fastECM.m') ) ...
		 ' -O betweenness_bin.m ' ...
		 'sites.google.com/site/bctnet/Home/functions/betweenness_bin.m ']); 	 
	  
	end % if exist
	betweenness (:, d) = betweenness_bin (vcurr_bin (:, :, d) );
	
	% apply bct measures to vcurr_bin: clustering (per node) 
	if (exist ('clustering_coef_bu') ~= 2);
	  
	  eval (['!wget -P ' fileparts (which ('fastECM.m') ) ...
		 ' -O clustering_coef_bu.m ' ...
		 ' sites.google.com/site/bctnet/Home/functions/clustering_coef_bu.m ']); 	 
	  
	end % if exist
	clustering (:, d) = clustering_coef_bu (vcurr_bin (:, :, d) );
	
	% apply bct measures to vcurr_bin: path length (per node) 
	if (exist ('distance_bin') ~= 2);
	  
	  eval (['!wget -P ' fileparts (which ('fastECM.m') ) ...
		 ' -O distance_bin.m ' ...
		 ' sites.google.com/site/bctnet/Home/functions/distance_bin.m ']); 	 
	  
	end % if exist
	pathlengthtmp = distance_bin (vcurr_bin (:, :, d) );
	pathlength (:, d) = mean (pathlengthtmp);
	
      end % if ~iter
      
    end % if wholemat 
    
    vcurr = vcurr_3/norm (vcurr_3, 2);  % normalise L2-norm
    
    if ( (~iter) & (degmap ~= 0) ) 
      
      dvcurr (:, d) = vcurr_2 (:);      % save 'node power' of this dynamic if requested
    end % if degmap
    
    iter = iter+1;                      % increase iteration counter
    dnorm = norm (vcurr-vprev, 2);      % L2-norm of difference prev-curr estimate
    cnorm = norm (vcurr, 2) *eps;       % L2-norm of current estimate
    fprintf ('dynamic %04d (%04d - %04d), iteration %03d, || v_i - v_ (i-1) || / || v_i * epsilon || = %0.16f / %0.16f\r', ...
	     d, d-1, ddiff+d-1, iter, dnorm, cnorm)   
    
  end % while
  
  if ( (rankmap ~= 0) | (normmap~= 0) )  
    
    rvcurr (:, d) = tiedrank (vcurr) ... % tied ranks: equal values lead to equal ranks
	/ (length (vcurr) +1);           % division: from uniform [1, N] to uniform ]0, 1[
    
    if ( normmap ~= 0) 
      
      mu = 0;                           % mean for N (0, 1) distributed centralities
      sig = 1;                          % standard deviation for N (0, 1) distribution
      
      % produce a map of gaussianised EC ranks if requested and write to normECM  
      nvcurr (:, d) = mu+sqrt (2) *sig ... % probit function (en.wikipedia.org/wiki/Probit) 
	  * erfinv ...                     % based on inverse error function 
	  (2*rvcurr (:, d) -1);            % uniform ]0, 1[ to N (0, 1) via inverse transform sampling
      
    end %if (normmap) 
    
  end % if (rankmap) 
  
  vcurr_out (:, d) = vcurr;
  fprintf ('\n');
  
end % for d



% write the eigenvactor centrality map to the file fastECM.nii (.gz) 
% and write orther centralities (node power, uniform, normal) if requested
write_map (inputfile, M, msk, atl, vcurr_out, 'fastECM', sprintf ('ECM [iterations: %d]', iter) );
if ( exist ('dvcurr') == 1 ) 
  
  write_map (inputfile, M, msk, atl, dvcurr, 'degCM', 'weighted degree centrality (node power) ');
  
end % if dvcurr

if ( (exist ('rvcurr') == 1) & rankmap ) 
  
  write_map (inputfile, M, msk, atl, rvcurr, 'rankECM', 'ECM [converted to ]0, 1[ ranks]');
  
end % if rvcurr

if ( exist ('nvcurr') == 1 ) 
  
  write_map (inputfile, M, msk, atl, nvcurr, 'normECM', 'ECM [converted to N (0, 1) values]');
  
end % if nvcurr



% if graph measures have been computed -> write them to a spreadsheet-readable
% XML-file (take care that native XML may take a long time to load even for small files) 
if (exist ('connmat_out') == 1) 
  
  % put matrices also in [vector_per_volume #volumes] format
  connmat_out = reshape (connmat_out, [prod (size (clus) ) dynamics]); 
  vcurr_bin = reshape (vcurr_bin, [prod (size (clus) ) dynamics]); 
  vcurr_mst = reshape (vcurr_mst, [prod (size (clus) ) dynamics]); 
  
  write_map (inputfile, M, ones (size (clus) ), 1, connmat_out, 'connections', 'connectivity matrix');
  write_map (inputfile, M, ones (size (clus) ), 1, vcurr_bin, 'backbone',  'binary backbone');
  write_map (inputfile, M, ones (size (clus) ), 1, vcurr_mst, 'min_span',  'minimal spanning tree');
  write_map (inputfile, M, msk,    atl, communities, 'communities', 'community_louvain');
  write_map (inputfile, M, msk,    atl, betweenness, 'betweenness', 'betweenness');
  write_map (inputfile, M, msk,    atl, clustering, 'clustering', 'clustering');
  write_map (inputfile, M, msk,    atl, pathlength, 'path_length', 'pathlength');
  
  % also make a readable file
  xmlfile = [fd filesep fn '_fastECMstats.xml'];
  fid = fopen (xmlfile, 'w');
  if (fid> (-1) ) 
    
    fprintf (fid, '<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>\n');
    fprintf (fid, ['<ss:Workbook xmlns:ss = "urn:schemas-microsoft-com:office:spreadsheet">\n']);
    writetosheet (fid, vcurr_out, 'fastECM');
    
    if ( exist ('dvcurr') == 1 ) 
      
      writetosheet (fid, dvcurr, 'degCM');
      
    end
    
    if ( (exist ('rvcurr') == 1) & rankmap ) 
      
      writetosheet (fid, rvcurr, 'rankECM');
      
    end
    
    if ( exist ('nvcurr') == 1 ) 
      
      writetosheet (fid, nvcurr, 'normECM');
      
    end
    
    %writetosheet (fid, connmat_out, 'connections'); % crashes spreadsheet program
    %writetosheet (fid, vcurr_bin, 'backbone'); % crashes spreadsheet program
    writetosheet (fid, communities, 'communities'); 
    writetosheet (fid, betweenness, 'betweenness');
    writetosheet (fid, clustering, 'clustering');
    writetosheet (fid, pathlength, 'path_length');
    fprintf (fid, '</ss:Workbook>\n');
    fclose (fid);
    
  end % if fid
  
end % if connmat_out

% if nifti support was added before, remove it again to not leave prints
if (usenifti) 
  
  rmpath (npath);                       % was own nifti added to path -- if yes then remove now
  
end % if usenifti 

return % fastECM

% tiedrank - Tied rank of each element in a set
% [r] = tiedrank (X, dim) 
% Computes tied rank of each element of X along given dimension 
% (default is to use the first non-singleton dimension) 
% "Tied ranking" is such that ex-aequo elements share the same (possibly half) rank: 
%  >> tiedrank (['DABBC']) %-> 6 1 2.5 2.5 4
% from https://github.com/kndiaye/matlab/blob/master/tiedrank.m
function [r] = tiedrank (X, dim) 

if nargin<2
  
  dim = min (find (size (X) >1) );
  
  if isempty (dim) 
    
    error ('X is empty!') 
    
  end % if isempty
  
end % if nargin

[ignore, Y] = sort (X, dim);
[ignore, r1] = sort (Y, dim);
[ignore, Y] = sort (-X, dim);
[ignore, r2] = sort (Y, dim);
r2 = size (X, dim) -r2+1;
r = (r1+r2) /2;

return

% function to write a variable as an XML sheet
% the option of writing flattened upper triangular matrices is implemented but never used!
function writetosheet (fid, variable, sheetname);

fprintf (fid, '<ss:Worksheet ss:Name = "%s">\n<ss:Table>\n', sheetname);

if (length (size (variable) ) == 3) 
  
  for m = 1:size (variable, 3)          % write matrix by flattening upper triangular matrix
    
    variable2d (:, m) = nonzeros (triu (squeeze (variable (:, :, m) ), 1) ); 
    
  end % for m
  variable = variable2d;
  clear variable2d;
  
end % if length

for c = 1:size (variable, 2) 
  
  fprintf (fid, '<ss:Row>\n');
  
  for r = 1:size (variable, 1) 
    
    fprintf (fid, '<ss:Cell><ss:Data ss:Type = "Number">%0.8f</ss:Data></ss:Cell>', variable (r, c) );
    
  end % for r
  
  fprintf (fid, '</ss:Row>\n');
  
end % for c

fprintf (fid, '</ss:Table>\n</ss:Worksheet>\n');

return



% function to write a map/matrix as a nifti file
function write_map (inputfile, M, msk, atl, vcurr, mapfile, descrip) 
%
% writes a map based on a template file
% - inputfile = filename of the template
% - M   = nifti record of the template
% - msk  = brain mask of the template
% - atl  = atlas volume (if used, otherwise 1) 
% - vcurr  = map values inside brain mask
% - mapfile = filename of the map file
% - descrip = description of contents
%
%
%
% this routine uses the following external code
%
% " Tools for Nifti / Analyze " by Jimmy Shen
% http://www.mathworks.com/matlabcentral/fileexchange/8797
% " quantile.m " by Anders Holtsberg
% www.spatial-econometrics.com/distrib/quantile.m
%



% if atlas used, make a map of regional centralities
if (prod (size (atl) ) ~= 1) 
  
  vcurr_2 = zeros (size (atl) );
  reg = unique (atl);     % region labels: nonzero values found inside atl
  for d = 1:size (vcurr, 2) 
    
    for r = 1:length (reg) 
      
      nvox = find (atl == reg (r) );
      vcurr_2 (nvox, d) = vcurr (r, d);%/ (length (nvox) );
      
    end % for r
    
  end % for d
  
  vcurr = vcurr_2;
  clear vcurr_2;
  
end % if prod

% create output array mout (and initialise 0) 
mout = zeros ([size(msk) size(vcurr, 2)]);
mout = nan*mout;                        % NaN outside mask

qnts = quantiles (vcurr ...
		  (vcurr>0), ...
		  [.01 .99]);           % quantiles to set contrast

% put vcurr in mout
smout = size (mout);
mout = reshape (mout, [prod(size(msk)) size(vcurr, 2)]);
mout (find (msk), :) = vcurr;
mout = reshape (mout, smout);

% write output nifti file based on mout
[fd, fn, fx2] = fileparts (inputfile);
[fd1, fn, fx] = fileparts (fn);         % base file name and 1st extension (2nd may be .gz) 
fx = [fx fx2];                          % in which case, concatenate by adding fx2
outputfile = [fd filesep fn '_' mapfile fx];
fprintf ('writing %s ...\n', outputfile);
M.hdr.hist.descrip = sprintf ('generated by fastECM - %s', descrip);
M.hdr.dime.dim = [ length(size(mout)) size(mout) ];  
M.hdr.dime.dim ( (end+1) :8) = 1;       % map is not 4D 
M.hdr.dime.pixdim ( (end+1) :8) = 1;    % neither are its voxels
if ~max ([32, 64] == M.hdr.dime.bitpix) % if data type not at least float
  M.hdr.dime.datatype = 64;             % make it double
end;       
M.hdr.dime.cal_min = qnts (1);          % min of range (for win/lev) 
M.hdr.dime.cal_max = qnts (2);          % max
M.hdr.dime.glmin = qnts (1);            % min of range
M.hdr.dime.glmax = qnts (2);            % max
M.hdr.dime.scl_slope = 1;               % the quantiles only make sense
M.hdr.dime.scl_inter = 0;               % with slope 1 and intercept 0
M.img = mout;                           % add voxel data to map 
save_untouch_nii (M, outputfile);       % write the file

return
