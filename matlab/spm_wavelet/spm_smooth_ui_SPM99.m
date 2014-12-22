function spm_smooth_ui_SPM99 
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

% Smoothing or convolving
%___________________________________________________________________________
%
% Convolves image files with an isotropic (in real space) Gaussian kernel 
% of a specified width.
%
% Uses:
%
% As a preprocessing step to suppress noise and effects due to residual 
% differences in functional and gyral anatomy during inter-subject 
% averaging.
%
% Inputs
%
% *.img conforming to SPM data format (see 'Data')
%
% Outputs
%
% The smoothed images are written to the same subdirectories as the 
% original *.img and are prefixed with a 's' (i.e. s*.img)
%
%__________________________________________________________________________
% @(#)spm_smooth_ui.m	2.7	99/09/20

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'smooth',index_of_Analysis};
%_______________________________________________________________________

global BCH;

% get filenames and kernel width
%----------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,'2.4');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Smooth');
spm_help('!ContextHelp','spm_smooth_ui.m');

% wavelet denoising or Gaussian smoothing
denoising = spm_input('Denoising method','1','m', ...
                      'Wavelet-based denoising|Gaussian smoothing', ...
		      [1 2],1);

if ((denoising==1) & (check_packages>2))
  kernel=0;
  dscheme = spm_input('Denoising Scheme','+1','m', ...
                      'HybridThresh|InvShrink|JamesStein|MADThresh|MinMaxTresh|SUREThresh(hard)|SUREThresh(soft)|VisuShrink(hard)|VisuShrink(soft)', ...
		      [1 2 3 4 5 6 7 8 9],5);
  M={'HYBRID' 'INV' 'JAMES' 'MAD' 'MINMAX' 'SURH' 'SURS' 'VISH' 'VISS'}; 
  M=M{dscheme};
  level   = spm_input('Decomposition level',       '+1', 'r', 4);
  degree  = spm_input('Degree of spline wavelets', '+1', 'r', 3);
else
  level=0;degree=0;
  M='SMOOTH';
  kernel = spm_input('Smoothing {FWHM in mm}','+1','e','0 0 0');
end;

% select files
P = spm_get(Inf,'.img','select scans');
n = size(P,1);

% select output prefix
op = spm_input('output filename prefix','+1','s','s');

% select output dir
[fpp fpf fpe]=fileparts(P(1,:));
[a b]=fileparts(fpp);
od = spm_get(-1,'*','output directory',a);

spatial_fMRI(M,P,level,degree,kernel,op,od);
		      
%
% Wavelet-based denoising of MR images
%
% (c) Alle Meije Wink 2004
%












