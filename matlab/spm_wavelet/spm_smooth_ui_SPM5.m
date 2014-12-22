function spm_smooth_ui_SPM5 
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
% Modified version of the original (spm99-compatible) program
% Conversion to SPM2 and SPM5, compatibility tests by Marko Wilke
%                                        University of Tuebingen
%                                                 22 - 04 - 2004
%

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

kernel=0;
if ((denoising==1) & (check_packages>2))
  dmethod = spm_input('Available Methods','+1','m', ...
		      'Based on N[0,1] noise (WaveLab)|Generalized Likelihood', ...
		      [1 2],1);
  
  if (dmethod==1)
    dscheme = spm_input('Denoising Scheme','+1','m', ...
			'HybridThresh|InvShrink|JamesStein|MADThresh|MinMaxTresh|SUREThresh(hard)|SUREThresh(soft)|VisuShrink(hard)|VisuShrink(soft)', ...
			[1 2 3 4 5 6 7 8 9],5);
    M={'HYBRID' 'INV' 'JAMES' 'MAD' 'MINMAX' 'SURH' 'SURS' 'VISH' 'VISS'}; 
    M=M{dscheme};
    level   = spm_input('Decomposition level',       '+1', 'r', 4);
    degree  = spm_input('Degree of spline wavelets', '+1', 'r', 3);
  else
    dscheme = spm_input('Denoising Scheme','+1','m', ...
			'gen.lik.(BOLD)|gen.lik.(MRI)', ...
			[1 2],2);
    M={'genlik_BOLD' 'genlik_MRI'};
    M=M{dscheme};
    level=0;degree=0;    
  end
else
  level=0;degree=0;
  M='SMOOTH';
  kernel = spm_input('Smoothing {FWHM in mm}','+1','e','0 0 0');
end;

% select files
P = spm_get5(Inf,'.img','select scans');
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

