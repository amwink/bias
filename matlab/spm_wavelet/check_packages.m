function TotalPack=check_packages;
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

% test if WaveLab is present
if (~exist('SUREThresh','file'))
  WaveLabPresent=0;
else
  WaveLabPresent=1;
end

% test if Fractional Spline Wavelets are present
if (~exist('FFTfractsplinefilters','file'))
  FSplinesPresent=0;
  
else
  FSplinesPresent=1;
end

% test if the FDR plugin for SPM  is present
if (~strcmp(fileparts(which('spm_getSPM')),fileparts(which('spm_P_FDR'))))
  FDRpluginPresent=0;
else
  FDRpluginPresent=1;
end

% test if the FDR plugin for SPM  is present
if (~exist('genlik_MRI','file'))
  SanjaPresent=0;
  fprintf('hoi\n');
  if (exist([fileparts(which('spatial_fMRI')) '/sanja']))
    addpath([fileparts(which('spatial_fMRI')) '/sanja']);
    if (exist('genlik_MRI','file'))
      SanjaPresent=1;
    end
  end    
else
  SanjaPresent=1;
end

TotalPack=WaveLabPresent+FSplinesPresent+FDRpluginPresent+SanjaPresent;

if (~SanjaPresent)
  fprintf('Please note that spm_wavelet supports Wavelet-based denoising based on the Generalized Likelihood Ratio (GLR).\n');
  fprintf('Please make sure that the matlab file genlik_MRI.m\n');
  fprintf('is in your matlab path.\n')
end;

if(TotalPack<3)
  ErrorTitle='Packages Missing';
  ErrorMessg={'Warning','', ...
	      'Apart  from SPM, the wavelet denoising package', ...
	      'relies on three other  packages: WaveLab, the FDR', ...
	      'plugin for SPM and fractional spline wavelets. ', ...
	      'These can be found at:', ...
	      '', ...
	      'http://www-stat.stanford.edu/~wavelab', ...
	      'http://www.sph.umich.edu/~nichols/FDR/spm_fdr.tar.gz', ...
	      'http://bigwww.epfl.ch/blu/fractsplinewavelets', ...
	      '','Please make sure that you have these packages:'};
  if (~WaveLabPresent)
    ErrorMessg{end+1}='* WaveLab';
  end
  if (~FSplinesPresent)
    ErrorMessg{end+1}='* Fractional spline wavelets';
  end;
  if (~FDRpluginPresent)
    ErrorMessg{end+1}='* FDR plugin for SPM';
  end;
  ErrorMessg{end+1}='';
  ErrorMessg{end+1}='and that they precede SPM99 in $MATLABPATH';
  ErrorMessg{end+1}='';
  ErrorMessg{end+1}='Continuing with Gaussian smoothing';
  warndlg(ErrorMessg,ErrorTitle)
end

return






