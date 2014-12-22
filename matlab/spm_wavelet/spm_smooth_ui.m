function spm_smooth_ui 
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

% determine the version of SPM and select the appropriate UI 
spmver=spm('ver');
fprintf('SPM version: %s\n',spmver);

switch(spmver)
  
 case 'SPM99' 
  spm_smooth_ui_SPM99  
 case 'SPM2'  
  spm_smooth_ui_SPM2  
 case 'SPM5'  
  spm_smooth_ui_SPM5  
 otherwise    
  spm_smooth_ui_SPM5  
  
end
