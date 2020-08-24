%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%
matlabbatch{1}.spm.spatial.coreg.estwrite.ref    = { '/usr/local/spm12/canonical/avg305T1.nii,1'         };
matlabbatch{1}.spm.spatial.coreg.estwrite.source = { [pwd filesep 'images' filesep 'T1_noneck.nii,1'   ] };
matlabbatch{1}.spm.spatial.coreg.estwrite.other  = { [pwd filesep 'images' filesep 'c1T1_noneck.nii,1' ] };
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'ncc';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'r';