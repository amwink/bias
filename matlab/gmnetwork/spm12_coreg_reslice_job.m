%-----------------------------------------------------------------------
% Job saved on 12-Mar-2020 16:30:21 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
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