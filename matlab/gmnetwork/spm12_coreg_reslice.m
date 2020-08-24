%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%
nrun = X; % enter the number of runs here
jobfile = { [ pwd filesep 'spm12_coreg_reslice_job.m' ] };
jobs    = repmat(jobfile, 1, nrun);
inputs  = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
