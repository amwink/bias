%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%

nrun     = 1; % enter the number of runs here
jobfile  = { [ pwd filesep 'spm12_segment_job.m' ] };
jobs     = repmat ( jobfile, 1, nrun );
inputs   = cell ( 0, nrun );
for crun = 1:nrun
end
spm        ( 'defaults', 'FMRI'          );
spm_jobman ( 'run',      jobs, inputs{:} );
