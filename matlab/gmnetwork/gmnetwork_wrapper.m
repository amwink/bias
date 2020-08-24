function gmnetwork_wrapper;
%
% Original function (in Amsterdam UMC):
% apply grey matter network scripts to the scans from
% /home/projects/AD_Niftis/work/amwink/upgradeScans 
% (on lnx-rad-0X - should be NII files in this dir).
%
% You can use this as an example to build your own script.
% (this is by no means a need-to-do-it-this-way template!)
%
% This script may be used to support the script 'gmnetworks.m' in 
% this directory; please read its header for credit and copyright.
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%

% from the Matlab prompt, run
% >> profile on; gmnet_wrapper; profile off

% some constants
niftiaddr = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/8797/versions/28/download/zip';
gmnetaddr = 'https://github.com/bettytijms/Single_Subject_Grey_Matter_Networks/archive/master.zip';

currdir   = pwd;
gmn_dir   = fileparts( which ( 'gmnetwork.m' ) );
spmdir    = '/usr/local/spm12';
spmdir    = '/usr/local/apps/spm/spm12_7487';
niftidir  = [ gmn_dir filesep 'tools4nifti' ];
newscrdir = [ gmn_dir filesep 'ScriptsFromEllen' ];

gmnstring  = 'Single_Subject_Grey_Matter_Networks';
gmnstring2 = [ gmnstring filesep 'Extract_individual_GM_networks' ];
gmnstring3 = [ 'batch_extract_networks' ];



% make sure we have SPM12
if ( exist ( 'spm' ) ~= 2 )
    fprintf ( '  SPM12 not found, adding to path ...\n' );
    addpath ( genpath ( spmdir ) );
    if ( exist ( 'spm' ) ~= 2 )
        error ( '  SPM still not found, try (re-)installing ...' );
    end
end
fprintf ( '  found spm: %s\n', which ( 'spm' ) );



% make sure we have Jimmy Shen's NIfTI tools'
if ( exist ( 'load_nii' ) ~= 2 )
  if ( exist ( niftidir ) ~= 7 )
    fprintf ( '  Jimmy Shens NIfTI toolbox not found, downloading ...\n' );
    s = strfind ( evalc ( [ '! /usr/bin/wget -O NIfTI-toolbox.zip --no-check-certificate ' niftiaddr ] ), 'saved' );
    if ( isempty ( s ) )
      error ( 'NIfTI_20XXXXXX.zip not found' );
    else
      [ a, b, c ] = mkdir ( niftidir );
      eval ( [ '! mv NIfTI*.zip ' niftidir ] );
      cd ( niftidir );
      evalc ( '! unzip NIfTI*.zip ' );
      evalc ( '! rm NIfTI*.zip*' );
      cd ( '..' );
    end; % if isempty
  end % if ~exist niftidir 
  fprintf ( '  adding %s to the path\n', [ niftidir ] );
  addpath ( genpath ( niftidir ) );
  if ( exist ( 'load_nii' ) ~= 2 )
    error ( 'toolbox still not installed, exiting ...' );
  end; % if exist loadnii (2)
end % if exist loadnii
fprintf ( '  found load_nii: %s\n', which ( 'load_nii' ) );



%%
%% this extra code -or something similar- can be used to generate
%% neck-removed T1s (cropped to the MNI bounding box) and c1 maps
%% (that is SPM speak for grey matter segmentations of a T1 scan)
%% ----> this code will probably not run outside lnx-rad-0X <----
%%
%%     % make NIfTI files of all T1 scans
%%     outputdirs={};
%%     for d=1:length (t1dirs)
%%         [ x y z ] = fileparts ( fileparts ( t1dirs { d } ) );
%%         outputdirs {d} = [ currdir filesep y ];
%%         mkdir ( outputdirs {d} );
%%         fprintf ( 'converting %s ... ', t1dirs { d } );
%%         c = evalc ( [ '! /usr/local/dcm2nii -o ' outputdirs{d} ' ' t1dirs{d} filesep '003-0001.img'] );
%%         fprintf ( 'done, output in %s\n', outputdirs {d} );
%%     end
%%     
%%     % use Hugo Vrenken's clipneck script to do just that
%%     for d=1:length (outputdirs)
%%         fprintf ( 'removing neck of %s ... ', outputdirs { d } );
%%         t1co = ls ( [ outputdirs{d} filesep 'co*1.nii.gz' ] );
%%         c = evalc ( [ '! clipneckscript.sh ' t1co ] );
%%         fprintf ( '\n' );
%%     end
%%     
%%     % run the SPM batch for skull stripping etc
%%     niis={};
%%     for d=1:length (outputdirs)
%%         niis {d} = ls ( [ outputdirs{d} filesep 'co*_noneck.nii.gz' ] );
%%         fprintf ( 'uncompressing %s ... ', niis { d } );
%%         c = evalc ( [ '! gzip -d ' niis{d} ] );
%%         niis {d} = strrep ( niis {d}, 'nii.gz', 'nii' );
%%     end
%%     
%%     spm12_segment; % SPM12 batch script
%%



% make sure we have Betty's code'
if ( exist ( [ gmnstring3 '_v20150902' ] ) ~= 2 )
    if ( exist ( [ gmn_dir filesep gmnstring ] ) ~= 7 )
        fprintf ( '  Bettys GM network scripts not found, downloading ...\n' );
        s = strfind ( evalc ( [ '! /usr/bin/wget --no-check-certificate ' gmnetaddr ] ), 'saved' );
        if ( isempty ( s ) )
            warning ( [ gmnstring '-master.zip not found' ] );
        else
            evalc ( [ '! unzip ' gmnstring '-master.zip' ] );
            evalc ( [ '! mv '    gmnstring '-master ' gmnstring ] );
        end
    else
        if ( ~exist ( [ gmnstring2 '_v20150902/' gmnstring3 '_v20150902.m' ] ) )
            evalc ( [ '! mv ' ...
                [ gmnstring2 '_v20150902/' gmnstring3 '_v20150902 ' ] ...
                [ gmnstring2 '_v20150902/' gmnstring3 '_v20150902.m'] ...
                ] )
        end
        addpath (genpath ( [ gmn_dir filesep gmnstring2 '_v20150902/' ] ) )
    end
end
fprintf ( '  found %s: %s\n', [ gmnstring3 '_v20150902' ], which ( [ gmnstring3 '_v20150902' ] ) );



% if we can find c1*_noneck.nii* images, we have GM masks and we can start
rootdir = currdir;                          % or change to suit your needs
c1files = evalc ( [ '!ls ' rootdir filesep '*' filesep 'c1*_noneck.nii*' ] );
c1files = textscan ( c1files, '%s' );
c1files = c1files {1};



% standard space template for resampling
MNIt1 = [ fileparts(which('spm')) filesep 'canonical/avg305T1.nii' ];



% make the networks with the function gmnetwork for each T1 / c1 pair
% 
% obs_net - observed grey matter correspondence
% ran_net - correspondence of permuted intensities
% 
for c1=1:length(c1files)  
  
  prefixes = { 'betty' 'allemeije' };
  for p=1:length(prefixes)
    
    [ obs_net, ran_net ] = gmnetwork ( c1files{c1}, MNIt1, prefixes{p} );  
  
  end; % for p
  
end; % for c1



return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


