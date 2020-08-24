function [ ori_corr, ran_corr ] = gmnetwork ( gmmap, template, prefix );
%
% compute the 'network' of correspondence of local cortex shape
% respresented in cubes of grey matter mini environments
%
% ori_corr - returns cube similarities of the observed GM map
% ran_corr - returns the similarities after permuting voxel locations
% [optional]
% prefix   - alternative name to look for already registered gmmap
%
% Please give credit where it's due: this is an implementation of the method described in:
%
%      Similarity-based extraction of individual networks from gray matter MRI scans 
%      Betty Marije Tijms, Peggy Seriès, David John Willshaw, Stephen MacGregor Lawrie
%      Cerebral Cortex 22 (7): pages 1530-41, July 2012.
%      https://doi.org/10.1093/cercor/bhr221
%
% (C) Alle Meije Wink, 2020
%     a.wink@amsterdamumc.nl
%



% we need SPM12 and nifti toolbox
if ( exist ( 'spm' ) ~= 2 )
  error ( '  SPM not found, try (re-)installing ...' );
end
if ( exist ( 'load_nii' ) ~= 2 )
  error ( '  Jimmy Shen''s NIfTI toolbox for .nii.gz support not found, try (re-)installing ...' );
end;



% set prefix if it wasn't given
if ( ~exist ( 'prefix', 'var' ) );
  prefix=[];
end;



% side of each cube (currently only 3 works)
cubesize = 3;



% copy t1 acd gmmap, make names of MNI-space gmmap & t1 files
% if c1 has not been rigid-registered to MNI, do that as well

wmmap  = strrep ( gmmap, 'c1', 'c2' );
csfmap = strrep ( gmmap, 'c1', 'c3' );
t1     = strrep ( gmmap, 'c1', ''   );

tmp_t1         = strrep (        t1,      'co',     'tmp_co'   );
tmp_gmmap      = strrep (     gmmap,      'c1',     'tmp_c1'   );
resliced_t1    = strrep (    tmp_t1,      'tmp_co', 'rtmp_co'  );
resliced_gmmap = strrep ( tmp_gmmap,      'tmp_c1',  'rtmp_c1' );
new_gmmap      = strrep ( resliced_gmmap, 'rtmp_c1', 'mni_c1'  );

gmnetname = strrep ( new_gmmap, '.nii',   '_gmnet.nii' );
gmnetname = strrep ( gmnetname, '.nii',   '.nii.gz'    );
gmnetname = strrep ( gmnetname, '.gz.gz', '.gz'        );

fprintf ( '  output will be saved to %s\n', gmnetname  );

if ( ~exist ( new_gmmap, 'file' ) )
  
  if ( ~isempty ( prefix ) )
    
    % This script assumes that files end in '_prefix.nii*' if a
    % prefix is given. This makes it possible to compare multiple
    % versions of preprocessed files (SPM, FSL, FreeSurfer, etc).
    %
    % If there is a coregistered gmmap with _prefix_noneck -> use that
    %
    prefixgmmap = strrep ( new_gmmap, '_noneck', [ '_' prefix '_noneck' ] );
    fprintf  ( '  using %s\n    for %s\n', prefixgmmap, new_gmmap );
    copyfile ( prefixgmmap, new_gmmap );
    
  else
    
    % do coregistration on copies
    copyfile ( t1,    tmp_t1    );
    copyfile ( gmmap, tmp_gmmap );

    % reslice t1 + gmmap to the MNI space
    spm_jobman('initcfg');
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.source            = {      tmp_t1 };
    matlabbatch{1}.spm.spatial.coreg.estwrite.other             = {   tmp_gmmap };
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref               = {    template };
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'ncc'; % normalised cross_correlation    -- fast + good enough for t1 -> t1
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 0;     % nearest neighbour interpolation -- does not mix intensities
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'r';
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
    fprintf ( '   co-registering %s to template ... ', tmp_t1 );
    tmp_capture = evalc ( 'spm_jobman(''run'', matlabbatch)' );
    fprintf ( 'done\n' );
    
    % remove tmp files, rename the C1 registered to MNI   
    delete   ( tmp_t1                   );
    delete   ( tmp_gmmap                );
    delete   ( resliced_t1              );
    movefile ( esliced_gmmap, new_gmmap );
    
  end % if isempty
  
end; % if ~exist



% load the voxel data of the registered gmmap
gm_mni = load_untouch_nii ( new_gmmap );
gm_img = double ( gm_mni.img );



% find the cube grid placement that has the maximum # of non-empty cubes
% 
% cube_nonzeros (scalar)             is the number of cubes with non_zero variance
% cube_offsets [ 27 cube_nonzeros ]  are the 1D offsets in the volume used by those cubes
% gm_incubes [ 27 cube_nonzeros ]    are the grey matter densities at those offsets
[ cube_nonzeros, cube_offsets, gm_incubes ] = cube_grid_position ( template, gm_img, cubesize );
fprintf ( '  found %d cubes\n', cube_nonzeros );

% save cube positions as an atlas
cubesname = strrep ( new_gmmap, '.nii',   '_cubes.nii' );
cubesname = strrep ( cubesname, '.nii',   '.nii.gz'    );
cubesname = strrep ( cubesname, '.gz.gz', '.gz'        );
cubesmap                  = zeros ( size(gm_img) );
cubesmap ( cube_offsets ) =  ones ( size( cube_offsets, 1 ), 1 ) * ( 1:size( cube_offsets, 2 ) );
gm_mni.hdr.hist.descrip           = 'grey matter cube indices';
gm_mni.hdr.dime.dim               = [ 3 size(cubesmap) ];
gm_mni.hdr.dime.datatype          =  8;
gm_mni.hdr.dime.cal_min           =  0;
gm_mni.hdr.dime.cal_max           =  size ( cube_offsets, 2 );
gm_mni.hdr.dime.glmin             = -1;
gm_mni.hdr.dime.glmax             =  size ( cube_offsets, 2 );
gm_mni.hdr.dime.dim ( (end+1):8 ) =  1;
gm_mni.hdr.dime.pixdim ( 1:8 )    =  1;
gm_mni.hdr.dime.scl_slope         =  1;
gm_mni.hdr.dime.scl_inter         =  0; 
gm_mni.img = cubesmap;
save_untouch_nii ( gm_mni, cubesname );



% randomise grey matter intensities for threshold detection
% check with Betty: should coefficients be exchanged between cubes???
% 1. exchange freely
gm_random = reshape ( gm_incubes ( randperm ( prod ( size ( gm_incubes ) ) ) ), size ( gm_incubes ) );
%
% 2. exchange whole cubes only
%%gm_random = gm_incubes ( :, randperm ( size ( gm_incubes, 2 ) ) );
%
% 3. exchange cubes and exchange coefficients inside each cube
%%gm_random = gm_incubes ( randperm ( size ( gm_incubes, 1 ) ), randperm ( size ( gm_incubes, 2 ) ) ); 



% compute cubes' correlations (taking into account cube rotations)
% add_diag = 0;   % only multiples of 90 degrees, OR
add_diag = 1;     % add multiples of 45 degrees 
[ ori_corr ran_corr ] = cube_cross_correlation ( gm_incubes, gm_random, cubesize, add_diag ); 



% save both to a file
ori_corr ( :, :, 2 ) = ran_corr;
gm_mni.hdr.hist.descrip = 'grey matter shape correspondence';
gm_mni.hdr.dime.dim               = [ 3 size(ori_corr) ];
gm_mni.hdr.dime.datatype          = 16;
gm_mni.hdr.dime.cal_min           = -1;
gm_mni.hdr.dime.cal_max           =  1;
gm_mni.hdr.dime.glmin             = -1;
gm_mni.hdr.dime.glmax             =  1;
gm_mni.hdr.dime.dim ( (end+1):8 ) =  1;
gm_mni.hdr.dime.pixdim ( 1:8 )    =  1;
gm_mni.hdr.dime.scl_slope         =  1;
gm_mni.hdr.dime.scl_inter         =  0; 
gm_mni.hdr.hist.qform_code        =  0;
gm_mni.hdr.hist.sform_code        =  0;
gm_mni.img = ori_corr;
save_untouch_nii ( gm_mni, gmnetname );



% if we copied a prefixed file as the coregistered GM map, delete the copy
if ( ~isempty ( prefix ) )
  copyfile ( gmnetname, strrep ( gmnetname, '_noneck', [ '_' prefix '_noneck' ] ) );
  copyfile ( cubesname, strrep ( cubesname, '_noneck', [ '_' prefix '_noneck' ] ) );
  delete   ( new_gmmap );
end



return;
