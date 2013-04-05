function gonglabels;
% builds the atlas used in the paper 
%   Gong, G. et al, Cerebral Cortex 19(3):524-36. 
%   doi:10.1093/cercor/bhn102
% by sampling 78 of the 116 AAL regions
%
% (C) 2013 Alle Meije Wink

% mapping of gong labels to AAL labels
% AAL label 27 --> Gong label 1
% AAL label 21 --> Gong label 2
%              ...
% AAL label 36 --> Gong label 77
% AAL label 30 --> Gong label 78

labels = [27 21  5 25  9 15  3  7 11 13 23 19 69 ... %  1 .. 13
	   1 17 57 59 61 63 65 67 49 51 53 43 45 ... % 14 .. 26
	  47 55 79 81 85 89 83 87 39 31 33 35 29 ... % 27 .. 39
	  28 22  6 26 10 16  4  8 12 14 24 20 70 ... % 40 .. 52
	   2 18 58 60 62 64 66 68 50 52 54 44 46 ... % 53 .. 65
	  48 56 80 82 86 90 84 88 40 32 34 36 30];   % 66 .. 78

% these are the names of the AAL atlases to be processed
% aal_MNI_V4.nii - the original atlas with a resolution of 2mm
% aal_MNI_V4.nii - the atlas, resampled at a resolution of 4mm

names = {'aal_MNI_V4.nii' 'aal_MNI_V4_4mm.nii'};

for lf=1:length(names)

  labelfile=names{lf};                       % select the AAL label file
  n=nifti(labelfile);                        % load it using SPM or niftimatlib
  
  nd=n.dat(:,:,:);                           % load 3D label data
  ng=zeros(size(nd));                        % 3D gong label data
  
  for ln=1:length(labels)
    
    lvox=find(nd==labels(ln));               % find AAL region label(ln) corresponding to gong label ln
    ng(lvox)=ln;                             % set value to ln instead of label(ln)
    
  end % for ln
  
  [a b]=fileparts(n.dat.fname);              % get directory of the AAL label file
  n.dat.fname=[pwd filesep b '_gong.nii'];   % make a 'neighbour' file in the same directory
  fprintf('writing %s\n',n.dat.fname);        
  n.dat(:,:,:)=ng;                           % assign the new label volume to the file
  create(n);                                 % create the new file

end % for lf

return