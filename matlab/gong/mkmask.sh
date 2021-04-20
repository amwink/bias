#!/bin/bash
#
# make brain mask exculding the ventricles
# using the data from FSL
#
# syntax: mkmask.sh <dim>
# <dim> = isotropic voxel size in mm
#
# (c) Alle Meije Wink 2015
# a.m.wink@gmail.com

# show or not show command before executing
function doit { echo $1; $1; }
#function doit { $1; }

# required parameter: (isotropic) voxel size
# voxel size N -> voxels will be NxNxN mm^3
if [[ ${1} == "" ]]; then
    echo "syntax:"
    echo "${0} <num>"
    echo "results in a MNI brain mask with NxNxN mm^3 voxels"
    exit 1
fi

# check for proper installation of FSL
# i.e., make sure that $FSLDIR is set so 
# that ${FSLDIR}/data/standard can be found
if [[ ${FSLDIR} == "" || ! -e ${FSLDIR} ]]; then
    echo "${0}"
    echo "requires FSL to be installed and set up properly." 
    echo "variable FSLDIR (or its directory) doesn't exist."
    exit
fi
TEMPLATE_BRAIN=${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz

# threshold for excluding inside (lower for less ventricle excluded)
HITHRESH=3000;
# threshold for lenient mask (should always be lower than HITHRESH)
LOTHRESH=2000;

# Do everything in the current directory
CD=${PWD}

# Make a mask that excludes at least the ventricles.
CM="flirt -in ${TEMPLATE_BRAIN} -ref ${TEMPLATE_BRAIN} -out ${CD}/standard -applyisoxfm ${1}";doit "${CM}"
CM="${FSLDIR}/bin/fslmaths ${CD}/standard -thr ${HITHRESH} -bin ${CD}/ventrmask";doit "${CM}"

# Use this mask to make a hull mask: purge 
# the outside, but leave the inside intact.
# The resulting image is called standardhull.
CM="python3 ${CD}/ventrmask.py";doit "${CM}"

# Now the first (ventricle mask) is inverted 
# and multiplied by the standard hull mask.
# This guarantees that ventricles and CSF are 1
# and the brain + voxels outside the hull are 0
CM="fslmaths ${CD}/ventrmask -mul -1 -add 1 -mul ${CD}/standardhull ${CD}/ventrmask";doit "${CM}"

# Now a more lenient threshold is applied to 
# make a bigger brain mask.
# The ventricle mask is subtracted and only
# the surviving voxels (value >0) are included.
CM="fslmaths ${CD}/standard -thr ${LOTHRESH} -add ${CD}/standardhull -bin -sub ventrmask -bin standardmask";doit "${CM}"

exit 0;

