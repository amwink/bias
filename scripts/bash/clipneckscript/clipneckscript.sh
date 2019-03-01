#!/bin/bash

# Script written by Alle Meije Wink based on a script from Veronica Popescu & Hugo Vrenken 
#
# The should have as argument the input file.
#
# Dit script maakt eerst een warp van standard space brain naar patient brain. Die matrix 
# toepassen op Hugo's (Heronica's?) 152 cut brain. Daarna de minimum Z berekenen van dit 
# brein en de nieuwe Z space. Daarna met fslroi het brein snijden dat het laagste puntje 
# van het standard-brain-naar-patient-space de ondergrens van het no-neck brain wordt. 

if [[ "$1" == "" ]]; then
	echo "usage: ${0} <anatomical MR scan> "
	exit 2
fi

INPUTFILE="${1}"
SCRIPTDIR=$(dirname $(which ${0}))

# clip .nii(.gz) extension
T1orig=${INPUTFILE%%.nii*} 

# map MNI head to input scan
MFILE=$(dirname {T1orig%})/avg152_to_$(basename ${T1orig%.nii*}).mat
if [[ ! -f ${MFILE} ]]; then
    CM="nice flirt -in   $FSLDIR/data/standard/avg152T1    -ref  $T1orig    -omat ${MFILE}"
    echo "$CM";$CM
fi

# map a box of 0s around an "inner MNI volume" of 1s to input scan
# made via this recipe
# fslmaths $FSLDIR/data/standard/MNI152_T1_2mm -add 1 -thr 0 -roi 1 89 1 107 1 89 0 -1 -bin MNI_allowed -odt char 
#
NFILE=${T1orig}_removeneck
if [[ `ls -d1 ${NFILE}.nii* 2> /dev/null|wc -l` -lt 1 ]]; then
    CM="nice flirt -in  ${SCRIPTDIR}/MNI_allowed.nii    -ref $T1orig    -interp nearestneighbour    -applyxfm    -init ${MFILE}    -out ${NFILE}"
    echo "$CM";$CM
fi

# make a clipped file
CFILE=${T1orig}_noneck
if [[ `ls -d1 ${CFILE}.nii* 2> /dev/null|wc -l` -lt 1 ]]; then
    
    # find the bounding box of the "inner MNI volume" warped to the input's space and the lowest slice in that box
    list=`fslstats ${NFILE} -w`
    zmin=`echo ${list}|awk '{print $5}'`

    echo "zmin = " $zmin
    
    # set new sizes
    xsize=`fslval $T1orig dim1`
    ysize=`fslval $T1orig dim2`
    zsize=`fslval $T1orig dim3`
    let newzsize=${zsize}-${zmin}

    CM="nice fslroi $T1orig ${T1orig}_noneck 0 ${xsize} 0 ${ysize} ${zmin} ${newzsize}"
    echo "$CM";$CM
fi

#rm -f ${MFILE} # maybe better not, this is the most time-consuming bit!
if [[ -f ${T1orig}_noneck.nii || -f ${T1orig}_noneck.nii.gz ]]; then
    echo rm -f ${NFILE}*    
else
    echo "file ${T1orig}_noneck.nii(.gz) not found -- keeping intermediates"
    exit 1
fi

exit 0
