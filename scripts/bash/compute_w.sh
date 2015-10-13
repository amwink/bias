#!/bin/bash
CD=${PWD}

# 
# syntax: compute_w.sh <rootdir> <reference> <study>
# 
# where:
#  <rootdir>   = main directory where all files can be found
#  <reference> = xls  file with reference subjects
#  <study>     = xls  file with study subjects
# 
#  format of xls files:
#   row 1: "filename"   "variate 1" "variate 2" etc
#   row 2: /path/to/nii     <value>     <value> ...
#   row 3: /path/to/nii     <value>     <value> ...
#   ...
# 
#  reference subjects: used to build GLM coefficients and SD
#  study subjects:     used to compute W-score maps based on reference
# 
#  (c) Alle Meije Wink 2015 (a.wink@vumc.nl)
# 

if [[ ${3} == "" ]]; then
   cat ${0}|grep '# '|head -19|sed -e 's|# |    |g'
   exit 1;
fi

RDR=${1} # root directory 
REF=${2} # XLS file describing reference group
STU=${3} # XLS file describing study group
	
# find xls2csv
printf "looking if xls2csv can be found                ... "
if [[ `which xls2csv 2> /dev/null` == "" ]]; then
    if [[ ! -f ${PWD}/bin/xls2csv ]]; then
	# build catdoc / xls2csv
	if [[ `which git` != "" ]]; then
	    git config --global http.sslverify false
	    git clone http://www.wagner.pp.ru/git/oss/catdoc.git;
	else
	    if [[ `which unzip` != "" ]]; then 
		wget -N http://www.cbv.ns.ca/owl/DOCS/tools/catdoc-0.94.2-win32.zip
		unzip -o catdoc-0.94.2-win32.zip 
	    else
		wget -N http://www.cbv.ns.ca/owl/DOCS/tools/catdoc-0.94.2.tar.gz
		tar zxvf catdoc-0.94.2.tar.gz
	    fi   
	    mv catdoc-0.94.2 catdoc
	fi
	cd catdoc; 
	./configure --prefix=${PWD%/*}; 
	make; 
	make install; 
	cd ..;

	# make config file
	echo source_charset=8859-1 > ~/.catdocrc
	echo target_charset=8859-1 >> ~/.catdocrc
	echo charset_path=${CD}/catdoc/charsets >> ~/.catdocrc
	echo map_path=${CD}/catdoc/charsets:${CD}/share/catdoc >> ~/.catdocrc
	echo format=ascii >> ~/.catdocrc
	echo unknown_char='?' >> ~/.catdocrc
	echo use_locale=yes >> ~/.catdocrc	
    fi
fi    
export PATH=${CD}/bin:${PATH//'${CD}/bin:'/}
printf "done\n"

################################################################################
######################### process reference dataset ############################
################################################################################

# Build reference GLM
printf "building reference GLM                         ... "
REFCSV=${REF%.xls*}.csv
CM="xls2csv ${REF}"
$CM > ${REFCSV} 2> /dev/null
nlines=`wc ${REFCSV} -l|awk '{print $1}'`;
REFGLM=${REF%.xls*}.glm
RFFGLM=${REF%.xls*}.glmfiles
rm -f ${CD}/tmp.txt
cat ${REFCSV}| \
    tail -${nlines}| \
    sed -e 's|"||g'| \
    sed -e 's|,| |g'| \
    sed -e 's|  | |g' \
    > ${CD}/tmp.txt
cat ${CD}/tmp.txt| \
    awk '{print $1}'| \
    grep 'nii' \
    > ${RFFGLM}
cat ${CD}/tmp.txt| \
    awk '{$1="1";print $0}'| \
    sed -e 's|^[ \t]*||'| \
    grep ' ' \
    > ${REFGLM}
rm -f ${CD}/tmp.txt
printf "done\n"

# read columns of reference GLM
REFHEADINGS=`head -1 ${REFCSV}|tr "," " "`

# For grouping directories: if first glm name not found, create grouping dirs (level: above sbject)
printf "checking nifti files for reference GLM         ... "
firstmap=`head -1 ${RFFGLM}|awk '{print $1}'`
if [[ ! -f ${firstmap} ]]; then
    echo $firstmap not found
    cat ${RFFGLM}|while read map; do
	if [[ ${map//nii/} != ${map} ]]; then
	    group1=`echo ${map//${RDR}/}|awk '{print $1}'`;
	    group=${group1%%/*};
	    if [[ -f ${group1//$group/} ]];then
		printf "MOVE ${group1//$group/} to new dir ${group}? (y/n)";
		read yn
		if [[ "$yn" == "y" ]]; then
		    mkdir -p ${RDR}/${group}
		    subjdir1=${group1//${group}\//}
		    subjdir=${subjdir1%%/*}		
		    mv ${RDR}/${subjdir} ${RDR}/${group}
		fi
	    fi
	fi
    done
fi
if [[ ! -f ${firstmap} ]]; then
    echo ${firstmap} still not found
    exit 1
fi
printf "done\n"

# make brain mask at right res
HAVEMASK=0
if [[ -f ${CD}/standardhull.nii || -f ${CD}/standardhull.nii.gz ]]; then    
    pixdim=`fslinfo ${firstmap} |grep pixdim1|awk '{print $2}'`
    pixdim2=`fslinfo ${CD}/standardhull |grep pixdim1|awk '{print $2}'`    
    if [[ `echo "$pixdim - $pixdim2"|bc` == 0 ]]; then
	HAVEMASK=1;
    fi      
fi

if [[ $HAVEMASK == 0 ]]; then
    pixdim=`fslinfo ${firstmap} |grep pixdim1|awk '{print $2}'`
    printf "generating brain mask                          ... "
    if [[ ! -f ${CD}/bias/matlab/gong/mkmask.sh ]]; then    
	git config --global http.sslverify false
	git clone https://github.com/amwink/bias.git
    fi
    ${CD}/bias/matlab/gong/mkmask.sh ${pixdim} > /dev/null
    flirt -in `locate grey.nii|grep spm8|head -1` -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -out ${CD}/standardgrey -applyisoxfm $pixdim
else
    printf "brain mask already computed                    ... "
fi
printf "done\n"

# compute coefficients
if [[ -f ${REFGLM%.glm*}_betas.nii.gz || -f ${REFGLM%.glm*}_betas.nii ]] && [[ -f ${REFGLM%.glm*}_sigsq.nii.gz || -f ${REFGLM%.glm*}_sigsq.nii ]]; then
    printf "GLM coefficients for this model already exist  ... "
else
    printf "computing GLM coefficients for reference group ... "
    R4D=${RFFGLM%.glm*}
    fslmerge -t ${R4D} `cat ${RFFGLM}`
    
    # make common mask (threshold smoothed maps at 10% robust range of nonzero vox)
    fslmaths ${R4D} -thrP 10 -bin -Tmean -thr 1 -bin ${REFGLM%.glm*}_commask
    #fslmaths ${CD}/standardgrey -thr 0.1 -bin -mas ${REFGLM%.glm*}_commask -mas ${CD}/standardhull ${REFGLM%.glm*}_commask  
    fslmaths ${REFGLM%.glm*}_commask -mas ${CD}/standardhull ${REFGLM%.glm*}_commask  
    
    # compute the glm for the reference group
    fsl_glm -i $R4D -d ${REFGLM} -m ${REFGLM%.glm*}_commask -o ${REFGLM%.glm*}_betas --out_sigsq=${REFGLM%.glm*}_sigsq
    
    # remove 4D reference file
    rm -f ${R4D}.nii*
    
    # restrict the mask for W-scores to where, inside the mask, betas are computed
    fslmaths ${REFGLM%.glm*}_betas -abs -bin -mas ${REFGLM%.glm*}_commask ${REFGLM%.glm*}_commask
fi
printf "done\n"

################################################################################
########################### process study dataset ##############################
################################################################################

printf "building study GLM                             ... "
# Build study GLM
STUCSV=${STU%.xls*}.csv
CM="xls2csv ${STU}"
$CM > ${STUCSV} 2> /dev/null
nlines=`wc ${STUCSV} -l|awk '{print $1}'`;
STUGLM=${STU%.xls*}.glm
STFGLM=${STU%.xls*}.glmfiles
rm -f ${CD}/tmp.txt
cat ${STUCSV}| \
    tail -${nlines}| \
    sed -e 's|"||g'| \
    sed -e 's|,| |g'| \
    sed -e 's|  | |g' \
    > ${CD}/tmp.txt
cat ${CD}/tmp.txt| \
    awk '{print $1}'| \
    grep 'nii' \
    > ${STFGLM}
cat ${CD}/tmp.txt| \
    awk '{$1="1";print $0}'| \
    sed -e 's|^[ \t]*||'| \
    grep ' ' \
    > ${STUGLM}
rm -f ${CD}/tmp.txt
printf "done\n"

# read columns of study GLM
STUHEADINGS=`head -1 ${STUCSV}|tr "," " "`

if [[ $STUHEADINGS == $REFHEADINGS ]]; then
    printf "study GLM matches reference GLM -> safe to proceed\r"
else
    printf "study GLM has headings:      ${STUHEADINGS}\nbut reference GLM has these: ${REFHEADINGS}.\nQuitting...\n"
    exit 1
fi

# For grouping directories: if first glm name not found, create grouping dirs (level: above sbject)
printf "checking nifti files for reference GLM         ... "
firstmap=`head -1 ${STFGLM}|awk '{print $1}'`
if [[ ! -f ${firstmap} ]]; then
    echo $firstmap not found
    cat ${STFGLM}|while read map; do
	if [[ ${map//nii/} != ${map} ]]; then
	    group1=`echo ${map//${RDR}/}|awk '{print $1}'`;
	    group=${group1%%/*};
	    if [[ -f ${group1//$group/} ]];then
		printf "MOVE ${group1//$group/} to new dir ${group}? (y/n)";
		read yn
		if [[ "$yn" == "y" ]]; then
		    mkdir -p ${RDR}/${group}
		    subjdir1=${group1//${group}\//}
		    subjdir=${subjdir1%%/*}
		    mv ${RDR}/${subjdir} ${RDR}/${group}
		fi
	    fi
	fi
    done
fi
if [[ ! -f ${firstmap} ]]; then
    echo ${firstmap} still not found
    exit 1
fi

# compute std dev as square root of _sigsq (see Rik's script #L76)
fslmaths ${REFGLM%.glm*}_sigsq -sqrt ${REFGLM%.glm*}_sig
printf "done\n"

#
# now the computation of W can be done via
#
#     real intensity - estimated intensity
# W = ------------------------------------
#               SD of reference
#
# where: real intensity = 
#             value in subject image
#        estimated intensity =
#             reference beta 0 (reference mean) +
#             subject regr 1 * beta 1           +
#             subject regr 2 * beta 2 ...
#        SD of reference =
#             sqt of residual mean squares
#             computed in the reference GLM
#

# number of betas (beta 0 is expected to be the mean)
num_beta=`fslinfo ${REFGLM%.glm*}_betas | grep dim4 | head -1 | awk '{print $2}'`
echo "Betas in GLM: $num_beta"

# number of subjects
num_subj=`cat ${STFGLM}|wc -l`
echo "Subjects in study: $num_subj"

rm -f ${CD}/beta_tmp.nii*
for (( subno=0; subno<num_subj; subno++ )); do

    printf "processing subj %04d:\n" ${subno}

    # get the lines from the file list and the GLM
    let subline=$subno+1;
    subnii=`head -${subline} ${STFGLM}|tail -1`
    subreg=`head -${subline} ${STUGLM}|tail -1`
    numreg=`echo ${subreg}|wc -w`

    # GLM used for reference should match GLM used for study
    if (( num_beta != numreg )); then
	echo "number of regressors (${numreg}) does not match number of betas (${num_beta})"
	exit 1
    fi

    # first beta should be the mean / intercept / estimated intensity w/o regressors
    regs=(${subreg#* });
    if [[ ${subreg%% *} != 1 ]];then
	echo "first beta is expected to be 1 for everyone (intercept)"
	exit 1
    fi
    printf "beta %04d\r" 0;   
    fslroi ${REFGLM%.glm*}_betas ${CD}/beta_tmp 0 1

    # add beta images weighted by the value for the current subject
    beta=1;
    for reg in "${regs[@]}"; do
	printf "beta %04d\r" $beta;
	CM="fslmaths ${REFGLM%.glm*}_betas -roi 0 -1 0 -1 0 -1 $beta 1 -Tmean -mul $num_beta -mul $reg -add ${CD}/beta_tmp ${CD}/beta_tmp"
	let beta=$beta+1;
	$CM
    done

    # beta_tmp now contains the expected intensity
    # w-score = (<subject intensity> - <expected intensity>) / (reference SD)
    subwnii=${subnii//.nii/_w.nii}
    CM="fslmaths ${subnii} -sub ${CD}/beta_tmp -div ${REFGLM%.glm*}_sig ${subwnii}"
    ${CM}
    printf "w-score map <rootdir>%s written\n" ${subwnii//${RDR}/}
    
    # remove estimated intensity for future consideration
    rm -f ${CD}/beta_tmp.nii*
    
done

