#
# Note: the most advanced implementation of the fastECM method is the MatLab one
#       https://github.com/amwink/bias/blob/master/matlab/fastECM/fastECM.m 
#
# It supports external masks, atlases for regional centrality, dynamic ECM and
# computation of graph metric from the connectivity matrix using the brain 
# connectivity toolbox BCT, see: http://www.brain-connectivity-toolbox.net
#
# There is a BCT for Python called BCTpy [ https://pypi.python.org/pypi/bctpy ]
# but of course getting this integrated into fastECM takes time - bear with me
#
# AMW
#

import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

def fastECM(inputfile = '',
            maskfile  = '',
            atlasfile = '',
            verbose   = True,
            dynamics  = 1,
            rankmap   = 0,
            normmap   = 0,
            degmap    = 0,
            maxiter   = 42,
            wholemat  = 0):
    """
    module fastECM:

    import fastECM
    fastECM.fastECM(inputfile, rankmap, normmap, degmap, maxiter, maskfile, atlasfile, wholemat, dynamics, verbose)
        - inputfile is the name of a 4D 
          fMRI image time series in the NifTI format
        - maskfile selects voxels inside a mask
        - atlasfile groups signals by pre-defined regions
        - verbose: print ouput if True
        - dynamics: number of ECM to extract from a 4D fMRI
        - rankmap ~=0 produces a file rankECM.nii
          with uniformly distributed [0,1] centrality (tied) ranks
        - normmap ~=0 produces a file normECM.nii
          with normally distributed N(0,1) centralities generated from ranks
        - degmap ~=0 produces a node power map (similar to degree)
          this is the column sum of connection strengths
        - maxiter limits the number of iterations of the algorithm
        - wholemat: do not use the 'fast' recipe, write out
          connectivity matrix and compute graph measures (costly!)

    (c) Alle Meije Wink -- 16/03/2012
        a.m.wink AT gmail.com

    If you use this method in your research, remember to cite this paper in the journal "Brain Connectivity":
    Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
                     "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI:
                     implementation, validation and interpretation."
                         Brain Connectivity 2012, Vol. 2 No. 5, pages 265-274
                         URL:    http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087  
                         or:     http://dare.ubvu.vu.nl/handle/1871/48750
    """   
    err=0   
   
    # load the input file   
    img            = nib.load(inputfile)        # open nifti file using nibabel
    m              = np.asarray(img.dataobj)    # quickest way to get to the data
    m[np.isnan(m)] = 0                          # set NaNs to 0
    
    # get the data to re-format to the right size
    msz = m.shape                               # get the (assumed 4) dimensions of the data
    tln = msz[3]                                # 4th dimension (3): time series length
    msz = msz[:3]                               # first 3 dimensions: size of volume
    msk = (m.min(3)>0)                          # mask: minimum over time > 0 -> 1, else 0
        
    # if a maskfile was supplied, use this as a mask as well
    if maskfile:
        img_mask=nib.load(maskfile)             # open nifti file
        mmsk=np.asarray(img_mask.dataobj)       # get voxel data
        msk=msk*(mmsk>0).astype(float)          # mask within available data

    # if an atlasfile was supplied, sample with (atlas & mask) and correlate labels
    if atlasfile:
        img_atl=nib.load(atlasfile)             # open nifti file
        matl=np.asarray(img_atl.dataobj)        # get voxel data
        msk=msk*(matl>0).astype(float)          # mask within available atlas data
        matl=matl*(msk>0).astype(float)         # only sample atlas within mask
    else:   
        matl=0

    # apply the final mask and count the number of voxels 
    msk = msk.reshape(np.prod(msz))             # reshape msk as a 1D vector
    m = m.reshape(np.prod(msz),tln).transpose() # reshape m as a sequence of vectors
    m = m[:,np.nonzero(msk)].squeeze()          # remove singleton dimensions
    msk = msk.reshape(msz)
    npt = m.shape[1]                            # count the number of voxels in the mask

    # when an atlas is used, signals are grouped within atlas regions
    if atlasfile:
        msk = msk.reshape(np.prod(msz))         # reshape msk as a 1D vector
        matl=matl.reshape(np.prod(msz))         # same 1D shape as mask above
        matl=matl[np.nonzero(msk)].squeeze()    # sample matl in mask (same size as m)
        regs=np.unique(matl[np.where(matl>0)])  # find the used region values
        rgts=np.zeros((tln,len(regs)))          # regional time series
        for t in range(len(regs)):            
            #print regs[t]
            rgts[:,t]=np.mean(m[:,np.where(matl==regs[t])], axis=2).squeeze()
        m=rgts                                  # now m contains regional mean time series
        np.savetxt(os.path.abspath(inputfile).replace(".nii.gz","_fastECM.txt"),m,fmt="%d")
        npt=len(regs)
        msk  = msk.reshape(msz)
    
    # Initialize eigenvector estimate    
    vprev = 0                                   # initialize previous ECM estimate
    vcurr = np.ones((npt,1))/np.sqrt(npt)       # initialize estimate with L2-norm == 1
    
    i = 0                                       # reset iteration counter
    dnorm = 1                                   # initial value for difference L2-norm
    cnorm = 0                                   # initial value to estimate L2-norm

    # Efficient power iteration
    while (i < maxiter) & (dnorm > cnorm):
        vprev   = vcurr                         # start with previous estimate
        prevsum = vprev.sum()                   # sum of estimate
        vcurr_1 = m.dot(vprev)                  # part one of M*v
        vcurr_2 = m.T.dot(vcurr_1)              # part two of M*v
        vcurr_3 = vcurr_2 + prevsum             # adding sum -- same effect as [M+1]*v
        vcurr   = vcurr_3/np.linalg.norm(vcurr_3,2)    # normalize L2-norm
          
        i += 1                                         # next iteration
        dnorm = np.linalg.norm(vcurr-vprev, 2)         # L2 norm of difference with previous result
        cnorm = np.linalg.norm(vcurr,2)*np.spacing(1)  # L2 norm of current result
        
        if verbose:
            print "iteration %02d, || v_i - v_(i-1) || / || v_i * epsilon || = %0.16f / %0.16f" % (i, dnorm, cnorm) #  some stats for the users
    
    if (i >= maxiter) & (dnorm > cnorm):        # test if we have converged (or reached max_iter)
        err = 128
        print "fast ECM error: power iteration algorithm did not converge"

    print vcurr
        
    # write the ECM
    write_map(inputfile,img,vcurr,msk,matl)
    
    # return 0 if no error   
    return err  
    
def write_map (inputfile='',
               img=0,
               vcurr=0,
               msk=0,
               atl=0, 
               mapfile='',
               descrip=''):
    """
    write_map(inputfile, img, msk, atl, mapfile, descrip)       
            - inputfile is the name of the 4D 
              fMRI image time series in the NifTI format
            - msk is the mask of in-brain voxels
            - atl groups signals by pre-defined regions
            - vcurr is the vector with coefficients
            - mapfile: file suffix for the mapfile
            - descrip: description of contents
            
    (c) Alle Meije Wink -- 16/03/2012
    a.m.wink AT gmail.com
            
    If you use this method in your research, remember to cite this paper in the journal "Brain Connectivity":
    Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
        "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI:
        implementation, validation and interpretation."
        Brain Connectivity 2012, Vol. 2 No. 5, pages 265-274
        URL:    http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087  
        or:     http://dare.ubvu.vu.nl/handle/1871/48750
    """
    
    if np.prod(atl):
        vcurr_2=np.zeros(atl.shape)             # allocate in-brain voxels
        regs=np.unique(atl)                     # find the used region values
        for t in range(0,len(regs)):
            all_t=np.where(atl==t)
            vcurr_2[all_t] = vcurr[t]#/len(all_t)
        vcurr=vcurr_2                           # now m contains regional mean time series
    
    outputfile=os.path.abspath(inputfile).replace(".nii","_fastECM.nii")
    msz=msk.shape                                # size of the volume
    msk=np.reshape(msk,np.prod(msz))             # 1D vector of the same size
    vout=msk.astype('float32')                   # output: zero outside mask                     
    vout[np.nonzero(msk)]=vcurr                  # put vcurr inside mask
    vout=np.reshape(vout,msz).astype('float32')  # resize to volume 
    print vout.dtype
    
    img_out=nib.Nifti1Image(vout,img.get_affine())
    img_out.to_filename(outputfile)   
