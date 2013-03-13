#!/bin/bash
#
# compiles the binaries for rwt(Plus)
# (C) Alle Meije Wink (a.wink@vumc.nl)
#

# rwt mex files -largeArrayDims
mex mdwt.c mdwt_r.c
mex midwt.c midwt_r.c
mex mrdwt.c mrdwt_r.c
mex mirdwt.c mirdwt_r.c

# Rice Wavelet Toolbox 1D extensions
mex multi1Drdwt.c
mex multi1Dirdwt.c

# Frequency domain shift-invariant wavelet transforms
mex rFWD.c polyphase.c
mex riFWD.c polyphase.c
mex fsidwt.c polyphase.c
mex fisidwt.c polyphase.c