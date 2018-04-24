/******************************************************************************
structures.h
Copyright (C) 2016  {fullname}

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Correspondence concerning RMSynth_GPU should be addressed to:
sarrvesh.ss@gmail.com

******************************************************************************/
#include "fitsio.h"
#include "hdf5.h"
#include "constants.h"

#ifdef __cplusplus
extern "C"
#endif

/* Structure to store the input options */
struct optionsList {
    char *qCubeName;
    char *uCubeName;
    char *freqFileName;
    char *outPrefix;

    int plotRMSF;
    double phiMin, dPhi;
    int nPhi;

    int nGPU;
    int fileFormat;
};

struct fits_header_parameters {
    // These variables are needed for the fits
    int maskAxisLen1, maskAxisLen2;
    float crval1, crval2, crval3;
    float crpix1, crpix2, crpix3;
    float cdelt1, cdelt2, cdelt3;
    char ctype1[CTYPE_LEN], ctype2[CTYPE_LEN], ctype3[CTYPE_LEN];
}

struct parameters {
    double phiMin, dPhi;
    int nPhi;

    int qAxisNum, uAxisNum;
    int qAxisLen1, qAxisLen2, qAxisLen3;
    int uAxisLen1, uAxisLen2, uAxisLen3;
    float lambda20;
    float K;
}

struct IOFileDescriptors {
    fitsfile *qFile, *uFile;
    fitsfile *qDirty, *uDirty, *pDirty;

    FILE *freq;

    hid_t qFileh5, uFileh5;
    hid_t qDirtyH5, uDirtyH5, pDirtyH5;

    hid_t qDataspace, uDataspace;
    hid_t qOutDataspace, uOutDataspace, pOutDataspace;
    hid_t qDataset, uDataset;
    hid_t qOutDataset, uOutDataset, pOutDataset;
    hid_t qMemspace, uMemspace;
    hid_t qOutMemspace, uOutMemspace, pOutMemspace;

};

/* Structure to store all information related to RM Synthesis */
struct parList {
    fitsfile *qFile, *uFile;
    fitsfile *qDirty, *uDirty, *pDirty;

    hid_t qFileh5, uFileh5;
    hid_t qDirtyH5, uDirtyH5, pDirtyH5;

    FILE *freq;

    int qAxisNum, uAxisNum;
    int qAxisLen1, qAxisLen2, qAxisLen3;
    int uAxisLen1, uAxisLen2, uAxisLen3;
    int maskAxisLen1, maskAxisLen2;
    float crval1, crval2, crval3;
    float crpix1, crpix2, crpix3;
    float cdelt1, cdelt2, cdelt3;
    char ctype1[CTYPE_LEN], ctype2[CTYPE_LEN], ctype3[CTYPE_LEN];

    float *freqList;
    float *lambda2;
    float lambda20;

    float *phiAxis;
    float *rmsf, *rmsfReal, *rmsfImag;
    float K;
};

struct DataArrays {
    float *freqList;
    int nFreq;
    float *lambda2;
    float *phiAxis;
    int nPhi;
    float *rmsf, *rmsfReal, *rmsfImag;
}

/* Structure to store useful GPU device information */
struct deviceInfoList {
    int deviceID;
    long unsigned int globalMem, constantMem, sharedMemPerBlock;
    int maxThreadPerMP;
    int maxThreadPerBlock;
    int threadBlockSize[N_DIMS];
    int warpSize;
    int nSM;
};

/* Structure to store timing information */
struct timeInfoList {
   clock_t startRead, stopRead;  /* Read time */
   float msRead;
   clock_t startWrite, stopWrite;/* Write time */
   float msWrite;
   clock_t startProc, stopProc;  /* Processing time */
   float msProc;
   clock_t startX, stopX;        /* Transfer time */
   float msX;
   clock_t startTime, endTime;   /* Total time */
   unsigned int cpuTime;
};
