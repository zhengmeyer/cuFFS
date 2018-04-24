params/******************************************************************************
fileaccess.c
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
#include "structures.h"
#include "constants.h"
#include "fileaccess.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#define BUNIT   "JY/BEAM"
#define RM      "PHI"

/*************************************************************
*
* Check Fitsio error and exit if required.
*
*************************************************************/
void checkFitsError(int status) {
    if(status) {
        printf("ERROR:");
        fits_report_error(stdout, status);
        printf("\n");
        exit(FAILURE);
    }
}

/*************************************************************
*
* Check of the input files are open-able
*
*************************************************************/
void checkInputFiles(struct optionsList *inOptions, struct IOFileDescriptors *descriptors) {
   int fitsStatus = SUCCESS;
   herr_t error;
   char buf[STRING_BUF_LEN];

   if(inOptions->fileFormat == FITS) {
      /* Check if all the input fits files are accessible */
      fits_open_file(&(descriptors->qFile), inOptions->qCubeName, READONLY, &fitsStatus);
      fits_open_file(&(descriptors->uFile), inOptions->uCubeName, READONLY, &fitsStatus);
      checkFitsError(fitsStatus);
   }
   else if(inOptions->fileFormat == HDF5) {
      /* Open HDF5 files */
      descriptors->qFileh5 = H5Fopen(inOptions->qCubeName, H5F_ACC_RDONLY, H5P_DEFAULT);
      descriptors->uFileh5 = H5Fopen(inOptions->uCubeName, H5F_ACC_RDONLY, H5P_DEFAULT);
      if(descriptors->qFileh5 < 0 || descriptors->uFileh5 < 0) {
         printf("Error: Unable to open the input HDF5 files\n\n");
         exit(FAILURE);
      }
      /* Check if the hdf5 files are compatible with HDFITS format */
      error = H5LTget_attribute_string(descriptors->qFileh5, "/", "CLASS", buf);
      if(error < 0) {
         printf("ERROR: Specified HDF5 file is not in HDFITS format\n\n");
         exit(FAILURE);
      }
   }
   else {}

   /* Check if you can open the frequency file */
   descriptors->freq = fopen(inOptions->freqFileName, FILE_READONLY);
   if(descriptors->freq == NULL) {
      printf("Error: Unable to open the frequency file\n\n");
      exit(FAILURE);
   }
}

/*************************************************************
*
* Read header information from the fits files
*
*************************************************************/
int getFitsHeader(struct optionsList *inOptions, struct fits_header_parameters *params, struct IOFileDescriptors *descriptors) {
    int fitsStatus = SUCCESS;
    char fitsComment[FLEN_COMMENT];

    /* Remember that the input fits images are rotated. */
    /* Frequency is the first axis */
    /* RA is the second */
    /* Dec is the third */

    /* Get the image dimensions from the Q cube */
    fits_read_key(descriptors->qFile, TINT, "NAXIS", &params->qAxisNum,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TINT, "NAXIS1", &params->qAxisLen3,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TINT, "NAXIS2", &params->qAxisLen1,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TINT, "NAXIS3", &params->qAxisLen2,
      fitsComment, &fitsStatus);
    /* Get the image dimensions from the U cube */
    fits_read_key(descriptors->uFile, TINT, "NAXIS", &params->uAxisNum,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->uFile, TINT, "NAXIS1", &params->uAxisLen3,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->uFile, TINT, "NAXIS2", &params->uAxisLen1,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->uFile, TINT, "NAXIS3", &params->uAxisLen2,
      fitsComment, &fitsStatus);
    /* Get WCS information */
    fits_read_key(descriptors->qFile, TFLOAT, "CRVAL1", &params->crval3,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CRVAL2", &params->crval1,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CRVAL3", &params->crval2,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CRPIX1", &params->crpix3,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CRPIX2", &params->crpix1,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CRPIX3", &params->crpix2,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CDELT1", &params->cdelt3,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CDELT2", &params->cdelt1,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TFLOAT, "CDELT3", &params->cdelt2,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TSTRING, "CTYPE1", &params->ctype3,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TSTRING, "CTYPE2", &params->ctype1,
      fitsComment, &fitsStatus);
    fits_read_key(descriptors->qFile, TSTRING, "CTYPE3", &params->ctype2,
      fitsComment, &fitsStatus);

    return(fitsStatus);
}

/*************************************************************
*
* Read header information from the HDF5 files
*
*************************************************************/
int getHDF5Header(struct optionsList *inOptions,
     struct fits_header_parameters *header_parameters,
     struct parameters *params,
     struct IOFileDescriptors *descriptors) {
    hsize_t tempArr[N_DIMS];
    herr_t error;

    /* Remember that the input fits images are NOT rotated. */
    /* RA is the first axis */
    /* Dec is the second */
    /* Frequency is the third */

    /* Get the dimensionality of the input datasets */
    error = H5LTget_dataset_ndims(descriptors->qFileh5, PRIMARYDATA, &(params->qAxisNum));
    error = H5LTget_dataset_ndims(descriptors->uFileh5, PRIMARYDATA, &(params->uAxisNum));
    /* Get the sizes of each dimension */
    error = H5LTget_dataset_info(descriptors->qFileh5, PRIMARYDATA, tempArr, NULL, NULL);
    params->qAxisLen1 = tempArr[1];
    params->qAxisLen2 = tempArr[2];
    params->qAxisLen3 = tempArr[0];
    error = H5LTget_dataset_info(descriptors->uFileh5, PRIMARYDATA, tempArr, NULL, NULL);
    params->uAxisLen1 = tempArr[1];
    params->uAxisLen2 = tempArr[2];
    params->uAxisLen3 = tempArr[0];
    /* Get WCS information */
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CRVAL1", &header_parameters->crval1);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CRVAL2", &header_parameters->crval2);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CRVAL3", &header_parameters->crval3);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CDELT1", &header_parameters->cdelt1);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CDELT2", &header_parameters->cdelt2);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CDELT3", &header_parameters->cdelt3);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CRPIX1", &header_parameters->crpix1);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CRPIX2", &header_parameters->crpix2);
    H5LTget_attribute_float(descriptors->qFileh5, PRIMARY, "CRPIX3", &header_parameters->crpix3);
    H5LTget_attribute_string(descriptors->qFileh5, PRIMARY, "CTYPE1", header_parameters->ctype1);
    H5LTget_attribute_string(descriptors->qFileh5, PRIMARY, "CTYPE2", header_parameters->ctype2);
    H5LTget_attribute_string(descriptors->qFileh5, PRIMARY, "CTYPE3", header_parameters->ctype3);

    return error;
}

/*************************************************************
*
* Create output FITS images
*
*************************************************************/
void makeOutputFitsImages(struct optionsList *inOptions,
    struct IOFileDescriptors *descriptors,
    struct fits_header_parameters *header_parameters,
    struct parameters *params) {
   int stat = SUCCESS;
   char filenamefull[FILENAME_LEN];
   long naxis[FITS_OUT_NAXIS];
   char fComment[FILENAME_LEN];
   float tempVar;

   /* Create the output Q, U, and P images */
   sprintf(filenamefull, "%s%s.fits", inOptions->outPrefix, Q_DIRTY);
   fits_create_file(&descriptors->qDirty, filenamefull, &stat);
   sprintf(filenamefull, "%s%s.fits", inOptions->outPrefix, U_DIRTY);
   fits_create_file(&descriptors->uDirty, filenamefull, &stat);
   sprintf(filenamefull, "%s%s.fits", inOptions->outPrefix, P_DIRTY);
   fits_create_file(&descriptors->pDirty, filenamefull, &stat);
   checkFitsError(stat);

   /* Assign empty string to fComment */
   sprintf(fComment, " ");

   /* What are the output cube sizes */
   naxis[0] = params->nPhi;
   naxis[1] = params->qAxisLen1;
   naxis[2] = params->qAxisLen2;

   /* Create the header for each output image */
   fits_create_img(descriptors->qDirty, FLOAT_IMG, FITS_OUT_NAXIS, naxis, &stat);
   fits_create_img(descriptors->uDirty, FLOAT_IMG, FITS_OUT_NAXIS, naxis, &stat);
   fits_create_img(descriptors->pDirty, FLOAT_IMG, FITS_OUT_NAXIS, naxis, &stat);
   checkFitsError(stat);

   /* Set the relevant keyvalues */
   fits_write_key(descriptors->qDirty, TSTRING, "BUNIT", BUNIT, fComment, &stat);
   fits_write_key(descriptors->uDirty, TSTRING, "BUNIT", BUNIT, fComment, &stat);
   fits_write_key(descriptors->pDirty, TSTRING, "BUNIT", BUNIT, fComment, &stat);

   fits_write_key(descriptors->qDirty, TDOUBLE, "CRVAL1", &params->phiMin, fComment, &stat);
   fits_write_key(descriptors->uDirty, TDOUBLE, "CRVAL1", &params->phiMin, fComment, &stat);
   fits_write_key(descriptors->pDirty, TDOUBLE, "CRVAL1", &params->phiMin, fComment, &stat);

   fits_write_key(descriptors->qDirty, TDOUBLE, "CDELT1", &params->dPhi, fComment, &stat);
   fits_write_key(descriptors->uDirty, TDOUBLE, "CDELT1", &params->dPhi, fComment, &stat);
   fits_write_key(descriptors->pDirty, TDOUBLE, "CDELT1", &params->dPhi, fComment, &stat);

   tempVar = 1;
   fits_write_key(descriptors->qDirty, TFLOAT, "CRPIX1", &tempVar, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CRPIX1", &tempVar, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CRPIX1", &tempVar, fComment, &stat);

   fits_write_key(descriptors->qDirty, TSTRING, "CTYPE1", RM, fComment, &stat);
   fits_write_key(descriptors->uDirty, TSTRING, "CTYPE1", RM, fComment, &stat);
   fits_write_key(descriptors->pDirty, TSTRING, "CTYPE1", RM, fComment, &stat);

   fits_write_key(descriptors->qDirty, TFLOAT, "CRVAL2", &header_parameters->crval1, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CRVAL2", &header_parameters->crval1, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CRVAL2", &header_parameters->crval1, fComment, &stat);

   fits_write_key(descriptors->qDirty, TFLOAT, "CDELT2", &header_parameters->cdelt1, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CDELT2", &header_parameters->cdelt1, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CDELT2", &header_parameters->cdelt1, fComment, &stat);

   fits_write_key(descriptors->qDirty, TFLOAT, "CRPIX2", &header_parameters->crpix1, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CRPIX2", &header_parameters->crpix1, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CRPIX2", &header_parameters->crpix1, fComment, &stat);

   fits_write_key(descriptors->qDirty, TSTRING, "CTYPE2", header_parameters->ctype1, fComment, &stat);
   fits_write_key(descriptors->uDirty, TSTRING, "CTYPE2", header_parameters->ctype1, fComment, &stat);
   fits_write_key(descriptors->pDirty, TSTRING, "CTYPE2", header_parameters->ctype1, fComment, &stat);

   fits_write_key(descriptors->qDirty, TFLOAT, "CRVAL3", &header_parameters->crval2, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CRVAL3", &header_parameters->crval2, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CRVAL3", &header_parameters->crval2, fComment, &stat);

   fits_write_key(descriptors->qDirty, TFLOAT, "CDELT3", &header_parameters->cdelt2, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CDELT3", &header_parameters->cdelt2, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CDELT3", &header_parameters->cdelt2, fComment, &stat);

   fits_write_key(descriptors->qDirty, TFLOAT, "CRPIX3", &header_parameters->crpix2, fComment, &stat);
   fits_write_key(descriptors->uDirty, TFLOAT, "CRPIX3", &header_parameters->crpix2, fComment, &stat);
   fits_write_key(descriptors->pDirty, TFLOAT, "CRPIX3", &header_parameters->crpix2, fComment, &stat);

   fits_write_key(descriptors->qDirty, TSTRING, "CTYPE3", header_parameters->ctype2, fComment, &stat);
   fits_write_key(descriptors->uDirty, TSTRING, "CTYPE3", header_parameters->ctype2, fComment, &stat);
   fits_write_key(descriptors->pDirty, TSTRING, "CTYPE3", header_parameters->ctype2, fComment, &stat);

   checkFitsError(stat);
}

/*************************************************************
*
* Create output HDF5 cubes
*
*************************************************************/
void makeOutputHDF5Images(struct optionsList *inOptions,
    struct IOFileDescriptors *descriptors,
    struct parameters *params,
    struct fits_header_parameters *header) {
   char filenamefull[FILENAME_LEN];
   float tempVar = 1;
   hsize_t dims[N_DIMS];
   herr_t qErr, uErr, pErr;
   hid_t qGrp, uGrp, pGrp;
   int positionID = POSITION_ID;

   /* Create the output Q, U, and P images */
   sprintf(filenamefull, "%s%s.h5", inOptions->outPrefix, Q_DIRTY);
   descriptors->qDirtyH5 = H5Fcreate(filenamefull, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
   sprintf(filenamefull, "%s%s.h5", inOptions->outPrefix, U_DIRTY);
   descriptors->uDirtyH5 = H5Fcreate(filenamefull, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
   sprintf(filenamefull, "%s%s.h5", inOptions->outPrefix, P_DIRTY);
   descriptors->pDirtyH5 = H5Fcreate(filenamefull, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

   /* Create the primary group */
   qGrp = H5Gcreate(descriptors->qDirtyH5, PRIMARY, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   uGrp = H5Gcreate(descriptors->uDirtyH5, PRIMARY, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   pGrp = H5Gcreate(descriptors->pDirtyH5, PRIMARY, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   if( qGrp<0 || uGrp<0 || pGrp<0) {
      printf("Error: Unable to create groups in output HDF5 files\n");
      exit(FAILURE);
   }
   H5Gclose(qGrp); H5Gclose(uGrp); H5Gclose(pGrp);

   /* Create the /PRIMNARY/DATA dataset */
   dims[0] = params->nPhi;
   dims[1] = params->qAxisLen1;
   dims[2] = params->qAxisLen2;
   qErr = H5LTmake_dataset_float(descriptors->qDirtyH5, PRIMARYDATA, N_DIMS, dims, NULL);
   uErr = H5LTmake_dataset_float(descriptors->uDirtyH5, PRIMARYDATA, N_DIMS, dims, NULL);
   pErr = H5LTmake_dataset_float(descriptors->pDirtyH5, PRIMARYDATA, N_DIMS, dims, NULL);
   if( qErr<0 || uErr<0 || pErr<0) {
      printf("Error: Unable to create output datasets in HDF5\n");
      exit(FAILURE);
   }

   /* CLASS attribute of ROOT should be set to HDFITS */
   H5LTset_attribute_string(descriptors->qDirtyH5, ROOT, "CLASS", HDFITS);
   H5LTset_attribute_string(descriptors->uDirtyH5, ROOT, "CLASS", HDFITS);
   H5LTset_attribute_string(descriptors->pDirtyH5, ROOT, "CLASS", HDFITS);

   /* Position attribute of /PRIMARY must be set to 1 */
   H5LTset_attribute_int(descriptors->qDirtyH5, PRIMARY, "POSITION", &positionID, sizeof(positionID));
   H5LTset_attribute_int(descriptors->uDirtyH5, PRIMARY, "POSITION", &positionID, sizeof(positionID));
   H5LTset_attribute_int(descriptors->pDirtyH5, PRIMARY, "POSITION", &positionID, sizeof(positionID));

   /* Create attributes for the /PRIMARY group */
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CRVAL1", &(header->crval1), sizeof(header->crval1));
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CRVAL2", &(header->crval2), sizeof(header->crval2));
   H5LTset_attribute_double(descriptors->qDirtyH5, PRIMARY, "CRVAL3", &(params->phiMin), sizeof(params->phiMin));
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CRPIX1", &(header->crpix1), sizeof(header->crpix1));
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CRPIX2", &(header->crpix2), sizeof(header->crpix2));
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CRPIX3", &tempVar, sizeof(tempVar));
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CDELT1", &(header->cdelt1), sizeof(header->cdelt1));
   H5LTset_attribute_float(descriptors->qDirtyH5, PRIMARY, "CDELT2", &(header->cdelt2), sizeof(header->cdelt2));
   H5LTset_attribute_double(descriptors->qDirtyH5, PRIMARY, "CDELT3", &(params->dPhi), sizeof(params->dPhi));
   H5LTset_attribute_string(descriptors->qDirtyH5, PRIMARY, "CTYPE1", header->ctype1);
   H5LTset_attribute_string(descriptors->qDirtyH5, PRIMARY, "CTYPE2", header->ctype2);
   H5LTset_attribute_string(descriptors->qDirtyH5, PRIMARY, "CTYPE3", RM);

   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CRVAL1", &(header->crval1), sizeof(header->crval1));
   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CRVAL2", &(header->crval2), sizeof(header->crval2));
   H5LTset_attribute_double(descriptors->uDirtyH5, PRIMARY, "CRVAL3", &(params->phiMin), sizeof(params->phiMin));
   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CRPIX1", &(header->crpix1), sizeof(header->crpix1));
   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CRPIX2", &(header->crpix2), sizeof(header->crpix2));
   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CRPIX3", &tempVar, sizeof(tempVar));
   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CDELT1", &(header->cdelt1), sizeof(header->cdelt1));
   H5LTset_attribute_float(descriptors->uDirtyH5, PRIMARY, "CDELT2", &(header->cdelt2), sizeof(header->cdelt2));
   H5LTset_attribute_double(descriptors->uDirtyH5, PRIMARY, "CDELT3", &(params->dPhi), sizeof(params->dPhi));
   H5LTset_attribute_string(descriptors->uDirtyH5, PRIMARY, "CTYPE1", header->ctype1);
   H5LTset_attribute_string(descriptors->uDirtyH5, PRIMARY, "CTYPE2", header->ctype2);
   H5LTset_attribute_string(descriptors->uDirtyH5, PRIMARY, "CTYPE3", RM);

   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CRVAL1", &(header->crval1), sizeof(header->crval1));
   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CRVAL2", &(header->crval2), sizeof(header->crval2));
   H5LTset_attribute_double(descriptors->pDirtyH5, PRIMARY, "CRVAL3", &(params->phiMin), sizeof(params->phiMin));
   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CRPIX1", &(header->crpix1), sizeof(header->crpix1));
   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CRPIX2", &(header->crpix2), sizeof(header->crpix2));
   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CRPIX3", &tempVar, sizeof(tempVar));
   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CDELT1", &(header->cdelt1), sizeof(header->cdelt1));
   H5LTset_attribute_float(descriptors->pDirtyH5, PRIMARY, "CDELT2", &(header->cdelt2), sizeof(header->cdelt2));
   H5LTset_attribute_double(descriptors->pDirtyH5, PRIMARY, "CDELT3", &(params->dPhi), sizeof(params->dPhi));
   H5LTset_attribute_string(descriptors->pDirtyH5, PRIMARY, "CTYPE1", header->ctype1);
   H5LTset_attribute_string(descriptors->pDirtyH5, PRIMARY, "CTYPE2", header->ctype2);
   H5LTset_attribute_string(descriptors->pDirtyH5, PRIMARY, "CTYPE3", RM);

   /* Create the /PRIMARY/DATA dataset */
   H5LTset_attribute_string(descriptors->qDirtyH5, PRIMARYDATA, "CLASS", H5IMAGE);
   H5LTset_attribute_string(descriptors->uDirtyH5, PRIMARYDATA, "CLASS", H5IMAGE);
   H5LTset_attribute_string(descriptors->pDirtyH5, PRIMARYDATA, "CLASS", H5IMAGE);
}

/*************************************************************
*
* Read the list of frequencies from the input freq file
*
*************************************************************/
int getFreqList(struct IOFileDescriptors *descriptors,
    struct parameters *params,
    struct DataArrays *data_array) {
    int i;
    float tempFloat;

    params->freqList = calloc(data_array->nFreqList, sizeof(data_array->freqList[0]));
    if(params->freqList == NULL) {
        printf("Error: Mem alloc failed while reading in frequency list\n\n");
        return(FAILURE);
    }
    for(i=0; i<params->nFreq; i++) {
        fscanf(descriptors->freq, "%f", &data_array->freqList[i]);
        if(feof(descriptors->freq)) {
            printf("Error: Frequency values and fits frames don't match\n");
            return(FAILURE);
        }
    }
    fscanf(descriptors->freq, "%f", &tempFloat);
    if(! feof(descriptors->freq)) {
        printf("Error: More frequency values present than fits frames\n\n");
        return(FAILURE);
    }

    /* Compute \lambda^2 from the list of generated frequencies */
    params->lambda2  = calloc(data_array->nFreq, sizeof(data_array->lambda2[0]));
    if(params->lambda2 == NULL) {
        printf("Error: Mem alloc failed while reading in frequency list\n\n");
        return(FAILURE);
    }
    // TODO refactor
    params->lambda20 = 0.0;
    for(i=0; i<data_array->nFreq; i++)
        params->lambda2[i] = (LIGHTSPEED / data_array->freqList[i]) *
                             (LIGHTSPEED / data_array->freqList[i]);

    return(SUCCESS);
}
