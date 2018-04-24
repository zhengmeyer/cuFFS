/******************************************************************************
rmsf.c
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
#include "structures.h"
#include<math.h>

#include "rmsf.h"

/*************************************************************
*
* Generate Rotation Measure Spread Function
*
*************************************************************/
int generateRMSF(struct optionsList *inOptions, struct DataArrays *data_arrays, struct parameters *params) {
    int i, j;

    data_arrays->nPhi = inOptions->nPhi;
    data_arrays->rmsf     = calloc(inOptions->nPhi, sizeof(data_arrays->rmsf));
    data_arrays->rmsfReal = calloc(inOptions->nPhi, sizeof(data_arrays->rmsfReal));
    data_arrays->rmsfImag = calloc(inOptions->nPhi, sizeof(data_arrays->rmsfImag));
    data_arrays->phiAxis  = calloc(inOptions->nPhi, sizeof(data_arrays->phiAxis));

    if(data_arrays->rmsf     == NULL || data_arrays->rmsfReal == NULL ||
       data_arrays->rmsfImag == NULL || data_arrays->phiAxis  == NULL)
        return(FAILURE);

    /* Get the normalization factor K */
    params->K = 1.0 / params->qAxisLen3;

    /* First generate the phi axis */
    for(i=0; i<inOptions->nPhi; i++) {
        data_arrays->phiAxis[i] = inOptions->phiMin + i * inOptions->dPhi;

        /* For each phi value, compute the corresponding RMSF */
        for(j=0; j<params->qAxisLen3; j++) {
            data_arrays->rmsfReal[i] += cos(2 * data_arrays->phiAxis[i] *
                                   (data_arrays->lambda2[j] - data_arrays->lambda20 ));
            data_arrays->rmsfImag[i] -= sin(2 * data_arrays->phiAxis[i] *
                                   (data_arrays->lambda2[j] - data_arrays->lambda20 ));
        }
        // Normalize with K
        data_arrays->rmsfReal[i] *= params->K;
        data_arrays->rmsfImag[i] *= params->K;
        data_arrays->rmsf[i] = sqrt( data_arrays->rmsfReal[i] * data_arrays->rmsfReal[i] +
                                data_arrays->rmsfImag[i] * data_arrays->rmsfImag[i] );
    }
    return(SUCCESS);
}

/*************************************************************
*
* Comparison function used by quick sort.
*
*************************************************************/
int compFunc(const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
}

/*************************************************************
*
* Find the median \lambda^2_0
*
*************************************************************/
void getMedianLambda20(struct parameters *params, struct DataArrays *data_arrays) {
    double *tempArray;
    int i;

    tempArray = calloc(params->qAxisLen3, sizeof(tempArray));
    for(i=0; i<params->qAxisLen3; i++)
        tempArray[i] = data_arrays->lambda2[i];

    /* Sort the list of lambda2 freq */
    qsort(tempArray, params->qAxisLen3, sizeof(tempArray), compFunc);

    /* Find the median value of the sorted list */
    params->lambda20 = tempArray[params->qAxisLen3/2];
    free(tempArray);
}

/*************************************************************
*
* Write RMSF to disk
*
*************************************************************/
int writeRMSF(struct optionsList inOptions, struct DataArrays data_arrays) {
    FILE *rmsf;
    char filename[FILENAME_LEN];
    int i;

    /* Open a text file */
    sprintf(filename, "%srmsf.txt", inOptions.outPrefix);
    printf("INFO: Writing RMSF to %s\n", filename);
    rmsf = fopen(filename, FILE_READWRITE);
    if(rmsf == NULL)
        return(FAILURE);

    for(i=0; i<inOptions.nPhi; i++)
        fprintf(rmsf, "%f\t%f\t%f\t%f\n", data_arrays.phiAxis[i], data_arrays.rmsfReal[i],
                data_arrays.rmsfImag[i], data_arrays.rmsf[i]);

    fclose(rmsf);
    return(SUCCESS);
}

#ifdef GNUPLOT_ENABLE
/*************************************************************
*
* Plot RMSF
*
*************************************************************/
int plotRMSF(struct optionsList inOptions) {
    FILE *gnuplotPipe;
    char commands[STRING_BUF_LEN];

    gnuplotPipe = popen("gnuplot -persist", FILE_READWRITE);
    if(gnuplotPipe == NULL)
        return(FAILURE);

    /* Plot the RMSF using the file that was written in writeRMSF() */
    sprintf(commands, "set title \"Rotation Measure Spread Function\"\n");
    sprintf(commands, "%sset xlabel \"Faraday Depth\"\n", commands);
    sprintf(commands, "%sset autoscale\n", commands);
    sprintf(commands,"%splot \"%srmsf.txt\" using 1:2 title 'RMSF' with lines,",
            commands, inOptions.outPrefix);
    sprintf(commands, "%s \"%srmsf.txt\" using 1:3 title 'Real' with lines,",
            commands, inOptions.outPrefix);
    sprintf(commands, "%s \"%srmsf.txt\" using 1:4 title 'Imag' with lines\n",
            commands, inOptions.outPrefix);
    fprintf(gnuplotPipe, "%s", commands);
    pclose(gnuplotPipe);

    return(SUCCESS);
}
#endif
