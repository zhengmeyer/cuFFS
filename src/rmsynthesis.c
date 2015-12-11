/******************************************************************************

A GPU Based implementation of RM Synthesis

Version: 0.1
Last edited: July 11, 2015

Version history:
================
v0.1    Assume FITS cubes have at least 3 frames with frequency being the 
         3rd axis. Also, implicitly assumes each freq channel has equal weight.
         

******************************************************************************/
#include"rmsynthesis.h"

/*************************************************************
*
* Parse the input file and extract the relevant keywords
*
*************************************************************/
struct optionsList parseInput(char *parsetFileName) {
    config_t cfg;
    struct optionsList inOptions;
    const char *str;
    
    /* Initialize configuration */
    config_init(&cfg);
    
    /* Read in the configuration file */
    if(!config_read_file(&cfg, parsetFileName)) {
        printf("\nError: Error reading parset file. Exiting with message: %s\n\n", 
               config_error_text(&cfg));
        config_destroy(&cfg);
        exit(FAILURE);
    }
    
    /* Get the names of fits files */
    if(config_lookup_string(&cfg, "qCubeName", &str)) {
        inOptions.qCubeName = malloc(strlen(str));
        strcpy(inOptions.qCubeName, str);
    }
    else {
        printf("\nError: 'qCubeName' undefined in parset");
        config_destroy(&cfg);
        exit(FAILURE);
    }
    if(config_lookup_string(&cfg, "uCubeName", &str)) {
        inOptions.uCubeName = malloc(strlen(str));
        strcpy(inOptions.uCubeName, str);
    }
    else {
        printf("\nError: 'uCubeName' undefined in parset");
        config_destroy(&cfg);
        exit(FAILURE);
    }
    
    /* Get the name of the frequency file */
    if(config_lookup_string(&cfg, "freqFileName", &str)) {
        inOptions.freqFileName = malloc(strlen(str));
        strcpy(inOptions.freqFileName, str);
    }
    else {
        printf("\nError: 'freqFileName' undefined in parset");
        config_destroy(&cfg);
        exit(FAILURE);
    }

    /* Check if an image mask is defined */
    if(! config_lookup_string(&cfg, "imageMask", &str)) {
        printf("\nINFO: Image mask not specified");
        inOptions.isImageMaskDefined = FALSE;
    }
    else {
        inOptions.imageMask = malloc(strlen(str));
        strcpy(inOptions.imageMask, str);
        inOptions.isImageMaskDefined = TRUE;
    }

    /* Get prefix for output files */
    if(config_lookup_string(&cfg, "outPrefix", &str)) {
        inOptions.outPrefix = malloc(strlen(str));
        strcpy(inOptions.outPrefix, str);
    }
    else {
        printf("\nINFO: 'outPrefix' is not defined. Defaulting to %s", DEFAULT_OUT_PREFIX);
        inOptions.outPrefix = malloc(strlen(DEFAULT_OUT_PREFIX));
        strcpy(inOptions.outPrefix, DEFAULT_OUT_PREFIX);
    }
    
    /* Get Faraday depth */
    if(! config_lookup_float(&cfg, "phiMin", &inOptions.phiMin)) {
        printf("\nError: 'phiMin' undefined in parset");
        config_destroy(&cfg);
        exit(FAILURE);
    }
    if(! config_lookup_float(&cfg, "dPhi", &inOptions.dPhi)) {
        printf("\nError: 'dPhi' undefined in parset");
        config_destroy(&cfg);
        exit(FAILURE);
    }
    if(! config_lookup_int(&cfg, "nPhi", &inOptions.nPhi)) {
        printf("\nError: 'nPhi' undefined in parset");
        config_destroy(&cfg);
        exit(FAILURE);
    }
    
    config_destroy(&cfg);
    return(inOptions);
}

/*************************************************************
*
* Print parsed input to screen
*
*************************************************************/
void printOptions(struct optionsList inOptions) {
    int i;
    
    printf("\n");
    for(i=0; i<SCREEN_WIDTH; i++) { printf("#"); }
    printf("\n");
    printf("\nQ Cube: %s", inOptions.qCubeName);
    printf("\nU Cube: %s", inOptions.uCubeName);
    printf("\n");
    if(inOptions.isImageMaskDefined == TRUE)
        printf("\nImage mask: %s\n", inOptions.imageMask);
    printf("\nphi min: %.2f", inOptions.phiMin);
    printf("\n# of phi planes: %d", inOptions.nPhi);
    printf("\ndelta phi: %.2lf", inOptions.dPhi);
    printf("\n\n");
    for(i=0; i<SCREEN_WIDTH; i++) { printf("#"); }
    printf("\n");
}

/*************************************************************
*
* Read header information from the fits files
*
*************************************************************/
int getFitsHeader(struct optionsList *inOptions, struct parList *params) {
    int fitsStatus = SUCCESS;
    char fitsComment[FLEN_COMMENT];
    
    /* Get the image dimensions from the Q cube */
    fits_read_key(params->qFile, TINT, "NAXIS", &params->qAxisNum, fitsComment, &fitsStatus);
    fits_read_key(params->qFile, TINT, "NAXIS1", &params->qAxisLen1, fitsComment, &fitsStatus);
    fits_read_key(params->qFile, TINT, "NAXIS2", &params->qAxisLen2, fitsComment, &fitsStatus);
    fits_read_key(params->qFile, TINT, "NAXIS3", &params->qAxisLen3, fitsComment, &fitsStatus);
    /* Get the image dimensions from the Q cube */
    fits_read_key(params->uFile, TINT, "NAXIS", &params->uAxisNum, fitsComment, &fitsStatus);
    fits_read_key(params->uFile, TINT, "NAXIS1", &params->uAxisLen1, fitsComment, &fitsStatus);
    fits_read_key(params->uFile, TINT, "NAXIS2", &params->uAxisLen2, fitsComment, &fitsStatus);
    fits_read_key(params->uFile, TINT, "NAXIS3", &params->uAxisLen3, fitsComment, &fitsStatus);
    
    return(fitsStatus);
}

/*************************************************************
*
* Read the list of frequencies from the input freq file
*
*************************************************************/
int getFreqList(struct optionsList *inOptions, struct parList *params) {
    int i;
    float tempFloat;
    
    params->freqList = calloc(params->qAxisLen3, sizeof(params->freqList));
    if(params->freqList == NULL) {
        printf("\nError: Mem alloc failed while reading in frequency list\n\n");
        return(FAILURE);
    }
    for(i=0; i<params->qAxisLen3; i++) {
        fscanf(params->freq, "%f", &params->freqList[i]);
        if(feof(params->freq)) {
            printf("\nError: Less frequency values present than fits frames\n\n");
            return(FAILURE);
        }
    }
    fscanf(params->freq, "%f", &tempFloat);
    if(! feof(params->freq)) {
        printf("\nError: More frequency values present than fits frames\n\n");
        return(FAILURE);
    }
    
    /* Compute \lambda^2 from the list of generated frequencies */
    params->lambda2  = calloc(params->qAxisLen3, sizeof(params->lambda2));
    if(params->lambda2 == NULL) {
        printf("\nError: Mem alloc failed while reading in frequency list\n\n");
        return(FAILURE);
    }
    params->lambda20 = 0.0;
    for(i=0; i<params->qAxisLen3; i++)
        params->lambda2[i] = (LIGHTSPEED / params->freqList[i]) * 
                             (LIGHTSPEED / params->freqList[i]);
    
    return(SUCCESS);
}

/*************************************************************
*
* Generate Rotation Measure Spread Function
*
*************************************************************/
int generateRMSF(struct optionsList *inOptions, struct parList *params) {
    int i, j;
    double K;
    
    params->rmsf     = calloc(inOptions->nPhi, sizeof(params->rmsf));
    params->rmsfReal = calloc(inOptions->nPhi, sizeof(params->rmsfReal));
    params->rmsfImag = calloc(inOptions->nPhi, sizeof(params->rmsfImag));
    params->phiAxis  = calloc(inOptions->nPhi, sizeof(params->phiAxis));
    
    if(params->rmsf     == NULL || params->rmsfReal == NULL ||
       params->rmsfImag == NULL || params->phiAxis  == NULL)
        return(FAILURE);
    
    /* Get the normalization factor K */
    K = 1.0 / params->qAxisLen3;
    
    /* First generate the phi axis */
    for(i=0; i<inOptions->nPhi; i++) {
        params->phiAxis[i] = inOptions->phiMin + i * inOptions->dPhi;
        
        /* For each phi value, compute the corresponding RMSF */
        for(j=0; j<params->qAxisLen3; j++) {
            params->rmsfReal[i] += cos(2 * params->phiAxis[i] *
                                   (params->lambda2[j] - params->lambda20 ));
            params->rmsfImag[i] -= sin(2 * params->phiAxis[i] *
                                   (params->lambda2[j] - params->lambda20 ));
        }
        // Normalize with K
        params->rmsfReal[i] *= K;
        params->rmsfImag[i] *= K;
        params->rmsf[i] = sqrt( params->rmsfReal[i] * params->rmsfReal[i] +
                                params->rmsfImag[i] * params->rmsfImag[i] );
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
void getMedianLambda20(struct parList *params) {
    double *tempArray;
    int i;
    
    tempArray = calloc(params->qAxisLen3, sizeof(tempArray));
    for(i=0; i<params->qAxisLen3; i++)
        tempArray[i] = params->lambda2[i];
    
    /* Sort the list of lambda2 freq */
    qsort(tempArray, params->qAxisLen3, sizeof(tempArray), compFunc);
    
    /* Find the median value of the sorted list */
    params->lambda20 = tempArray[params->qAxisLen3/2];
}

/*************************************************************
*
* Write RMSF to disk
*
*************************************************************/
int writeRMSF(struct optionsList inOptions, struct parList params) {
    FILE *rmsf;
    char filename[FILENAME_LEN];
    int i;
    
    /* Open a text file */
    sprintf(filename, "%srmsf.txt", inOptions.outPrefix);
    printf("\nINFO: Writing RMSF to %s", filename);
    rmsf = fopen(filename, FILE_READWRITE);
    if(rmsf == NULL)
        return(FAILURE);
    
    for(i=0; i<inOptions.nPhi; i++)
        fprintf(rmsf, "%f\t%f\t%f\t%f\n", params.phiAxis[i], params.rmsfReal[i],
                params.rmsfImag[i], params.rmsf[i]);
    
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
    sprintf(commands, "%splot \"%srmsf.txt\" using 1:2 title 'RMSF' with lines,", 
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

/*************************************************************
*
* Read in the Stokes-Q and -U images
*
*************************************************************/
int getImageData(struct optionsList *inOptions, struct parList *params) {
    long *fPixel;
    int i;
    LONGLONG nElements;
    int status = SUCCESS;
    
    /* Set the starting pixel to read from the FITS file */
    fPixel = calloc(params->qAxisNum, sizeof(fPixel));
    for(i=1; i<=params->qAxisNum; i++)
        fPixel[i-1] = 1;
    
    /* Allocate memory to store the Q and U image arrays */
    nElements = params->qAxisLen1 * params->qAxisLen2 * params->qAxisLen3;
    params->qImageArray = calloc(nElements, sizeof(params->qImageArray));
    params->uImageArray = calloc(nElements, sizeof(params->uImageArray));
    
    /* Read pixel values */
    fits_read_pix(params->qFile, TFLOAT, fPixel, nElements, NULL, params->qImageArray, NULL, &status);
    fits_read_pix(params->uFile, TFLOAT, fPixel, nElements, NULL, params->uImageArray, NULL, &status);
    fits_close_file(params->qFile, &status);
    fits_close_file(params->uFile, &status);
    if(status) {
        fits_report_error(stdout, status);
        return(FAILURE);
    }
    return(SUCCESS);
}

/*************************************************************
*
* Read in the input mask image
*
*************************************************************/
int getImageMask(struct optionsList *inOptions, struct parList *params) {
    int fitsStatus = SUCCESS;
    char fitsComment[FLEN_COMMENT];

    fits_read_key(params->maskFile, TINT, "NAXIS1", &params->maskAxisLen1, fitsComment, &fitsStatus);
    fits_read_key(params->maskFile, TINT, "NAXIS2", &params->maskAxisLen2, fitsComment, &fitsStatus);

    /*fits_read_key(params->qFile, TINT, "NAXIS", &params->qAxisNum, fitsComment, &fitsStatus);*/
    return(SUCCESS);
}

/*************************************************************
*
* Check for valid CUDA supported devices. If detected, 
*  print useful device information
*
*************************************************************/
struct deviceInfoList * getDeviceInformation(int *nDevices) {
    int dev;
    cudaError_t errorID;
    int deviceCount = NO_DEVICE;
    struct cudaDeviceProp deviceProp;
    struct deviceInfoList *gpuList;
    
    /* Check for valid devices */
    cudaDeviceReset();
    errorID = cudaGetDeviceCount(&deviceCount);
    if(errorID != cudaSuccess) {
        printf("cudaGetDeviceCount returned %d\n%s", (int)errorID, 
               cudaGetErrorString(errorID));
        exit(FAILURE);
    }
    if(deviceCount == NO_DEVICE) {
        printf("\nError: Could not detect CUDA supported GPU(s)\n\n");
        exit(FAILURE);
    }
    printf("\nINFO: Detected %d CUDA-supported GPU(s)\n", deviceCount);
    *nDevices = deviceCount;

    /* Store useful information about each GPU in a structure array */
    gpuList = malloc(deviceCount * sizeof(struct deviceInfoList));
    for(dev=0; dev < deviceCount; dev++) {
        cudaSetDevice(dev);
        cudaGetDeviceProperties(&deviceProp, dev);
        printf("\nDevice %d: %s", dev, deviceProp.name);
        gpuList[dev].deviceID    = dev;
        gpuList[dev].globalMem   = deviceProp.totalGlobalMem;
        gpuList[dev].constantMem = deviceProp.totalConstMem;
        gpuList[dev].sharedMemPerBlock = deviceProp.sharedMemPerBlock;
        gpuList[dev].maxThreadPerMP = deviceProp.maxThreadsPerMultiProcessor;
        gpuList[dev].maxThreadPerBlock = deviceProp.maxThreadsPerBlock;
        gpuList[dev].threadBlockSize[0] = deviceProp.maxThreadsDim[0];
        gpuList[dev].threadBlockSize[1] = deviceProp.maxThreadsDim[1];
        gpuList[dev].threadBlockSize[2] = deviceProp.maxThreadsDim[2];
    }
    printf("\n");
    return(gpuList);
}

/*************************************************************
*
* Main code
*
*************************************************************/
int main(int argc, char *argv[]) {
    /* Variable declaration */
    char *parsetFileName = argv[1];
    struct optionsList inOptions;
    struct parList params;
    int fitsStatus;
    int status;
    int nDevices, deviceID;
    struct deviceInfoList *gpuList;
    int i;

    fitsfile *x;
    
    printf("\nRM Synthesis v%s", VERSION_STR);
    printf("\nWritten by Sarrvesh S. Sridhar\n");
    
    /* Verify command line input */
    if(argc!=2) {
        printf("\nERROR: Invalid command line input. Terminating Execution!");
        printf("\nUsage: %s <parset filename>\n\n", argv[0]);
        return(FAILURE);
    } 
    if(strcmp(parsetFileName, "-h") == 0) {
        /* Print help and exit */
        printf("\nUsage: %s <parset filename>\n\n", argv[0]);
        return(SUCCESS);
    }
    
    /* Parse the input file */
    printf("\nINFO: Parsing input file %s", parsetFileName);
    inOptions = parseInput(parsetFileName);
    
    /* Print parset input options to screen */
    printOptions(inOptions);
    
    /* Retreive information about all connected GPU devices */
    gpuList = getDeviceInformation(&nDevices);
    
    /* Open the input files */
    printf("\nINFO: Accessing the input files");
    printf("\nWARN: Assuming the 3rd axis in the fits files is the frequency axis");
    fitsStatus = SUCCESS;
    fits_open_file(&params.qFile, inOptions.qCubeName, READONLY, &fitsStatus);
    fits_open_file(&params.uFile, inOptions.uCubeName, READONLY, &fitsStatus);
    if(fitsStatus != SUCCESS) {
        fits_report_error(stdout, fitsStatus);
        return(fitsStatus);
    }
    params.freq = fopen(inOptions.freqFileName, FILE_READONLY);
    if(params.freq == NULL) {
        printf("\nError: Unable to open the frequency file\n\n");
        return(FAILURE);
    }
    if(inOptions.isImageMaskDefined == TRUE) {
        printf("\nINFO: Accessing the input image mask %s", inOptions.imageMask);
        fitsStatus = SUCCESS;
        fits_open_file(&params.maskFile, inOptions.imageMask, READONLY, &fitsStatus);
        if(fitsStatus != SUCCESS) {
            fits_report_error(stdout, fitsStatus);
            return(fitsStatus);
        }
    }
    
    /* Gather information from input image fits header */
    status = getFitsHeader(&inOptions, &params);
    if(status) {
        fits_report_error(stdout, status);
        return(FAILURE);
    }
    
    /* Read frequency list */
    if(getFreqList(&inOptions, &params))
        return(FAILURE);
    
    /* Find median lambda20 */
    getMedianLambda20(&params);
    
    /* Generate RMSF */
    printf("\nINFO: Computing RMSF");
    if(generateRMSF(&inOptions, &params)) {
        printf("\nError: Mem alloc failed while generating RMSF");
        return(FAILURE);
    }    
    /* Write RMSF to disk */
    if(writeRMSF(inOptions, params)) {
        printf("\nError: Unable to write RMSF to disk\n\n");
        return(FAILURE);
    }
    /* Plot RMSF */
    #ifdef GNUPLOT_ENABLE
    printf("\nINFO: Plotting RMSF with gnuplot");
    if(plotRMSF(inOptions)) {
        printf("\nError: Unable to plot RMSF\n\n");
        return(FAILURE);
    }
    #endif
    
    /* Read image planes from the Q and U cubes */
    printf("\nINFO: Reading in FITS images");
    if(getImageData(&inOptions, &params))
        return(FAILURE);

    /* Read image mask */
    printf("\nINFO: Reading in input image mask");
    if(inOptions.isImageMaskDefined == TRUE) {
        printf("\nImage mask is defined. Use the image mask.");
        if(getImageMask(&inOptions, &params))
            return(FAILURE);
    }
    else
        printf("\nImage mask is not defined. Making a mask to include all pixels.");

    printf("\n\n");
    return(SUCCESS);
}
