/*************************************************************
* Allocates host memory
*
*************************************************************/
int allocateHostMemoryForComputation(float **lambdaDiff2Pointer, float **qImageArrayPointer, float **uImageArrayPointer,
									 float **qPhi, float **uPhi, float **pPhi,
									 int nInFrequencies, int nRA, int nPhi,  int nInElements, int nOutElements){
	/* Allocate memory on the host */
	t->startProc = clock();

	*lambdaDiff2Pointer = (float *)calloc(nInFrequencies, sizeof(**lambdaDiff2));
	*qImageArrayPointer = (float *)calloc(nInElements, sizeof(**qImageArray));
	*uImageArrayPointer = (float *)calloc(nInElements, sizeof(**uImageArray));
	*qPhiPointer = (float *)calloc(nOutElements, sizeof(**qPhi));
	*uPhiPointer = (float *)calloc(nOutElements, sizeof(**uPhi));
	*pPhiPointer = (float *)calloc(nOutElements, sizeof(**pPhi));
	if(*lambdaDiff2Pointer == NULL || *qImageArrayPointer == NULL ||
			*uImageArrayPointer == NULL || *qPhiPointer == NULL ||
			*uPhiPointer == NULL || *pPhiPointer == NULL) {
		printf("ERROR: Unable to allocate memory on host\n");
		return FAILURE;
	}

}

int computeLambdaSquareDifference(float lambdaDiff2, float *lambda2, float lambda20, size_t size){
	/* Compute \lambda^2 - \lambda^2_0 once. Common for all threads */
	for(i=0;i<size;i++) lambdaDiff2[i] = 2.0*(lambda2[i]-lambda20);
}

int allocateDeviceMemoryForComputation(int cudaDeviceId, float **lambdaDiff2, float **phiAxis,
		float **qImageArray, float **uImageArray,
		float **qPhi, float **uPhi, float **pPhi,
		int nInElements, int nOutElements, int nFrequencies, int nPhi){
	int currentCudaDevice;
	// reads the current selected device and switch back after
	cudaGetDevice(&currentCudaDevice);

	// set the desired device
	cudaSetDevice(deviceId);

	/* Allocate memory on the device */
	cudaMalloc(&lambdaDiff2, sizeof(**lambdaDiff2)*nFrequencies);
	cudaMalloc(&phiAxis, sizeof(**phiAxis)*inOptions->nPhi);
	cudaMalloc(&qImageArray, nInElements*sizeof(**qImageArray));
	cudaMalloc(&uImageArray, nInElements*sizeof(**uImageArray));
	cudaMalloc(&qPhi, nOutElements*sizeof(**qPhi));
	cudaMalloc(&uPhi, nOutElements*sizeof(**uPhi));
	cudaMalloc(&pPhi, nOutElements*sizeof(**pPhi));
	checkCudaError();
	// switch back to the previous device
	cudaSetDevice(currentCudaDevice);
}

int copyLambdaDifferenceToDevice(int deviceId, float *lambdaDiff2, float *device_lambdaDiff2, size_t length){
	int currentCudaDevice;
	// reads the current selected device and switch back after
	cudaGetDevice(&currentCudaDevice);
	// set the desired device
	cudaSetDevice(deviceId);

	/* Transfer \lambda^2 - \lambda^2_0 to device */
	cudaMemcpy(device_lambdaDiff2, lambdaDiff2, sizeof(*lambdaDiff)*length, cudaMemcpyHostToDevice);
	checkCudaError();

	// switch back to the previous device
	cudaSetDevice(currentCudaDevice);
}

int copyPhiToDevice(int deviceId, float *phi, float *device_phi, size_t length){
	int currentCudaDevice;
	// reads the current selected device and switch back after
	cudaGetDevice(&currentCudaDevice);
	// set the desired device
	cudaSetDevice(deviceId);

	/* Allocate and transfer phi axis info. Common for all threads */
	cudaMemcpy(device_phi, phi, sizeof(*phi)*length, cudaMemcpyHostToDevice);
	checkCudaError();

	// switch back to the previous device
	cudaSetDevice(currentCudaDevice);
}

int copyStepToDevice(int deviceId, float *qImageArray, float *d_qImageArray float *uImageArray, float *d_uImageArray, size_t length){
	int currentCudaDevice;
	// reads the current selected device and switch back after
	cudaGetDevice(&currentCudaDevice);
	// set the desired device
	cudaSetDevice(deviceId);

	cudaMemcpy(d_qImageArray, qImageArray,
	                  nInElements*sizeof(*qImageArray),
	                  cudaMemcpyHostToDevice);
	cudaMemcpy(d_uImageArray, uImageArray,
	                  nInElements*sizeof(*qImageArray),
	                  cudaMemcpyHostToDevice);
	checkCudaError();

	// switch back to the previous device
	cudaSetDevice(currentCudaDevice);
}

int copyStepToDevice(int deviceId, float *qImageArray, float *d_qImageArray float *uImageArray, float *d_uImageArray, size_t length){
	int currentCudaDevice;
	// reads the current selected device and switch back after
	cudaGetDevice(&currentCudaDevice);
	// set the desired device
	cudaSetDevice(deviceId);

	cudaMemcpy(d_qImageArray, qImageArray,
	                  nInElements*sizeof(*qImageArray),
	                  cudaMemcpyHostToDevice);
	cudaMemcpy(d_uImageArray, uImageArray,
	                  nInElements*sizeof(*qImageArray),
	                  cudaMemcpyHostToDevice);
	checkCudaError();

	// switch back to the previous device
	cudaSetDevice(currentCudaDevice);
}


/*************************************************************
*
* GPU accelerated RM Synthesis function
*
*************************************************************/
extern "C"
int doRMSynthesis(struct optionsList *inOptions,
    struct IOFileDescriptors *descriptors,
    struct parameters *params,
    strcut DataArrays *arrays,
    struct deviceInfoList selectedDeviceInfo,
    struct timeInfoList *t) {

    int i, j;

    float *lambdaDiff2, *d_lambdaDiff2;

    float *qImageArray, *uImageArray;
    float *d_qImageArray, *d_uImageArray;

    float *qPhi, *uPhi, *pPhi;
    float *d_qPhi, *d_uPhi, *d_pPhi;

    float *d_phiAxis;

    // Dimension sizes
    int nFrequencies, int nRa, int nDec, int nPhi;
    long nInElements, nOutElements;

    dim3 calcThreadSize, calcBlockSize;
    long *fPixel;
    int fitsStatus = 0;

    herr_t h5ErrorQ, h5ErrorU, h5ErrorP;
    herr_t qerror, uerror, perror;
    hsize_t offsetIn[N_DIMS], countIn[N_DIMS], dimIn;
    hsize_t offsetOut[N_DIMS], countOut[N_DIMS], dimOut;

    /* Set mode-specific configuration */
    switch(inOptions->fileFormat) {
       case FITS:
          /* For FITS, set some pixel access limits */
          fPixel = (long *)calloc(params->qAxisNum, sizeof(*fPixel));
          fPixel[0] = 1; fPixel[1] = 1;
          /* Determine what the appropriate block and grid sizes are */
          calcThreadSize.x = selectedDeviceInfo.warpSize;
          calcBlockSize.x  = params->qAxisLen1; // Number of RA or LOS in this frame
          calcBlockSize.y  = inOptions->nPhi/calcThreadSize.x + 1;
          printf("INFO: Launching %dx%d blocks each with %d threads\n",
                  calcBlockSize.x, calcBlockSize.y, calcThreadSize.x);
          break;
       case HDF5:
          /* For HDF5, set up the hyperslab and data subset for input */
          dimIn = params->qAxisLen2 * params->qAxisLen3;
          descriptors.qDataset   = H5Dopen2(descriptors->qFileh5, PRIMARYDATA, H5P_DEFAULT);
          descriptors.qDataspace = H5Dget_space(qDataset);
          descriptors.uDataset   = H5Dopen2(descriptors->uFileh5, PRIMARYDATA, H5P_DEFAULT);
          descriptors.uDataspace = H5Dget_space(uDataset);
          countIn[0] = params->qAxisLen3;
          countIn[1] = 1; countIn[2] = params->qAxisLen2;
          offsetIn[0] = 0; offsetIn[1] = 0; offsetIn[2] = 0;
          descriptors.qMemspace = H5Screate_simple(1, &dimIn, NULL);
          descriptors.uMemspace = H5Screate_simple(1, &dimIn, NULL);
          if( qDataset<0 || uDataset<0 || qDataspace<0 || uDataspace<0 || qMemspace<0 || uMemspace<0 )
          { printf("\nError: HDF5 allocation failed\n"); }

          /* Set up the hyperslab and data subset for output */
          dimOut = params->qAxisLen2 * inOptions->nPhi;
          descriptors.qOutDataset   = H5Dopen2(descriptors->qDirtyH5, PRIMARYDATA, H5P_DEFAULT);
          descriptors.qOutDataspace = H5Dget_space(qOutDataset);
          descriptors.uOutDataset   = H5Dopen2(descriptors->uDirtyH5, PRIMARYDATA, H5P_DEFAULT);
          descriptors.uOutDataspace = H5Dget_space(uOutDataset);
          descriptors.pOutDataset   = H5Dopen2(descriptors->pDirtyH5, PRIMARYDATA, H5P_DEFAULT);
          descriptors.pOutDataspace = H5Dget_space(uOutDataset);
          countOut[0] = inOptions->nPhi;
          countOut[1] = 1; countOut[2] = descriptors->qAxisLen2;
          offsetOut[0] = 0; offsetOut[1] = 0; offsetOut[2] = 0;
          descriptors.qOutMemspace = H5Screate_simple(1, &dimOut, NULL);
          descriptors.uOutMemspace = H5Screate_simple(1, &dimOut, NULL);
          descriptors.pOutMemspace = H5Screate_simple(1, &dimOut, NULL);
          if( descriptors.qOutDataset<0 || descriptors.uOutDataset<0 || descriptors.pOutDataset<0 ||
              descriptors.qOutDataspace<0 || descriptors.uOutDataspace<0 || descriptors.pOutDataspace<0 ||
              descriptors.qOutMemspace<0 || descriptors.uOutMemspace<0 || descriptors.pOutMemspace<0 ) {
             printf("\nError: HDF5 output allocation failed\n");
             exit(FAILURE);
          }

          /* Determine what the appropriate block and grid sizes are */
          calcThreadSize.x = selectedDeviceInfo.warpSize;
          calcBlockSize.y  = params->qAxisLen2;
          calcBlockSize.x  = inOptions->nPhi/calcThreadSize.x + 1;
          printf("INFO: Launching %dx%d blocks each with %d threads\n",
                 calcBlockSize.x, calcBlockSize.y, calcThreadSize.x);
          break;
    }
    t->startProc = clock();

    // Computes the dimension of the computation
    nRa = params->qAxisLen1;
    nDec = params->qAxisLen2;
    nFrequencies = params->qAxisLen3;
    nPhi = inOptions->nPhi;
    nInElements = nInFrequencies * nRA;
    nOutElements= nPhi * nRa;


    /* Allocate memory on the host */
    if (allocateHostMemoryForComputation(&lambdaDiff2, &qImageArray, &uImageArray, // input arrays
    		&qPhi, &uPhi, &pPhi,                                                   // output arrays
			nFrequencies, nRa, nPhi,				                               // dimensions of the axes
			nInElements, nOutElements) == FAILURE) exit(FAILURE);


    /* Allocate memory on the device */
    allocateDeviceMemoryForComputation(selectedDeviceInfo.deviceID, &d_lambdaDiff2, &d_phiAxis,
    		&d_qImageArray, &d_uImageArray,
    		&d_qPhi, &d_uPhi, &d_pPhi,
			nInElements, nOutElements, nRa, nPhi){){

    /* Compute \lambda^2 - \lambda^2_0 once. Common for all threads */
	computeLambdaSquareDifference(lambdaDiff2, data_arrays->lambda2, data_arrays->lambda20, nFrequencies);


    t->stopProc = clock();
    t->msProc += ((float)(t->stopProc - t->startProc))/CLOCKS_PER_SEC;

    /* Transfer \lambda^2 - \lambda^2_0 to device */
    t->startX = clock();

    copyLambdaDifferenceToDevice(selectedDeviceInfo.deviceID, lambdaDiff2, d_lambdaDiff2, nFrequencies);


    /* Allocate and transfer phi axis info. Common for all threads */
    copyPhiToDevice(selectedDeviceInfo.deviceID, data_arrays->phiAxis, d_phiAxis, nPhi);

    t->stopX = clock();
    t->msX += ((float)(t->stopX - t->startX))/CLOCKS_PER_SEC;

    /* Process each line of sight individually */
    //cudaEventRecord(totStart);
    for(j=1; j<=params->qAxisLen2; j++) {
       /* Read one frame at a time. In the original cube, this is
          all sightlines in one DEC row */
       t->startRead = clock();
       switch(inOptions->fileFormat) {
          case FITS:
             fPixel[2] = j;
             fits_read_pix(descriptors->qFile, TFLOAT, fPixel, nInElements, NULL,
                           qImageArray, NULL, &fitsStatus);
             fits_read_pix(descriptors->uFile, TFLOAT, fPixel, nInElements, NULL,
                           uImageArray, NULL, &fitsStatus);
             checkFitsError(fitsStatus);
             break;
          case HDF5:
             offsetIn[1] = j-1;
             qerror = H5Sselect_hyperslab(descriptors.qDataspace, H5S_SELECT_SET, offsetIn,
                                    NULL, countIn, NULL);
             uerror = H5Sselect_hyperslab(descriptors.uDataspace, H5S_SELECT_SET, offsetIn,
                                    NULL, countIn, NULL);
             h5ErrorQ = H5Dread(descriptors.qDataset, H5T_NATIVE_FLOAT, qMemspace,
                                   descriptors.qDataspace, H5P_DEFAULT, qImageArray);
             h5ErrorU = H5Dread(descriptors.uDataset, H5T_NATIVE_FLOAT, uMemspace,
                                   descriptors.uDataspace, H5P_DEFAULT, uImageArray);
             if(h5ErrorQ<0 || h5ErrorU<0 || qerror<0 || uerror<0 ) {
                printf("\nError: Unable to read input data cubes\n\n");
                exit(FAILURE);
             }
             break;
       }
       t->stopRead = clock();
       t->msRead += ((float)(t->stopRead - t->startRead))/CLOCKS_PER_SEC;

       /* Transfer input images to device */
       t->startX = clock();

       copyStepToDevice(selectedDeviceInfo.deviceID, qImageArray, d_qImageArray uImageArray, d_uImageArray, nInElements);

       t->stopX = clock();
       t->msX += ((float)(t->stopX - t->startX))/CLOCKS_PER_SEC;

       /* Launch kernels to compute Q(\phi), U(\phi), and P(\phi) */
       t->startProc = clock();
       switch(inOptions->fileFormat) {
       case FITS:
          computeQUP_fits<<<calcBlockSize, calcThreadSize>>>(d_qImageArray,
                   d_uImageArray, params->qAxisLen3, params->nPhi,
                   params->K, d_qPhi, d_uPhi, d_pPhi, d_phiAxis, d_lambdaDiff2);
          break;
       case HDF5:
          computeQUP_hdf5<<<calcBlockSize, calcThreadSize>>>(d_qImageArray, d_uImageArray,
                         params->qAxisLen2, params->qAxisLen3, params->K, d_qPhi,
                         d_uPhi, d_pPhi, d_phiAxis, inOptions->nPhi, d_lambdaDiff2);
          break;
       }
       cudaThreadSynchronize();
       t->stopProc = clock();
       t->msProc += ((float)(t->stopProc - t->startProc))/CLOCKS_PER_SEC;

       /* Move Q(\phi), U(\phi) and P(\phi) to host */
       t->startX = clock();
       cudaMemcpy(qPhi, d_qPhi, nOutElements*sizeof(*qPhi), cudaMemcpyDeviceToHost);
       cudaMemcpy(uPhi, d_uPhi, nOutElements*sizeof(*qPhi), cudaMemcpyDeviceToHost);
       cudaMemcpy(pPhi, d_pPhi, nOutElements*sizeof(*qPhi), cudaMemcpyDeviceToHost);
       t->stopX = clock();
       t->msX += ((float)(t->stopX - t->startX))/CLOCKS_PER_SEC;

       /* Write the output cubes to disk */
       t->startWrite = clock();
       switch(inOptions->fileFormat) {
          case FITS:
             fits_write_pix(descriptors->qDirty, TFLOAT, fPixel, nOutElements, qPhi, &fitsStatus);
             fits_write_pix(descriptors->uDirty, TFLOAT, fPixel, nOutElements, uPhi, &fitsStatus);
             fits_write_pix(descriptors->pDirty, TFLOAT, fPixel, nOutElements, pPhi, &fitsStatus);
             checkFitsError(fitsStatus);
             break;
          case HDF5:
             offsetOut[1] = j-1;
             qerror = H5Sselect_hyperslab(descriptors.qOutDataspace, H5S_SELECT_SET,
                                          offsetOut, NULL, countOut, NULL);
             uerror = H5Sselect_hyperslab(descriptors.uOutDataspace, H5S_SELECT_SET,
                                          offsetOut, NULL, countOut, NULL);
             perror = H5Sselect_hyperslab(descriptors.pOutDataspace, H5S_SELECT_SET,
                                          offsetOut, NULL, countOut, NULL);
             h5ErrorQ = H5Dwrite(descriptors.qOutDataset, H5T_NATIVE_FLOAT, descriptors.qOutMemspace,
                                 descriptors.qOutDataspace, H5P_DEFAULT, qPhi);
             h5ErrorU = H5Dwrite(descriptors.uOutDataset, H5T_NATIVE_FLOAT, descriptors.uOutMemspace,
                                 descriptors.uOutDataspace, H5P_DEFAULT, uPhi);
             h5ErrorP = H5Dwrite(descriptors.pOutDataset, H5T_NATIVE_FLOAT, descriptors.pOutMemspace,
                                 descriptors.pOutDataspace, H5P_DEFAULT, pPhi);
             if(h5ErrorQ<0 || h5ErrorU || h5ErrorP<0 ||
                qerror<0 || uerror<0 || perror<0 ) {
                printf("\nError: Unable to write output data cubes\n\n");
                exit(FAILURE);
             }
             break;
       }
       t->stopWrite = clock();
       t->msWrite += ((float)(t->stopWrite - t->startWrite))/CLOCKS_PER_SEC;
    }

    /* Free all the allocated memory */
    free(qImageArray); free(uImageArray);
    cudaFree(d_qImageArray); cudaFree(d_uImageArray);
    free(qPhi); free(uPhi); free(pPhi);
    cudaFree(d_qPhi); cudaFree(d_uPhi); cudaFree(d_pPhi);
    free(lambdaDiff2); cudaFree(d_lambdaDiff2);
    cudaFree(d_phiAxis);
    switch(inOptions->fileFormat) {
    case FITS:
       free(fPixel);
       break;
    case HDF5:
       H5Sclose(descriptors.qMemspace);  H5Sclose(descriptors.uMemspace);
       H5Sclose(descriptors.qDataspace); H5Sclose(descriptors.uDataspace);
       H5Dclose(descriptors.qDataset);   H5Dclose(descriptors.uDataset);
       H5Sclose(descriptors.qOutMemspace);  H5Sclose(descriptors.uOutMemspace);
       H5Sclose(descriptors.qOutDataspace); H5Sclose(descriptors.uOutDataspace);
       H5Dclose(descriptors.qOutDataset);   H5Dclose(descriptors.uOutDataset);
       H5Sclose(descriptors.pOutMemspace); H5Sclose(descriptors.pOutDataspace);
       H5Dclose(descriptors.pOutDataset);
       H5Fclose(descriptors->qDirtyH5);
       H5Fclose(descriptors->uDirtyH5);
       H5Fclose(descriptors->pDirtyH5);
       break;
    }

    return(SUCCESS);
}
