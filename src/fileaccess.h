/******************************************************************************
fileaccess.h
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
#ifndef FILEACCESS_H
#define FILEACCESS_H

#ifdef __cplusplus
extern "C"
#endif

void checkFitsError(int status);
void checkInputFiles(struct optionsList *inOptions, struct IOFileDescriptors *descriptors);

int getFitsHeader(struct optionsList *inOptions, struct parameters *params, struct fits_header_parameters *header, struct IOFileDescriptors *descriptors);
int getHDF5Header(struct optionsList *inOptions, struct fits_header_parameters *header_parameters, struct parameters *params, struct IOFileDescriptors *descriptors);

void makeOutputFitsImages(struct optionsList *inOptions, struct IOFileDescriptors *descriptors, struct fits_header_parameters *header_parameters, struct parameters *params);
void makeOutputHDF5Images(struct optionsList *inOptions, struct IOFileDescriptors *descriptors, struct parameters *params, struct fits_header_parameters *header);

int getFreqList(struct IOFileDescriptors *descriptors, struct parameters *params, struct DataArrays *data_array);

/* Define the output file names here */
#define DIRTY_P "dirtyP.fits"
#define DIRTY_Q "dirtyQ.fits"
#define DIRTY_U "dirtyU.fits"

#endif
