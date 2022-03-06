/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
 * Copyright by The HDF Group.                                                   !
 * All rights reserved.                                                          !
 *                                                                               !
 * The full HDF5 copyright notice, including terms governing use, modification,  !
 * and redistribution, is contained in the file COPYING.  COPYING can be found   !
 * at the root of the source code distribution tree.                             !
 * For questions contact the help desk at help@hdfgroup.org                      !
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * created by:  Isayah Reed <ireed@hdfgroup.org>
 *              Wednesday, May  4, 2011
 *
 * Purpose:
 *  This is 1 of 4 simple HDF5 programs that demonstrates how to add CF attributes
 *  to an h5 file. This example shows how to use chunking and compression. A file 
 *  is created with 3 datasets: lat, lon, temp. Lat contains the CF attributes: 
 *  units and long_name. Lon has the same CF attributes as the latitude dataset. 
 *  Temp contains the CF attributes: units, long_name, _FillValue, coordinates. 
 *  It is has a chunk size is 900x1800. The deflate compression is used with a 
 *  compression level of 1. Outputs data to chunk_compress.h5
 *
 * Compile:
 *  gcc hdf5_chunk_compress.c -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lsz -lz -lm
 *
*/



#include "hdf5.h"
#include "string.h"
#include <stdlib.h>

#define H5FILE_NAME "chunk_compress.h5"
#define NX1     1800
#define NY1     3600
#define STRINGLISTSIZE  2
#define TEMP "temp"
#define LAT "lat"
#define LON "lon"

/* attributes */
#define UNITS "units"
#define FILLVALUE "_FillValue"
#define LONGNAME "long_name"
#define COORDINATES "coordinates"

int main (void)
{
    hid_t       file, dataset, att;         /* file, dataset, attribute handles */
    hid_t       dataprop;    /* data properties */
    hid_t       floatType, stringType, arrayType; /* datatypes */
    hid_t       floatSpace, stringSpace, arraySpace; /* dataspaces */
    hsize_t     dimsa[2], dimsa3[3], dimsf[1], dimsc[2];  /* dataset dimensions */
    herr_t      status;
    hvl_t       hvl[STRINGLISTSIZE]; /* varialbe length string list */
    int         i, j, k;
    char  *degrees_east= "degrees_east", *degrees_north= "degrees_north", 
            *kelvin= "kelvin", *latitude= "latitude", 
            *longitude= "longitude", *temperature= "temperature"; 
    char *longname, *units;
    char *coorlist[STRINGLISTSIZE]= {"lat", "lon"};
    float **temp_array, lat_array[NX1], lon_array[NY1], fillvalue;



   file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

/* PART 1: temp */
  
    temp_array= malloc(NX1*sizeof(float*));
    for(i=0; i<NX1; i++)
        temp_array[i]= malloc(NY1*sizeof(float));

   /* values between a[60][*] and a[120][*] is around 300.0 */ 
    for(i=60; i<=120; i++)
        for(j=0; j<NY1; j++)
            temp_array[i][j]= 300.0;
   
   /* values between a[0][*] and a[59][*], a[121][*] and a[179][*] 
        is around 280.0 */
    for(i=0; i<60; i++)
        for(j=0; j<NY1; j++)
            temp_array[i][j]= 280.0;
            
    for(i=121; i<NX1; i++)
        for(j=0; j<NY1; j++)
            temp_array[i][j]= 280.0;
            
    /* a few values should be -999.0 */ 
   


    /* initialize handles  */
    floatType= H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_precision(floatType, 32);
    dimsa[0]= NX1;
    dimsa[1]= NY1;
    floatSpace= H5Screate_simple(2, dimsa, NULL);
    stringType= H5Tcopy(H5T_C_S1);
    stringSpace= H5Screate(H5S_SCALAR);

    /* set chunking properties; chunk size = 900x1800 */
    dimsc[0]= 900;
    dimsc[1]= 1800;
    dataprop= H5Pcreate(H5P_DATASET_CREATE);
    status= H5Pset_chunk(dataprop,2,dimsc);
    /* use gzip (deflate) compression with a level 1 compression filter */
    H5Pset_deflate(dataprop,1);
   
    /* write dataset with chunking properties enabled */
    dataset= H5Dcreate(file, TEMP, floatType, floatSpace,
        H5P_DEFAULT, dataprop, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, temp_array);
            
   
   
    /* add 1st attribute */
    units= kelvin;
    status= H5Tset_size(stringType, (hsize_t)strlen(units));
    att= H5Acreate(dataset, UNITS, stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,units);
    H5Aclose(att);
   
    /* add 2nd attribute */
    dimsf[0]= 1;
    fillvalue= -999.0;
    floatSpace= H5Screate_simple(1, dimsf, NULL);
    att= H5Acreate(dataset, FILLVALUE, floatType, floatSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,floatType,&fillvalue);
    H5Aclose(att);
   
   /* add 3rd attribute */
    longname= temperature;
    status= H5Tset_size(stringType, (hsize_t)strlen(longname));
    att= H5Acreate(dataset, LONGNAME, stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,longname);
    H5Aclose(att);


   
   /* add 4th attribute */
    dimsa[0]= STRINGLISTSIZE;
    arraySpace= H5Screate_simple(1, &dimsa[0], NULL);
    arrayType= H5Tcopy(H5T_C_S1);
    status= H5Tset_size(arrayType, H5T_VARIABLE);
    att= H5Acreate(dataset, COORDINATES, arrayType, arraySpace,
        H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,arrayType,coorlist);
    H5Aclose(att);
   

    H5Dclose(dataset);
   




/* PART 2: lat */

    /*initialize: arithmetic sequence from -900.0 - 900.0 */
    lat_array[0]= -90.0;
    for(i = 1; i < NX1; i++)
        lat_array[i] = lat_array[i-1]+0.1;

    dimsf[0] = NX1;
    floatSpace = H5Screate_simple(1, dimsf, NULL);


    /* set chunking properties; chunk size = 900 */
    dimsc[0]= 900;
    dataprop= H5Pcreate(H5P_DATASET_CREATE);
    status= H5Pset_chunk(dataprop,1,dimsc);
    /* use gzip (deflate) compression with a level 1 compression filter */
    H5Pset_deflate(dataprop,1);

    /* write dataset with chunking properties enabled */
    dataset= H5Dcreate(file, LAT, floatType, floatSpace,
        H5P_DEFAULT, dataprop, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, lat_array);


    /* add 1st attribute */
    longname= latitude;
    status= H5Tset_size(stringType, (hsize_t)strlen(longname));
    att= H5Acreate(dataset, LONGNAME, stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,longname);
    H5Aclose(att);

    /* add 2nd attribute */
    units= degrees_north;
    status= H5Tset_size(stringType, (hsize_t)strlen(units));
    att= H5Acreate(dataset, UNITS, stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,units);
    H5Aclose(att);

    H5Dclose(dataset);
   
   
/* Part 3: lon  */

    /*initialize: arithmetic sequence from -1800.0 - 1800.0 */
    lon_array[0]= -180.0;
    for(i = 1; i < NY1; i++)
        lon_array[i] = lon_array[i-1]+0.1;

    dimsf[0] = NY1;
    floatSpace = H5Screate_simple(1, dimsf, NULL);

    /* set chunking properties; chunk size = 1800 */
    dimsc[0]= 1800;
    dataprop= H5Pcreate(H5P_DATASET_CREATE);
    status= H5Pset_chunk(dataprop,1,dimsc);
    /* use gzip (deflate) compression with a level 1 compression filter */
    H5Pset_deflate(dataprop,1);

    /* write dataset with chunking properties enabled */
    dataset= H5Dcreate(file, LON, floatType, floatSpace,
        H5P_DEFAULT, dataprop, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, lon_array);

    /* add 1st attribute */
    longname= longitude;
    status= H5Tset_size(stringType, (hsize_t)strlen(longname));
    att= H5Acreate(dataset, LONGNAME, stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,longname);
    H5Aclose(att);

    /* add 2nd attribute */
    units= degrees_east;
    status= H5Tset_size(stringType, (hsize_t)strlen(units));
    att= H5Acreate(dataset, UNITS, stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,units);
    H5Aclose(att);

    H5Dclose(dataset);

    
    /* Close/release resources. */
    for(i=0; i<NX1; i++)
        free(temp_array[i]);
    free(temp_array);
    H5Pclose(dataprop);
    H5Sclose(floatSpace);
    H5Sclose(stringSpace);
    H5Sclose(arraySpace);
    H5Tclose(floatType);
    H5Tclose(stringType);
    H5Tclose(arrayType);
    H5Fclose(file);

    return 0;
}

