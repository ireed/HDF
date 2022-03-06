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
 *              May  4, 2011
 *
 * Purpose:
 *  This is 1 of 4 simple HDF5 programs that demonstrates how to add CF attributes
 *  to an h5 file. This example demonstrates how to use the HDF5 dimension scale 
 *  API. A file is created with 3 datasets: lat, lon, temp. Lat contains the CF 
 *  attributes: units and long_name. Lon has the same CF attributes as the latitude 
 *  dataset. Temp contains the CF attributes: units, long_name, _FillValue, and 
 *  coordinates. Outputs data to ds_ex1.h5
 *
 * Compile:
 *  gcc dim_scale.c -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lsz -lz -lm
 *
*/


#include "hdf5.h"
#include "string.h"

#define H5FILE_NAME "ds_ex1.h5"
#define TEMP "temp"
#define LAT "lat"
#define LON "lon"
#define NX1 180
#define NY1 360
#define STRINGLISTSIZE  2

#define UNITS "Units"
#define FILLVALUE "_FillValue"
#define LONGNAME "long_name"
#define COORDINATES "coordinates"

int
main (void)
{
    hid_t       file, dataset[3], att;         /* file, dataset, attribute handles */
    hid_t       floatType, stringType, arrayType; /* datatypes */
    hid_t       floatSpace, stringSpace, arraySpace; /* dataspaces */
    hsize_t     dimsa[2], dimsa3[3], dimsf[1], dimsc[2];  /* dataset dimensions */
    herr_t      status;
    int         i, j, k;
    float temp_array[NX1][NY1], lat_array[NX1], lon_array[NY1], fillvalue;
    char *longname, *units;
    char  *degrees_east= "degrees_east", *degrees_north= "degrees_north",
            *kelvin= "kelvin", *latitude= "latitude",
            *longitude= "longitude", *temperature= "temperature";
    char *coorlist[STRINGLISTSIZE]= {"lat", "lon"};


    /* create h5 file */
   file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);



   /* initialize handles */
   floatType= H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_precision(floatType, 32);
   dimsa[0]= NX1;
   dimsa[1]= NY1;
   floatSpace= H5Screate_simple(2, dimsa, NULL);
   stringType= H5Tcopy(H5T_C_S1);
    stringSpace= H5Screate(H5S_SCALAR);
   


/* PART 1: temp */


    /* initialize temperature array  */
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
            


	dataset[0]= H5Dcreate(file, TEMP, floatType, floatSpace,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status= H5Dwrite(dataset[0], H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
         H5P_DEFAULT, temp_array);
			
   
   
   /* add 1st attributei: set units= kelvin */
	units= kelvin;
    status= H5Tset_size(stringType, (hsize_t)strlen(units));
	att= H5Acreate(dataset[0], UNITS, stringType, stringSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
	status= H5Awrite(att,stringType,units);
    H5Aclose(att);
   
   /* add 2nd attribute: set _fillvalue= -999.0 */
   dimsf[0]= 1;
   fillvalue= -999.0;
   floatSpace= H5Screate_simple(1, dimsf, NULL);
   att= H5Acreate(dataset[0], FILLVALUE, floatType, floatSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
   status= H5Awrite(att,floatType,&fillvalue);
    H5Aclose(att);
   
   /* add 3rd attribute: set longname= temperature */
	longname= temperature;
   status= H5Tset_size(stringType, (hsize_t)strlen(longname));
   att= H5Acreate(dataset[0], LONGNAME, stringType, stringSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
   status= H5Awrite(att,stringType,longname);
    H5Aclose(att);


   
   /* add 4th attribute: set coordinates= lat lon */

   dimsa[0]= STRINGLISTSIZE;
    arraySpace= H5Screate_simple(1, &dimsa[0], NULL);
    arrayType= H5Tcopy(H5T_C_S1);
   status= H5Tset_size(arrayType, H5T_VARIABLE);
   att= H5Acreate(dataset[0], COORDINATES, arrayType, arraySpace, 
        H5P_DEFAULT, H5P_DEFAULT);
   status= H5Awrite(att,arrayType,coorlist);
    H5Aclose(att);
   

   




/* PART 2: lat */

   /*initialize: arithmetic sequence from -90.0 - 89.0 */
   lat_array[0]= -90.0;
   for(i = 1; i < NX1; i++)
       lat_array[i] = lat_array[i-1]+1.0;

    dimsf[0] = NX1;
    floatSpace = H5Screate_simple(1, dimsf, NULL);


    dataset[1] = H5Dcreate(file, LAT, floatType, floatSpace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status= H5Dwrite(dataset[1], H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
         H5P_DEFAULT, lat_array);

    /* add 1st attribute */
	longname= latitude;
    status= H5Tset_size(stringType, (hsize_t)strlen(longname));
    att= H5Acreate(dataset[1], LONGNAME, stringType, stringSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,longname);
    H5Aclose(att);

    /* add 2nd attribute */
	units= degrees_north;
    status= H5Tset_size(stringType, (hsize_t)strlen(units));
    att= H5Acreate(dataset[1], UNITS, stringType, stringSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,units);
    H5Aclose(att);

    
    /* add dimension scale extensions */
    /* add CLASS attribute, set NAME= "lat" */
    H5DSset_scale(dataset[1], LAT);

    /* add REFERENCE_LIST attribute */
    H5DSattach_scale(dataset[0], dataset[1], 0);


   
   
/* Part 3: lon */

    /*initialize: arithmetic sequence from -180.0 - 179.0 */
   lon_array[0]= -180.0;
   for(i = 1; i < NY1; i++)
       lon_array[i] = lon_array[i-1]+1.0;

    dimsf[0] = NY1;
    floatSpace = H5Screate_simple(1, dimsf, NULL);

    dataset[2] = H5Dcreate(file, LON, floatType, floatSpace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Dwrite(dataset[2], H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
         H5P_DEFAULT, lon_array);

    /* add 1st attribute */
	longname= longitude;
    status= H5Tset_size(stringType, (hsize_t)strlen(longname));
    att= H5Acreate(dataset[2], LONGNAME, stringType, stringSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,longname);
    H5Aclose(att);

    /* add 2nd attribute */
	units= degrees_east;
    status= H5Tset_size(stringType, (hsize_t)strlen(units));
    att= H5Acreate(dataset[2], UNITS, stringType, stringSpace, 
        H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,units);
    H5Aclose(att);

    /* add dimension scale extensions */
    /* add CLASS attribute, set NAME= "lon" */
    H5DSset_scale(dataset[2], LON);
    /* add REFERENCE_LIST attribute */
    H5DSattach_scale(dataset[0], dataset[2], 1);


	

    /* Close/release resources. */
    for(i=0; i<3; i++)
        H5Dclose(dataset[i]);
    H5Sclose(floatSpace);
    H5Sclose(stringSpace);
    H5Sclose(arraySpace);
    H5Tclose(floatType);
    H5Tclose(stringType);
    H5Tclose(arrayType);
    H5Fclose(file);

    return 0;
}

