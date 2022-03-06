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
 *  to an h5 file. This example shows how to use 3-D datasets. A file is created 
 *  with 4 datasets: radiation, latitude, longitude, pressure. Radiation is a 3-D
 *  dataset that contains the CF attributes: units, fillvalue, long_name, and 
 *  coordinates. Latitude, longitude, and pressure contain the CF attributes units 
 *  and long_name. Data is written to 3D.h5
 *
 * Compile:
 *  gcc hdf5_3D.c -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lz -lm
 *
*/


#include "hdf5.h"
#include "string.h"
#include <stdlib.h>

#define H5FILE_NAME "3D.h5"
#define RADIATION "Radiation"
#define LATITUDE "latitude"
#define LONGITUDE "longitude"
#define PRESSURE "pressure"


int main (void)
{
    hid_t       file, dataset, att;         /* file, dataset, attribute handles */
    hid_t       floatType, stringType, arrayType; /* datatypes */
    hid_t       floatSpace, stringSpace, arraySpace; /* dataspaces */
    hsize_t     dimsa[2], dimsa3[3], dimsf[1], dimsc[2];  /* dataset dimensions */
    herr_t      status;
    int         i, j, k;
    char *longname, *units;
    char *coorlist[3]= {"latitude","longitude","pressure"};
    float *rad_array, lat_array[180],
            lon_array[360], press_array[20], fillvalue;


   file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

/* PART 1: radiation */
  
    /* allocate the radiation array */
    /* rad_array[180][360][120] */
    rad_array = (float *) malloc(180*360*120*sizeof(float));
/*    rad_array= (float***)malloc(180*sizeof(float**));
    for(i=0; i<180; i++)
    {
        rad_array[i]= malloc(360*sizeof(float*));
        for(j=0; j<360; j++)
            rad_array[i][j]= malloc(120*sizeof(float));
    }

*/
    /* initialize array */
    for(i=0; i<180*360*20; i++)
           	*(rad_array+i)= 280.0;


   /* initialize handles */
    floatType= H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_precision(floatType, 32);
    stringType= H5Tcopy(H5T_C_S1);
    stringSpace= H5Screate(H5S_SCALAR);

    dimsa3[0]= 180;
    dimsa3[1]= 360;
    dimsa3[2]= 20;
    floatSpace= H5Screate_simple(3, dimsa3, NULL);
    dataset= H5Dcreate(file, RADIATION, floatType, floatSpace,
    	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
    	H5P_DEFAULT, rad_array);
			
   
   
   /* add 1st attribute */
    status= H5Tset_size(stringType, (hsize_t)strlen("Watts/(m^2)"));
    att= H5Acreate(dataset, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"Watts/(m^2)");
    H5Aclose(att);
   
    /* add 2nd attribute */
    dimsf[0]= 1;
    fillvalue= -999.0;
    floatSpace= H5Screate(H5S_SCALAR);
    att= H5Acreate(dataset, "_FillValue", floatType, floatSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,floatType,&fillvalue);
    H5Aclose(att);
   
    /* add 3rd attribute */
    status= H5Tset_size(stringType, (hsize_t)strlen("outgoing long-wave radiation"));
    att= H5Acreate(dataset, "long_name", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"outgoing long-wave radiation");
    H5Aclose(att);


   
   /* add 4th attribute */
    dimsa[0]= 3;
    arraySpace= H5Screate_simple(1, &dimsa[0], NULL);
    arrayType= H5Tcopy(H5T_C_S1);
    status= H5Tset_size(arrayType, H5T_VARIABLE);
    att= H5Acreate(dataset, "coordinates", arrayType, arraySpace, H5P_DEFAULT,
		H5P_DEFAULT);
    status= H5Awrite(att,arrayType,coorlist);
    H5Aclose(att);


    H5Dclose(dataset);
   




/* PART 2: latitude  */

    /*initialize: arithmetic sequence from -90.0 - 89.0 */
    lat_array[0]= -90.0;
    for(i = 1; i < 180; i++)
        lat_array[i] = lat_array[i-1]+1.0;

    dimsf[0] = 180;
    floatSpace = H5Screate_simple(1, dimsf, NULL);


    dataset = H5Dcreate(file, LATITUDE, floatType, floatSpace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, lat_array);

    /* add 1st attribute */
    status= H5Tset_size(stringType, (hsize_t)strlen("latitude"));
    att= H5Acreate(dataset, "long_name", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"latitude");
    H5Aclose(att);

    /* add 2nd attribute */
    status= H5Tset_size(stringType, (hsize_t)strlen("degrees_north"));
    att= H5Acreate(dataset, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"degrees_north");
    H5Aclose(att);

    
   
   
/* Part 3: lon  */

    /*initialize: arithmetic sequence from -180.0 - 179.0 */
   lon_array[0]= -180.0;
   for(i = 1; i < 360; i++)
       lon_array[i] = lon_array[i-1]+1.0;

    dimsf[0] = 360;
    floatSpace = H5Screate_simple(1, dimsf, NULL);

    dataset = H5Dcreate(file, LONGITUDE, floatType, floatSpace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
         H5P_DEFAULT, lon_array);

    /* add 1st attribute */
    status= H5Tset_size(stringType, (hsize_t)strlen("longitude"));
    att= H5Acreate(dataset, "long_name", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"longitude");

    /* add 2nd attribute */
    status= H5Tset_size(stringType, (hsize_t)strlen("degrees_east"));
    att= H5Acreate(dataset, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"degrees_east");


	
/* Part 4: pressure  */


    /* initialize: arithmetic sequence from 1000 to 250 */
    press_array[0]= 1000.0;
    for(i=1; i<20; i++)
	press_array[i]= press_array[i-1]-37.5;
			
    dimsf[0] = 20;
    floatSpace = H5Screate_simple(1, dimsf, NULL);

    dataset = H5Dcreate(file, PRESSURE, floatType, floatSpace,
	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	    H5P_DEFAULT, press_array);
			 
    /* add 1st attribute */
    status= H5Tset_size(stringType, strlen("pressure"));
    att= H5Acreate(dataset, "long_name", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"pressure");

    /* add 2nd attribute */
    status= H5Tset_size(stringType, strlen("hpa"));
    att= H5Acreate(dataset, "units", stringType, stringSpace, H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(att,stringType,"hpa");
		
	
    

    /* Close/release resources. */
/*    for(i=0; i<180; i++)
        for(j=0; j<360; j++)
            free(rad_array[i][j]);
    for(i=0; i<180; i++)
            free(rad_array[i]);
*/

     free(rad_array);
    H5Sclose(floatSpace);
    H5Sclose(stringSpace);
    H5Sclose(arraySpace);
    H5Tclose(floatType);
    H5Tclose(stringType);
    H5Tclose(arrayType);
    H5Fclose(file);

    return 0;
}

