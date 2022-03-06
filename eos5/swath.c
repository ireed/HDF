/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at      *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * created by:  Isayah Reed <ireed@hdfgroup.org>
 *		May  4, 2011
 * 
 *  This is a simple program that demonstrates how to use HDF-EOS5 to add CF 
 *  attributes to an h5 file. Creates a swath called Swath1 and adds the 
 *  temperature dataset with 2 dimensions: latitude and longitude. The
 *  temperature dataset has the CF attributes: units, long_name, _FillValue,
 *  coordinates, valid_min, valid_max, valid_range, scale_factor, add_offset. 
 *  Latitude and longitude have the CF attributes units and long_name. Outputs 
 *  data to general.h5
 *
 * Compile:
 *  $(HDF5_DIR)/bin/h5cc swath.c -I$(EOS_DIR)/include -L$(EOS_DIR)/lib -lhe5_hdfeos -lGctp
 *
*/


#include "HE5_HdfEosDef.h"

#define NX1     180
#define NY1     360

#define TEMP "temp"
#define LAT "lat"
#define LON "lon"

#define UNITS "Units"
#define FILLVALUE "_FillValue"
#define LONGNAME "long_name"
#define COORDINATES "coordinates"
#define STDNAME "standard_name"
#define ADDOFFSET "add_offset"
#define SCALEFACTOR "scale_factor"
#define VALIDMIN "valid_min"
#define VALIDMAX "valid_max"
#define VALIDRANGE "valid_range"

int main(void)
{
    herr_t  status= FAIL;
    hid_t   swid, file, dataset, att;  /* file, dataset, attribute handles */
    hid_t   dtype, dspace;
    hsize_t size[2]= {NX1, NY1}, temp_size;
    char  *degrees_east= "degrees_east", *degrees_north= "degrees_north",
            *kelvin= "kelvin", *latitude= "latitude",
            *longitude= "longitude", *temperature= "temperature",
             *coorlist[2]= {"lat", "lon"};
    float    temp_array[NX1][NY1];
    float  lat_array[NX1][NY1];
    float   lon_array[NX1][NY1];
    float value[1];
    int     i, j;

/* PART 1: temp */

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

    /* a few random values should be -999.0; ~5% */
/*
    for(i=0; i<NX1*.05; i++)
        for(j=0; j<NY1*.05; j++)
            temp_array[i][j]= rand() %
*/
    for(i=0; i<NX1; i++)
        lat_array[i][0]= -90.0;
    for(i=0; i<NX1; i++)
        for(j=1; j<NY1; j++)
            lat_array[i][j]= lat_array[i][j-1]+1.0;

    for(j=0; j<NY1; j++)
        lon_array[0][j]= -180.0;
    for(j=1; j<NY1; j++)
        for(i=1; i<NX1; i++)
            lon_array[i][j]= lon_array[i-1][j]+1.0;

    /* Open a new HDF-EOS swath file, "Swath.he5" */ 
    file = HE5_SWopen("swath.he5", H5F_ACC_TRUNC); 

    /* Create the swath, "Swath1", within the file */ 
    swid = HE5_SWcreate(file, "Swath1"); 

    /* Define dimensions and specify their sizes */ 
    status = HE5_SWdefdim(swid, "GeoXtrack", NX1); 
    status = HE5_SWdefdim(swid, "GeoTrack", NY1); 

    /* fill value MUST be called before defining the data field */
    value[0]= -999.0;
    HE5_SWsetfillvalue(swid, TEMP, H5T_NATIVE_FLOAT, &value[0]);

    /* Define data fields (datasets) */
    HE5_SWdefdatafield(swid, TEMP, "GeoTrack,GeoXtrack", NULL, 
            H5T_NATIVE_FLOAT, 0);
 
    status = HE5_SWwritefield(swid, TEMP, NULL, NULL, NULL, temp_array);

/*************   add cf attributes **********/

/* part 1: temperature */
    /* add units */
    size[0]= strlen(kelvin);
    status= HE5_SWwritelocattr(swid, TEMP, UNITS, H5T_C_S1, 
            &size[0], (void*)kelvin);

    /* add long_name */
    size[0]= strlen(temperature);
    status= HE5_SWwritelocattr(swid, TEMP, LONGNAME, H5T_C_S1, 
            &size[0], (void*)temperature);

    /* add coordinates attribute */
    size[0]= 2;
    dtype= H5Tcopy(H5T_C_S1);
    status= H5Tset_size(dtype, H5T_VARIABLE);
    status= HE5_SWwritelocattr(swid, TEMP, COORDINATES, dtype,
            &size[0], coorlist);

    /* add valid_range */
    size[0]= 2;
    value[0]= 0;
    value[1]= 400;
    status= HE5_SWwritelocattr(swid, TEMP, VALIDRANGE, H5T_NATIVE_FLOAT, 
            &size[0], value);
    
    /* add valid_min */
    size[0]= 1;
    status= HE5_SWwritelocattr(swid, TEMP, VALIDMIN, H5T_NATIVE_FLOAT, 
            &size[0], &value[0]);

    /* add valid_min */
    status= HE5_SWwritelocattr(swid, TEMP, VALIDMAX, H5T_NATIVE_FLOAT, 
            &size[0], &value[1]);

    /* add add_offset */
    value[0]= 10;
    status= HE5_SWwritelocattr(swid, TEMP, ADDOFFSET, H5T_NATIVE_FLOAT, 
            &size[0], &value[0]);

    /* add scale_factor */
    value[0]= 1.01;
    status= HE5_SWwritelocattr(swid, TEMP, SCALEFACTOR, H5T_NATIVE_FLOAT, 
            &size[0], &value[0]);




    /* define dimensions */
    status = HE5_SWdefgeofield(swid, "lon", "GeoTrack,GeoXtrack", 
            NULL, H5T_NATIVE_FLOAT, 0);
    status = HE5_SWdefgeofield(swid, "lat", "GeoTrack,GeoXtrack", 
            NULL, H5T_NATIVE_FLOAT, 0);

/* part 2: latitude */

    /* write latitude dataset */
    status= HE5_SWwritefield(swid, LAT, NULL, NULL, 
            NULL, lat_array);

    /* add units */
    size[0]= strlen(degrees_north);
    status= HE5_SWwritelocattr(swid, LAT, UNITS, H5T_C_S1,
            &size[0], (void*)degrees_north);

    /* add long_name */
    size[0]= strlen(latitude);
    status= HE5_SWwritelocattr(swid, LAT, LONGNAME, H5T_C_S1,
            &size[0], (void*)latitude);

    /* add standard_name */
    status= HE5_SWwritelocattr(swid, LAT, STDNAME, H5T_C_S1,
            &size[0], (void*)latitude);

/* part 3: longitude */

    status= HE5_SWwritefield(swid, LON, NULL, NULL, 
            NULL, lon_array);

    /* add units */
    size[0]= strlen(degrees_east);
    status= HE5_SWwritelocattr(swid, LON, UNITS, H5T_C_S1,
            &size[0], (void*)degrees_east);

    /* add long_name */
    size[0]= strlen(longitude);
    status= HE5_SWwritelocattr(swid, LON, LONGNAME, H5T_C_S1,
            &size[0], (void*)longitude);

    /* add standard_name */
    status= HE5_SWwritelocattr(swid, LON, STDNAME, H5T_C_S1,
            &size[0], (void*)longitude);


    /* Close the swath interface */ 
    status = HE5_SWdetach(swid); 
    /* Close the swath file */ 
    status = HE5_SWclose(file); 

return 0;
}


