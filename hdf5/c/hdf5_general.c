/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
 * Copyright by The HDF Group.                                                   !
 * All rights reserved.                                                          !
 *                                                                               !
 * This code is provided as open source but it NOT licensed, has no limitations, !
 * and has no expressed or implied warranties.
 * For questions contact the help desk at help@hdfgroup.org                      !
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * created by:  Isayah Reed
 *              May  4, 2011
 *
 * Purpose: 
 *  This is 1 of 4 simple HDF5 programs that demonstrates how to add CF attributes 
 *  to an h5 file. Creates a file with 3 datasets: lat, lon, temp. Lat contains 
 *  the CF attributes: units, long_name, and standard_name. Lon has the same CF 
 *  attributes as the latitude dataset. Temp contains the CF attributes: units, 
 *  long_name, _FillValue, coordinates, valid_min, valid_max, valid_range, 
 *  scale_factor, add_offset. Outputs data to general.h5 
 *
 * Compile:
 *  gcc hdf5_general.c -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lz -lm
 *
*/


#include "hdf5.h"
#include "string.h"
#include <time.h>

#define H5FILE_NAME "general.h5"
/* dataset names */
#define TEMP "temp"
#define LAT "lat"
#define LON "lon"
/* CF attributes */
#define UNITS "units"
#define FILLVALUE "_FillValue"
#define LONGNAME "long_name"
#define COORDINATES "coordinates"
#define VALIDMIN "valid_min"
#define VALIDMAX "valid_max"
#define VALIDRANGE "valid_range"
#define SCALEFACTOR "scale_factor"
#define ADDOFFSET "add_offset"
#define STDNAME "standard_name"



int main(void)
{
    hid_t file, dataset, attr;	/* file, dataset, attribute handles */
    hid_t floatType, stringType, arrayType;	/* datatypes */
    hid_t floatSpace, stringSpace, arraySpace;	/* dataspaces */
    hsize_t dimsa[2], dimsf;	/* dataset dimensions */
    herr_t status;
    int i, j;
    float temp_array[180][360], lat_array[180], lon_array[360], value[2];
    char *coorlist[2] = { "lat", "lon" };

    int xrand[30], yrand[30]; /* array to store the random fillvalues */
    srand(time(NULL));


    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

/* PART 1: TEMPERATURE */

    /* values between a[60][*] and a[120][*] is around 300.0 */
    for(i = 60; i <= 120; i++)
	    for(j = 0; j < 360; j++)
	        temp_array[i][j] = 300.0;

    /* values between a[0][*] and a[59][*], a[121][*] and a[179][*] 
       is around 280.0 */
    for(i = 0; i < 60; i++)
	    for(j = 0; j < 360; j++)
	        temp_array[i][j] = 280.0;

    for(i = 121; i < 180; i++)
	    for(j = 0; j < 360; j++)
	        temp_array[i][j] = 280.0;

/* UNCOMMENT TO ADD RANDOM FILLVALUES TO THE TEMPERATURE DATASET */
    /* a few random values should be -999.0  */
/*
    for(i=0; i<30; i++)
        xrand[i]= rand()%180;
    for(i=0; i<30; i++)
        yrand[i]= rand()%360;
*/
    /* insert the fillvalue (-999) in random parts of the temp array */
/*
    for(i=0; i<30; i++)
        temp_array[xrand[i]][yrand[i]]= -999.0;
*/


    /* initialize handles */
    floatType = H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_precision(floatType, 32); /* set as 32bit */
    dimsa[0] = 180;
    dimsa[1] = 360;
    floatSpace = H5Screate_simple(2, dimsa, NULL);
    stringType = H5Tcopy(H5T_C_S1);
    stringSpace = H5Screate(H5S_SCALAR);

    /* create the temperature dataset  */
    dataset = H5Dcreate(file, TEMP, floatType, floatSpace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, temp_array);


    /* add UNITS attribute */
    /* UNITS: a string that represents the quantity of measurement */
    status = H5Tset_size(stringType, (hsize_t) strlen("kelvin"));
    attr= H5Acreate(dataset, UNITS, stringType, stringSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, stringType, "kelvin");
    H5Aclose(attr);

    /* add FILLVALUE attribute */
    /* FILLVALUE: represents a missing or undefined value */
    dimsf = 1;
    value[0] = -999.0;
    floatSpace = H5Screate_simple(1, &dimsf, NULL);
    attr= H5Acreate(dataset, FILLVALUE, floatType, floatSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, floatType, &value[0]);
    H5Aclose(attr);

    /* add LONGNAME attribute */
    /* LONGNAME: a string representing a long descriptive name for the data */
    status = H5Tset_size(stringType, (hsize_t) strlen("temperature"));
    attr= H5Acreate(dataset, LONGNAME, stringType, stringSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, stringType, "temperature");
    H5Aclose(attr);

    /* add COORDINATES attribute */
    /* COORDINATES: a string list of the associated coordinate variable 
          names of the variable */
    dimsa[0] = 2;
    arraySpace = H5Screate_simple(1, &dimsa[0], NULL);
    arrayType = H5Tcopy(H5T_C_S1);
    /* set the size of each string entry to be unlimited length */
    status = H5Tset_size(arrayType, H5T_VARIABLE);
    attr= H5Acreate(dataset, COORDINATES, arrayType, arraySpace,
		    H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, arrayType, coorlist);
    H5Aclose(attr);

    /* add VALID_MIN attribute */
    /* VALID_MIN: smallest valid value of a variable */
    value[0]= 0.0;
    floatSpace= H5Screate(H5S_SCALAR);
    attr= H5Acreate(dataset, VALIDMIN, floatType, floatSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,floatType,&value[0]);
    H5Sclose(floatSpace);
    H5Aclose(attr);

    /* add VALID_MAX attribute */
    /* VALID_MAX: Smallest valid value of a variable */
    value[0]= 400.0;
    floatSpace= H5Screate(H5S_SCALAR);
    attr= H5Acreate(dataset, VALIDMAX, floatType, floatSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,floatType,&value[0]);
    H5Sclose(floatSpace);
    H5Aclose(attr);


    /*NOTE: a dataset should not contain both VALID_RANGE and VALID_MIN/MAX */


    /* add VALID_RANGE attribute */
    /* VALID_RANGE: smallest and largest valid values of a variable */
    dimsf= 2;
    arraySpace= H5Screate_simple(1, &dimsf, NULL);
    value[0]= 275;
    value[1]= 305;
    attr= H5Acreate(dataset, VALIDRANGE, floatType, arraySpace,
        H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,floatType,value);
    H5Sclose(arraySpace);
    H5Aclose(attr);

    /* add ADD_OFFSET attribute */
    /* ADD_OFFSET: number added to the data after it is read by an application */
    value[0]= 1000.0;
    floatSpace= H5Screate(H5S_SCALAR);
    attr= H5Acreate(dataset, ADDOFFSET, floatType, floatSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,floatType,&value[0]);
    H5Sclose(floatSpace);
    H5Aclose(attr);

    /* add SCALE_FACTOR attribute */
    /* SCALE_FACTOR: factor to multiply after data is read */
    value[0]= 5.0;
    floatSpace= H5Screate(H5S_SCALAR);
    attr= H5Acreate(dataset, SCALEFACTOR, floatType, floatSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,floatType,&value[0]);
    H5Sclose(floatSpace);
    H5Aclose(attr);


    H5Dclose(dataset);



/**** PART 2: LATITUDE ****/

    /*initialize data array: arithmetic sequence from -90.0 - 89.0 */
    lat_array[0] = -90.0;
    for (i = 1; i < 180; i++)
	lat_array[i] = lat_array[i - 1] + 1.0;

    dimsf = 180;
    floatSpace = H5Screate_simple(1, &dimsf, NULL);


    dataset = H5Dcreate(file, LAT, floatType, floatSpace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, lat_array);

    /* add LONG_NAME attribute */
    /*stringSpace= H5Screate(H5S_SCALAR); */
    status = H5Tset_size(stringType, (hsize_t) strlen("latitude"));
    attr= H5Acreate(dataset, LONGNAME, stringType, stringSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, stringType, "latitude");
    H5Aclose(attr);

    /* add UNITS attribute */
    status = H5Tset_size(stringType, (hsize_t) strlen("degrees_north"));
    attr= H5Acreate(dataset, UNITS, stringType, stringSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, stringType, "degrees_north");
    H5Aclose(attr);

    /* add STANDARD_NAME attribute */
    /* STANDARD_NAME: identifies the variable is the coordinate variable */
    status= H5Tset_size(stringType, (hsize_t)strlen("latitude"));
    attr= H5Acreate(dataset, STDNAME, stringType, stringSpace,
            H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,stringType,"latitude");
    H5Aclose(attr);

    H5Dclose(dataset);




/**** Part 3: LONGITUDE ****/

    /*initialize: arithmetic sequence from -180.0 - 179.0 */
    lon_array[0] = -180.0;
    for (i = 1; i < 360; i++)
	lon_array[i] = lon_array[i - 1] + 1.0;

    dimsf = 360;
    floatSpace = H5Screate_simple(1, &dimsf, NULL);

    dataset = H5Dcreate(file, LON, floatType, floatSpace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, lon_array);

    /* add LONGNAME attribute */
    status = H5Tset_size(stringType, (hsize_t) strlen("longitude"));
    attr= H5Acreate(dataset, LONGNAME, stringType, stringSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, stringType, "longitude");
    H5Aclose(attr);

    /* add UNITS attribute */
    status = H5Tset_size(stringType, (hsize_t) strlen("degrees_east"));
    attr= H5Acreate(dataset, UNITS, stringType, stringSpace, 
            H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, stringType, "degrees_east");
    H5Aclose(attr);

    /* add STANDARD_NAME attribute */
    /* STANDARD_NAME: identifies the variable is the coordinate variable */
    status= H5Tset_size(stringType, (hsize_t)strlen("longitude"));
    attr= H5Acreate(dataset, STDNAME, stringType, stringSpace,
            H5P_DEFAULT, H5P_DEFAULT);
    status= H5Awrite(attr,stringType,"longitude");
    H5Aclose(attr);
    H5Dclose(dataset);

    /* Close/release resources. */
    H5Sclose(floatSpace);
    H5Sclose(stringSpace);
    H5Sclose(arraySpace);
    H5Tclose(floatType);
    H5Tclose(stringType);
    H5Tclose(arrayType);
    H5Fclose(file);

    return 0;
}
