
/*
 * Programmer:  Isayah Reed <help@hdfgroup.org> 
 *              May  4, 2011
 *
 * Purpose:
 *  This is a simple netCDF4 program that demonstrates how to add CF attributes
 *  to an h5 file. Creates a file with 3 variables: lat, lon, temp. Lat contains
 *  the CF attributes: units, long_name, and standard_name. Lon has the same CF
 *  attributes as the latitude variable. Temp contains the CF attributes: units,
 *  long_name, _FillValue, and coordinates. Outputs data to netcdf_ex1.h5
 *
 * Compile:
 *  gcc netcdf_ex1.c -I$NETCDF4_DIR/include -I$HDF5_DIR/include -L$NETCDF4_DIR/lib -L$HDF5_DIR/lib $NETCDF4_DIR/lib/libnetcdf.a $HDF5_DIR/lib/libhdf5_hl.a $HDF5_DIR/lib/libhdf5.a -lsz -lz -lm
 *
*/


    #include <stdlib.h>
    #include <stdio.h>
    #include "netcdf.h"
    #include "hdf5.h"
    #include "hdf5.h"
    #include "string.h"
     
     /* This is the name of the data file we will create. */
    #define FILE_NAME "netcdf_ex1.nc"
     
     /* We are writing 2D data, a 180 x 360 grid. */
    #define NDIMS 2
    #define NX1     180      /* dataset dimensions */
    #define NY1     360
    #define STRINGLISTSIZE  2
    /* datasets */
    #define TEMP "temp"
    #define LAT "lat"
    #define LON "lon"
    /* attributes */
    #define UNITS "units"
    #define FILLVALUE "_FillValue"
    #define LONGNAME "long_name"
    #define COORDINATES "coordinates"
     
     /* Handle errors by printing an error message and exiting with a
      * non-zero status. */
     #define ERRCODE 2
     #define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

     int main()
     {
        int ncid, varid[3];
        int dimsf[1], dimsa[2], dimst; /* float, array, text dimensions */
        float temp_array[NX1][NY1], lat_array[NX1], 
                lon_array[NY1], fillvalue= -999.0;
        int retval, i, j;
        char  *degrees_east= "degrees_east", *degrees_north= "degrees_north",
            *kelvin= "kelvin", *latitude= "latitude",
            *longitude= "longitude", *temperature= "temperature";
        char *coorlist[STRINGLISTSIZE]= {"lat\0", "lon\0"};        
     
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

        lat_array[0]= -90.0;
        for(i=1; i<NX1; i++)
            lat_array[i]= lat_array[i-1]+1.0;

        lon_array[0]= -180;
        for(i=1; i<NY1; i++)
            lon_array[i]= lon_array[i-1]+1.0;
     
        /* Always check the return code of every netCDF function call. In
         * this example program, any retval which is not equal to NC_NOERR
         * (0) will cause the program to print an error message and exit
         * with a non-zero return code. */
     
        /* Create the file. The NC_CLOBBER parameter tells netCDF to
         * overwrite this file, if it already exists.*/
        if ((retval = nc_create(FILE_NAME, NC_NETCDF4|NC_CLOBBER, &ncid)))
           ERR(retval);
     


        /* Define the dimensions. NetCDF will hand back an ID for each. */
        if ((retval = nc_def_dim(ncid, LAT, NX1, &dimsa[0])))
           ERR(retval);
        if ((retval = nc_def_dim(ncid, LON, NY1, &dimsa[1])))
           ERR(retval);
     
     
        /* Define the variable. The type of the variable in this case is
         * NC_INT (4-byte integer). */
        if ((retval = nc_def_var(ncid, TEMP, NC_FLOAT, 2,
                    dimsa, &varid[0])))
           ERR(retval);
        if ((retval = nc_def_var(ncid, LAT, NC_FLOAT, 1,
                    &dimsa[0], &varid[1])))
           ERR(retval);
        if ((retval = nc_def_var(ncid, LON, NC_FLOAT, 1,
                    &dimsa[1], &varid[2])))
           ERR(retval);

     
        /* define fillvalue, set no_fill=0 to write fillvalue */
        if ((retval= nc_def_var_fill(ncid, varid[0], 0, &fillvalue)))
            ERR(retval);


        /* End define mode. This tells netCDF we are done defining
         * metadata. */
        if ((retval = nc_enddef(ncid)))
           ERR(retval);


        /* Write the coordinate variable data. This will put the latitudes
           and longitudes into the netCDF file. */
        if ((retval = nc_put_var_float(ncid, varid[0], 
                &temp_array[0][0])))
           ERR(retval);

        if ((retval = nc_put_var_float(ncid, varid[1], &lat_array[0])))
           ERR(retval);

        if ((retval = nc_put_var_float(ncid, varid[2], &lon_array[0])))
           ERR(retval);

        if ((retval = nc_put_att_text(ncid, varid[0], UNITS, strlen(kelvin),
                kelvin)))
            ERR(retval);

        if ((retval = nc_put_att_text(ncid, varid[1], UNITS, 
                strlen(degrees_north), degrees_north)))
            ERR(retval);

        if ((retval = nc_put_att_text(ncid, varid[2], UNITS, 
                strlen(degrees_east), degrees_east)))
            ERR(retval);
     
        if ((retval = nc_put_att_text(ncid, varid[0], LONGNAME, 
                strlen(temperature), temperature)))
            ERR(retval);

        if ((retval = nc_put_att_text(ncid, varid[1], LONGNAME, 
                strlen(latitude), latitude)))
            ERR(retval);

        if ((retval = nc_put_att_text(ncid, varid[2], LONGNAME, 
                strlen(longitude), longitude)))
            ERR(retval);

        /* use nc_put_att_string for string lists */
        if ((retval = nc_put_att_string(ncid, varid[0], COORDINATES, 
                STRINGLISTSIZE, (const char**)&coorlist)))
            ERR(retval);

        /* Close the file. This frees up any internal netCDF resources
         * associated with the file, and flushes any buffers. */
        if ((retval = nc_close(ncid)))
           ERR(retval);
     
        return 0;
     }

