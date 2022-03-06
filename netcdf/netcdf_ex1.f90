 ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 ! Copyright by The HDF Group.                                               *
 ! Copyright by the Board of Trustees of the University of Illinois.         *
 ! All rights reserved.                                                      *
 !                                                                           *
 ! This file is part of HDF5.  The full HDF5 copyright notice, including     *
 ! terms governing use, modification, and redistribution, is contained in    *
 ! the files COPYING and Copyright.html.  COPYING can be found at the root   *
 ! of the source code distribution tree; Copyright.html can be found at      *
 ! http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 ! access to either file, you may request a copy from help@hdfgroup.org.     *
 ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

 ! Programmer:  Isayah Reed <ireed@hdfgroup.org>
 !              May  4, 2011
 !
 ! Purpose:
 !  This is a simple netCDF4 program that demonstrates how to add CF attributes
 !  to an h5 file. Creates a file with 3 variables: lat, lon, temp. Lat contains
 !  the CF attributes: units, long_name, and standard_name. Lon has the same CF
 !  attributes as the latitude variable. Temp contains the CF attributes: units,
 !  long_name, _FillValue, and coordinates. Outputs data to netcdf_ex1.h5
 !  
 ! Compile:
 !  ifort netcdf_ex1.f90 -I$(NETCDF4_DIR)/include -I$(HDF5_DIR)/include -L$(NETCDF4_DIR)/lib -L$(HDF5_DIR)/lib $(NETCDF4_DIR)/lib/libnetcdff.a $(NETCDF4_DIR)/lib/libnetcdf.a  $(HDF5_DIR)/lib/libhdf5_hl.a $(HDF5_DIR)/lib/libhdf5.a -lsz -lz -lm
 !  


PROGRAM main

  USE NETCDF

  IMPLICIT NONE

  CHARACTER(LEN=13), PARAMETER :: filename = "netcdf_ex1.nc"
  CHARACTER(LEN=4) , PARAMETER :: TEMP= "temp"
  CHARACTER(LEN=3) , PARAMETER :: LAT= "lat"
  CHARACTER(LEN=3) , PARAMETER :: LON= "lon"
  CHARACTER(LEN=9) , PARAMETER :: LONGNAME = "long_name"
  CHARACTER(LEN=5) , PARAMETER :: UNITS = "units"
  CHARACTER(LEN=10) , PARAMETER :: FILLVALUE = "_FillValue"
  CHARACTER(LEN=11) , PARAMETER :: COORDINATES = "coordinates"

  INTEGER          , PARAMETER :: nx1     = 180
  INTEGER          , PARAMETER :: ny1     = 360

  INTEGER :: status
  
  INTEGER :: ncid, varid0, varid1, varid2 ! handles
  CHARACTER(LEN=20), DIMENSION(2) ::  attr_data  ! Attribute string data
  ! switch dimensions to adhere to C standards
  INTEGER, DIMENSION(1:2) :: temp_dims = (/ny1, nx1/) ! size read/write buffer
  INTEGER :: lat_dims != (/nx1/) 
  INTEGER :: lon_dims != (/ny1/) 
  INTEGER, DIMENSION(1) :: dims1
  INTEGER, DIMENSION(2) :: dims2
  DOUBLE PRECISION, DIMENSION(1:nx1,1:ny1) :: temp_data
  DOUBLE PRECISION, DIMENSION(nx1) :: lat_data
  DOUBLE PRECISION, DIMENSION(ny1) :: lon_data
  DOUBLE PRECISION, DIMENSION(1) :: fill_data
  INTEGER :: i, j, strlen
  DOUBLE PRECISION :: fill_value= -999.0
  CHARACTER, DIMENSION(2) ::  hvl_data

  !
  ! Initialize FORTRAN interface.
  !
  ! CALL h5open_f(hdferr)
  !
  ! Initialize temperature data.
  !
  DO i = 1, 60
     DO j = 1, ny1 
        temp_data(i,j) = 280.0
     ENDDO
  ENDDO
  DO i = 61, 120
     DO j = 1, ny1 
        temp_data(i,j) = 300.0
     ENDDO
  ENDDO
  DO i = 121, nx1
     DO j = 1, ny1 
        temp_data(i,j) = 280.0
     ENDDO
  ENDDO


  lat_data(1)= -90.0
  DO i = 2, nx1
    lat_data(i) = lat_data(i-1)+1.0
  ENDDO


  lon_data(1)= -180.0
  DO i = 2, ny1
    lon_data(i) = lon_data(i-1)+1.0
  ENDDO

  !
  ! Create a new file using the default properties.
  ! 
  status=  nf90_create(path=filename, cmode=IOR(NF90_CLOBBER,NF90_HDF5), ncid=ncid)

  ! define dimensions
  status= nf90_def_dim(ncid, LAT, NX1, lat_dims)
  status= nf90_def_dim(ncid, LON, NY1, lon_dims)

  ! Define the variable. The type of the variable in this case is
  ! NC_FLOAT
  status= nf90_def_var(ncid, TEMP, NF90_FLOAT, (/lat_dims, lon_dims/), varid0)
  status= nf90_def_var(ncid, LAT, NF90_FLOAT, lat_dims, varid1)
  status= nf90_def_var(ncid, LON, NF90_FLOAT, lon_dims, varid2)

  ! define fillvalue, set no_fill=0 to write fillvalue
  status= nf90_def_var_fill(ncid, varid0, 0, -999)

  ! End define mode. This tells netCDF we are done defining metadata.
  status= nf90_enddef(ncid)


  ! don't need to adjust array for fortran-C compatibility
  status= nf90_put_var(ncid, varid0, temp_data) 
  status= nf90_put_var(ncid, varid1, lat_data) 
  status= nf90_put_var(ncid, varid2, lon_data) 

  status= nf90_put_att(ncid, varid0, UNITS, "kelvin")
  status= nf90_put_att(ncid, varid1, UNITS, "degrees_north")
  status= nf90_put_att(ncid, varid2, UNITS, "degrees_east")
  
  status= nf90_put_att(ncid, varid0, LONGNAME, "temperature")
  status= nf90_put_att(ncid, varid1, LONGNAME, "latitude")
  status= nf90_put_att(ncid, varid2, LONGNAME, "longitude")

  status= nf90_put_att(ncid, varid0, FILLVALUE, -999.0)

  ! COORDINATES
  status= nf90_put_att(ncid, varid0, COORDINATES, LAT)
  status= nf90_put_att(ncid, varid0, COORDINATES, LON)

  status= nf90_close(ncid)
  
END PROGRAM main
