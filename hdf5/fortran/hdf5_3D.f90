 ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
 ! Copyright by The HDF Group.                                                   !
 ! All rights reserved.                                                          !
 !                                                                               !
 ! The full HDF5 copyright notice, including terms governing use, modification,  !
 ! and redistribution, is contained in the file COPYING.  COPYING can be found   !
 ! at the root of the source code distribution tree.                             !
 ! For questions contact the help desk at help@hdfgroup.org                      !
 ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
 !
 ! Programmer:  Isayah Reed
 !              May  4, 2011
 !
 ! Purpose:
 !  This is 1 of 3 simple HDF5 programs that demonstrates how to add CF attributes
 !  to an h5 file. This example shows how to use 3-D datasets. A file is created
 !  with 4 datasets: radiation, latitude, longitude, pressure. Radiation is a 3-D
 !  dataset that contains the CF attributes: units, fillvalue, long_name, and
 !  coordinates. Latitude, longitude, and pressure contain the CF attributes units
 !  and long_name. Data is written to 3D.h5
 !  
 ! Compile:
 !  $(HDF5_DIR)/bin/h5fc hdf5_3D.f90 -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl  -lhdf5hl_fortran -lhdf5_fortran -lz -lm
 !

PROGRAM main

  USE HDF5

  IMPLICIT NONE

  CHARACTER(LEN=6), PARAMETER :: filename = "3D.h5"
  CHARACTER(LEN=9) , PARAMETER :: RADIATION= "Radiation"
  CHARACTER(LEN=8) , PARAMETER :: LATITUDE= "latitude"
  CHARACTER(LEN=9) , PARAMETER :: LONGITUDE= "longitude"
  CHARACTER(LEN=8) , PARAMETER :: PRESSURE= "pressure"
  CHARACTER(LEN=9) , PARAMETER :: LONGNAME = "long_name"
  CHARACTER(LEN=5) , PARAMETER :: UNITS = "units"
  CHARACTER(LEN=10) , PARAMETER :: FILLVALUE = "_FillValue"
  CHARACTER(LEN=11) , PARAMETER :: COORDINATES = "coordinates"

  INTEGER          , PARAMETER :: nx1     = 180
  INTEGER          , PARAMETER :: ny1     = 360
  INTEGER          , PARAMETER :: nz1     = 20

  INTEGER :: hdferr
  INTEGER(HID_T) :: file, space, dset, attr_id, attr_space ! handles
  INTEGER(HID_T) :: atype_id
  CHARACTER(LEN=10), DIMENSION(3) ::  coorlist  ! Attribute string data
  ! switch dimensions to adhere to C standards
  INTEGER(HSIZE_T), DIMENSION(1:3) :: temp_dims = (/nz1,ny1,nx1/) ! size read/write buffer
  INTEGER(HSIZE_T), DIMENSION(1:1) :: lat_dims = (/nx1/) 
  INTEGER(HSIZE_T), DIMENSION(1:1) :: lon_dims = (/ny1/) 
  INTEGER(HSIZE_T), DIMENSION(1:1) :: press_dims = (/nz1/) 
  INTEGER(HSIZE_T), DIMENSION(1) :: dimsf ! dimension for scalars
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsa ! dimension for arrays
  DOUBLE PRECISION, DIMENSION(1:nx1,1:ny1,1:nz1) :: temp_data
  DOUBLE PRECISION, DIMENSION(nx1) :: lat_data
  DOUBLE PRECISION, DIMENSION(ny1) :: lon_data
  DOUBLE PRECISION, DIMENSION(nz1) :: press_data
  DOUBLE PRECISION, DIMENSION(1) :: fill_data
  CHARACTER(LEN=30) :: attr_data
  INTEGER :: i, j, k, strlen


  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  !
  ! Initialize temperature data.
  !
  DO i = 1, 60
     DO j = 1, ny1 
        DO k = 1, nz1
            temp_data(i,j,k) = 200.0
        ENDDO
     ENDDO
  ENDDO
  DO i = 61, 120
     DO j = 1, ny1 
        DO k = 1, nz1
            temp_data(i,j,k) = 250.0
        ENDDO
     ENDDO
  ENDDO
  DO i = 121, nx1
     DO j = 1, ny1 
        DO k = 1, nz1
            temp_data(i,j,k) = 200.0
        ENDDO
     ENDDO
  ENDDO
  !
  ! Create a new file using the default properties.
  !
  
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)

  ! PART 1: temperature

  !
  ! Create dataspace.  Setting size to be the current size.
  !
  CALL h5screate_simple_f(3, temp_dims, space, hdferr)
  !
  ! Create the dataset.  We will use all default properties for this
  ! example.
  !
  CALL h5dcreate_f(file, RADIATION, h5t_ieee_f32le, space, dset, hdferr)
  !
  ! Write the data to the dataset.
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, temp_data, temp_dims, hdferr)
  !
  !
  ! add 1st attribute: units
  !
  dimsf = 1
  attr_data= "Watts/(m^2)"
  strlen= 11
  !
  ! Create scalar data space for the attribute.
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  ! Create datatype for the attribute.
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  !
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  !Create dataset attribute.
  !
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, hdferr)
  !
  ! Write the attribute data.
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  !
  ! Close the dataspace and attribute.
  !
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: fillvalue
  !
  dimsf = 1
  fill_data(1)= -999.0
  CALL h5screate_simple_f(1, dimsf, space, hdferr)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, hdferr)
  CALL h5acreate_f(dset, FILLVALUE, atype_id, space, &
                      attr_id, hdferr)
  dimsf = 1
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, fill_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 3rd attribute: long name
  !
  dimsf = 1
  attr_data= "outgoing long-wave radiation"
  strlen= 28
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, hdferr)
  dimsf = 3
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 4th attribute: coordinates
  !
  dimsf = 3
  coorlist(1)= "latitude"
  coorlist(2)= "longitude"
  coorlist(3)= "pressure"
  strlen= 10
  CALL h5screate_simple_f(1, dimsf, space, hdferr)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  CALL h5acreate_f(dset, COORDINATES, atype_id, space, &
                      attr_id, hdferr)
  dimsf = 1
  CALL h5awrite_f(attr_id, atype_id, coorlist, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)


  CALL h5dclose_f(dset , hdferr)



  ! PART 2: latitude

  ! initialize data
  lat_data(1)= -90.0
  DO i = 2, nx1
    lat_data(i) = lat_data(i-1)+1.0
  ENDDO

  !
  CALL h5screate_simple_f(1, lat_dims, space, hdferr)
  !
  CALL h5dcreate_f(file, LATITUDE, h5t_ieee_f32le, space, dset, hdferr)
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, lat_data, lat_dims, hdferr)
  !
  ! add 1st attribute: long_name
  !
  dimsf = 1
  attr_data= "latitude"
  strlen= 8
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, hdferr)
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: units
  !
  dimsf = 1
  attr_data= "degrees_north"
  strlen= 13
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  !
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, hdferr)
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)


  CALL h5dclose_f(dset , hdferr)



  ! PART 3: longitude

  ! initialize data
  lon_data(1)= -180.0
  DO i = 2, ny1
    lon_data(i) = lon_data(i-1)+1.0
  ENDDO

  !
  CALL h5screate_simple_f(1, lon_dims, space, hdferr)
  !
  CALL h5dcreate_f(file, LONGITUDE, h5t_ieee_f32le, space, dset, hdferr)
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, lon_data, lon_dims, hdferr)

  !
  ! add 1st attribute: long_name
  !
  dimsf = 1
  attr_data= "longitude"
  strlen= 9
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  !
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, hdferr)
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: units
  !
  dimsf = 1
  attr_data= "degrees_east"
  strlen= 12
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  !
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, hdferr)
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)


  CALL h5dclose_f(dset , hdferr)





  ! PART 4: pressure

  ! initialize data
  press_data(1)= 1000.0
  DO i = 2, nz1
    press_data(i) = press_data(i-1)-37.5
  ENDDO

  !
  CALL h5screate_simple_f(1, press_dims, space, hdferr)
  !
  CALL h5dcreate_f(file, PRESSURE, h5t_ieee_f32le, space, dset, hdferr)
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, press_data, press_dims, hdferr)

  !
  ! add 1st attribute: long_name
  !
  dimsf = 1
  attr_data= "pressure"
  strlen= 8
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  !
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, hdferr)
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: units
  !
  dimsf = 1
  attr_data= "hpa"
  strlen= 3
  !
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  !
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  !
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  !
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, hdferr)
  !
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dimsf, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  CALL h5dclose_f(dset , hdferr)

  !
  ! Close and release resources.
  !
  CALL h5fclose_f(file , hdferr)
  !
END PROGRAM main
