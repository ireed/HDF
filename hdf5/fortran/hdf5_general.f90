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
 ! created by:  Isayah Reed
 !              May  4, 2011
 !
 ! Purpose:
 !  This is 1 of 3 simple HDF5 programs that demonstrates how to add CF attributes
 !  to an h5 file. Creates a file with 3 datasets: lat, lon, temp. Lat contains
 !  the CF attributes: units, long_name, and standard_name. Lon has the same CF
 !  attributes as the latitude dataset. Temp contains the CF attributes: units,
 !  long_name, _FillValue, coordinates, valid_min, valid_max, valid_range,
 !  scale_factor, add_offset. Outputs data to general.h5
 !  
 ! Compile:
 !  $(HDF5_DIR)/bin/h5fc hdf5_general.f90 -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl  -lhdf5hl_fortran -lhdf5_fortran -lz -lm
 !

PROGRAM main

  USE HDF5

  IMPLICIT NONE

  CHARACTER(LEN=10), PARAMETER :: filename = "general.h5"
  CHARACTER(LEN=4) , PARAMETER :: TEMP= "temp"
  CHARACTER(LEN=3) , PARAMETER :: LAT= "lat"
  CHARACTER(LEN=3) , PARAMETER :: LON= "lon"
  CHARACTER(LEN=9) , PARAMETER :: LONGNAME = "long_name"
  CHARACTER(LEN=5) , PARAMETER :: UNITS = "units"
  CHARACTER(LEN=10) , PARAMETER :: FILLVALUE = "_FillValue"
  CHARACTER(LEN=11) , PARAMETER :: COORDINATES = "coordinates"
  CHARACTER(LEN=11) , PARAMETER :: VALIDRANGE= "valid_range"
  CHARACTER(LEN=9) , PARAMETER :: VALIDMAX = "valid_max"
  CHARACTER(LEN=9) , PARAMETER :: VALIDMIN = "valid_min"
  CHARACTER(LEN=12) , PARAMETER :: SCALEFACTOR = "scale_factor"
  CHARACTER(LEN=10) , PARAMETER :: ADDOFFSET = "add_offset"
  CHARACTER(LEN=13) , PARAMETER :: STDNAME = "standard_name"

  INTEGER          , PARAMETER :: nx1     = 180
  INTEGER          , PARAMETER :: ny1     = 360

  INTEGER :: i, j, strlen, status
  INTEGER(HID_T) :: file, space, dset, attr_id, attr_space ! handles
  INTEGER(HID_T) :: atype_id, vtype_id
  CHARACTER(LEN=3), DIMENSION(1:2) ::  coorlist  ! coordinates {'lat', 'lon'}
  INTEGER(HSIZE_T), DIMENSION(1:2) :: temp_dims = (/ny1, nx1/) ! size read/write buffer
  INTEGER(HSIZE_T), DIMENSION(1:1) :: lat_dims = (/nx1/) 
  INTEGER(HSIZE_T), DIMENSION(1:1) :: lon_dims = (/ny1/) 
  INTEGER(HSIZE_T), DIMENSION(1)  :: dimsf ! dimension for scalars
  INTEGER(HSIZE_T), DIMENSION(2) :: dimsa   ! dimensions for arrays
  DOUBLE PRECISION, DIMENSION(1:nx1,1:ny1) :: temp_data
  DOUBLE PRECISION, DIMENSION(nx1) :: lat_data
  DOUBLE PRECISION, DIMENSION(ny1) :: lon_data
  DOUBLE PRECISION :: value
  DOUBLE PRECISION, DIMENSION(1:2) :: arrdata 
  CHARACTER(LEN=20) :: strdata
  CHARACTER, DIMENSION(2) ::  hvl_data


  ! Initialize FORTRAN interface.
  CALL h5open_f(status)
  
  ! Initialize temperature data.
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
  
  ! Create a new file using the default properties.
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, status)


  ! PART 1: temperature

  ! Create dataspace.  Setting size to be the current size of the 
  ! temperature array
  CALL h5screate_simple_f(2, temp_dims, space, status)
  
  ! Create the dataset.  We will use all default properties for this
  ! example.
  CALL h5dcreate_f(file, TEMP, h5t_ieee_f32le, space, dset, status)
  
  ! Write the data to the dataset.
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, temp_data, temp_dims, status)
  CALL h5sclose_f(space, status)
  
  ! add 1st attribute: units
  strdata= "kelvin"
  strlen= 6
  ! Create scalar data space for the attribute.
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  ! Create datatype for the attribute.
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  ! set the size of the datatype to be the size of the string "kelvin"
  CALL h5tset_size_f(atype_id, 6, status)
  ! set string as NULLTERM to imitate C-style strings
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  !Create dataset attribute.
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, status)
  ! Write the attribute data.
  dimsf = 1
  CALL h5awrite_f(attr_id, atype_id, "kelvin", dimsf, status)
  ! Close the dataspace
  CALL h5sclose_f(space, status)
  ! Close the attribute.
  CALL h5aclose_f(attr_id, status)

  ! add 2nd attribute: fillvalue
  dimsf = 1
  value= -999.0
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, status)
  CALL h5acreate_f(dset, FILLVALUE, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add 3rd attribute: long name
  dimsf = 1
  strdata= "temperature"
  strlen= 11
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add 4th attribute: coordinates
  dimsf = 2
  coorlist(1)= "lat"
  coorlist(2)= "lon"
  strlen= 3
  CALL h5screate_simple_f(1, dimsf, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, COORDINATES, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, atype_id, coorlist, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add attribute: add_offset
  dimsf = 1
  value= 5
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, status)
  CALL h5acreate_f(dset, ADDOFFSET, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add attribute: SCALEFACTOR
  dimsf = 1
  value= 1.01
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, status)
  CALL h5acreate_f(dset, SCALEFACTOR, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)


  ! add attribute: valid_min
  dimsf = 1
  value= 0
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, status)
  CALL h5acreate_f(dset, VALIDMIN, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add attribute: valid_max
  dimsf = 1
  value= 400
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, status)
  CALL h5acreate_f(dset, VALIDMAX, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add attribute: valid_range
  dimsf = 2
  arrdata(1)= 9
  arrdata(2)= 400
  CALL h5screate_simple_f(1, dimsf, space, status)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, status)
  !CALL h5tset_size_f(atype_id, H5T_VLEN_F, status)
  CALL h5acreate_f(dset, VALIDRANGE, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, atype_id, arrdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  CALL h5dclose_f(dset, status)
  

  ! PART 2: latitude

  ! initialize data
  lat_data(1)= -90.0
  DO i = 2, nx1
    lat_data(i) = lat_data(i-1)+1.0
  ENDDO

  CALL h5screate_simple_f(1, lat_dims, space, status)
  CALL h5dcreate_f(file, LAT, h5t_ieee_f32le, space, dset, status)
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, lat_data, lat_dims, status)
  CALL h5sclose_f(space, status)

  ! add 1st attribute: long_name
  dimsf = 1
  strdata= "latitude"
  strlen= 8
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, status)
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add 2nd attribute: units
  dimsf = 1
  strdata= "degrees_north"
  strlen= 13
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, status)
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add standard name
  dimsf = 1
  strdata= "latitude"
  strlen= 8
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, STDNAME, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  CALL h5dclose_f(dset, status)

  ! PART 3: longitude

  ! initialize data
  lon_data(1)= -180.0
  DO i = 2, ny1
    lon_data(i) = lon_data(i-1)+1.0
  ENDDO

  CALL h5screate_simple_f(1, lon_dims, space, status)
  CALL h5dcreate_f(file, LON, h5t_ieee_f32le, space, dset, status)
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, lon_data, lon_dims, status)
  CALL h5sclose_f(space, status)

  ! add 1st attribute: long_name
  dimsf = 1
  strdata= "longitude"
  strlen= 9
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, status)
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add 2nd attribute: units
  dimsf = 1
  strdata= "degrees_east"
  strlen= 12
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, UNITS, atype_id, space, &
                      attr_id, status)
  dimsf = 2
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  ! add standard name
  dimsf = 1
  strdata= "longitude"
  strlen= 9
  CALL h5screate_f(H5S_SCALAR_F, space, status)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, status)
  CALL h5tset_size_f(atype_id, strlen, status)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, status)
  CALL h5acreate_f(dset, STDNAME, atype_id, space, &
                      attr_id, status)
  CALL h5awrite_f(attr_id, atype_id, strdata, dimsf, status)
  CALL h5sclose_f(space, status)
  CALL h5aclose_f(attr_id, status)

  CALL h5dclose_f(dset , status)

  ! Close and release resources.
  CALL h5fclose_f(file , status)
  !
END PROGRAM main
