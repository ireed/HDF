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
 !  to an h5 file. This example shows how to use chunking and compression. A file
 !  is created with 3 datasets: lat, lon, temp. Lat contains the CF attributes:
 !  units and long_name. Lon has the same CF attributes as the latitude dataset.
 !  Temp contains the CF attributes: units, long_name, _FillValue, coordinates.
 !  It is has a chunk size is 900x1800. The deflate compression is used with a
 !  compression level of 1. Outputs data to chunk_compress.h5
 !  
 ! Compile:
 !  $(HDF5_DIR)/bin/h5fc hdf5_chunk_compress.f90 -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl  -lhdf5hl_fortran -lhdf5_fortran -lsz -lz -lm
 !

PROGRAM main

  USE HDF5

  IMPLICIT NONE

  CHARACTER(LEN=17), PARAMETER :: filename = "chunk_compress.h5"
  CHARACTER(LEN=4) , PARAMETER :: TEMP= "temp"
  CHARACTER(LEN=3) , PARAMETER :: LAT= "lat"
  CHARACTER(LEN=3) , PARAMETER :: LON= "lon"
  CHARACTER(LEN=9) , PARAMETER :: LONGNAME = "long_name"
  CHARACTER(LEN=5) , PARAMETER :: UNITS = "units"
  CHARACTER(LEN=10) , PARAMETER :: FILLVALUE = "_FillValue"
  CHARACTER(LEN=11) , PARAMETER :: COORDINATES = "coordinates"

  INTEGER          , PARAMETER :: nx1     = 1800
  INTEGER          , PARAMETER :: ny1     = 3600

  INTEGER :: hdferr
  INTEGER(HID_T) :: file, space, dset, attr_id, attr_space ! handles
  INTEGER(HID_T) :: atype_id, prop_id
  CHARACTER(LEN=20), DIMENSION(2) ::  attr_data  ! Attribute string data
  ! switch dimensions to adhere to C standards
  INTEGER(HSIZE_T), DIMENSION(1:2) :: chunk_dims = (/900, 900/)
  INTEGER(HSIZE_T), DIMENSION(1:2) :: temp_dims = (/ny1, nx1/)
  INTEGER(HSIZE_T), DIMENSION(1:1) :: lat_dims = (/nx1/) 
  INTEGER(HSIZE_T), DIMENSION(1:1) :: lon_dims = (/ny1/) 
  INTEGER(HSIZE_T), DIMENSION(1) :: dims1
  INTEGER(HSIZE_T), DIMENSION(2) :: dims2
  DOUBLE PRECISION, DIMENSION(1:nx1,1:ny1) :: temp_data
  DOUBLE PRECISION, DIMENSION(nx1) :: lat_data
  DOUBLE PRECISION, DIMENSION(ny1) :: lon_data
  DOUBLE PRECISION, DIMENSION(1) :: fill_data
  INTEGER :: i, j, strlen







  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  !
  ! Initialize temperature data.
  !
  !DO i = 1, nx1
  !   DO j = 1, ny1 
  !      wdata(i,j) = (i-1)*(j-1)-(j-1)
  !   ENDDO
  !ENDDO
  DO i = 1, 600
     DO j = 1, ny1 
        temp_data(i,j) = 280.0
     ENDDO
  ENDDO
  DO i = 601, 1200
     DO j = 1, ny1 
        temp_data(i,j) = 300.0
     ENDDO
  ENDDO
  DO i = 1201, nx1
     DO j = 1, ny1 
        temp_data(i,j) = 280.0
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
  CALL h5screate_simple_f(2, temp_dims, space, hdferr)
  !
  ! Set chunking property parameters
  ! chunk size = 900x1800 
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
  CALL h5pset_chunk_f(prop_id, 2, chunk_dims, hdferr)
  ! set compression
  CALL h5pset_deflate_f(prop_id, 1, hdferr)
  !
  ! Create the dataset.
  !
  CALL h5dcreate_f(file, TEMP, h5t_ieee_f32le, space, dset, hdferr, &
            prop_id)
  !
  ! Write the data to the dataset.
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, temp_data, temp_dims, hdferr)
  !
  CALL h5sclose_f(space, hdferr)
  !
  ! add 1st attribute: units
  !
  dims1(1) = 1
  attr_data(1)= "kelvin"
  strlen= 6
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
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  !
  ! Close/release dataspace and attribute resources
  !
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: fillvalue
  !
  dims1(1) = 1
  fill_data(1)= -999.0
  CALL h5screate_simple_f(1, dims1, space, hdferr)
  CALL h5tcopy_f(h5t_ieee_f32le, atype_id, hdferr)
  CALL h5acreate_f(dset, FILLVALUE, atype_id, space, &
                      attr_id, hdferr)
  dims1(1) = 1
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, fill_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 3rd attribute: long name
  !
  dims1(1) = 1
  attr_data(1)= "temperature"
  strlen= 11
  CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  CALL h5acreate_f(dset, LONGNAME, atype_id, space, &
                      attr_id, hdferr)
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 4th attribute: coordinates
  !
  dims1(1) = 2
  attr_data(1)= "lat"
  attr_data(2)= "lon"
  strlen= 20
  CALL h5screate_simple_f(1, dims1, space, hdferr)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
  CALL h5tset_size_f(atype_id, strlen, hdferr)
  CALL h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)
  CALL h5acreate_f(dset, COORDINATES, atype_id, space, &
                      attr_id, hdferr)
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! Close the dataset
  !
  CALL h5dclose_f(dset , hdferr)




  ! PART 2: latitude

  ! initialize data
  lat_data(1)= -90.0
  DO i = 2, nx1
    lat_data(i) = lat_data(i-1)+0.1
  ENDDO

  !
  CALL h5screate_simple_f(1, lat_dims, space, hdferr)
  !
  CALL h5dcreate_f(file, LAT, h5t_ieee_f32le, space, dset, hdferr)
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, lat_data, lat_dims, hdferr)
  !
  ! add 1st attribute: long_name
  !
  dims1(1) = 1
  attr_data(1)= "latitude"
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
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: units
  !
  dims1(1) = 1
  attr_data(1)= "degrees_north"
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
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  CALL h5dclose_f(dset , hdferr)




  ! PART 3: longitude

  ! initialize data
  lon_data(1)= -180.0
  DO i = 2, ny1
    lon_data(i) = lon_data(i-1)+0.1
  ENDDO

  !
  CALL h5screate_simple_f(1, lon_dims, space, hdferr)
  !
  CALL h5dcreate_f(file, LON, h5t_ieee_f32le, space, dset, hdferr)
  !
  CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, lon_data, lon_dims, hdferr)
  !
  CALL h5sclose_f(space, hdferr)

  !
  ! add 1st attribute: long_name
  !
  dims1(1) = 1
  attr_data(1)= "longitude"
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
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)

  !
  ! add 2nd attribute: units
  !
  dims1(1) = 1
  attr_data(1)= "degrees_east"
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
  dims1(1) = 2
  CALL h5awrite_f(attr_id, atype_id, attr_data, dims1, hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5aclose_f(attr_id, hdferr)


  CALL h5dclose_f(dset , hdferr)


  !
  ! Close and release resources.
  !
  CALL h5fclose_f(file , hdferr)
  !
END PROGRAM main
