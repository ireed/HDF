c created by:  Isayah Reed <ireed@hdfgroup.org>
c              May  4, 2011
c  
c  This is a simple program that demonstrates how to use HDF-EOS5 to add CF
c  attributes to an h5 file. Creates a swath called Swath1 and adds the
c  temperature dataset with 2 dimensions: latitude and longitude. The
c  temperature dataset has the CF attributes: units, long_name, _FillValue,
c  coordinates, valid_min, valid_max, valid_range, scale_factor, add_offset.
c  Latitude and longitude have the CF attributes units and long_name. Outputs
c  data to general.h5
c  
c Compile:
c  $(HDF5_DIR)/bin/h5fc swath.f -I$(EOS_DIR)/include -L$(EOS_DIR)/lib -lhe5_hdfeos -lGctp

      program     swath

      implicit    none
      include     'hdfeos5.inc'

      integer     status
      integer     he5_swopen
      integer     he5_swcreate
      integer     he5_swsetfill
      integer     he5_swdefgfld
      integer     he5_swdefdfld
      integer     he5_swdefdim
      integer     he5_swwrfld
      integer     he5_swwrlattr
      integer     he5_swdetach
      integer     he5_swclose
      integer     swfid, swid

      integer*4   dtrack, strlen(2), value, NX1, NY1
      integer*4   start(2), stride(2), edge(2)
      parameter   (NX1=180, NY1=360)
      real   lon_array(NX1,NY1), lat_array(NX1,NY1)
      real   temp_data(NY1,NX1)
      character*50 attr(2)
      integer     i, j

      integer     FAIL
      parameter   (FAIL=-1)

c     Open the HDF-EOS file, "swath.he5" using "READ/WRITE" access code
c     -----------------------------------------------------------------
      swfid = he5_swopen('swath.he5',HE5F_ACC_TRUNC)

      swid = he5_swcreate(swfid, "Swath1")

c   Initialize temperature data.
      do i = 1, 60
          do j = 1, NY1
              temp_data(i,j) = 280.0
          enddo
      enddo
      do i = 61, 120
          do j = 1, NY1
              temp_data(i,j) = 300.0
          enddo
      enddo
      do i = 121, NX1
          do j = 1, NY1
              temp_data(i,j) = 280.0
          enddo
      enddo

c   Initialize latitude data.
        do i=1, NX1
            lat_array(i,1)= -90.0;
        enddo
        do i=1, NX1
            do j=2, NY1
                lat_array(i,j)= lat_array(i,j-1)+1.0;
            enddo
        enddo

c   Initialize longitude data.
        do j=1, NY1
            lon_array(1,j)= -180.0;
        enddo
        do j=2, NY1
            do i=2, NX1
                lon_array(i,j)= lon_array(i-1,j)+1.0;
            enddo
        enddo

c   define latitude and longitude dimensions
      dtrack = NX1
      status = he5_swdefdim(swid, "GeoXtrack", dtrack)

      dtrack = NY1
      status = he5_swdefdim(swid, "GeoTrack", dtrack)


c     Define Geolocation and Data fields
c     ----------------------------------
c     ---------------------------------------------------------------
c     We define six fields.  The first three, "Time", "Longitude"
c	  and "Latitude" are geolocation fields and thus we use the
c	  geolocation dimensions "GeoTrack" and "GeoXtrack" in the field
c	  definitions.  We also must specify the data type using the
c	  standard HDF data type codes.  In this example the geolocation
c	  are 4-byte (32 bit) floating point numbers.
c     
c	  The next three fields are data fields.  Note that either
c	  geolocation or data dimensions can be used. 
c     ---------------------------------------------------------------

c   set _FillValue
        value= -999
        status = he5_swsetfill(swid, "temp", HE5T_NATIVE_FLOAT, value)            

c   define data fields
        status = he5_swdefdfld(swid, "temp", "GeoTrack,GeoXtrack",
     1      " ", HE5T_NATIVE_FLOAT, 0)     

c   write the temperature dataset
        start(1)= 0
        start(2)= 0
        stride(1)= 1
        stride(2)= 1
        edge(1)= NY1
        edge(2)= NX1
        status= he5_swwrfld(swid, "temp", start, stride, edge, 
     1          temp_data)

c   add long_name attribute
        status = he5_swwrlattr(swid,"temp","long_name",
     1      HE5T_NATIVE_CHAR,11,"temperature")

c   add units attribute
        strlen(1)= 6
        status = he5_swwrlattr(swid,"temp","units",
     1      HE5T_NATIVE_CHAR,strlen,"kelvin")

c   add valid_min attribute
        strlen(1)= 1
        status = he5_swwrlattr(swid,"temp","valid_min",
     1      HE5T_NATIVE_FLOAT,strlen,0.0)

c   add valid_max attribute
        status = he5_swwrlattr(swid,"temp","valid_max",
     1      HE5T_NATIVE_FLOAT,strlen,400.0)

c   add valid_range attribute
        strlen(1)= 2
        edge(1)= 275
        edge(2)= 305
        status = he5_swwrlattr(swid,"temp","valid_range",
     1      HE5T_NATIVE_FLOAT,strlen,edge)

c   add add_offset attribute
        status = he5_swwrlattr(swid,"temp","add_offset",
     1      HE5T_NATIVE_FLOAT,1,1)

c   add scale_factor
        status = he5_swwrlattr(swid,"temp","scale_factor",
     1      HE5T_NATIVE_FLOAT,1,1.01)

c    add coordinates attribute
        strlen(1)= 3
        strlen(2)= 3
        attr(1)= "lat"
        attr(2)= "lon"
        status = he5_swwrlattr(swid,"temp","coordinates",
     1      HE5T_NATIVE_CHAR,strlen,attr)



c   define latitude dimension
        status = he5_swdefgfld(swid, "lat", "GeoTrack,GeoXtrack",
     1      " ", HE5T_NATIVE_FLOAT, 0)
c   define longitude dimension
        status = he5_swdefgfld(swid, "lon", "GeoTrack,GeoXtrack",
     1      " ", HE5T_NATIVE_FLOAT, 0)




c   write latitude dataset
        status= he5_swwrfld(swid, "lat", start, stride, edge, 
     1      lat_array)

c   add long_name attribute
        strlen(1)= 8
        attr(1)= "latitude"
        status = he5_swwrlattr(swid,"lat","long_name",
     1      HE5T_NATIVE_CHAR,strlen,attr)

c   add units attribute
        status = he5_swwrlattr(swid,"lat","units",
     1      HE5T_NATIVE_CHAR,13,"degrees_north")

        status= he5_swwrfld(swid, "lon", start, stride, edge, 
     1      lon_array)
            




c   write longitude dataset
        edge(1)= NX1
        edge(2)= NY1
        status= he5_swwrfld(swid, "lon", start, stride, edge,
     1      lon_array)

c   add long_name attribute
        strlen(1)= 9
        attr(1)= "longitude"
        status = he5_swwrlattr(swid,"lon","long_name",
     1      HE5T_NATIVE_CHAR,strlen,attr)

c   add units attribute
        status = he5_swwrlattr(swid,"lon","units",
     1      HE5T_NATIVE_CHAR,12,"degrees_east")


      
c     Detach from the swath
c     ---------------------      
      status = he5_swdetach(swid)
      
c     Close the file
c     --------------      
      status = he5_swclose(swfid)

      stop
      end




