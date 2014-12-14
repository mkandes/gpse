! ==========================================================================
! NAME
!
!     io [ io ] - io module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     io is a Fortran module ... 
!
! OPTIONS
!
! SEE ALSO
!
! BUGS
!
! HISTORY
!
! AUTHOR(S)
!
!     Marty Kandes
!     Computational Science Research Center
!     San Diego State University
!     5500 Campanile Drive
!     San Diego, California 92182
!
! COPYRIGHT
!     
!     Copyright (c) 2014 Martin Charles Kandes
!
! LAST UPDATED
!
!     Friday, December 12th, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      INTEGER, PUBLIC :: ioByteSizeInt   = -1 
      INTEGER, PUBLIC :: ioByteSizeReal  = -1
      INTEGER, PUBLIC :: ioByteSizeCmplx = -1

      PUBLIC :: io_vtk_write_header
      PUBLIC :: io_vtk_write_title
      PUBLIC :: io_vtk_write_type
      PUBLIC :: io_vtk_write_geometry
      PUBLIC :: io_vtk_write_attributes


      PUBLIC :: io_get_byte_sizes

      PUBLIC :: io_read_real3
      PUBLIC :: io_read_cmplx3
      PUBLIC :: io_write_real3
      PUBLIC :: io_write_cmplx3

      PUBLIC :: io_read_psi3
      PUBLIC :: io_read_vex3
      PUBLIC :: io_write_psi3
      PUBLIC :: io_write_vex3

       

      PUBLIC :: io_write_vtk_header
      PUBLIC :: io_write_vtk_xcoordinates
      PUBLIC :: io_write_vtk_ycoordinates
      PUBLIC :: io_write_vtk_zcoordinates
      PUBLIC :: io_write_vtk_repsi
      PUBLIC :: io_write_vtk_impsi

      PUBLIC :: io_write_vtk_psi3
      PUBLIC :: io_write_vtk_vex3

      PUBLIC :: io_read_bin_real3
      PUBLIC :: io_write_bin_real3
      PUBLIC :: io_read_bin_cmplx3
      PUBLIC :: io_write_bin_cmplx3
      PUBLIC :: io_read_vtk_real3
      PUBLIC :: io_write_vtk_real3
      PUBLIC :: io_read_vtk_cmplx3
      PUBLIC :: io_write_vtk_cmplx3
                   
      CONTAINS

         SUBROUTINE io_vtk_write_header ( fileHeader , fileName , fileUnit , filePosition )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileHeader
            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition

            OPEN    ( UNIT = fileUnit , FILE = fileName , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
            
            WRITE   ( UNIT = fileUnit , FMT = * ) fileHeader
            INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE   ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_vtk_write_title ( fileTitle , fileName , fileUnit , filePosition ) ! Add exception: 256 characters maximum

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileTitle
            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition

            OPEN    ( UNIT = fileUnit , FILE = fileName , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
                 
            WRITE   ( UNIT = fileUnit , FMT = * ) fileTitle
            INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE   ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_vtk_write_type ( fileType , fileName , fileUnit , filePosition ) ! Add exception: only 'ASCII' or 'BINARY' is allowed
         
            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileType
            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition

            OPEN    ( UNIT = fileUnit , FILE = fileName , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

            WRITE   ( UNIT = fileUnit , FMT = * ) fileType
            INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE   ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_vtk_write_geometry ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE io_vtk_write_attributes ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE io_get_byte_sizes ( ioByteSizeInt, ioByteSizeReal , ioByteSizeCmplx )

            IMPLICIT NONE

            INTEGER, INTENT ( INOUT ) :: ioByteSizeInt
            INTEGER, INTENT ( INOUT ) :: ioByteSizeReal
            INTEGER, INTENT ( INOUT ) :: ioByteSizeCmplx

            ioByteSizeInt   = SIZEOF ( 0 )
            ioByteSizeReal  = SIZEOF ( 0.0 )
            ioByteSizeCmplx = SIZEOF ( CMPLX ( 0.0 , 0.0 ) )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_read_real3 ( filePrefix , fileUnit , filePos , Real3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: filePrefix

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePos

            REAL, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Real3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit

            OPEN ( UNIT = fileUnit , FILE = TRIM ( filePrefix//'-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'READ' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            READ ( UNIT = fileUnit , POS = filePos ) Real3
            INQUIRE ( UNIT = fileUnit , POS = filePos )
            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_read_cmplx3 ( filePrefix , fileUnit , filePos , Cmplx3 ) 

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: filePrefix

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePos

            COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Cmplx3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit

            OPEN ( UNIT = fileUnit , FILE = TRIM ( filePrefix//'-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'READ' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            READ ( UNIT = fileUnit , POS = filePos ) Cmplx3
            INQUIRE ( UNIT = fileUnit , POS = filePos )
            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_real3 ( filePrefix , fileUnit , filePos , Real3 ) 

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: filePrefix

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePos

            REAL, DIMENSION ( : , : , : ), INTENT ( IN ) :: Real3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit

            OPEN ( UNIT = fileUnit , FILE = TRIM ( filePrefix//'-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'WRITE' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            WRITE ( UNIT = fileUnit , POS = filePos ) Real3
            INQUIRE ( UNIT = fileUnit , POS = filePos ) 
            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_cmplx3 ( filePrefix , fileUnit , filePos , Cmplx3 ) 

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: filePrefix

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePos

            COMPLEX, DIMENSION ( : , : , : ), INTENT ( IN ) :: Cmplx3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit

            OPEN ( UNIT = fileUnit , FILE = TRIM ( filePrefix//'-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'WRITE' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            WRITE ( UNIT = fileUnit , POS = filePos ) Cmplx3
            INQUIRE ( UNIT = fileUnit , POS = filePos ) 
            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_read_psi3 ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nXa 
            INTEGER, INTENT ( IN    ) :: nXb 
            INTEGER, INTENT ( IN    ) :: nXbc 
            INTEGER, INTENT ( IN    ) :: nYa 
            INTEGER, INTENT ( IN    ) :: nYb 
            INTEGER, INTENT ( IN    ) :: nYbc 
            INTEGER, INTENT ( IN    ) :: nZa 
            INTEGER, INTENT ( IN    ) :: nZb 
            INTEGER, INTENT ( IN    ) :: nZbc 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'READ' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        READ ( UNIT = fileUnit , POS = filePosition ) Psi3 ( j , k , l )
                        INQUIRE ( UNIT = fileUnit , POS = filePosition )

                     END DO

                  END DO

               END DO

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_read_vex3 ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nXa 
            INTEGER, INTENT ( IN    ) :: nXb 
            INTEGER, INTENT ( IN    ) :: nXbc 
            INTEGER, INTENT ( IN    ) :: nYa 
            INTEGER, INTENT ( IN    ) :: nYb 
            INTEGER, INTENT ( IN    ) :: nYbc 
            INTEGER, INTENT ( IN    ) :: nZa 
            INTEGER, INTENT ( IN    ) :: nZb 
            INTEGER, INTENT ( IN    ) :: nZbc 

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l 

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'READ' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )

               DO l = nZa , nZb 

                  DO k = nYa , nYb 

                     DO j = nXa , nXb 

                        READ ( UNIT = fileUnit , POS = filePosition ) Vex3 ( j , k , l ) 
                        INQUIRE ( UNIT = fileUnit , POS = filePosition )

                     END DO

                  END DO

               END DO

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_psi3 ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nXa
            INTEGER, INTENT ( IN    ) :: nXb
            INTEGER, INTENT ( IN    ) :: nXbc 
            INTEGER, INTENT ( IN    ) :: nYa
            INTEGER, INTENT ( IN    ) :: nYb
            INTEGER, INTENT ( IN    ) :: nYbc 
            INTEGER, INTENT ( IN    ) :: nZa
            INTEGER, INTENT ( IN    ) :: nZb
            INTEGER, INTENT ( IN    ) :: nZbc 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'WRITE' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        WRITE ( UNIT = fileUnit , POS = filePosition ) Psi3 ( j , k , l )
                        INQUIRE ( UNIT = fileUnit , POS = filePosition )

                     END DO

                  END DO

               END DO

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE 

         SUBROUTINE io_write_vex3 ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nXa 
            INTEGER, INTENT ( IN    ) :: nXb 
            INTEGER, INTENT ( IN    ) :: nXbc 
            INTEGER, INTENT ( IN    ) :: nYa 
            INTEGER, INTENT ( IN    ) :: nYb 
            INTEGER, INTENT ( IN    ) :: nYbc 
            INTEGER, INTENT ( IN    ) :: nZa 
            INTEGER, INTENT ( IN    ) :: nZb 
            INTEGER, INTENT ( IN    ) :: nZbc 

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l 

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN ( UNIT = fileUnit , FILE = TRIM ( 'vex-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'WRITE' ,FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
    
               DO l = nZa , nZb 

                  DO k = nYa , nYb 

                     DO j = nXa , nXb

                        WRITE ( UNIT = fileUnit , POS = filePosition ) Vex3 ( j , k , l )
                        INQUIRE ( UNIT = fileUnit , POS = filePosition )

                     END DO

                  END DO

               END DO

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

!         SUBROUTINE io_write_( characterString , fileName , fileUnit , filePosition )
!
!            IMPLICIT NONE
!
!            CHARACTER ( LEN = * ), INTENT ( IN ) :: characterString
!            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
!
!            INTEGER              , INTENT ( IN ) :: fileUnit
!            INTEGER              , INTENT ( IN ) :: filePosition
!
!            OPEN ( UNIT = fileUnit , FILE = fileName , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
!
!            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )            
!            
!
!            RETURN
!
!         END SUBROUTINE

         SUBROUTINE io_write_vtk_header ( fileName , fileUnit , filePosition , nX , nY , nZ )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nX
            INTEGER, INTENT ( IN    ) :: nY
            INTEGER, INTENT ( IN    ) :: nZ

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = fileUnit , FMT = '(A26)' ) '# vtk DataFile Version 3.0' ! VTK file identifier and version number
               WRITE ( UNIT = fileUnit , FMT = '(A26)' ) 'STANDARD LEGACY VTK FORMAT' ! VTK file header; 256 characters maximum
               WRITE ( UNIT = fileUnit , FMT = '(A5)'  ) 'ASCII'                      ! VTK file format: ASCII or BINARY
               WRITE ( UNIT = fileUnit , FMT = '(A24)' ) 'DATASET RECTILINEAR_GRID'             ! VTK dataset structure
               WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX , nY , nZ

               INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_xcoordinates ( fileName , fileUnit , filePosition , nX , nXa , nXb , nXbc , X )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN    ) :: fileName

            INTEGER              , INTENT ( IN    ) :: fileUnit
            INTEGER              , INTENT ( INOUT ) :: filePosition
            INTEGER              , INTENT ( IN    ) :: nX
            INTEGER              , INTENT ( IN    ) :: nXa
            INTEGER              , INTENT ( IN    ) :: nXb
            INTEGER              , INTENT ( IN    ) :: nXbc
     
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A13,1X,I4.1,1X,A6)' ) 'X_COORDINATES' , nX , 'double'
               DO j = nXa , nXb

                  WRITE ( UNIT = fileUnit , FMT = * ) X ( j )

               END DO

               INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_ycoordinates ( fileName , fileUnit , filePosition , nY , nYa , nYb , nYbc , Y )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nY
            INTEGER, INTENT ( IN    ) :: nYa
            INTEGER, INTENT ( IN    ) :: nYb
            INTEGER, INTENT ( IN    ) :: nYbc

            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: k

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY , 'double'

               DO k = nYa , nYb

                  WRITE ( UNIT = fileUnit , FMT = * ) Y ( k )

               END DO

               INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_zcoordinates ( fileName , fileUnit , filePosition , mpiSource , nZ , nZa , nZb , nZbc , Z )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN ) :: mpiSource
            INTEGER, INTENT ( IN ) :: nZ
            INTEGER, INTENT ( IN ) :: nZa 
            INTEGER, INTENT ( IN ) :: nZb 
            INTEGER, INTENT ( IN ) :: nZbc

            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: l

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

               IF ( mpiSource == 0 ) THEN ! i.e., if mpiSource == MPI_MASTER

                  WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , nZ , 'double'
                  DO l = nZa , nZb 

                     WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

                  END DO

               ELSE

                  WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(F23.15)' ) Z ( nZa )

                  DO l = nZa + 1 , nZb

                     WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

                  END DO

               END IF

               INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_repsi ( fileName , fileUnit , filePosition , mpiSource , nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: mpiSource
            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN    ) :: nXa
            INTEGER, INTENT ( IN    ) :: nXb
            INTEGER, INTENT ( IN    ) :: nXbc
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN    ) :: nYa
            INTEGER, INTENT ( IN    ) :: nYb
            INTEGER, INTENT ( IN    ) :: nYbc
            INTEGER, INTENT ( IN ) :: nZ
            INTEGER, INTENT ( IN    ) :: nZa
            INTEGER, INTENT ( IN    ) :: nZb
            INTEGER, INTENT ( IN    ) :: nZbc

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
 
            IF ( mpiSource == 0 ) THEN

               WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , nX * nY * nZ
               WRITE ( UNIT = fileUnit , FMT = '(A24)'          ) 'SCALARS Re(Psi) double 1'
               WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'LOOKUP_TABLE default'

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        WRITE ( UNIT = fileUnit , FMT = * ) REAL ( Psi3 ( j , k , l ) )

                     END DO

                  END DO

               END DO

            ELSE

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb
                        
                        WRITE ( UNIT = fileUnit , POS = filePosition , FMT = * ) REAL ( Psi3 ( j , k , l ) )
                        INQUIRE ( UNIT = fileUnit , POS = filePosition )

                     END DO

                  END DO

               END DO

            END IF

            INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_impsi ( fileName , fileUnit , filePosition , mpiSource , nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3)

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: mpiSource
            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN    ) :: nXa
            INTEGER, INTENT ( IN    ) :: nXb
            INTEGER, INTENT ( IN    ) :: nXbc
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN    ) :: nYa
            INTEGER, INTENT ( IN    ) :: nYb
            INTEGER, INTENT ( IN    ) :: nYbc
            INTEGER, INTENT ( IN ) :: nZ
            INTEGER, INTENT ( IN    ) :: nZa
            INTEGER, INTENT ( IN    ) :: nZb
            INTEGER, INTENT ( IN    ) :: nZbc

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

            IF ( mpiSource == 0 ) THEN

               WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A24)' ) 'SCALARS Im(Psi) double 1'
               WRITE ( UNIT = fileUnit , FMT = '(A20)' ) 'LOOKUP_TABLE default'

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        WRITE ( UNIT = fileUnit , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )

                     END DO

                  END DO

               END DO

            ELSE

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        WRITE ( UNIT = fileUnit , POS = filePosition , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )
                        INQUIRE ( UNIT = fileUnit , POS = filePosition )

                     END DO

                  END DO

               END DO

            END IF

            INQUIRE ( UNIT = fileUnit , POS = filePosition )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_psi3 ( fileUnit , filePosition , nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , X , Y , Z , Psi3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: fileUnit
            INTEGER, INTENT ( INOUT ) :: filePosition
            INTEGER, INTENT ( IN    ) :: nX
            INTEGER, INTENT ( IN    ) :: nXa
            INTEGER, INTENT ( IN    ) :: nXb
            INTEGER, INTENT ( IN    ) :: nXbc
            INTEGER, INTENT ( IN    ) :: nY
            INTEGER, INTENT ( IN    ) :: nYa
            INTEGER, INTENT ( IN    ) :: nYb
            INTEGER, INTENT ( IN    ) :: nYbc
            INTEGER, INTENT ( IN    ) :: nZ
            INTEGER, INTENT ( IN    ) :: nZa
            INTEGER, INTENT ( IN    ) :: nZb
            INTEGER, INTENT ( IN    ) :: nZbc

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

            CHARACTER ( LEN = 4 ) :: fileUnitChar

            INTEGER :: j , k , l

            WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
            OPEN  ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = fileUnit , FMT = '(A26)'                         ) '# vtk DataFile Version 3.0'
               WRITE ( UNIT = fileUnit , FMT = '(A31)'                         ) 'STANDARD LEGACY VTK FORMAT'
               WRITE ( UNIT = fileUnit , FMT = '(A5)'                          ) 'ASCII'
               WRITE ( UNIT = fileUnit , FMT = '(A24)'                         ) 'DATASET RECTILINEAR_GRID'
               WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX , nY , nZ
               WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)'           ) 'X_COORDINATES' , nX , 'double'
               INQUIRE ( UNIT = fileUnit , POS = filePosition )
               DO j = nXa , nXb

                  WRITE ( UNIT = fileUnit , POS = filePosition , FMT = * ) X ( j ) ! j + nX * [ ( k - 1 ) + nY * ( l - 1 ) ]

               END DO
               WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY , 'double'
               DO k = nYa , nYb

                  WRITE ( UNIT = fileUnit , FMT = * ) Y ( k ) 

               END DO
!               WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , nX , 'double'
!               DO l = nZa , nZb
!
!                  WRITE ( UNIT = fileUnit , FMT = * ) Z ( l ) 
!
!               END DO
!               WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , nX * nY * nZ
!               WRITE ( UNIT = fileUnit , FMT = '(A28)'          ) 'SCALARS Re(Cmplx3) double 1'
!               WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'LOOKUP_TABLE default'
!               DO l = nZa , nZb
!
!                  DO k = nYa , nYb
!
!                     DO j = nXa , nXb
!
!                        WRITE ( UNIT = fileUnit , FMT = * ) REAL ( Psi3 ( j , k , l ) )
!
!                     END DO
!
!                  END DO
!
!               END DO
!               WRITE ( UNIT = fileUnit , FMT = '(A28)' ) 'SCALARS Im(Cmplx3) double 1'
!               WRITE ( UNIT = fileUnit , FMT = '(A20)' ) 'LOOKUP_TABLE default'
!               DO l = nZa , nZb
!
!                  DO k =  nYa , nYb
!
!                     DO j = nXa , nXb
!
!                        WRITE ( UNIT = fileUnit , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )
!
!                     END DO
!
!                  END DO
!
!               END DO
!
            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )
            

            RETURN

         END SUBROUTINE 

         SUBROUTINE io_write_vtk_vex3 ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE io_read_bin_real3 ( fileName , fileUnit , Real3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN ) :: fileUnit

         REAL, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Real3

         OPEN  ( UNIT = fileUnit , FILE = trim(fileName//'.bin') , ACTION = 'READ' , FORM = 'UNFORMATTED' , STATUS = 'OLD' )
         READ  ( UNIT = fileUnit ) Real3
         CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

         RETURN

         END SUBROUTINE

         SUBROUTINE io_write_bin_real3 ( fileName , fileUnit , Real3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN ) :: fileUnit

         REAL, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Real3

         OPEN  ( UNIT = fileUnit , FILE = trim(fileName//'.bin') , ACTION = 'WRITE' , FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
         WRITE ( UNIT = fileUnit ) Real3
         CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

         RETURN

         END SUBROUTINE

         SUBROUTINE io_read_bin_cmplx3 ( fileName , fileUnit , Cmplx3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
     
         INTEGER, INTENT ( IN ) :: fileUnit

         COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Cmplx3

         OPEN  ( UNIT = fileUnit , FILE = trim(fileName//'.bin') , ACTION = 'READ' , FORM = 'UNFORMATTED' , STATUS = 'OLD' )
         READ  ( UNIT = fileUnit ) Cmplx3
         CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

         RETURN

         END SUBROUTINE

         SUBROUTINE io_write_bin_cmplx3 ( fileName , fileUnit , Cmplx3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN ) :: fileUnit

         COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Cmplx3

         OPEN  ( UNIT = fileUnit , FILE = trim(fileName//'.bin') , ACTION = 'WRITE' , FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
         WRITE ( UNIT = fileUnit ) Cmplx3
         CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

         RETURN

         END SUBROUTINE

         SUBROUTINE io_read_vtk_real3 ( fileName , fileUnit , nX , nY , nZ , X , Y , Z , Real3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN    ) :: fileUnit
         INTEGER, INTENT ( INOUT ) :: nX
         INTEGER, INTENT ( INOUT ) :: nY
         INTEGER, INTENT ( INOUT ) :: nZ

         REAL, DIMENSION ( :         ), INTENT ( INOUT ) :: X
         REAL, DIMENSION ( :         ), INTENT ( INOUT ) :: Y
         REAL, DIMENSION ( :         ), INTENT ( INOUT ) :: Z
         REAL, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Real3

         RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_real3 ( fileName , fileUnit , nX , nY , nZ , X , Y , Z , Real3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN ) :: fileUnit
         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nY
         INTEGER, INTENT ( IN ) :: nZ

         REAL, DIMENSION ( :         ), INTENT ( IN ) :: X
         REAL, DIMENSION ( :         ), INTENT ( IN ) :: Y
         REAL, DIMENSION ( :         ), INTENT ( IN ) :: Z
         REAL, DIMENSION ( : , : , : ), INTENT ( IN ) :: Real3

         CHARACTER ( LEN = 4 ) :: fileUnitChar

         INTEGER :: j , k , l

         WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
         OPEN  ( UNIT = fileUnit , FILE = trim(fileName//fileUnitChar//'.vtk') , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
         WRITE ( UNIT = fileUnit , FMT = '(A26)'                         ) '# vtk DataFile Version 3.0'
         WRITE ( UNIT = fileUnit , FMT = '(A31)'                         ) 'STANDARD LEGACY VTK FORMAT'
         WRITE ( UNIT = fileUnit , FMT = '(A5)'                          ) 'ASCII'
         WRITE ( UNIT = fileUnit , FMT = '(A24)'                         ) 'DATASET RECTILINEAR_GRID'
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX , nY , nZ
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)'           ) 'X_COORDINATES' , nX , 'double'
         DO j = 1 , nX

            WRITE ( UNIT = fileUnit , FMT = * ) X ( j )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY , 'double'
         DO k = 1 , nY

            WRITE ( UNIT = fileUnit , FMT = * ) Y ( k )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , nX , 'double'
         DO l = 1 , nZ

            WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , nX * nY * nZ
         WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'SCALARS Re3d double 1'
         WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'LOOKUP_TABLE default'
         DO l = 1 , nZ

            DO k = 1 , nY

               DO j = 1 , nX

                  WRITE ( UNIT = fileUnit , FMT = * ) Real3 ( j , k , l )

               END DO

            END DO

         END DO
         CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

         RETURN

         END SUBROUTINE

         SUBROUTINE io_read_vtk_cmplx3 ( fileName , fileUnit , nX , nY , nZ , X , Y , Z , Cmplx3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN    ) :: fileUnit
         INTEGER, INTENT ( INOUT ) :: nX
         INTEGER, INTENT ( INOUT ) :: nY
         INTEGER, INTENT ( INOUT ) :: nZ
 
         REAL, DIMENSION ( : ), INTENT ( INOUT ) :: X
         REAL, DIMENSION ( : ), INTENT ( INOUT ) :: Y
         REAL, DIMENSION ( : ), INTENT ( INOUT ) :: Z

         COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Cmplx3

         RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk_cmplx3 ( fileName , fileUnit , nX , nY , nZ , X , Y , Z , Cmplx3 )

         IMPLICIT NONE

         CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

         INTEGER, INTENT ( IN ) :: fileUnit
         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nY
         INTEGER, INTENT ( IN ) :: nZ
 
         REAL, DIMENSION ( : ), INTENT ( IN ) :: X
         REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
         REAL, DIMENSION ( : ), INTENT ( IN ) :: Z

         COMPLEX, DIMENSION ( : , : , : ), INTENT ( IN ) :: Cmplx3

         CHARACTER ( LEN = 4 ) :: fileUnitChar

         INTEGER :: j , k , l

         WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
         OPEN  ( UNIT = fileUnit , FILE = trim(fileName//fileUnitChar//'.vtk') , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
         WRITE ( UNIT = fileUnit , FMT = '(A26)'                         ) '# vtk DataFile Version 3.0'
         WRITE ( UNIT = fileUnit , FMT = '(A31)'                         ) 'STANDARD LEGACY VTK FORMAT'
         WRITE ( UNIT = fileUnit , FMT = '(A5)'                          ) 'ASCII'
         WRITE ( UNIT = fileUnit , FMT = '(A24)'                         ) 'DATASET RECTILINEAR_GRID'
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX , nY , nZ
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)'           ) 'X_COORDINATES' , nX , 'double'
         DO j = 1 , nX

            WRITE ( UNIT = fileUnit , FMT = * ) X ( j )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY , 'double'
         DO k = 1 , nY

            WRITE ( UNIT = fileUnit , FMT = * ) Y ( k )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , nX , 'double'
         DO l = 1 , nZ

            WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , nX * nY * nZ
         WRITE ( UNIT = fileUnit , FMT = '(A28)'          ) 'SCALARS Re(Cmplx3) double 1'
         WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'LOOKUP_TABLE default'
         DO l = 1 , nZ

            DO k = 1 , nY

               DO j = 1 , nX

                  WRITE ( UNIT = fileUnit , FMT = * ) REAL ( Cmplx3 ( j , k , l ) )

               END DO

            END DO

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A28)' ) 'SCALARS Im(Cmplx3) double 1'
         WRITE ( UNIT = fileUnit , FMT = '(A20)' ) 'LOOKUP_TABLE default'
         DO l = 1 , nZ

            DO k = 1 , nY

               DO j = 1 , nX

                  WRITE ( UNIT = fileUnit , FMT = * ) AIMAG ( Cmplx3 ( j , k , l ) )

               END DO

            END DO

         END DO
         CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

         RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
