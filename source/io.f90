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
!     Sunday, December 21st, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: io_read_bin_psi
      PUBLIC :: io_read_bin_vex
      PUBLIC :: io_write_bin_psi
      PUBLIC :: io_write_bin_vex
      PUBLIC :: io_write_vtk_header
      PUBLIC :: io_write_vtk_xcoordinates
      PUBLIC :: io_write_vtk_ycoordinates
      PUBLIC :: io_write_vtk_zcoordinates
      PUBLIC :: io_write_vtk_repsi
      PUBLIC :: io_write_vtk_impsi

      CONTAINS

         SUBROUTINE io_read_bin_psi ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

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

         SUBROUTINE io_read_bin_vex ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 )

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

         SUBROUTINE io_write_bin_psi ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

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

         SUBROUTINE io_write_bin_vex ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 )

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
               WRITE ( UNIT = fileUnit , FMT = '(A24)'          ) 'SCALARS RePsi double 1'
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

               WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A24)' ) 'SCALARS ImPsi double 1'
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

      END MODULE

! =========================================================================
