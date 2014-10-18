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
!     Tuesday, October 14th, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: io_read_bin
      PUBLIC :: io_write_bin
      PUBLIC :: io_write_vtk

      CONTAINS

         SUBROUTINE io_read_bin ( fileName , Psi3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Psi3

            OPEN  ( UNIT = 500 , FILE = trim(fileName//'.bin') , ACTION = 'READ' , FORM = 'UNFORMATTED' , STATUS = 'OLD' )

               READ ( UNIT = 500 ) Psi3

            CLOSE ( UNIT = 500 , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_bin ( fileName , Psi3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            COMPLEX, DIMENSION ( : , : , : ), INTENT ( IN ) :: Psi3

            OPEN  ( UNIT = 600 , FILE = trim(fileName//'.bin') , ACTION = 'WRITE' , FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = 600 ) Psi3

            CLOSE ( UNIT = 600 , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE io_write_vtk ( fileName , fileNumber , nX , nY , nZ , X , Y , Z , Vex3 , Psi3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
           
            INTEGER, INTENT ( IN ) :: fileNumber
            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN ) :: nZ

            REAL, DIMENSION ( : ), INTENT ( IN ) :: X
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( : , : , : ), INTENT ( IN ) :: Vex3
            
            COMPLEX, DIMENSION ( : , : , : ), INTENT ( IN ) :: Psi3

            CHARACTER ( LEN = 4 ) :: fileNumberString

            INTEGER :: j , k , l

            WRITE ( UNIT = fileNumberString , FMT = '(I4.4)' ) fileNumber
            OPEN  ( UNIT = fileNumber , FILE = trim(fileName//fileNumberString//'.vtk') , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = fileNumber , FMT = '(A26)' ) '# vtk DataFile Version 3.0'
               WRITE ( UNIT = fileNumber , FMT = '(A31)' ) 'GPSE STANDARD LEGACY VTK FORMAT'
               WRITE ( UNIT = fileNumber , FMT = '(A5)' ) 'ASCII'
               WRITE ( UNIT = fileNumber , FMT = '(A24)' ) 'DATASET RECTILINEAR_GRID'
               WRITE ( UNIT = fileNumber , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX , nY , nZ
               WRITE ( UNIT = fileNumber , FMT = '(A13,1X,I4.1,1X,A6)' ) 'X_COORDINATES' , nX , 'double'
               DO j = 1 , nX

                  WRITE ( UNIT = fileNumber , FMT = * ) X ( j )

               END DO
               WRITE ( UNIT = fileNumber , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY , 'double'
               DO k = 1 , nY

                  WRITE ( UNIT = fileNumber , FMT = * ) Y ( k )

               END DO
               WRITE ( UNIT = fileNumber , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , nX , 'double'
               DO l = 1 , nZ

                  WRITE ( UNIT = fileNumber , FMT = * ) Z ( l )

               END DO
               WRITE ( UNIT = fileNumber , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , nX * nY * nZ
               WRITE ( UNIT = fileNumber , FMT = '(A20)' ) 'SCALARS Vex double 1'
               WRITE ( UNIT = fileNumber , FMT = '(A20)' ) 'LOOKUP_TABLE default'
               DO l = 1 , nZ

                  DO k = 1 , nY

                     DO j = 1 , nX

                        WRITE ( UNIT = fileNumber , FMT = * ) Vex3 ( j , k , l ) 

                     END DO

                  END DO

               END DO
               WRITE ( UNIT = fileNumber , FMT = '(A22)' ) 'SCALARS RePsi double 1'
               WRITE ( UNIT = fileNumber , FMT = '(A20)' ) 'LOOKUP_TABLE default'
               DO l = 1 , nZ

                  DO k = 1 , nY

                     DO j = 1 , nX

                        WRITE ( UNIT = fileNumber , FMT = * ) REAL ( Psi3 ( j , k , l ) )

                     END DO

                  END DO

               END DO
               WRITE ( UNIT = fileNumber , FMT = '(A22)' ) 'SCALARS ImPsi double 1'
               WRITE ( UNIT = fileNumber , FMT = '(A20)' ) 'LOOKUP_TABLE default'
               DO l = 1 , nZ

                  DO k = 1 , nY

                     DO j = 1 , nX

                        WRITE ( UNIT = fileNumber , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )

                     END DO

                  END DO

               END DO

            CLOSE ( UNIT = fileNumber , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
