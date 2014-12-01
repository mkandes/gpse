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
!     Saturday, November 29th, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: io_read_bin_real3
      PUBLIC :: io_write_bin_real3
      PUBLIC :: io_read_bin_cmplx3
      PUBLIC :: io_write_bin_cmplx3
      PUBLIC :: io_read_vtk_real3
      PUBLIC :: io_write_vtk_real3
      PUBLIC :: io_read_vtk_cmplx3
      PUBLIC :: io_write_vtk_cmplx3

      CONTAINS

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
