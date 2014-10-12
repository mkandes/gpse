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
!     Tuesday, October 7th, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV
!      USE            :: LIB_VTK_IO
      USE            :: GRID, ONLY: nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc

      IMPLICIT NONE
      PRIVATE

      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_SIMAG     = '.vti'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_SPOLY     = '.vtp'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_SRECT     = '.vtr'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_SSTRUCT   = '.vts'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_SUNSTRUCT = '.vtu'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_PIMAG     = '.pvti'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_PPOLY     = '.pvtp'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_PRECT     = '.pvtr'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_PSTRUCT   = '.pvts'
      CHARACTER ( LEN = * ), PARAMETER, PUBLIC :: VTK_XML_PUNSTRUCT = '.pvtu'

      INTEGER, PARAMETER :: MAX_LEN = 80

      CHARACTER ( LEN = MAX_LEN ), PUBLIC  :: vtkFileName
      CHARACTER ( LEN = MAX_LEN ), PUBLIC  :: vtkOutputFmt           = 'ASCII'
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkMeshTopology        = 'RECTILINEAR_GRID'
      CHARACTER ( LEN = MAX_LEN ), PUBLIC  :: vtkMode                = 'LEGACY'
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkTitle 
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkVariableBlockAction = 'CLOSE'
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkVariableLocation
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkVariableName

      INTEGER, PUBLIC :: vtkError   = -1

      PUBLIC :: read_bin
      PUBLIC :: write_bin
!      PUBLIC :: read_gpi
!      PUBLIC :: write_gpi
      PUBLIC :: write_vtk

      CONTAINS

         SUBROUTINE read_bin ( fileName , Psi3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Psi3

            OPEN  ( UNIT = 500 , FILE = trim(fileName//'.bin') , ACTION = 'READ' , FORM = 'UNFORMATTED' , STATUS = 'OLD' )

               READ ( UNIT = 500 ) Psi3

            CLOSE ( UNIT = 500 , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE write_bin ( fileName , Psi3 )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Psi3

            OPEN  ( UNIT = 600 , FILE = trim(fileName//'.bin') , ACTION = 'WRITE' , FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )

               WRITE ( UNIT = 600 ) Psi3

            CLOSE ( UNIT = 600 , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

!         SUBROUTINE read_gpi ( fileName , fileNumber , X , Y , Z , Vex , Psi )
!
!            IMPLICIT NONE
!
!            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
!
!            INTEGER, INTENT ( IN ) :: fileNumber
!
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Vex 
!
!            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( IN ) :: Psi 
!
!            CHARACTER ( LEN = 4 ) :: fileNumberStr
!            CHARACTER ( LEN = 1 ) :: spaceStr
!
!            INTEGER :: j , k , l
!
!            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber
!           OPEN  ( UNIT = fileNumber , FILE = trim(fileName//fileNumberStr//'.gpi') , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
!            DO l = 1 , nZ
!
!               DO k = 1 , nY
!
!                  DO j = 1 , nX
!
!                     READ ( UNIT = fileNumber , FMT = '( 1X , 6 ( F23.15 ) )' ) X ( j ) , Y ( k ) , Z ( l ) , Vex ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) , REAL ( Psi ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) ) , AIMAG ( Psi ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) )
!
!                  END DO
!
!               END DO
!               READ ( UNIT = fileNumber , FMT = * ) spaceStr
!
!            END DO
!            CLOSE ( UNIT = fileNumber , STATUS = 'KEEP' )
!
!            RETURN
!
!         END SUBROUTINE 

!         SUBROUTINE write_gpi ( fileName , fileNumber , X , Y , Z , Vex , Psi )
!
!            IMPLICIT NONE
!
!            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
!
!            INTEGER, INTENT ( IN ) :: fileNumber
!
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Vex 
!
!            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( IN ) :: Psi 
!
!            CHARACTER ( LEN = 4 ) :: fileNumberStr
!
!            INTEGER :: j , k , l
!
!            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber
!            OPEN  ( UNIT = fileNumber , FILE = trim(fileName//fileNumberStr//'.gpi') , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
!            DO l = 1 , nZ
!            
!               DO k = 1 , nY
!            
!                  DO j = 1 , nX
!
!                     WRITE ( UNIT = fileNumber , FMT = '( 1X , 6 ( F23.15 ) )' ) X ( j ) , Y ( k ) , Z ( l ) , Vex ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) , REAL ( Psi ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) ) , AIMAG ( Psi ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) )
!
!                  END DO
!
!               END DO
!               WRITE ( UNIT = fileNumber , FMT = * )
!
!            END DO
!            CLOSE ( UNIT = fileNumber , STATUS = 'KEEP' )
!
!            RETURN
!
!         END SUBROUTINE

         SUBROUTINE write_vtk ( fileName , fileNumber , nX , nY , nZ , X , Y , Z , Vex3 , Psi3 )

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

!         SUBROUTINE write_vtk ( fileName , fileNumber , X , Y , Z , Vex , Psi )
!
!            IMPLICIT NONE
!
!            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
!
!            INTEGER, INTENT ( IN ) :: fileNumber
!
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( : ), INTENT ( IN ) :: Vex
!
!            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( IN ) :: Psi
!
!            CHARACTER ( LEN = 4 ) :: fileNumberStr
!
!            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber
!
!            IF ( vtkMode == 'LEGACY' ) THEN
!
!               vtkFileName = fileName//fileNumberStr//'.vtk'
!               vtkError = vtk_ini ( 'ASCII' , trim(vtkFileName) , 'gpse-vex-psi' , 'RECTILINEAR_GRID' )
!               vtkError = vtk_geo ( nX , nY , nZ , X , Y , Z )
!               vtkError = vtk_dat ( nX * nY * nZ , 'NODE' )
!               vtkError = vtk_var ( nX * nY * nZ , 'Vex3'   , Vex )
!               vtkError = vtk_var ( nX * nY * nZ , 'Re(Psi3)' , REAL  ( Psi ) )
!               vtkError = vtk_var ( nX * nY * nZ , 'Im(Psi3)' , AIMAG ( Psi ) )
!               vtkError = vtk_end ( )
!
!            ELSE IF ( vtkMode == 'XML' ) THEN
!
!            ELSE
!
!            END IF
!
!            RETURN
!
!         END SUBROUTINE

      END MODULE

! =========================================================================
