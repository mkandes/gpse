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
!     Monday, September 8th, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: LIB_VTK_IO
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
      PUBLIC :: write_gpi
      PUBLIC :: write_vtk

      CONTAINS

         SUBROUTINE read_bin ( fileName , fileNumber , Psi )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileNumber

            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( INOUT ) :: Psi

            CHARACTER ( LEN = 4 ) :: fileNumberStr

            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber
            OPEN  ( UNIT = fileNumber , FILE = trim(fileName//fileNumberStr//'.bin') , ACTION = 'READ' , FORM = 'UNFORMATTED' , STATUS = 'OLD' )
            READ ( UNIT = fileNumber ) Psi
            CLOSE ( UNIT = fileNumber , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE write_bin ( fileName , fileNumber , Psi )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileNumber

            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( IN ) :: Psi 

            CHARACTER ( LEN = 4 ) :: fileNumberStr

            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber
            OPEN  ( UNIT = fileNumber , FILE = trim(fileName//fileNumberStr//'.bin') , ACTION = 'WRITE' , FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            WRITE ( UNIT = fileNumber ) Psi
            CLOSE ( UNIT = fileNumber , STATUS = 'KEEP' )

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

         SUBROUTINE write_gpi ( fileName , fileNumber , X , Y , Z , Vex , Psi )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileNumber

            REAL, DIMENSION ( : ), INTENT ( IN ) :: X
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Vex 

            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( IN ) :: Psi 

            CHARACTER ( LEN = 4 ) :: fileNumberStr

            INTEGER :: j , k , l

            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber
            OPEN  ( UNIT = fileNumber , FILE = trim(fileName//fileNumberStr//'.gpi') , ACTION = 'WRITE' , FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
            DO l = 1 , nZ
            
               DO k = 1 , nY
            
                  DO j = 1 , nX

                     WRITE ( UNIT = fileNumber , FMT = '( 1X , 6 ( F23.15 ) )' ) X ( j ) , Y ( k ) , Z ( l ) , Vex ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) , REAL ( Psi ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) ) , AIMAG ( Psi ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) )

                  END DO

               END DO
               WRITE ( UNIT = fileNumber , FMT = * )

            END DO
            CLOSE ( UNIT = fileNumber , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE write_vtk ( fileName , fileNumber , X , Y , Z , Vex , Psi )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

            INTEGER, INTENT ( IN ) :: fileNumber

            REAL, DIMENSION ( : ), INTENT ( IN ) :: X
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( : ), INTENT ( IN ) :: Vex

            COMPLEX, ALLOCATABLE, DIMENSION ( : ), INTENT ( IN ) :: Psi

            CHARACTER ( LEN = 4 ) :: fileNumberStr

            WRITE ( UNIT = fileNumberStr , FMT = '(I4.4)' ) fileNumber

            IF ( vtkMode == 'LEGACY' ) THEN

               vtkFileName = fileName//fileNumberStr//'.vtk'
               vtkError = vtk_ini ( 'ASCII' , trim(vtkFileName) , 'gpse-vex-psi' , 'RECTILINEAR_GRID' )
               vtkError = vtk_geo ( nX , nY , nZ , X , Y , Z )
               vtkError = vtk_dat ( nX * nY * nZ , 'NODE' )
               vtkError = vtk_var ( nX * nY * nZ , 'Vex3'   , Vex )
               vtkError = vtk_var ( nX * nY * nZ , 'Re(Psi3)' , REAL  ( Psi ) )
               vtkError = vtk_var ( nX * nY * nZ , 'Im(Psi3)' , AIMAG ( Psi ) )
               vtkError = vtk_end ( )

            ELSE IF ( vtkMode == 'XML' ) THEN

            ELSE

            END IF

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
