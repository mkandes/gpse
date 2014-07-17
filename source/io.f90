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
!     Thursday, July 17th, 2014
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

      INTEGER, PUBLIC :: vtkFileNum = -1
      INTEGER, PUBLIC :: vtkError   = -1

      PUBLIC :: write_vtk

      CONTAINS

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
