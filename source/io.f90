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
!     Thursday, March 13th, 2014
!
! -------------------------------------------------------------------------

      MODULE IO

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: LIB_VTK_IO
      USE            :: GRID

      IMPLICIT NONE
      PRIVATE

      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_SIMAG = '.vti'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_SPOLY = '.vtp'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_SRECT = '.vtr'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_SSTRUCT = '.vts'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_SUNSTRUCT = '.vtu'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_PIMAG = '.pvti'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_PPOLY = '.pvtp'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_PRECT = '.pvtr'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_PSTRUCT = '.pvts'
      CHARACTER ( LEN = * ), PARAMETER :: VTK_XML_PUNSTRUCT = '.pvtu'

      INTEGER, PARAMETER :: MAX_LEN = 80 

      CHARACTER ( LEN = MAX_LEN ), PUBLIC  :: vtkFilename
      CHARACTER ( LEN = MAX_LEN ), PUBLIC  :: vtkOutputFmt = 'ASCII'
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkMeshTopology = 'RECTILINEAR_GRID'
      CHARACTER ( LEN = MAX_LEN ), PUBLIC  :: vtkMode = 'LEGACY'
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkTitle
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkVariableBlockAction = 'CLOSE'
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkVariableLocation
      CHARACTER ( LEN = MAX_LEN ), PRIVATE :: vtkVariableName

      INTEGER, PUBLIC :: vtkError = 0

      PUBLIC :: write_vtk

      CONTAINS

         SUBROUTINE write_vtk ( t , Vex , RePsi , ImPsi )

         IMPLICIT NONE

         REAL, INTENT ( IN ) :: t
         
         REAL, DIMENSION ( : ), INTENT ( IN ) :: Vex
         REAL, DIMENSION ( : ), INTENT ( IN ) :: RePsi
         REAL, DIMENSION ( : ), INTENT ( IN ) :: ImPsi

         IF ( vtkMode == 'LEGACY' ) THEN

            vtkError = vtk_ini ( 'ASCII' , 'TEST.vtk' , 'TEST' , 'RECTILINEAR_GRID' )
            vtkError = vtk_geo ( nX , nY , nZ , X , Y , Z )
            vtkError = vtk_dat ( nX * nY * nZ , 'NODE' )
            vtkError = vtk_var ( nX * nY * nZ , 'Vex' , Vex )
            vtkError = vtk_var ( nX * nY * nZ , 'RePsi' , RePsi )
            vtkError = vtk_var ( nX * nY * nZ , 'ImPsi' , ImPsi )
            vtkError = vtk_end ( )

         ELSE IF ( vtkMode == 'XML' ) THEN

         ELSE

         END IF

         RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
