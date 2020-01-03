! ==================================================================================================================================
! NAME
!
!     rot [ rot ] - Rotation Module
!
! SYNOPSIS
!
!     USE :: ROT
!
! DESCRIPTION  
!
!     ROT is a custom Fortran module written to compute the elemental rotation matrices of a Cartesian coordinate system.
!
! OPTIONS
!
! SEE ALSO
!
! BUGS
!
! HISTORY
!
! AUTHOR
!
!     Marty Kandes, Ph.D.
!     Computational & Data Science Research Specialist
!     High-Performance Computing User Services Group
!     San Diego Supercomputer Center
!     University of California, San Diego
!
! COPYRIGHT
!     
!     Copyright (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020 Martin Charles Kandes
!
! LAST UPDATED
!
!     Tuesday, April 2nd, 2019
!
! ----------------------------------------------------------------------------------------------------------------------------------

      MODULE rot

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      PRIVATE

! --- VARIABLE DECLARATIONS --------------------------------------------------------------------------------------------------------

      REAL, PUBLIC :: thetaX
      REAL, PUBLIC :: thetaY
      REAL, PUBLIC :: thetaZ

! --- VARIABLE DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     thetaX is a PUBLIC, REAL-valued variable that stores the rotation angle about the x-axis.
!
!     thetaY is a PUBLIC, REAL-valued variable that stores the rotation angle about the y-axis.
!
!     thetaZ is a PUBLIC, REAL-valued variable that stores the rotation angle about the z-axis.
!      
! --- ARRAY DECLARATIONS -----------------------------------------------------------------------------------------------------------

      REAL, DIMENSION ( 3 , 3 ), PUBLIC :: R

! --- ARRAY DEFINITIONS ------------------------------------------------------------------------------------------------------------
!
!     R is a PUBLIC, REAL-valued rank-two array that stores the 3-by-3 rotation matrix.
!
! --- FUNCTION DECLARATIONS --------------------------------------------------------------------------------------------------------

      PUBLIC :: rot_rx
      PUBLIC :: rot_ry
      PUBLIC :: rot_rz

! --- FUNCTION DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     rot_rx is a PUBLIC FUNCTION that computes the elemental rotation matrix about the x-axis of a Cartesian coordinate system.
!
!     rot_ry is a PUBLIC FUNCTION that computes the elemental rotation matrix about the y-axis of a Cartesian coordinate system.
!
!     rot_rz is a PUBLIC FUNCTION that computes the elemental rotation matrix about the z-axis of a Cartesian coordinate system.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      FUNCTION rot_rx ( theta )

      IMPLICIT NONE

      REAL, INTENT ( IN ) :: theta

      REAL, DIMENSION ( 3 , 3 ) :: rot_rx

      rot_rx ( 1 , 1 ) = 1.0
      rot_rx ( 2 , 1 ) = 0.0
      rot_rx ( 3 , 1 ) = 0.0
      rot_rx ( 1 , 2 ) = 0.0
      rot_rx ( 2 , 2 ) = COS ( theta )
      rot_rx ( 3 , 2 ) = SIN ( theta )
      rot_rx ( 1 , 3 ) = 0.0 
      rot_rx ( 2 , 3 ) = -SIN ( theta )
      rot_rx ( 3 , 3 ) = COS ( theta )

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      FUNCTION rot_ry ( theta )

      IMPLICIT NONE

      REAL, INTENT ( IN ) :: theta

      REAL, DIMENSION ( 3, 3 ) :: rot_ry

      rot_ry ( 1 , 1 ) = COS ( theta ) 
      rot_ry ( 2 , 1 ) = 0.0
      rot_ry ( 3 , 1 ) = -SIN ( theta )
      rot_ry ( 1 , 2 ) = 0.0
      rot_ry ( 2 , 2 ) = 1.0
      rot_ry ( 3 , 2 ) = 0.0
      rot_ry ( 1 , 3 ) = SIN ( theta )
      rot_ry ( 2 , 3 ) = 0.0
      rot_ry ( 3 , 3 ) = COS ( theta )

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      FUNCTION rot_rz ( theta )

      IMPLICIT NONE

      REAL, INTENT ( IN ) :: theta

      REAL, DIMENSION ( 3 , 3 ) :: rot_rz

      rot_rz ( 1 , 1 ) = COS ( theta )
      rot_rz ( 2 , 1 ) = SIN ( theta )
      rot_rz ( 3 , 1 ) = 0.0
      rot_rz ( 1 , 2 ) = -SIN ( theta )
      rot_rz ( 2 , 2 ) = COS ( theta )
      rot_rz ( 3 , 2 ) = 0.0
      rot_rz ( 1 , 3 ) = 0.0
      rot_rz ( 2 , 3 ) = 0.0
      rot_rz ( 3 , 3 ) = 1.0

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      END MODULE

! ==================================================================================================================================
