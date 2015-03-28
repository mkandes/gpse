! ==========================================================================
! NAME
!
!     rot [ rot ] - rot module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     rot is a Fortran module ... 
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
!     Copyright (c) 2014, 2015 Martin Charles Kandes
!
! LAST UPDATED
!
!     Monday, March 23rd, 2015
!
! -------------------------------------------------------------------------

      MODULE rot

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      REAL, PUBLIC :: thetaX
      REAL, PUBLIC :: thetaY
      REAL, PUBLIC :: thetaZ

      REAL, DIMENSION ( 3 , 3 ), PUBLIC :: R

      PUBLIC :: rot_rx
      PUBLIC :: rot_ry
      PUBLIC :: rot_rz

      CONTAINS

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

      END MODULE

! =========================================================================
