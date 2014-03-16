! ==========================================================================
! NAME
!
!     vex [ veks ] - VEX Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     VEX is a Fortran module ... 
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
!     Wednesday, March 5th, 2014
!
! -------------------------------------------------------------------------

      MODULE VEX

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         IMPLICIT NONE
         PRIVATE

         PUBLIC :: vex_3d_lin
         PUBLIC :: vex_3d_sho
         PUBLIC :: vex_3d_shor

         CONTAINS

            REAL FUNCTION vex_3d_lin ( x , y , z , xO , yO , zO , fX , fY , fZ )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: x
               REAL, INTENT ( IN ) :: y
               REAL, INTENT ( IN ) :: z
               REAL, INTENT ( IN ) :: fX
               REAL, INTENT ( IN ) :: fY
               REAL, INTENT ( IN ) :: fZ
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO

               vex_3d_lin = fX * ( x - xO ) + fY * ( y - yO ) + fZ * ( z - zO ) 

               RETURN

            END FUNCTION

            REAL FUNCTION vex_3d_sho ( x , y , z , xO , yO , zO , wX , wY , wZ )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: x 
               REAL, INTENT ( IN ) :: y
               REAL, INTENT ( IN ) :: z
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: wX
               REAL, INTENT ( IN ) :: wY
               REAL, INTENT ( IN ) :: wZ
               
               vex_3d_sho = 0.5 * ( ( wX * ( x - xO ) )**2 + ( wY * ( y - yO ) )**2 + ( wZ * ( z - zO ) )**2 )

               RETURN

            END FUNCTION

            REAL FUNCTION vex_3d_shor ( x , y , z , xO , yO , zO , rO , wR , wZ )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: x
               REAL, INTENT ( IN ) :: y
               REAL, INTENT ( IN ) :: z
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: rO
               REAL, INTENT ( IN ) :: wR
               REAL, INTENT ( IN ) :: wZ

               vex_3d_shor = 0.5 * ( wR * ( SQRT ( ( x - xO )**2 + ( y - yO )**2 ) - rO )**2 + ( wZ * ( z - zO ) )**2 )

               RETURN

            END FUNCTION

      END MODULE

! =========================================================================
