! ==========================================================================
! NAME
!
!     psi [ (p)sī ] - Psi Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     psi is a Fortran module ... 
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
!     Wednesday, April 9th, 2014
!
! -------------------------------------------------------------------------

      MODULE PSI

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MATH

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: psi_3d_se_sho_ani
      PUBLIC :: psi_3d_se_sho_axi
      !PUBLIC :: psi_3d_se_sho_iso

      CONTAINS

         COMPLEX FUNCTION psi_3d_se_sho_ani ( nX , nY , nZ , xO , yO , zO , wX , wY , wZ , x , y , z )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nY
         INTEGER, INTENT ( IN ) :: nZ
          
         REAL, INTENT ( IN ) :: xO
         REAL, INTENT ( IN ) :: yO
         REAL, INTENT ( IN ) :: zO
         REAL, INTENT ( IN ) :: wX
         REAL, INTENT ( IN ) :: wY
         REAL, INTENT ( IN ) :: wZ
         REAL, INTENT ( IN ) :: x
         REAL, INTENT ( IN ) :: y
         REAL, INTENT ( IN ) :: z

         psi_3d_se_sho_ani = CMPLX ( &
            & ( 1.0 / SQRT ( REAL ( 2**nX * factorial ( nX ) ) ) ) * SQRT ( SQRT ( wX / PI ) ) * &
            & ( 1.0 / SQRT ( REAL ( 2**nY * factorial ( nY ) ) ) ) * SQRT ( SQRT ( wY / PI ) ) * &
            & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
            & hermite ( nX , SQRT ( wX ) * ( x - xO ) ) * EXP ( -0.5 * wX * ( x - xO )**2 ) * &
            & hermite ( nY , SQRT ( wY ) * ( y - yO ) ) * EXP ( -0.5 * wY * ( y - yO )**2 ) * &
            & hermite ( nZ , SQRT ( wZ ) * ( z - zO ) ) * EXP ( -0.5 * wZ * ( z - zO )**2 ) , 0.0 )

         RETURN

         END FUNCTION

         COMPLEX FUNCTION psi_3d_se_sho_axi ( nR , mL , nZ , xO , yO , zO , wR , wZ , x , y , z )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nR
         INTEGER, INTENT ( IN ) :: mL
         INTEGER, INTENT ( IN ) :: nZ

         REAL, INTENT ( IN ) :: xO
         REAL, INTENT ( IN ) :: yO
         REAL, INTENT ( IN ) :: zO
         REAL, INTENT ( IN ) :: wR
         REAL, INTENT ( IN ) :: wZ
         REAL, INTENT ( IN ) :: x
         REAL, INTENT ( IN ) :: y
         REAL, INTENT ( IN ) :: z

         psi_3d_se_sho_axi = CMPLX ( &
            & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / &
            & ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * & 
            & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
            & SQRT ( ( x - xO )**2 + ( y - yO )**2 )**ABS ( mL ) * &
            & alaguerre ( nR , ABS ( mL ) , wR * ( ( x - xO )**2 + ( y - yO )**2 ) ) * &
            & EXP ( -0.5 * wR * ( ( x - xO )**2 + ( y - yO )**2 ) ) * &
            & hermite ( nZ , SQRT ( wZ ) * ( z - zO ) ) * EXP ( -0.5 * wZ * ( z - zO )**2 ) , 0.0 ) * &
            & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( y - yO , x - xO ) ) )

         RETURN

         END FUNCTION

      END MODULE

! =========================================================================
