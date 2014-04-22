! ==========================================================================
! NAME
!
!     grid [ grid ] - Grid Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     grid is a Fortran module ... 
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
!     Monday, April 21st, 2014
!
! -------------------------------------------------------------------------

      MODULE GRID

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: regular_grid_axis

      CONTAINS

         SUBROUTINE regular_grid_axis ( nQ , nQa , nQb , nQbc , qO , dQ , Q )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: nQ
            INTEGER, INTENT ( IN ) :: nQa
            INTEGER, INTENT ( IN ) :: nQb
            INTEGER, INTENT ( IN ) :: nQbc

            REAL, INTENT ( IN ) :: qO
            REAL, INTENT ( IN ) :: dQ
 
            REAL, DIMENSION ( nQa - nQbc : nQb + nQbc ), INTENT ( INOUT ) :: Q

            INTEGER :: j

            IF ( MODULO ( nQ , 2 ) == 0 ) THEN ! nQ is even.

               DO j = nQa - nQbc , nQb + nQbc

                  Q ( j ) = qO + ( REAL ( j - nQ / 2 ) - 0.5 ) * dQ

               END DO

            ELSE ! nQ is odd. 

               DO j = nQa - nQbc , nQb + nQbc

                  Q ( j ) = qO + REAL ( j - ( nQ + 1 ) / 2 ) * dQ

               END DO

            END IF

            RETURN

         END SUBROUTINE

      END MODULE
! =========================================================================
