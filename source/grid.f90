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
!     Friday, March 21st, 2014
!
! -------------------------------------------------------------------------

      MODULE GRID

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: regular_grid_axis

      CONTAINS

         SUBROUTINE regular_grid_axis ( nQ , nQa , nQb , qO , dQ , Q )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nQ
         INTEGER, INTENT ( IN ) :: nQa
         INTEGER, INTENT ( IN ) :: nQb

         REAL, INTENT ( IN ) :: qO
         REAL, INTENT ( IN ) :: dQ
 
         REAL, DIMENSION ( : ), INTENT ( INOUT ) :: Q

         INTEGER :: j

         IF ( nQ > 0 ) THEN

            IF ( dQ > 0.0 ) THEN

               IF ( SIZE ( Q ) >= nQ ) THEN ! maximum value of loop index is less than array upper bound. 

                  IF ( MODULO ( nQ , 2 ) == 0 ) THEN ! nQ is even.

                     DO j = nQa , nQb

                        Q ( j ) = qO + ( REAL ( j - nQ / 2 ) - 0.5 ) * dQ

                     END DO

                  ELSE ! nQ is odd. 

                     DO j = nQa , nQb

                        Q ( j ) = qO + REAL ( j - ( nQ + 1 ) / 2 ) * dQ

                     END DO

                  END IF

               ELSE ! range check failed. size of coordinate axis array must be greater than or (preferably) equal to the number of grid points, otherwise, loop index will run out-of-bounds.

               END IF

            ELSE ! range check failed. distance between grid points along coordinate axis must be greater than 0.0.

            END IF

         ELSE ! range check failed. number of grid points along coordinate axis must be greater than 0.

         END IF

         RETURN

         END SUBROUTINE

      END MODULE
! =========================================================================
