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
!     Wednesday, March 5th, 2014
!
! -------------------------------------------------------------------------

      MODULE GRID

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         IMPLICIT NONE
         PRIVATE

         INTEGER, PUBLIC :: nX = 0
         INTEGER, PUBLIC :: nY = 0
         INTEGER, PUBLIC :: nZ = 0
          
         REAL, PRIVATE :: xO = 0.0
         REAL, PRIVATE :: yO = 0.0
         REAL, PRIVATE :: zO = 0.0
         REAL, PUBLIC :: dX = 0.0
         REAL, PUBLIC :: dY = 0.0
         REAL, PUBLIC :: dZ = 0.0

         REAL, ALLOCATABLE, DIMENSION ( : ), PUBLIC :: X
         REAL, ALLOCATABLE, DIMENSION ( : ), PUBLIC :: Y
         REAL, ALLOCATABLE, DIMENSION ( : ), PUBLIC :: Z

         CONTAINS

            SUBROUTINE regular_grid_axis ( nQ , qO , dQ , Q )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: nQ

               REAL, INTENT ( IN ) :: qO
               REAL, INTENT ( IN ) :: dQ
 
               REAL, DIMENSION ( : ), INTENT ( INOUT ) :: Q

               INTEGER :: i

               IF ( nQ > 0 ) THEN

                  IF ( dQ > 0.0 ) THEN

                     IF ( SIZE ( Q ) >= nQ ) THEN ! maximum value of loop index is less than array upper bound. 

                        IF ( MODULO ( nQ , 2 ) == 0 ) THEN ! nQ is even.

                           DO i = 1 , nQ

                              Q ( i ) = qO + ( REAL ( i - nQ / 2 ) - 0.5 ) * dQ

                           END DO

                        ELSE ! nQ is odd. 

                           DO i = 1 , nQ

                              Q ( i ) = qO + REAL ( i - ( nQ + 1 ) / 2 ) * dQ

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
