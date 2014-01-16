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
!     Thursday, January 16th, 2014
!
! -------------------------------------------------------------------------

      MODULE GRID

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         IMPLICIT NONE
         PRIVATE

         CHARACTER ( LEN = * ), PARAMETER :: VERSION_NUMBER = '0.0.5'

         CONTAINS

            SUBROUTINE regular_grid_axis ( nX , xO , dX , X )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: nX

               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: dX
 
               REAL, DIMENSION ( : ), INTENT ( INOUT ) :: X

               INTEGER :: j

               IF ( nX > 0 ) THEN

                  IF ( dX > 0.0 ) THEN

                     IF ( SIZE ( X ) >= nX ) THEN ! maximum value of loop index is less than array upper bound. 

                        IF ( MODULO ( nX , 2 ) == 0 ) THEN ! nX is even.

                           DO j = 1 , nX

                              X ( j ) = xO + ( REAL ( j - nX / 2 ) - 0.5 ) * dX

                           END DO

                        ELSE ! nX is odd. 

                           DO j = 1 , nX

                              X ( j ) = xO + REAL ( j - ( nX + 1 ) / 2 ) * dX

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
