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
!     Thursday, October 16th, 2014
!
! -------------------------------------------------------------------------

      MODULE GRID

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: grid_bound_cond_size
      PUBLIC :: grid_regular
      PUBLIC :: grid_regular_axis

      CONTAINS

         SUBROUTINE grid_bound_cond_size ( fdOrder , nXbc , nYbc , nZbc )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: fdOrder
            INTEGER, INTENT ( INOUT ) :: nXbc 
            INTEGER, INTENT ( INOUT ) :: nYbc 
            INTEGER, INTENT ( INOUT ) :: nZbc 

            IF ( fdOrder == 2 ) THEN

               nXbc = 1
               nYbc = 1
               nZbc = 1

            ELSE IF ( fdOrder == 4 ) THEN

               nXbc = 2
               nYbc = 2
               nZbc = 2

            ELSE IF ( fdOrder == 6 ) THEN

               nXbc = 3
               nYbc = 3
               nZbc = 3

            ELSE IF ( fdOrder == 8 ) THEN

               nXbc = 4
               nYbc = 4
               nZbc = 4

            ELSE

               WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : grid : grid_bound_cond_size : ERROR - fdOrder not supported.'

            END IF

            RETURN

         END SUBROUTINE         

         SUBROUTINE grid_regular ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , X , Y , Z )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN ) :: nXa
            INTEGER, INTENT ( IN ) :: nXb
            INTEGER, INTENT ( IN ) :: nXbc
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN ) :: nYa
            INTEGER, INTENT ( IN ) :: nYb
            INTEGER, INTENT ( IN ) :: nYbc
            INTEGER, INTENT ( IN ) :: nZ
            INTEGER, INTENT ( IN ) :: nZa
            INTEGER, INTENT ( IN ) :: nZb
            INTEGER, INTENT ( IN ) :: nZbc 

            REAL, INTENT ( IN ) :: xO
            REAL, INTENT ( IN ) :: yO
            REAL, INTENT ( IN ) :: zO
            REAL, INTENT ( IN ) :: dX
            REAL, INTENT ( IN ) :: dY
            REAL, INTENT ( IN ) :: dZ

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( INOUT ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( INOUT ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Z

            CALL grid_regular_axis ( nX , nXa , nXb , nXbc , xO , dX , X )
            CALL grid_regular_axis ( nY , nYa , nYb , nYbc , yO , dY , Y )
            CALL grid_regular_axis ( nZ , nZa , nZb , nZbc , zO , dZ , Z )

            RETURN

         END SUBROUTINE         

         SUBROUTINE grid_regular_axis ( nQ , nQa , nQb , nQbc , qO , dQ , Q )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: nQ
            INTEGER, INTENT ( IN ) :: nQa
            INTEGER, INTENT ( IN ) :: nQb
            INTEGER, INTENT ( IN ) :: nQbc

            REAL, INTENT ( IN ) :: qO
            REAL, INTENT ( IN ) :: dQ

            REAL, DIMENSION ( nQa - nQbc : nQb + nQbc ), INTENT ( INOUT ) :: Q

            INTEGER :: j

            IF ( MODULO ( nQ , 2 ) == 0 ) THEN ! nQ is even

               DO j = nQa - nQbc , nQb + nQbc

                  Q ( j ) = qO + ( REAL ( j - nQ / 2 ) - 0.5 ) * dQ

               END DO

            ELSE ! nQ is odd 

               DO j = nQa - nQbc , nQb + nQbc

                  Q ( j ) = qO + REAL ( j - ( nQ + 1 ) / 2 ) * dQ

               END DO

            END IF

            RETURN

         END SUBROUTINE

      END MODULE
! =========================================================================
