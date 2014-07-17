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
!     Wednesday, July 16th, 2014
!
! -------------------------------------------------------------------------

      MODULE GRID

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      INTEGER, PARAMETER, PRIVATE :: unitGridIn = 501

      INTEGER, PUBLIC :: nX   = 0 
      INTEGER, PUBLIC :: nXa  = 0 
      INTEGER, PUBLIC :: nXb  = 0 
      INTEGER, PUBLIC :: nXbc = 0 
      INTEGER, PUBLIC :: nY   = 0
      INTEGER, PUBLIC :: nYa  = 0 
      INTEGER, PUBLIC :: nYb  = 0 
      INTEGER, PUBLIC :: nYbc = 0 
      INTEGER, PUBLIC :: nZ   = 0 
      INTEGER, PUBLIC :: nZa  = 0 
      INTEGER, PUBLIC :: nZb  = 0 
      INTEGER, PUBLIC :: nZbc = 0

      REAL, PUBLIC :: xO = 0.0
      REAL, PUBLIC :: yO = 0.0
      REAL, PUBLIC :: zO = 0.0
      REAL, PUBLIC :: dX = 0.0
      REAL, PUBLIC :: dY = 0.0
      REAL, PUBLIC :: dZ = 0.0 

      PUBLIC :: grid_read_inputs
      PUBLIC :: grid_bound_cond_size
      PUBLIC :: grid_regular
      PUBLIC :: grid_regular_axis

      NAMELIST /nmlGridIn/ nX , nY , nZ , xO , yO , zO , dX , dY , dZ

      CONTAINS

         SUBROUTINE grid_read_inputs ( )

            IMPLICIT NONE

            OPEN ( UNIT = unitGridIn, FILE = 'grid.in' , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
               READ ( UNIT = unitGridIn , NML = nmlGridIn )
            CLOSE ( UNIT = unitGridIn , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE grid_bound_cond_size ( fdOrder )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: fdOrder

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

               ! Error: fdOrder not supported.

            END IF

            RETURN

         END SUBROUTINE         

         SUBROUTINE grid_regular ( X , Y , Z )

            IMPLICIT NONE

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
