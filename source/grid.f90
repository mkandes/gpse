! ======================================================================
! NAME
!
!     grid [ grid ] - Grid Module
!
! SYNOPSIS
!
!     USE :: GRID
!
! DESCRIPTION  
!
!     GRID is a custom Fortran module written to build the regular 
!     Cartesian grid of points that defines the spatial computational 
!     domain on which external potentials, wave functions and other 
!     physical quantities may be computed. It also provides some utility 
!     functions and subroutines to determine information about any 
!     grid it is used to define. For example, see the recently added
!     grid_nearest_point function, which is used to find the 
!     index of the nearest grid point given a spatial coordinate. The 
!     module only has dependency on the ISO_FORTRAN_ENV module.
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
!     Thursday, August 15th, 2019
!
! ----------------------------------------------------------------------

      MODULE GRID

! --- MODULE DECLARATIONS ----------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV

! --- MODULE DEFINITIONS -----------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides 
!        information about the run-time environment.
!
! ----------------------------------------------------------------------

      IMPLICIT NONE
      PRIVATE

! --- SUBROUTINE DECLARATIONS ------------------------------------------

      PUBLIC :: grid_boundary_condition_size
      PUBLIC :: grid_regular
      PUBLIC :: grid_regular_axis
      PUBLIC :: grid_nearest_point

! --- SUBROUTINE DEFINITIONS -------------------------------------------
!
!     grid_boundary_condition_size is a PUBLIC SUBROUTINE that sets the
!        width of the grid boundaries (in number of grid points) by 
!        that required for the order-of-accuracy of the finite 
!        difference operators in use. 
!
!     grid_regular is a PUBLIC SUBROUTINE that builds a 
!        three-dimensional, regular Cartesian grid.
!
!     grid_regular_axis is a PUBLIC SUBROUTINE that sets the position 
!        of each grid point along a linear coordinate axis with a fixed
!        distance between each pair of points.  
!
!     grid_nearest_point is a PUBLIC SUBROUTINE that returns the 
!        indicies of the nearest grid point on the grid to a given 
!        coordinate position.
!
! --- FUNCTION DECLARATIONS --------------------------------------------

      PUBLIC :: grid_linear_search
      PRIVATE :: grid_binary_search

! --- FUNCTION DEFINITIONS ---------------------------------------------
!
!     grid_linear_search is a PRIVATE FUNCTION that performs a linear 
!         search for the index of the nearest grid point on a linear 
!         coordinate axis given a coordinate position.
!
!     grid_binary_search is a PRIVATE FUNCTION that performs a binary
!         search for the index of the nearest grid point on a linear 
!         coordinate axis given a coordinate position. NOTE: Not yet 
!         implemented.
!
! ----------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------

      SUBROUTINE grid_boundary_condition_size(fdOrder, nXbc, nYbc, nZbc)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: fdOrder
      INTEGER, INTENT(INOUT) :: nXbc, nYbc, nZbc 

      IF (fdOrder == 2) THEN

         nXbc = 1
         nYbc = 1
         nZbc = 1

      ELSE IF (fdOrder == 4) THEN

         nXbc = 2
         nYbc = 2
         nZbc = 2

      ELSE IF (fdOrder == 6) THEN

         nXbc = 3
         nYbc = 3
         nZbc = 3

      ELSE IF (fdOrder == 8) THEN

         nXbc = 4
         nYbc = 4
         nZbc = 4

      ELSE

         WRITE(UNIT=ERROR_UNIT, FMT=*) 'gpse: &
            & grid_boundary_condition_size: ERROR - fdOrder not &
            & supported.'
         STOP

      END IF

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE grid_regular(nX, nXa, nXb, nXbc, nY, nYa, nYb, nYbc, &
         & nZ, nZa, nZb, nZbc, xO, yO, zO, dX, dY, dZ, X, Y, Z)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nX, nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nY, nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZ, nZa, nZb, nZbc 

      REAL, INTENT(IN) :: xO, yO, zO
      REAL, INTENT(IN) :: dX, dY, dZ

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(INOUT) :: X
      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(INOUT) :: Y
      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: Z

      CALL grid_regular_axis(nX, nXa, nXb, nXbc, xO, dX, X)
      CALL grid_regular_axis(nY, nYa, nYb, nYbc, yO, dY, Y)
      CALL grid_regular_axis(nZ, nZa, nZb, nZbc, zO, dZ, Z)

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE grid_regular_axis(nQ, nQa, nQb, nQbc, qO, dQ, Q)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nQ, nQa, nQb, nQbc

      REAL, INTENT(IN) :: qO, dQ

      REAL, DIMENSION(nQa - nQbc : nQb + nQbc), INTENT(INOUT) :: Q

      INTEGER :: j

      IF (MODULO(nQ, 2) == 0) THEN ! nQ is even
         DO j = nQa - nQbc, nQb + nQbc
            Q(j) = qO + (REAL(j - nQ / 2) - 0.5) * dQ
         END DO
      ELSE ! nQ is odd 
         DO j = nQa - nQbc, nQb + nQbc
            Q(j) = qO + REAL(j - (nQ + 1) / 2) * dQ
         END DO
      END IF

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE grid_nearest_point(nXa, nXb, nXbc, nYa, nYb, nYbc, &
         & nZa, nZb, nZbc, nXo, nYo, nZo, xO, yO, zO, X, Y, Z)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc
      INTEGER, INTENT(INOUT) :: nXo, nYo, nZo

      REAL, INTENT(IN) :: xO, yO, zO

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(IN) :: X
      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(IN) :: Y
      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(IN) :: Z

      nXo = grid_linear_search(nXa, nXb, nXbc, xO, X)
      nYo = grid_linear_search(nYa, nYb, nYbc, yO, Y)
      nZo = grid_linear_search(nZa, nZb, nZbc, zO, Z)

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      INTEGER FUNCTION grid_linear_search(nQa, nQb, nQbc, qO, Q) &
         & RESULT(nQo)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nQa, nQb, nQbc

      REAL, INTENT(IN) :: qO

      REAL, DIMENSION(nQa - nQbc : nQb + nQbc), INTENT(IN) :: Q

      INTEGER :: j

      nQo = 0 ! must be zero for MPI_ALLREDUCE scheme used in pmca_current_3d_rect to work

      DO j = nQa, nQb
         IF ((qO >= Q(j)).AND.(qO < Q(j+1))) THEN
            ! add code to determine if q0 is closer to Q(j) or Q(j+1)
            nQo = j
         END IF
      END DO

      RETURN
      END FUNCTION

! ----------------------------------------------------------------------

      INTEGER RECURSIVE FUNCTION grid_binary_search()
      IMPLICIT NONE

      RETURN
      END FUNCTION

! ----------------------------------------------------------------------

      END MODULE

! ======================================================================
