! ======================================================================
! NAME
!
!     pmca [ pmca ] - Probability and Mass Current Analysis Module
!
! SYNOPSIS
!
!     USE :: PMCA
!
! DESCRIPTION  
!
!     PMCA is a custom Fortran module written to compute the probability 
!     and mass currents from a single-particle and/or single-component 
!     Bose-Einstein condensate wave function.
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
!     Copyright (c) 2014, 2015, 2016, 2017, 2018, 2019 Martin Charles Kandes
!
! LAST UPDATED
!
!     Thursday, July 4th, 2019
!
! ----------------------------------------------------------------------

      MODULE pmca

! --- MODULE DECLARATIONS ----------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE :: GRID
      USE :: MPI        

! --- MODULE DEFINITIONS -----------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides 
!        information about the run-time environment.
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

!      INCLUDE 'mpif.h' ! if MPI module is not supported, use header file.

      PRIVATE

! --- VARIABLE DECLARATIONS --------------------------------------------

      REAL, PRIVATE :: pmcaV12 = 0.0
      REAL, PRIVATE :: pmcaV34 = 0.0
      REAL, PRIVATE :: pmcaV41 = 0.0
      REAL, PRIVATE :: pmcaV23 = 0.0

      REAL, PRIVATE :: pmcaI12 = 0.0
      REAL, PRIVATE :: pmcaI34 = 0.0
      REAL, PRIVATE :: pmcaI41 = 0.0
      REAL, PRIVATE :: pmcaI23 = 0.0

      REAL, PRIVATE :: pmcaI1 = 0.0
      REAL, PRIVATE :: pmcaI2 = 0.0
      REAL, PRIVATE :: pmcaI3 = 0.0
      REAL, PRIVATE :: pmcaI4 = 0.0

! --- VARIABLE DEFINITIONS ---------------------------------------------
!      
! --- SUBROUTINE DECLARATIONS ------------------------------------------

      PUBLIC :: pmca_compute_velocities
      PUBLIC :: pmca_write_velocities

      PUBLIC :: pmca_compute_currents
      PUBLIC :: pmca_write_currents

      PUBLIC :: pmca_density
      PUBLIC :: pmca_phase
      PUBLIC :: pmca_velocity
      PUBLIC :: pmca_current_density
      PUBLIC :: pmca_quantum_potential
      PUBLIC :: pmca_vorticity

      PRIVATE :: pmca_grad_f_real_3d_rect_cd2
      PRIVATE :: pmca_div_f_cmplx_3d_rect_cd2
      PRIVATE :: pmca_curl_f_cmplx_3d_rect_cd2
      PRIVATE :: pmca_laplacian_f_cmplx_3d_rect_cd2

! --- SUBROUTINE DEFINITIONS -------------------------------------------
!
! --- FUNCTION DECLARATIONS --------------------------------------------

      PRIVATE :: pmca_current_3d_rect

! --- FUNCTION DEFINITIONS ---------------------------------------------
!
! ----------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------

      SUBROUTINE pmca_compute_velocities()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_write_velocities()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_compute_currents(mpiMaster, mpiReal, mpiError, &
         & pmcaQuadRule, nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, nZb, &
         & nZbc, dX, dY, dZ, X, Y, Z, J3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: mpiMaster, mpiReal
      INTEGER, INTENT(INOUT) :: mpiError
      INTEGER, INTENT(IN) :: pmcaQuadRule
      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc

      REAL, INTENT (IN) :: dX, dY, dZ

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(IN) :: X
      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(IN) :: Y
      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(IN) :: Z

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc) :: J3

      REAL :: temp

      IF (pmcaQuadRule == 1) THEN

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, 0.0, 0.0, 0.0, Y(nYb), Z(nZa), Z(nZb), &
            & dX, dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI12, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, 0.0, 0.0, Y(nYa), 0.0, Z(nZa), Z(nZb), &
            & dX, dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI34, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, 0.0, X(nXb), 0.0, 0.0, Z(nZa), Z(nZb), &
            & dX, dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI41, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, X(nXa), 0.0, 0.0, 0.0, Z(nZa), Z(nZb), &
            & dX, dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI23, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, 0.0, X(nXb), 0.0, Y(nYb), 0.0, 0.0, dX,&
            & dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI1, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, X(nXa), 0.0, 0.0, Y(nYb), 0.0, 0.0, dX,&
            & dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI2, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, X(nXa), 0.0, Y(nYa), 0.0, 0.0, 0.0, dX,&
            & dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI3, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         temp = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, 0.0, X(nXb), Y(nYa), 0.0, 0.0, 0.0, dX,&
            & dY, dZ, X, Y, Z, J3)
         CALL MPI_REDUCE(temp, pmcaI4, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

      ELSE

         WRITE(UNIT=ERROR_UNIT, FMT = *) 'gpse: pmca_compute_currents: &
            & ERROR - pmcaQuadRule is not supported.'
         STOP

      END IF

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_write_currents()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_density(nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, nZb,&
         & nZbc, Psi3, Rho3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa
      INTEGER, INTENT(IN) :: nXb
      INTEGER, INTENT(IN) :: nXbc
      INTEGER, INTENT(IN) :: nYa
      INTEGER, INTENT(IN) :: nYb
      INTEGER, INTENT(IN) :: nYbc
      INTEGER, INTENT(IN) :: nZa
      INTEGER, INTENT(IN) :: nZb
      INTEGER, INTENT(IN) :: nZbc

      COMPLEX, DIMENSION(nXa - nXbc : nXb + nXbc, &
                       & nYa - nYbc : nYb + nYbc, &
                       & nZa - nZbc : nZb + nZbc), INTENT(IN) :: Psi3

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: Rho3

      INTEGER :: j, k, l

      Rho3 = 0.0

!$OMP PARALLEL DO IF (nZa /= nZb) DEFAULT(SHARED) SCHEDULE(STATIC)
      DO l = nZa, nZb
!$OMP    PARALLEL DO IF (nZa == nZb) DEFAULT(SHARED) SCHEDULE(STATIC)
         DO k = nYa, nYb
            DO j = nXa, nXb
               Rho3(j,k,l) = ABS(Psi3(j,k,l))**2
            END DO
         END DO
!$OMP    END PARALLEL DO
      END DO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_phase(nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, nZb,&
         & nZbc, Psi3, Phi3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa
      INTEGER, INTENT(IN) :: nXb
      INTEGER, INTENT(IN) :: nXbc
      INTEGER, INTENT(IN) :: nYa
      INTEGER, INTENT(IN) :: nYb
      INTEGER, INTENT(IN) :: nYbc
      INTEGER, INTENT(IN) :: nZa
      INTEGER, INTENT(IN) :: nZb
      INTEGER, INTENT(IN) :: nZbc

      COMPLEX, DIMENSION(nXa - nXbc : nXb + nXbc, &
                       & nYa - nYbc : nYb + nYbc, &
                       & nZa - nZbc : nZb + nZbc), INTENT(IN) :: Psi3

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: Phi3

      INTEGER :: j, k, l

      Phi3 = 0.0

!$OMP PARALLEL DO IF (nZa /= nZb) DEFAULT(SHARED) SCHEDULE(STATIC)
      DO l = nZa, nZb
!$OMP    PARALLEL DO IF (nZa == nZb) DEFAULT(SHARED) SCHEDULE(STATIC)
         DO k = nYa, nYb
            DO j = nXa, nXb
               Phi3(j,k,l) = ATAN2(AIMAG(Psi3(j,k,l)), REAL(Psi3(j,k,l)))
            END DO
         END DO
!$OMP    END PARALLEL DO
      END DO
!$OMP END PARALLEL DO

      RETURN
      END 

! ----------------------------------------------------------------------

      SUBROUTINE pmca_velocity(nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, &
         & nZb, nZbc, dX, dY, dZ, Phi3, V3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa
      INTEGER, INTENT(IN) :: nXb
      INTEGER, INTENT(IN) :: nXbc
      INTEGER, INTENT(IN) :: nYa
      INTEGER, INTENT(IN) :: nYb
      INTEGER, INTENT(IN) :: nYbc
      INTEGER, INTENT(IN) :: nZa
      INTEGER, INTENT(IN) :: nZb
      INTEGER, INTENT(IN) :: nZbc

      REAL, INTENT(IN) :: dX
      REAL, INTENT(IN) :: dY
      REAL, INTENT(IN) :: dZ

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: Phi3

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: V3

      CALL pmca_grad_f_real_3d_rect_cd2(nXa, nXb, nXbc, nYa, nYb, &
         & nYbc, nZa, nZb, nZbc, dX, dY, dZ, Phi3, V3)

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_current_density(nXa, nXb, nXbc, nYa, nYb, nYbc, &
         & nZa, nZb, nZbc, Rho3, V3, J3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa
      INTEGER, INTENT(IN) :: nXb
      INTEGER, INTENT(IN) :: nXbc
      INTEGER, INTENT(IN) :: nYa
      INTEGER, INTENT(IN) :: nYb
      INTEGER, INTENT(IN) :: nYbc
      INTEGER, INTENT(IN) :: nZa
      INTEGER, INTENT(IN) :: nZb
      INTEGER, INTENT(IN) :: nZbc

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: Rho3

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: V3

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: J3

      INTEGER :: j, k, l

      J3 = 0.0

!$OMP PARALLEL DO IF (nZa /= nZb) DEFAULT(SHARED) SCHEDULE(STATIC)
      DO l = nZa , nZb
!$OMP    PARALLEL DO IF (nZa == nZb) DEFAULT(SHARED) SCHEDULE(STATIC)
         DO k = nYa , nYb
            DO j = nXa , nXb
               J3(1,j,k,l) = Rho3(j,k,l) * V3(1,j,k,l)
               J3(2,j,k,l) = Rho3(j,k,l) * V3(2,j,k,l)
               J3(3,j,k,l) = Rho3(j,k,l) * V3(3,j,k,l)
            END DO
         END DO
!$OMP    END PARALLEL DO
      END DO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_quantum_potential()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_vorticity()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_grad_f_real_3d_rect_cd2(nXa, nXb, nXbc, nYa, &
         & nYb, nYbc, nZa, nZb, nZbc, dX, dY, dZ, F3, GradF3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc

      REAL, INTENT(IN) :: dX, dY, dZ

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: F3

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: GradF3

      INTEGER :: j, k, l

      GradF3 = 0.0

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb
!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb
            DO j = nXa , nXb
               GradF3(1,j,k,l) = CMPLX(0.5/dX, 0.0) * (F3(j+1,k,l) - F3(j-1,k,l))
               GradF3(2,j,k,l) = CMPLX(0.5/dY, 0.0) * (F3(j,k+1,l) - F3(j,k-1,l))
               GradF3(3,j,k,l) = CMPLX(0.5/dZ, 0.0) * (F3(j,k,l+1) - F3(j,k,l-1))
            END DO
         END DO
!$OMP    END PARALLEL DO
      END DO
!$OMP END PARALLEL DO 

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_div_f_cmplx_3d_rect_cd2()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_curl_f_cmplx_3d_rect_cd2()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_laplacian_f_cmplx_3d_rect_cd2()
      IMPLICIT NONE

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      REAL FUNCTION pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, &
         & nYbc, nZa, nZb, nZbc, xI, xF, yI, yF, zI, zF, dX, dY, dZ, &
         & X, Y, Z, J3) RESULT(pmcaI)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc

      REAL, INTENT(IN) :: xI, xF
      REAL, INTENT(IN) :: yI, yF
      REAL, INTENT(IN) :: zI, zF

      REAL, INTENT(IN) :: dX, dY, dZ

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(IN) :: X
      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(IN) :: Y
      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(IN) :: Z

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: J3

      INTEGER :: j, k, l

      INTEGER :: nXi = -1, nXf = -1
      INTEGER :: nYi = -1, nYf = -1
      INTEGER :: nZi = -1, nZf = -1

      pmcaI = 0.0

      CALL grid_nearest_point(nXa, nXb, nXbc, nYa, nYb, nYbc, &
         & nZa, nZb, nZbc, nXi, nYi, nZi, xI, yI, zI, X, Y, Z)
      CALL grid_nearest_point(nXa, nXb, nXbc, nYa, nYb, nYbc, &
         & nZa, nZb, nZbc, nXf, nYf, nZf, xF, yF, zF, X, Y, Z)

      RETURN
      END FUNCTION

! ----------------------------------------------------------------------
!
!      REAL FUNCTION pmca_current_3d_rect(mpiRank, nXa, nXb, nXbc, nYa, &
!         & nYb, nYbc, nZa, nZb, nZbc, xI, xF, yI, yF, zI, zF, dX, dY, &
!         & dZ, X, Y, Z, J3) RESULT(pmcaI)
!      IMPLICIT NONE
!
!      INTEGER, INTENT(IN) :: mpiRank
!      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
!      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
!      INTEGER, INTENT(IN) :: nZa, nZb, nZbc
!
!      REAL, INTENT(IN) :: xI, xF
!      REAL, INTENT(IN) :: yI, yF
!      REAL, INTENT(IN) :: zI, zF
!
!      REAL, INTENT(IN) :: dX, dY, dZ
!
!      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(IN) :: X
!      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(IN) :: Y
!      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(IN) :: Z
!
!      REAL, DIMENSION(3, &
!                    & nXa - nXbc : nXb + nXbc, &
!                    & nYa - nYbc : nYb + nYbc, &
!                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: J3
!
!      INTEGER :: j, k, l
!
!      INTEGER :: nXi = -1, nXf = -1
!      INTEGER :: nYi = -1, nYf = -1
!      INTEGER :: nZi = -1, nZf = -1
!
!      pmcaI = 0.0
!
!      IF ((xI == xF).AND.(yI < yF).AND.(zI < zF)) THEN
!       
!         DO j = nXa, nXb
!            IF ((xI >= X(j-1)).AND.(xI <= X(j))) THEN
!               nXi = j
!               nXf = nXi
!            END IF 
!         END DO  
!
!         DO k = nYa, nYb
!            IF ((yI >= Y(k-1)).AND.(yI <= Y(k))) THEN
!               nYi = k
!            END IF
!         END DO
!
!         DO k = nYa, nYb
!            IF ((yF >= Y(k-1)).AND.(yF <= Y(k))) THEN
!               nYf = k
!            END IF
!         END DO
!
!         IF ((zI <= Z(nZa)).AND.(zF >= Z(nZb))) THEN
!
!            nZi = nZa
!            nZf = nZb
!
!         ELSE IF ((zI <= Z(nZa)) .AND. &
!               & ((zF >= Z(nZa)).AND.(zF <= Z(nZb)))) THEN
!
!            nZi = nZa
!
!            DO l = nZa, nZb
!               IF ((zF >= Z(l-1)).AND.(zF <= Z(l))) THEN
!                  nZf = l
!               END IF
!            END DO
!
!         ELSE IF (((zI >= Z(nZa)).AND.(zI <= Z(nZb))) .AND. &
!                 & (zF >= Z(nZb))) THEN
!
!            DO l = nZa, nZb
!               IF ((zI >= Z(l-1)).AND.(zI <= Z(l))) THEN
!                  nZi = l
!               END IF
!            END DO
!
!            nZf = nZb
!
!         ELSE IF (((zI >= Z(nZa)).AND.(zI <= Z(nZb))) .AND. &
!           & ((zF >= Z(nZa)).AND.(zF <= Z(nZb)))) THEN
!
!            DO l = nZa, nZb
!               IF ((zI >= Z(l-1)).AND.(zI <= Z(l))) THEN
!                  nZi = l
!               END IF
!            END DO
!
!            DO l = nZa, nZb
!               IF ((zF >= Z(l-1)).AND.(zF <= Z(l))) THEN
!                  nZf = l
!               END IF
!            END DO
!
!         ELSE
!
!            ! z-axis integration limits do not contain this slab; do nothing 
!
!         END IF
!
!      ELSE IF ((xI < xF).AND.(yI == yF).AND.(zI < zF)) THEN
!
!         DO j = nXa, nXb
!            IF ((xI >= X(j-1)).AND.(xI <= X(j))) THEN
!               nXi = j
!            END IF
!         END DO
!
!         DO j = nXa, nXb
!            IF ((xF >= X(j-1)).AND.(xF <= X(j))) THEN
!               nXf = j
!            END IF
!         END DO
!
!         DO k = nYa, nYb
!            IF ((yI >= Y(k-1)).AND.(yI <= Y(k))) THEN
!               nYi = k
!               nYf = nYi
!            END IF
!         END DO
!
!         IF ((zI <= Z(nZa)).AND.(zF >= Z(nZb))) THEN
!
!            nZi = nZa
!            nZf = nZb
!
!         ELSE IF ((zI <= Z(nZa)) .AND. &
!               & ((zF >= Z(nZa)).AND.(zF <= Z(nZb)))) THEN
!
!            nZi = nZa
!
!            DO l = nZa, nZb
!               IF ((zF >= Z(l-1)).AND.(zF <= Z(l))) THEN
!                  nZf = l
!               END IF
!            END DO
!
!         ELSE IF (((zI >= Z(nZa)).AND.(zI <= Z(nZb))) .AND. &
!                 & (zF >= Z(nZb))) THEN
!
!            DO l = nZa, nZb
!               IF ((zI >= Z(l-1)).AND.(zI <= Z(l))) THEN
!                  nZi = l
!               END IF
!            END DO
!
!            nZf = nZb
!
!          ELSE IF (((zI >= Z(nZa)).AND.(zI <= Z(nZb))) .AND. &
!           & ((zF >= Z(nZa)).AND.(zF <= Z(nZb)))) THEN
!
!            DO l = nZa, nZb
!               IF ((zI >= Z(l-1)).AND.(zI <= Z(l))) THEN
!                  nZi = l
!               END IF
!            END DO
!
!            DO l = nZa, nZb
!               IF ((zF >= Z(l-1)).AND.(zF <= Z(l))) THEN
!                  nZf = l
!               END IF
!            END DO
!
!         ELSE
!
!            ! z-axis integration limits do not contain this slap; do nothing
!
!         END IF
!
!      ELSE IF ((xI < xF).AND.(yI < yF).AND.(zI == zF)) THEN
!
!         DO j = nXa, nXb
!            IF ((xI >= X(j-1)).AND.(xI <= X(j))) THEN
!               nXi = j
!            END IF
!         END DO
!
!         DO j = nXa, nXb
!            IF ((xF >= X(j-1)).AND.(xF <= X(j))) THEN
!               nXf = j
!            END IF
!         END DO
!
!         DO k = nYa, nYb
!            IF ((yI >= Y(k-1)).AND.(yI <= Y(k))) THEN
!               nYi = k
!            END IF
!         END DO
!
!         DO k = nYa, nYb
!            IF ((yF >= Y(k-1)).AND.(yF <= Y(k))) THEN
!               nYf = k
!            END IF
!         END DO
!
!         DO l = nZa, nZb
!            IF ((zI >= Z(l-1)).AND.(zI <= Z(l))) THEN
!               nZi = l
!               nZf = nZi
!            END IF
!         END DO
!
!      ELSE
!
!         WRITE(UNIT=ERROR_UNIT, FMT=*) 'gpse: ERROR - Integration &
!            & limits are ill-defined.'
!
!      END IF
!
!      WRITE(UNIT=OUTPUT_UNIT, FMT=*) mpiRank, xI, xF, yI, yF, zI, zF
!      WRITE(UNIT=OUTPUT_UNIT, FMT=*) mpiRank, X(nXa), X(nXb), Y(nYa), Y(nYb), Z(nZa), Z(nZb)
!      WRITE(UNIT=OUTPUT_UNIT, FMT=*) mpiRank, nXa, nXb, nYa, nYb, nZa, nZb
!      WRITE(UNIT=OUTPUT_UNIT, FMT=*) mpiRank, nXi, nXf, nYi, nYf, nZi, nZf
!
!      ! check if X, Y, Z, are in bounds of integration limits; only compute integral partitions if they are; otherwise wait for all other processes to compute integral portions
!
!      ! it could be even worse; you'd need to limit the integration limits by the loop counters .... ah; i.e., nXa, nXb might only cover a portion of the domain to integrate on the surface that is avaiable for this individual MPI process
!
!      ! *******  write a new function or subroutine to return the index limits of some integration bounds on the grid; probably best way to do this? use some sort of counterpropagating bisection method to find index bounds? *******
!
!      RETURN
!      END FUNCTION
!
! ----------------------------------------------------------------------

      END MODULE

! ======================================================================
