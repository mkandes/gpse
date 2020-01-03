!  ======================================================================
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
!     Created to study the quantum Compton generator data, where, in 
!     particular, the more simple angular momentum expectation values 
!     used initially led to large disagreements with the theory due to 
!     the precense of many vortex-like structures and other complicating
!     factors in the raw simulation data.
!
! NOTES
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
!     Sunday, November 24th, 2019
!
! ----------------------------------------------------------------------

      MODULE pmca

! --- MODULE DECLARATIONS ----------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE :: MPI
      USE :: GRID

! --- MODULE DEFINITIONS -----------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides 
!        information about the run-time environment.
!
!     MPI is the standard Message Passing Interface module.
!
!     GRID is a custom Fortran module written to build the regular
!        Cartesian grid of points that defines the spatial computational
!        domain on which external potentials, wave functions and other
!        physical quantities may be computed.
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
      PUBLIC :: pmca_velocity2
      PUBLIC :: pmca_current_density
      PUBLIC :: pmca_current_density2
      PUBLIC :: pmca_quantum_potential
      PUBLIC :: pmca_vorticity

      PRIVATE :: pmca_grad_f_real_3d_rect_cd2
      PRIVATE :: pmca_grad_f_cmplx_3d_rect_cd2
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

      SUBROUTINE pmca_compute_currents(mpiRank, mpiMaster, mpiReal, & 
         & mpiError, pmcaQuadRule, nXa, nXb, nXbc, nYa, nYb, nYbc, nZa,&
         & nZb, nZbc, nZ, dX, dY, dZ, Xa, Ya, Za, Zc, J3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: mpiRank, mpiMaster, mpiReal
      INTEGER, INTENT(INOUT) :: mpiError
      INTEGER, INTENT(IN) :: pmcaQuadRule
      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc, nZ

      REAL, INTENT (IN) :: dX, dY, dZ

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(IN) :: Xa
      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(IN) :: Ya
      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(IN) :: Za
      REAL, DIMENSION(  1 - nZbc : nZ  + nZbc), INTENT(IN) :: Zc

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc) :: J3

      REAL :: pmcaI = 0.0 ! temporary, local slab current variable

      IF (pmcaQuadRule == 1) THEN

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc,&
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, 0.0, 0.0, 0.0, Ya(nYb), &
            & Zc(1), Zc(nZ), Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI12, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, Xa(nXa), 0.0, 0.0, 0.0, &
            & Zc(1), Zc(nZ), Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI23, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, 0.0, 0.0, Ya(nYa), 0.0, &
            & Zc(1), Zc(nZ), Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI34, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, 0.0, Xa(nXb), 0.0, 0.0, &
            & Zc(1), Zc(nZ), Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI41, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, 0.0, Xa(nXb), 0.0, &
            & Ya(nYb), 0.0, 0.0, Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI1, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, Xa(nXa), 0.0, 0.0, &
            & Ya(nYb), 0.0, 0.0, Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI2, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, Xa(nXa), 0.0, Ya(nYa), &
            & 0.0, 0.0, 0.0, Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI3, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

         pmcaI = pmca_current_3d_rect(nXa, nXb, nXbc, nYa, nYb, nYbc, &
            & nZa, nZb, nZbc, nZ, dX, dY, dZ, 0.0, Xa(nXb), Ya(nYa), &
            & 0.0, 0.0, 0.0, Xa, Ya, Za, Zc, J3)
         CALL MPI_REDUCE(pmcaI, pmcaI4, 1, mpiReal, MPI_SUM, &
            & mpiMaster, MPI_COMM_WORLD, mpiError)

      ELSE

         WRITE(UNIT=ERROR_UNIT, FMT = *) 'gpse: pmca_compute_currents: &
            & ERROR - pmcaQuadRule is not supported.'
         STOP

      END IF

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_write_currents(mpiRank, mpiMaster, tN)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: mpiRank, mpiMaster

      REAL, INTENT (IN) :: tN

!     Write computed currents to file from mpiMaster
      IF ( mpiRank == mpiMaster ) THEN

         OPEN(UNIT=990, FILE='pmca.output', ACCESS='APPEND', &
            & ACTION='WRITE', FORM='FORMATTED', STATUS='UNKNOWN')
         WRITE(UNIT=990, FMT='(9(F23.15))') tN, pmcaI1, &
            & pmcaI12, pmcaI2, pmcaI23, pmcaI3, pmcaI34, pmcaI4, &
            & pmcaI41
         CLOSE(UNIT=990, STATUS='KEEP')

      END IF

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
         & nZb, nZbc, dX, dY, dZ, Rho3, J3, V3)
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
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: Rho3

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: J3

      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: V3

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_velocity2(nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, &
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
         & nZa, nZb, nZbc, dX, dY, dZ, Psi3, GradPsi3, J3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc

      REAL, INTENT(IN) :: dX, dY, dZ

      COMPLEX, DIMENSION(nXa - nXbc : nXb + nXbc, &
                       & nYa - nYbc : nYb + nYbc, &
                       & nZa - nZbc : nZb + nZbc), INTENT(IN) :: Psi3

      COMPLEX, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: GradPsi3

      COMPLEX, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: J3

      INTEGER :: j, k, l

      J3 = CMPLX(0.0,0.0)

      CALL pmca_grad_f_cmplx_3d_rect_cd2(nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, nZb, nZbc, dX, dY, dZ, Psi3, GradPsi3)
      J3(1,:,:,:) = CONJG(Psi3) * GradPsi3(1,:,:,:)
      J3(2,:,:,:) = CONJG(Psi3) * GradPsi3(2,:,:,:)
      J3(3,:,:,:) = CONJG(Psi3) * GradPsi3(3,:,:,:)
      CALL pmca_grad_f_cmplx_3d_rect_cd2(nXa, nXb, nXbc, nYa, nYb, nYbc, nZa, nZb, nZbc, dX, dY, dZ, CONJG(Psi3), GradPsi3)
      J3(1,:,:,:) = J3(1,:,:,:) - Psi3 * GradPsi3(1,:,:,:)
      J3(2,:,:,:) = J3(2,:,:,:) - Psi3 * GradPsi3(2,:,:,:)
      J3(3,:,:,:) = J3(3,:,:,:) - Psi3 * GradPsi3(3,:,:,:)
      J3 = CMPLX(0.0,0.5) * J3

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_current_density2(nXa, nXb, nXbc, nYa, nYb, nYbc, &
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
               GradF3(1,j,k,l) = (0.5/dX) * (F3(j+1,k,l) - F3(j-1,k,l))
               GradF3(2,j,k,l) = (0.5/dY) * (F3(j,k+1,l) - F3(j,k-1,l))
               GradF3(3,j,k,l) = (0.5/dZ) * (F3(j,k,l+1) - F3(j,k,l-1))
            END DO
         END DO
!$OMP    END PARALLEL DO
      END DO
!$OMP END PARALLEL DO 

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE pmca_grad_f_cmplx_3d_rect_cd2(nXa, nXb, nXbc, nYa, &
         & nYb, nYbc, nZa, nZb, nZbc, dX, dY, dZ, F3, GradF3)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc

      REAL, INTENT(IN) :: dX, dY, dZ

      COMPLEX, DIMENSION(nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: F3

      COMPLEX, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(INOUT) :: GradF3

      INTEGER :: j, k, l

      GradF3 = CMPLX(0.0,0.0)

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb
!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb
            DO j = nXa , nXb
               GradF3(1,j,k,l) = CMPLX(0.5/dX,0.0) * (F3(j+1,k,l) - F3(j-1,k,l))
               GradF3(2,j,k,l) = CMPLX(0.5/dY,0.0) * (F3(j,k+1,l) - F3(j,k-1,l))
               GradF3(3,j,k,l) = CMPLX(0.5/dZ,0.0) * (F3(j,k,l+1) - F3(j,k,l-1))
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
         & nYbc, nZa, nZb, nZbc, nZ, dX, dY, dZ, xI, xF, yI, yF, zI, zF,&
         & Xa, Ya, Za, Zc, J3) RESULT(pmcaI)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nXa, nXb, nXbc
      INTEGER, INTENT(IN) :: nYa, nYb, nYbc
      INTEGER, INTENT(IN) :: nZa, nZb, nZbc, nZ

      REAL, INTENT(IN) :: dX, dY, dZ

      REAL, INTENT(IN) :: xI, xF
      REAL, INTENT(IN) :: yI, yF
      REAL, INTENT(IN) :: zI, zF

      REAL, DIMENSION(nXa - nXbc : nXb + nXbc), INTENT(IN) :: Xa
      REAL, DIMENSION(nYa - nYbc : nYb + nYbc), INTENT(IN) :: Ya
      REAL, DIMENSION(nZa - nZbc : nZb + nZbc), INTENT(IN) :: Za
      REAL, DIMENSION(  1 - nZbc : nZ  + nZbc), INTENT(IN) :: Zc


      REAL, DIMENSION(3, &
                    & nXa - nXbc : nXb + nXbc, &
                    & nYa - nYbc : nYb + nYbc, &
                    & nZa - nZbc : nZb + nZbc), INTENT(IN) :: J3

      INTEGER :: j, k, l

      INTEGER :: nXiG = -1, nXfG = -1
      INTEGER :: nYiG = -1, nYfG = -1
      INTEGER :: nZiG = -1, nZfG = -1
      INTEGER :: nZiL = -1, nZfL = -1

      pmcaI = 0.0

!     Search for global location of integration limits in x, y, and z
      nXiG = grid_linear_search(nXa, nXb, nXbc, xI, Xa)
      nXfG = grid_linear_search(nXa, nXb, nXbc, xF, Xa)

      nYiG = grid_linear_search(nYa, nYb, nYbc, yI, Ya)
      nYfG = grid_linear_search(nYa, nYb, nYbc, yF, Ya)

      nZiG = grid_linear_search(1, nZ, nZbc, zI, Zc)
      nZfG = grid_linear_search(1, nZ, nZbc, zF, Zc)

!     1                        nZiG           nZfG                     nZ
!     |---------------|--------*-------|------*-------|----------------|
!     nZa1        nZb1 nZa2        nZb2 nZa3      nZb3 nZa4         nZa4
!     |------*-----*--|----------------|--------------|----------------|
!            nZiG  nZfG
!     |------*--------|----------------|--------------|-----*----------|
!            nZiG                                           nZfG

!     Search for slab-local location of integration limits in z
      IF ((nZa <= nZiG).AND.(nZiG <= nZb)) THEN
         nZiL = nZiG
      ELSE IF ((nZa > nZiG).AND.(nZa < nZfG)) THEN
         nZiL = nZa
      END IF

      IF ((nZa <= nZfG).AND.(nZfG <= nZb)) THEN
         nZfL = nZfG
      ELSE IF ((nZb > nZiG).AND.(nZb < nZfG)) then
         nZfL = nZb
      END IF

      IF (nXiG == nXfG) THEN

         DO l = nZiL, nZfL
            DO k = nYiG, nYfG
               pmcaI = pmcaI + J3(1,nXiG,k,l)
            END DO
         END DO
         pmcaI = pmcaI * dY * dZ

      ELSE IF (nYiG == nYfG) THEN

         DO l = nZiL, nZfL
            DO j = nXiG, nXfG
               pmcaI = pmcaI + J3(2,j,nYiG,l)
            END DO
         END DO
         pmcaI = pmcaI * dX * dZ

      ELSE IF (nZiG == nZfG) THEN 

         IF ((nZiG == nZiL).AND.(nZfG == nZfL)) THEN

            DO k = nYiG, nYfG
               DO j = nXiG, nXfG
                  pmcaI = pmcaI + J3(3,j,k,nZiG)
               END DO
            END DO
            pmcaI = pmcaI * dX * dY

         END IF

      ELSE

         WRITE(UNIT=ERROR_UNIT, FMT = *) 'gpse: pmca_current_3d_rect: &
            & ERROR - Surface integral must currently lie in one of &
            & the Cartesian planes.'
         STOP

      END IF

      RETURN
      END FUNCTION

! ----------------------------------------------------------------------

      END MODULE

! ======================================================================
