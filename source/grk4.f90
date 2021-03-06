! ==================================================================================================================================
! NAME
!
!     grk4 [ ] - Generalized 4th-Order Runge-Kutta Module
!
! SYNOPSIS
!
!     USE :: GRK4
!
! DESCRIPTION  
!
!     GRK4 is a custom Fortran module written to numerically approximate solutions of the time-dependent, three-dimensional 
!        Gross-Pitaveskii equation for a single component Bose-Einstein condensate in a rotating reference frame using the GRK4 + 
!        CDX algorithm.
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
!     Tuesday, April 2nd, 2019
!
! ----------------------------------------------------------------------------------------------------------------------------------

      MODULE GRK4

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      PRIVATE

! --- SUBROUTINE DECLARATIONS ------------------------------------------------------------------------------------------------------

      PUBLIC  :: grk4_y_3d_stgx
      PUBLIC  :: grk4_f_gp_3d_rrf_cdx

      PRIVATE :: grk4_y_3d_stg1
      PRIVATE :: grk4_y_3d_stg2
      PRIVATE :: grk4_y_3d_stg3
      PRIVATE :: grk4_y_3d_stg4
      PRIVATE :: grk4_f_gp_3d_rrf_cd2
      PRIVATE :: grk4_f_gp_3d_rrf_cd4
      PRIVATE :: grk4_f_gp_3d_rrf_cd6
      PRIVATE :: grk4_f_gp_3d_rrf_cd8

! --- SUBROUTINE DEFINITIONS -------------------------------------------------------------------------------------------------------
!
!     grk4_y_3d_stgx is a PUBLIC SUBROUTINE that determines which stage of the GRK4 + CDX algorithm is called.
!
!     grk4_f_gp_3d_rrf_cdx is a PUBLIC SUBROUTINE that determines the order-of-accuracy for the central finite differences used to 
!        approximate the spatial derivaties of the three-dimensional Gross-Pitaevskii equation in a rotating reference frame. 
!
!     grk4_y_3d_stg1 is a PRIVATE SUBROUTINE that computes the intermediate wave function for the 2nd stage of the GRK4 + CDX 
!        algorithm using the result computed for the increment of the 1st stage. 
!
!     grk4_y_3d_stg2 is a PRIVATE SUBROUTINE that computes the intermediate wave function for the 3rd stage of the GRK4 + CDX 
!        algorithm using the result computed for the increments of the 1st and 2nd stages.
!
!     grk4_y_3d_stg3 is a PRIVATE SUBROUTINE that computes the intermediate wave function for the 4th stage of the GRK4 + CDX 
!        algorithm using the result computed for the increments of the 2nd and 3rd stages.
!
!     grk4_y_3d_stg4 is a PRIVATE SUBROUTINE that computes the wave function at the next time step using all four of the increments 
!        computed for GRK4 + CDX algorithm.
!        
!     grk4_f_gp_3d_rrf_cd2 is a PRIVATE SUBROUTINE that computes the GRK4 increment for the three-dimensional Gross-Pitaevskii 
!        equation in a rotating reference frame at any stage using 2nd-order central finite differences.
!
!     grk4_f_gp_3d_rrf_cd4 is a PRIVATE SUBROUTINE that computes the GRK4 increment for the three-dimensional Gross-Pitaevskii 
!        equation in a rotating reference frame at any stage using 4th-order central finite differences.
!
!     grk4_f_gp_3d_rrf_cd6 is a PRIVATE SUBROUTINE that computes the GRK4 increment for the three-dimensional Gross-Pitaevskii 
!        equation in a rotating reference frame at any stage using 6th-order central finite differences. ** NOTE: THIS 
!        ORDER-OF-ACCURACY HAS NOT BEEN WELL-TESTED.
!
!     grk4_f_gp_3d_rrf_cd8 is a PRIVATE SUBROUTINE that computes the GRK4 increment for the three-dimensional Gross-Pitaevskii 
!        equation in a rotating reference frame at any stage using 8th-order central finite differences. ** NOTE: THIS 
!        ORDER-OF-ACCURACY HAS NOT BEEN WELL-TESTED.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE grk4_y_3d_stgx ( grk4Stage , grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 ,&
         K3 , K4 , Psi3a , Psi3b )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: grk4Stage
      INTEGER, INTENT ( IN ) :: grk4Lambda
      INTEGER, INTENT ( IN ) :: nXa 
      INTEGER, INTENT ( IN ) :: nXb 
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa 
      INTEGER, INTENT ( IN ) :: nYb 
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa 
      INTEGER, INTENT ( IN ) :: nZb 
      INTEGER, INTENT ( IN ) :: nZbc

      COMPLEX, INTENT ( IN ) :: dTz 

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K1
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K2
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K3
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K4
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi3a
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3b

      IF ( grk4Stage == 1 ) THEN
 
         CALL grk4_y_3d_stg1 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , Psi3a , Psi3b )

      ELSE IF ( grk4Stage == 2 ) THEN

         CALL grk4_y_3d_stg2 ( grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , Psi3a , Psi3b )

      ELSE IF ( grk4Stage == 3 ) THEN

         CALL grk4_y_3d_stg3 ( grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K2 , K3 , Psi3a , Psi3b )

      ELSE IF ( grk4Stage == 4 ) THEN

         CALL grk4_y_3d_stg4 ( grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , &
            & Psi3a , Psi3b )

      ELSE

         WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : grk4_y_3d_stgx : ERROR - grk4Stage not recognized.'
         STOP

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY ,&
         & dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: fdOrder
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc 
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc 
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO 
      REAL, INTENT ( IN ) :: dX
      REAL, INTENT ( IN ) :: dY
      REAL, INTENT ( IN ) :: dZ
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: gS

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F

      IF ( fdOrder == 2 ) THEN

         CALL grk4_f_gp_3d_rrf_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , wX , &
            & wY , wZ , gS , X , Y , Z , Vex , Psi , F )

      ELSE IF ( fdOrder == 4 ) THEN

         CALL grk4_f_gp_3d_rrf_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , wX , &
            & wY , wZ , gS , X , Y , Z , Vex , Psi , F )

      ELSE IF ( fdOrder == 6 ) THEN

         CALL grk4_f_gp_3d_rrf_cd6 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , &
            & X , Y , Z , Vex , Psi , F )

      ELSE IF ( fdOrder == 8 ) THEN

         CALL grk4_f_gp_3d_rrf_cd8 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , &
            & X , Y , Z , Vex , Psi , F )

      ELSE

         WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : grk4_f_gp_3d_rrf_cdx : ERROR - fdOrder not supported.'
         STOP

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

!     y_n + 0.5 * dT * k_1

      SUBROUTINE grk4_y_3d_stg1 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , Psi3a , Psi3b )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      COMPLEX, INTENT ( IN ) :: dTz

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K1
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi3a
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3b

      INTEGER :: j , k , l 

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) 
         DO k = nYa , nYb

            DO j = nXa , nXb

               Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + CMPLX ( 0.5 , 0.0 ) * dTz * K1 ( j , k , l )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

!     y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2

      SUBROUTINE grk4_y_3d_stg2 ( grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , Psi3a , &
         & Psi3b )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: grk4Lambda
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc
                     
      COMPLEX, INTENT ( IN ) :: dTz

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K1
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K2
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi3a
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3b

      INTEGER :: j , k , l 

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + dTz * ( CMPLX ( 0.5 - 1.0 / REAL ( grk4Lambda ) , 0.0 ) * &
                  & K1 ( j , k , l ) + CMPLX ( 1.0 / REAL ( grk4Lambda ) , 0.0 ) * K2 ( j , k , l ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

!     y_n + ( 1 - lambda / 2 ) * dT * k_2 + ( lambda / 2 ) * dT * k_3

      SUBROUTINE grk4_y_3d_stg3 ( grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K2 , K3 , Psi3a , &
         & Psi3b )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: grk4Lambda
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      COMPLEX, INTENT ( IN ) :: dTz

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K2
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K3
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi3a
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3b

      INTEGER :: j , k , l 

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + dTz * ( CMPLX ( 1.0 - 0.5 * REAL ( grk4Lambda ) , 0.0 ) * &
                  & K2 ( j , k , l ) + CMPLX ( 0.5 * REAL ( grk4Lambda ) , 0.0 ) * K3 ( j , k , l ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

!     y_{ n + 1 } = y_n + ( dT / 6 ) * [ k_1 + ( 4 - lambda ) * k_2 + lambda * k_3 + k_4 ]

      SUBROUTINE grk4_y_3d_stg4 ( grk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , &
         & Psi3a , Psi3b )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: grk4Lambda
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      COMPLEX, INTENT ( IN ) :: dTz

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K1
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K2
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K3
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: K4
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi3a
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3b

      INTEGER :: j , k , l 

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + CMPLX ( 1.0 / 6.0 , 0.0 ) * dTz * ( K1 ( j , k , l ) + &
                  & CMPLX ( 4.0 - REAL ( grk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( REAL ( grk4Lambda ) , 0.0 ) * &
                  & K3 ( j , k , l ) + K4 ( j , k , l ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE grk4_f_gp_3d_rrf_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , wX ,&
         & wY , wZ , gS , X , Y , Z , Vex , Psi , F )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: dX
      REAL, INTENT ( IN ) :: dY
      REAL, INTENT ( IN ) :: dZ
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: gS

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               F ( j , k , l ) = &
                  & CMPLX ( 0.5 * ( wY * ( X ( j ) - xO ) - wX * ( Y ( k ) - yO ) ) / dZ , 0.5 / dZ**2 ) * Psi ( j , k , l - 1 ) + &
                  & CMPLX ( 0.5 * ( wX * ( Z ( l ) - zO ) - wZ * ( X ( j ) - xO ) ) / dY , 0.5 / dY**2 ) * Psi ( j , k - 1 , l ) + &
                  & CMPLX ( 0.5 * ( wZ * ( Y ( k ) - yO ) - wY * ( Z ( l ) - zO ) ) / dX , 0.5 / dX**2 ) * Psi ( j - 1 , k , l ) - &
                  & CMPLX ( 0.0 , 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 + Vex ( j , k , l ) + &
                  &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                  & CMPLX ( 0.5 * ( wY * ( Z ( l ) - zO ) - wZ * ( Y ( k ) - yO ) ) / dX , 0.5 / dX**2 ) * Psi ( j + 1 , k , l ) + &
                  & CMPLX ( 0.5 * ( wZ * ( X ( j ) - xO ) - wX * ( Z ( l ) - zO ) ) / dY , 0.5 / dY**2 ) * Psi ( j , k + 1 , l ) + &
                  & CMPLX ( 0.5 * ( wX * ( Y ( k ) - yO ) - wY * ( X ( j ) - xO ) ) / dZ , 0.5 / dZ**2 ) * Psi ( j , k , l + 1 )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO 

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE grk4_f_gp_3d_rrf_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , wX ,&
         & wY , wZ , gS , X , Y , Z , Vex , Psi , F ) 

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: dX
      REAL, INTENT ( IN ) :: dY
      REAL, INTENT ( IN ) :: dZ
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: gS

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               F ( j , k , l ) = &
                  & CMPLX ( ( wX * ( Y ( k ) - yO ) - wY * ( X ( j ) - xO ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 2 ) + &
                  & CMPLX ( 0.75 * ( wY * ( X ( j ) - xO ) - wX * ( Y ( k ) - yO ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 1 ) + &
                  & CMPLX ( ( wZ * ( X ( j ) - xO ) - wX * ( Z ( l ) - zO ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * &
                  &    Psi ( j , k - 2 , l ) + &
                  & CMPLX ( 0.75 * ( wX * ( Z ( l ) - zO ) - wZ * ( X ( j ) - xO ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * &
                  &    Psi ( j , k - 1 , l ) + &
                  & CMPLX ( ( wY * ( Z ( l ) - zO ) - wZ * ( Y ( k ) - yO ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * &
                  &    Psi ( j - 2 , k , l ) + &
                  & CMPLX ( 0.75 * ( wZ * ( Y ( k ) - yO ) - wY * ( Z ( l ) - zO ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * &
                  &    Psi ( j - 1 , k , l ) - &
                  & CMPLX ( 0.0 , 1.25 * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex ( j , k , l ) + &
                  &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                  & CMPLX ( 0.75 * ( wY * ( Z ( l ) - zO ) - wZ * ( Y ( k ) - yO ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * &
                  &    Psi ( j + 1 , k , l ) + &
                  & CMPLX ( ( wZ * ( Y ( k ) - yO ) - wY * ( Z ( l ) - zO ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * &
                  &    Psi ( j + 2 , k , l ) + &
                  & CMPLX ( 0.75 * ( wZ * ( X ( j ) - xO ) - wX * ( Z ( l ) - zO ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 1 , l ) + &
                  & CMPLX ( ( wX * ( Z ( l ) - zO ) - wZ * ( X ( j ) - xO ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 2 , l ) + &
                  & CMPLX ( 0.75 * ( wX * ( Y ( k ) - yO ) - wY * ( X ( j ) - xO ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l + 1 ) + &
                  & CMPLX ( ( wY * ( X ( j ) - xO ) - wX * ( Y ( k ) - yO ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l + 2 )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE grk4_f_gp_3d_rrf_cd6 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS ,&
         & X , Y , Z , Vex , Psi , F ) 

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      REAL, INTENT ( IN ) :: dX
      REAL, INTENT ( IN ) :: dY
      REAL, INTENT ( IN ) :: dZ
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: gS

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               F ( j , k , l ) = &
                  & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi ( j , k , l - 3 ) + &
                  & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 2 ) + &
                  & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 1 ) + &
                  & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi ( j , k - 3 , l ) + &
                  & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * &
                  &    Psi ( j , k - 2 , l ) + &
                  & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * &
                  &    Psi ( j , k - 1 , l ) + &
                  & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi ( j - 3 , k , l ) + &
                  & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * &
                  &    Psi ( j - 2 , k , l ) + &
                  & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * &
                  &    Psi ( j - 1 , k , l ) - &
                  & CMPLX ( 0.0 , ( 49.0 / 36.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex ( j , k , l ) + &
                  &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                  & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * &
                  &    Psi ( j + 1 , k , l ) + &
                  & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * &
                  &    Psi ( j + 2 , k , l ) + &
                  & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi ( j + 3 , k , l ) + &
                  & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 1 , l ) + &
                  & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 2 , l ) + &
                  & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi ( j , k + 3 , l ) + &
                  & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l + 1 ) + &
                  & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l + 2 ) + &
                  & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi ( j , k , l + 3 )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE grk4_f_gp_3d_rrf_cd8 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS ,&
         & X , Y , Z , Vex , Psi , F ) 

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc

      REAL, INTENT ( IN ) :: dX
      REAL, INTENT ( IN ) :: dY
      REAL, INTENT ( IN ) :: dZ
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: gS

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN    ) :: Psi
      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               F ( j , k , l ) = &
                  & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 4 ) + &
                  & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 105.0 * dZ ) ,  4.0 / ( 315.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 3 ) + &
                  & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dZ ) , -1.0 / ( 10.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 2 ) + &
                  & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dZ ) ,  4.0 / ( 5.0 * dZ**2 ) ) * &
                  &    Psi ( j , k , l - 1 ) + &
                  & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 280.0 * dY ) , -1.0 / ( 1120 * dY**2 ) ) * Psi ( j , k - 4 , l ) + &
                  & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * &
                  &    Psi ( j , k - 3 , l ) + &
                  & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi ( j , k - 2 , l ) + &
                  & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * &
                  &    Psi ( j , k - 1 , l ) + &
                  & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * &
                  &    Psi ( j - 4 , k , l ) + &
                  & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * &
                  &    Psi ( j - 3 , k , l ) + &
                  & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi ( j - 2 , k , l ) + &
                  & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * &
                  &    Psi ( j - 1 , k , l ) - &
                  & CMPLX ( 0.0 , ( 205.0 / 144.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex ( j , k , l ) + &
                  &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                  & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * &
                  &    Psi ( j + 1 , k , l ) + &
                  & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi ( j + 2 , k , l ) + &
                  & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * &
                  &    Psi ( j + 3 , k , l ) + &
                  & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * &
                  &    Psi ( j + 4 , k , l ) + &
                  & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 1 , l ) + &
                  & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi ( j , k + 2 , l ) + &
                  & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 3 , l ) + &
                  & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 280.0 * dY ) , -1.0 / ( 1120.0 * dY**2 ) ) * &
                  &    Psi ( j , k + 4 , l ) + &
                  & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * &
                  &    Psi ( j , k , l + 1 ) + &
                  & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi ( j , k , l + 2 ) + &
                  & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dY**2 ) ) * &
                  &    Psi ( j , k , l + 3 ) + &
                  & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi ( j , k , l + 4 )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      END MODULE

! ==================================================================================================================================
