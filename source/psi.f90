! ==========================================================================
! NAME
!
!     psi [ (p)sī ] - Psi Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     psi is a Fortran module ... 
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
!     Saturday, June 7th, 2014
!
! -------------------------------------------------------------------------

      MODULE PSI

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MATH

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: psi_3d_se_sho_ani
      PUBLIC :: psi_3d_se_sho_axi
      !PUBLIC :: psi_3d_se_sho_iso

      CONTAINS

         SUBROUTINE psi_3d_se_sho_ani ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , psiNx , psiNy , psiNz , psiXo , &
            & psiYo , psiZo , psiWx , psiWy , psiWz , X , Y , Z , Psi3 )
 
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
            INTEGER, INTENT ( IN ) :: psiNx
            INTEGER, INTENT ( IN ) :: psiNy
            INTEGER, INTENT ( IN ) :: psiNz

            REAL, INTENT ( IN ) :: psiXo
            REAL, INTENT ( IN ) :: psiYo
            REAL, INTENT ( IN ) :: psiZo
            REAL, INTENT ( IN ) :: psiWx
            REAL, INTENT ( IN ) :: psiWy
            REAL, INTENT ( IN ) :: psiWz
 
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Psi3 ( j , k , l )  = CMPLX ( & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNx * factorial ( psiNx ) ) ) ) * SQRT ( SQRT ( psiWx / PI ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNy * factorial ( psiNy ) ) ) ) * SQRT ( SQRT ( psiWy / PI ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNz * factorial ( psiNz ) ) ) ) * SQRT ( SQRT ( psiWz / PI ) ) * & 
                        & hermite ( psiNx , SQRT ( psiWx ) * ( X ( j ) - psiXo ) ) * EXP ( -0.5 * psiWx * ( X ( j ) - psiXo )**2 ) * & 
                        & hermite ( psiNy , SQRT ( psiWy ) * ( Y ( k ) - psiYo ) ) * EXP ( -0.5 * psiWy * ( Y ( k ) - psiYo )**2 ) * & 
                        & hermite ( psiNz , SQRT ( psiWz ) * ( Z ( l ) - psiZo ) ) * EXP ( -0.5 * psiWz * ( Z ( l ) - psiZo )**2 ) , 0.0 )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

!         COMPLEX FUNCTION psi_3d_se_sho_ani ( nX , nY , nZ , xO , yO , zO , wX , wY , wZ , x , y , z )
!
!         IMPLICIT NONE
!
!         INTEGER, INTENT ( IN ) :: nX
!         INTEGER, INTENT ( IN ) :: nY
!         INTEGER, INTENT ( IN ) :: nZ
!          
!         REAL, INTENT ( IN ) :: xO
!         REAL, INTENT ( IN ) :: yO
!         REAL, INTENT ( IN ) :: zO
!         REAL, INTENT ( IN ) :: wX
!         REAL, INTENT ( IN ) :: wY
!         REAL, INTENT ( IN ) :: wZ
!         REAL, INTENT ( IN ) :: x
!         REAL, INTENT ( IN ) :: y
!         REAL, INTENT ( IN ) :: z
!
!         psi_3d_se_sho_ani = CMPLX ( &
!            & ( 1.0 / SQRT ( REAL ( 2**nX * factorial ( nX ) ) ) ) * SQRT ( SQRT ( wX / PI ) ) * &
!            & ( 1.0 / SQRT ( REAL ( 2**nY * factorial ( nY ) ) ) ) * SQRT ( SQRT ( wY / PI ) ) * &
!            & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
!            & hermite ( nX , SQRT ( wX ) * ( x - xO ) ) * EXP ( -0.5 * wX * ( x - xO )**2 ) * &
!            & hermite ( nY , SQRT ( wY ) * ( y - yO ) ) * EXP ( -0.5 * wY * ( y - yO )**2 ) * &
!            & hermite ( nZ , SQRT ( wZ ) * ( z - zO ) ) * EXP ( -0.5 * wZ * ( z - zO )**2 ) , 0.0 )
!
!         RETURN
!
!         END FUNCTION

         SUBROUTINE psi_3d_se_sho_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , psiNr , psiMl , psiNz , psiXo , & 
            & psiYo , psiZo , psiWr , psiWz , X , Y , Z , Psi3 )

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
            INTEGER, INTENT ( IN ) :: psiNr
            INTEGER, INTENT ( IN ) :: psiMl
            INTEGER, INTENT ( IN ) :: psiNz

            REAL, INTENT ( IN ) :: psiXo
            REAL, INTENT ( IN ) :: psiYo
            REAL, INTENT ( IN ) :: psiZo
            REAL, INTENT ( IN ) :: psiWr
            REAL, INTENT ( IN ) :: psiWz
 
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb  

                     Psi3 ( j , k , l ) = CMPLX ( &
                        & SQRT ( ( psiWr**( ABS ( psiMl ) + 1 ) * REAL ( factorial ( psiNr ) ) ) / & 
                        & ( PI * REAL ( factorial ( psiNr + ABS ( psiMl ) ) ) ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNz * factorial ( psiNz ) ) ) ) * SQRT ( SQRT ( psiWz / PI ) ) * & 
                        & SQRT ( ( X ( j ) - psiXo )**2 + ( Y ( k ) - psiYo )**2 )**ABS ( psiMl ) * & 
                        & alaguerre ( psiNr , ABS ( psiMl ) , psiWr * ( ( X ( j ) - psiXo )**2 + ( Y ( k ) - psiYo )**2 ) ) * & 
                        & EXP ( -0.5 * psiWr * ( ( X ( j ) - psiXo )**2 + ( Y ( k ) - psiYo )**2 ) ) * & 
                        & hermite ( psiNz , SQRT ( psiWz ) * ( Z ( l ) - psiZo ) ) * EXP ( -0.5 * psiWz * ( Z ( l ) - psiZo )**2 ) , 0.0 ) * & 
                        & EXP ( CMPLX ( 0.0 , REAL ( psiMl ) * ATAN2 ( Y ( k ) - psiYo , X ( j ) - psiXo ) ) ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

!         COMPLEX FUNCTION psi_3d_se_sho_axi ( nR , mL , nZ , xO , yO , zO , wR , wZ , x , y , z )
!
!         IMPLICIT NONE
!
!         INTEGER, INTENT ( IN ) :: nR
!         INTEGER, INTENT ( IN ) :: mL
!         INTEGER, INTENT ( IN ) :: nZ
!
!         REAL, INTENT ( IN ) :: xO
!         REAL, INTENT ( IN ) :: yO
!         REAL, INTENT ( IN ) :: zO
!         REAL, INTENT ( IN ) :: wR
!         REAL, INTENT ( IN ) :: wZ
!         REAL, INTENT ( IN ) :: x
!         REAL, INTENT ( IN ) :: y
!         REAL, INTENT ( IN ) :: z
!
!         psi_3d_se_sho_axi = CMPLX ( &
!            & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / &
!            & ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * & 
!            & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
!            & SQRT ( ( x - xO )**2 + ( y - yO )**2 )**ABS ( mL ) * &
!            & alaguerre ( nR , ABS ( mL ) , wR * ( ( x - xO )**2 + ( y - yO )**2 ) ) * &
!            & EXP ( -0.5 * wR * ( ( x - xO )**2 + ( y - yO )**2 ) ) * &
!            & hermite ( nZ , SQRT ( wZ ) * ( z - zO ) ) * EXP ( -0.5 * wZ * ( z - zO )**2 ) , 0.0 ) * &
!            & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( y - yO , x - xO ) ) )
!
!         RETURN
!
!         END FUNCTION

      END MODULE

! =========================================================================
