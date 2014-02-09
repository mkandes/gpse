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
!     Sunday, February 9th, 2014
!
! -------------------------------------------------------------------------

      MODULE PSI

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         USE MATH

         IMPLICIT NONE
         PRIVATE

         PUBLIC :: psi_3d_se_sho_ani
         PUBLIC :: psi_3d_se_sho_axi
         PUBLIC :: psi_se_sho_3d_ani
         PUBLIC :: psi_se_sho_3d_axi
         !PUBLIC :: psi_se_sho_3d_iso

         CHARACTER ( LEN = * ), PARAMETER :: VERSION_NUMBER = '0.0.7'

         CONTAINS

            COMPLEX FUNCTION psi_3d_se_sho_ani ( nX , nY , nZ , x , y , z , xO , yO , zO , wX , wY , wZ )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: nX
               INTEGER, INTENT ( IN ) :: nY
               INTEGER, INTENT ( IN ) :: nZ
          
               REAL, INTENT ( IN ) :: x
               REAL, INTENT ( IN ) :: y
               REAL, INTENT ( IN ) :: z
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: wX
               REAL, INTENT ( IN ) :: wY
               REAL, INTENT ( IN ) :: wZ

               psi_3d_sho_ani = CMPLX ( &
                  & ( 1.0 / SQRT ( REAL ( 2**nX * factorial ( nX ) ) ) ) * SQRT ( SQRT ( wX / PI ) ) * &
                  & ( 1.0 / SQRT ( REAL ( 2**nY * factorial ( nY ) ) ) ) * SQRT ( SQRT ( wY / PI ) ) * &
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
                  & hermite ( nX , SQRT ( wX ) * ( x - xO ) ) * EXP ( -0.5 * wX * ( x - xO )**2 ) * &
                  & hermite ( nY , SQRT ( wY ) * ( y - yO ) ) * EXP ( -0.5 * wY * ( y - yO )**2 ) * &
                  & hermite ( nZ , SQRT ( wZ ) * ( z - zO ) ) * EXP ( -0.5 * wZ * ( z - zO )**2 ) , 0.0 )

               RETURN

            END FUNCTION

            COMPLEX FUNCTION psi_3d_se_sho_axi ( nR , mL , nZ , x , y , z , xO , yO , zO , wR , wZ )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: nR
               INTEGER, INTENT ( IN ) :: mL
               INTEGER, INTENT ( IN ) :: nZ

               REAL, INTENT ( IN ) :: x
               REAL, INTENT ( IN ) :: y
               REAL, INTENT ( IN ) :: z
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: wR
               REAL, INTENT ( IN ) :: wZ

               psi_3d_sho_axi = CMPLX ( &
                  & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / &
                  & ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * & 
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
                  & SQRT ( ( x - xO )**2 + ( y - yO )**2 )**ABS ( mL ) * &
                  & alaguerre ( nR , ABS ( mL ) , wR * ( ( x - xO )**2 + ( y - yO )**2 ) ) * &
                  & EXP ( -0.5 * wR * ( ( x - xO )**2 + ( y - yO )**2 ) ) * &
                  & hermite ( nZ , SQRT ( wZ ) * ( z - zO ) ) * EXP ( -0.5 * wZ * ( z - zO )**2 ) , 0.0 ) * &
                  & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( y - yO , x - xO ) ) )

               RETURN

            END FUNCTION

            SUBROUTINE psi_se_sho_3d_ani ( nX , nY , nZ , wX , wY , wZ , xO , yO , zO , X , Y , Z , Psi3 )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: nX
               INTEGER, INTENT ( IN ) :: nY
               INTEGER, INTENT ( IN ) :: nZ

               REAL, INTENT ( IN ) :: wX
               REAL, INTENT ( IN ) :: wY
               REAL, INTENT ( IN ) :: wZ
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO

               REAL, DIMENSION ( : ), INTENT ( IN ) :: X
               REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( : ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Psi3

               INTEGER :: j , k , l

               DO l = 1 , SIZE ( Z ) 

                  DO k = 1 , SIZE ( Y ) 

                     DO j = 1 , SIZE ( X ) 

                        Psi3 ( j , k , l ) = CMPLX ( & 
                           & hermite ( nX , SQRT ( wX ) * ( X ( j ) - xO ) ) * EXP ( -0.5 * wX * ( X ( j ) - xO )**2 ) * &
                           & hermite ( nY , SQRT ( wY ) * ( Y ( k ) - yO ) ) * EXP ( -0.5 * wY * ( Y ( k ) - yO )**2 ) * & 
                           & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 ) 

                     END DO

                  END DO

               END DO

               Psi3 = CMPLX ( & 
                  & ( 1.0 / SQRT ( REAL ( 2**nX * factorial ( nX ) ) ) ) * SQRT ( SQRT ( wX / PI ) ) * &
                  & ( 1.0 / SQRT ( REAL ( 2**nY * factorial ( nY ) ) ) ) * SQRT ( SQRT ( wY / PI ) ) * &
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) , 0.0 ) * Psi3

               RETURN

            END SUBROUTINE

            SUBROUTINE psi_se_sho_3d_axi ( nR , mL , nZ , wR , wZ , xO , yO , zO , X , Y , Z , Psi3 )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: nR
               INTEGER, INTENT ( IN ) :: mL
               INTEGER, INTENT ( IN ) :: nZ
       
               REAL, INTENT ( IN ) :: wR
               REAL, INTENT ( IN ) :: wZ
               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
          
               REAL, DIMENSION ( : ), INTENT ( IN ) :: X
               REAL, DIMENSION ( : ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( : ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( : , : , : ), INTENT ( INOUT ) :: Psi3

               INTEGER :: j , k , l

               DO l = 1 , SIZE ( Z )

                  DO k = 1 , SIZE ( Y )

                     DO j = 1 , SIZE ( X )

                        Psi3 ( j , k , l ) = CMPLX ( &
                           & SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 )**ABS ( mL ) * & 
                           & alaguerre ( nR , ABS ( mL ) , wR * ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) ) * &
                           & EXP ( -0.5 * wR * ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) ) * &
                           & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 ) * &
                           & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( Y ( k ) - yO , X ( j ) - xO ) ) ) 

                     END DO

                  END DO

               END DO

               Psi3 = CMPLX ( &
                  & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * &
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) , 0.0 ) * Psi3

               RETURN

            END SUBROUTINE

      END MODULE
! =========================================================================
