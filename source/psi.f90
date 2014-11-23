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
!     Saturday, October 18th, 2014
!
! -------------------------------------------------------------------------

      MODULE PSI

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MATH

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: psi_init
      PUBLIC :: psi_normalize
      PUBLIC :: psi_boost

      PRIVATE :: psi_3d_se_sho_ani
      PRIVATE :: psi_3d_se_sho_axi
      PRIVATE :: psi_3d_se_sho_iso

      CONTAINS

         SUBROUTINE psi_init ( initPsi , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nX , nY , nZ , nR , mL , xO , yO , zO , wX , wY , wZ , wR , X , Y , Z , Psi3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: initPsi
            INTEGER, INTENT ( IN ) :: nXa 
            INTEGER, INTENT ( IN ) :: nXb 
            INTEGER, INTENT ( IN ) :: nXbc
            INTEGER, INTENT ( IN ) :: nYa 
            INTEGER, INTENT ( IN ) :: nYb 
            INTEGER, INTENT ( IN ) :: nYbc
            INTEGER, INTENT ( IN ) :: nZa 
            INTEGER, INTENT ( IN ) :: nZb 
            INTEGER, INTENT ( IN ) :: nZbc
            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN ) :: nZ
            INTEGER, INTENT ( IN ) :: nR
            INTEGER, INTENT ( IN ) :: mL

            REAL, INTENT ( IN ) :: xO
            REAL, INTENT ( IN ) :: yO
            REAL, INTENT ( IN ) :: zO
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: wR 

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            IF ( initPsi == 0 ) THEN

               ! ?

            ELSE IF ( initPsi == 1 ) THEN 

               CALL psi_3d_se_sho_ani ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nX , nY , nZ , xO , yO , zO , wX , wY , wZ , X , Y , Z , Psi3 )

            ELSE IF ( initPsi == 2 ) THEN 

               CALL psi_3d_se_sho_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nR , mL , nZ , xO , yO , zO , wR , wZ , X , Y , Z , Psi3 )

            ELSE 

               WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : psi : psi_init : ERROR - initPsi not supported.'

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_3d_se_sho_ani ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nX , nY , nZ , xO , yO , zO , wX , wY , wZ , X , Y , Z , Psi3 )
 
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
            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN ) :: nZ

            REAL, INTENT ( IN ) :: xO
            REAL, INTENT ( IN ) :: yO
            REAL, INTENT ( IN ) :: zO
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ

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
                        & ( 1.0 / SQRT ( REAL ( 2**nX * factorial ( nX ) ) ) ) * SQRT ( SQRT ( wX / PI ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**nY * factorial ( nY ) ) ) ) * SQRT ( SQRT ( wY / PI ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * & 
                        & hermite ( nX , SQRT ( wX ) * ( X ( j ) - xO ) ) * EXP ( -0.5 * wX * ( X ( j ) - xO )**2 ) * & 
                        & hermite ( nY , SQRT ( wY ) * ( Y ( k ) - yO ) ) * EXP ( -0.5 * wY * ( Y ( k ) - yO )**2 ) * & 
                        & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_3d_se_sho_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nR , mL , nZ , xO , yO , zO , wR , wZ , X , Y , Z , Psi3 )

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
            INTEGER, INTENT ( IN ) :: nR
            INTEGER, INTENT ( IN ) :: mL
            INTEGER, INTENT ( IN ) :: nZ

            REAL, INTENT ( IN ) :: xO
            REAL, INTENT ( IN ) :: yO
            REAL, INTENT ( IN ) :: zO
            REAL, INTENT ( IN ) :: wR
            REAL, INTENT ( IN ) :: wZ

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
                        & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / & 
                        & ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * & 
                        & SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 )**ABS ( mL ) * & 
                        & alaguerre ( nR , ABS ( mL ) , wR * ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) ) * & 
                        & EXP ( -0.5 * wR * ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) ) * & 
                        & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 ) * & 
                        & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( Y ( k ) - yO , X ( j ) - xO ) ) ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_3d_se_sho_iso ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_normalize ( ) 

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_boost ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , pX , pY , pZ , X , Y , Z , Psi3 )

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
            REAL, INTENT ( IN ) :: pX
            REAL, INTENT ( IN ) :: pY
            REAL, INTENT ( IN ) :: pZ

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

                     Psi3 ( j , k , l ) = Psi3 ( j , k , l ) * EXP ( CMPLX ( 0.0 , pX * ( X ( j ) - xO ) + pY * ( Y ( k ) - yO ) + pZ * ( Z ( l ) - zO ) ) )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
