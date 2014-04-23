! ==========================================================================
! NAME
!
!     vex [ veks ] - VEX Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     VEX is a Fortran module ... 
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
!     Wednesday, April 23rd, 2014
!
! -------------------------------------------------------------------------

      MODULE VEX

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: vex_3d_lin
      PUBLIC :: vex_3d_sho
      PUBLIC :: vex_3d_shor

      CONTAINS

         SUBROUTINE vex_3d_lin ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , fX , fY , fZ , X , Y , Z , Vex3 )

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
            REAL, INTENT ( IN ) :: fX
            REAL, INTENT ( IN ) :: fY
            REAL, INTENT ( IN ) :: fZ

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + fX * ( X ( j ) - xO ) + fY * ( Y ( k ) - yO ) + fZ * ( Z ( l ) - zO )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

!         REAL FUNCTION vex_3d_lin ( xO , yO , zO , fX , fY , fZ , x , y , z )
!
!         IMPLICIT NONE
!
!         REAL, INTENT ( IN ) :: xO
!         REAL, INTENT ( IN ) :: yO
!         REAL, INTENT ( IN ) :: zO
!         REAL, INTENT ( IN ) :: fX
!         REAL, INTENT ( IN ) :: fY
!         REAL, INTENT ( IN ) :: fZ
!         REAL, INTENT ( IN ) :: x
!         REAL, INTENT ( IN ) :: y
!         REAL, INTENT ( IN ) :: z
!
!         vex_3d_lin = fX * ( x - xO ) + fY * ( y - yO ) + fZ * ( z - zO ) 
!
!         RETURN
!
!         END FUNCTION

         SUBROUTINE vex_3d_sho ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , wX , wY , wZ , X , Y , Z , Vex3 )

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
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( ( wX * ( X ( j ) - xO ) )**2 + ( wY * ( Y ( k ) - yO ) )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

!         REAL FUNCTION vex_3d_sho ( xO , yO , zO , wX , wY , wZ , x , y , z )
!
!         IMPLICIT NONE
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
!         vex_3d_sho = 0.5 * ( ( wX * ( x - xO ) )**2 + ( wY * ( y - yO ) )**2 + ( wZ * ( z - zO ) )**2 )
!
!         RETURN
!
!         END FUNCTION

         SUBROUTINE vex_3d_shor ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , rO , wR , wZ , X , Y , Z , Vex3 )

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
            REAL, INTENT ( IN ) :: rO
            REAL, INTENT ( IN ) :: wR
            REAL, INTENT ( IN ) :: wZ

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l 

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb 

               DO k = nYa , nYb 

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( wR * ( SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) - rO )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

                  END DO
                 
               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

!         REAL FUNCTION vex_3d_shor ( xO , yO , zO , rO , wR , wZ , x , y , z )
!
!         IMPLICIT NONE
!
!         REAL, INTENT ( IN ) :: xO
!         REAL, INTENT ( IN ) :: yO
!         REAL, INTENT ( IN ) :: zO
!         REAL, INTENT ( IN ) :: rO
!         REAL, INTENT ( IN ) :: wR
!         REAL, INTENT ( IN ) :: wZ
!         REAL, INTENT ( IN ) :: x
!         REAL, INTENT ( IN ) :: y
!         REAL, INTENT ( IN ) :: z
!
!         vex_3d_shor = 0.5 * ( wR * ( SQRT ( ( x - xO )**2 + ( y - yO )**2 ) - rO )**2 + ( wZ * ( z - zO ) )**2 )
!
!         RETURN
!
!         END FUNCTION

      END MODULE

! =========================================================================
