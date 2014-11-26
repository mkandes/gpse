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
!     Saturday, November 22nd, 2014
!
! -------------------------------------------------------------------------

      MODULE VEX

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: vex_init

      PUBLIC :: vex_3d_lin
      PUBLIC :: vex_3d_sho
      PUBLIC :: vex_3d_shor

      CONTAINS

         SUBROUTINE vex_init ( initVex , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , rO , fX , fY , fZ , wX , wY , wZ , wR , X , Y , Z , Vex3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: initVex
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
            REAl, INTENT ( IN ) :: rO
            REAL, INTENT ( IN ) :: fX
            REAL, INTENT ( IN ) :: fY
            REAL, INTENT ( IN ) :: fZ
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: wR 

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            IF ( initVex == 0 ) THEN 

               CALL vex_3d_lin ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , fX , fY , fZ , X , Y , Z , Vex3 )

            ELSE IF ( initVex == 1 ) THEN 

               CALL vex_3d_sho ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , wX , wY , wZ , X , Y , Z , Vex3 )

            ELSE IF ( initVex == 2 ) THEN 

               CALL vex_3d_shor ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , rO , wR , wZ , X , Y , Z , Vex3 )

            ELSE 

               WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : vex : vex_init : ERROR - initVex not supported.'

            END IF

            RETURN

         END SUBROUTINE

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

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + fX * ( X ( j ) - xO ) + fY * ( Y ( k ) - yO ) + fZ * ( Z ( l ) - zO )

                  END DO

               END DO
!$OMP          END PARALLEL DO

            END DO
!$OMP       END PARALLEL DO

            RETURN

         END SUBROUTINE

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

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( ( wX * ( X ( j ) - xO ) )**2 + ( wY * ( Y ( k ) - yO ) )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

                  END DO

               END DO
!$OMP          END PARALLEL DO

            END DO
!$OMP       END PARALLEL DO

            RETURN

         END SUBROUTINE

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

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb 

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb 

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( wR * ( SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) - rO )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

                  END DO
                 
               END DO
!$OMP          END PARALLEL DO

            END DO
!$OMP       END PARALLEL DO

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
