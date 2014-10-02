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
!     Tuesday, September 9th, 2014
!
! -------------------------------------------------------------------------

      MODULE VEX

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: GRID, ONLY: nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc

      IMPLICIT NONE
      PRIVATE

      INTEGER, PARAMETER, PRIVATE :: vexUnitIn = 503
      INTEGER, PARAMETER, PRIVATE :: vexUnitInit = 998

      LOGICAL, PUBLIC :: vexRead   = .FALSE.

      INTEGER, PUBLIC :: vexFmtIn  = -1      ! 0 = Unformatted ( Binary ) ; 1 = Formatted ( GPI ) 
      INTEGER, PUBLIC :: vexInit = -1 ! 0 = Linear ; 1 = Simple Harmonic Oscillator ; 2 = Simple Harmonic Oscillator Ring

      REAL, PUBLIC :: vexXo = 0.0
      REAL, PUBLIC :: vexYo = 0.0
      REAL, PUBLIC :: vexZo = 0.0
      REAL, PUBLIC :: vexRo = 0.0
      REAL, PUBLIC :: vexFx = 0.0  
      REAL, PUBLIC :: vexFy = 0.0  
      REAL, PUBLIC :: vexFz = 0.0  
      REAL, PUBLIC :: vexWx = 0.0  
      REAL, PUBLIC :: vexWy = 0.0  
      REAL, PUBLIC :: vexWz = 0.0
      REAL, PUBLIC :: vexWr = 0.0

      PUBLIC :: vex_read_inputs
      PUBLIC :: vex_read_init
      PUBLIC :: vex_compute_init

      PUBLIC :: vex_3d_lin
      PUBLIC :: vex_3d_sho
!      PUBLIC :: vex_3d_sho_rot
      PUBLIC :: vex_3d_shor

      NAMELIST /nmlVexIn/ vexRead , vexFmtIn , vexInit , vexXo , vexYo , vexZo , vexRo , vexFx , vexFy , vexFz , vexWx , vexWy , vexWz , vexWr 

      CONTAINS

         SUBROUTINE vex_read_inputs ( )

            IMPLICIT NONE

            OPEN ( UNIT = vexUnitIn, FILE = 'vex.in' , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
               READ ( UNIT = vexUnitIn , NML = nmlVexIn )
            CLOSE ( UNIT = vexUnitIn , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_read_init ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_compute_init ( X , Y , Z , Vex3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            IF ( vexInit == 0 ) THEN 

               CALL vex_3d_lin ( X , Y , Z , Vex3 )

            ELSE IF ( vexInit == 1 ) THEN 

               CALL vex_3d_sho ( X , Y , Z , Vex3 )

            ELSE IF ( vexInit == 2 ) THEN 

               CALL vex_3d_shor ( X , Y , Z , Vex3 )

            ELSE 

               ! ErvexRor: vexInit not defined.

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_3d_lin ( X , Y , Z , Vex3 )

            IMPLICIT NONE

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

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + vexFx * ( X ( j ) - vexXo ) + vexFy * ( Y ( k ) - vexYo ) + vexFz * ( Z ( l ) - vexZo )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_3d_sho ( X , Y , Z , Vex3 )

            IMPLICIT NONE

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

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( ( vexWx * ( X ( j ) - vexXo ) )**2 + ( vexWy * ( Y ( k ) - vexYo ) )**2 + ( vexWz * ( Z ( l ) - vexZo ) )**2 )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_3d_shor ( X , Y , Z , Vex3 )

            IMPLICIT NONE

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

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( vexWr * ( SQRT ( ( X ( j ) - vexXo )**2 + ( Y ( k ) - vexYo )**2 ) - vexRo )**2 + ( vexWz * ( Z ( l ) - vexZo ) )**2 )

                  END DO
                 
               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
