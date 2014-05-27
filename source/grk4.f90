! ==========================================================================
! NAME
!
!     grk4 [ ] - 
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     physcon is a Fortran module ... 
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
!     Monday, May 26th, 2014
!
! -------------------------------------------------------------------------

      MODULE GRK4

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         IMPLICIT NONE
         PRIVATE

         PUBLIC :: f_gp_3d_rrf_mol_grk4_cd2
         PUBLIC :: f_gp_3d_rrf_mol_grk4_cd4

         CONTAINS

            SUBROUTINE f_gp_3d_rrf_mol_grk4_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

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

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi
               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

               INTEGER :: j , k , l

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        F ( j , k , l ) = &
                           & CMPLX ( 0.5 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 0.5 / dZ**2 ) * Psi ( j , k , l - 1 ) + &
                           & CMPLX ( 0.5 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 0.5 / dY**2 ) * Psi ( j , k - 1 , l ) + &
                           & CMPLX ( 0.5 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 0.5 / dX**2 ) * Psi ( j - 1 , k , l ) - &
                           & CMPLX ( 0.0 , -1.0 / dX**2 - 1.0 / dY**2 - 1.0 / dZ**2 - Vex ( j , k , l ) - &
                           &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                           & CMPLX ( 0.5 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 0.5 / dX**2 ) * Psi ( j + 1 , k , l ) + &
                           & CMPLX ( 0.5 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 0.5 / dY**2 ) * Psi ( j , k + 1 , l ) + &
                           & CMPLX ( 0.5 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 0.5 / dZ**2 ) * Psi ( j , k , l + 1 )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               RETURN

            END SUBROUTINE

            SUBROUTINE f_gp_3d_rrf_mol_grk4_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F)

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

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi
               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

               INTEGER :: j , k , l

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        F ( j , k , l ) = &
                           & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi ( j , k , l - 2 ) + &
                           & CMPLX ( 0.75 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi ( j , k , l - 1 ) + &
                           & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi ( j , k - 2 , l ) + &
                           & CMPLX ( 0.75 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi ( j , k - 1 , l ) + &
                           & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi ( j - 2 , k , l ) + &
                           & CMPLX ( 0.75 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi ( j - 1 , k , l ) + &
                           & CMPLX ( 0.0 , -1.25 * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) - Vex ( j , k , l ) - &
                           &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                           & CMPLX ( 0.75 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi ( j + 1 , k , l ) + &
                           & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi ( j + 2 , k , l ) + &
                           & CMPLX ( 0.75 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi ( j , k + 1 , l ) + &
                           & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi ( j , k + 2 , l ) + &
                           & CMPLX ( 0.75 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi ( j , k , l + 1 ) + &
                           & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi ( j , k , l + 2 )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL


               RETURN

            END SUBROUTINE

      END MODULE
! =========================================================================
