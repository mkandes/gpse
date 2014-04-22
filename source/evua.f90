! ==========================================================================
! NAME
!
!     evua [ evua ] - evua module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     evua is a Fortran module ... 
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
!     Tuesday, April 22nd, 2014
!
! -------------------------------------------------------------------------

      MODULE EVUA

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         IMPLICIT NONE
         PRIVATE

         PUBLIC :: l2_norm_3d_rect
         PUBLIC :: x_3d_rect
         PUBLIC :: y_3d_rect
         PUBLIC :: z_3d_rect
!         PUBLIC :: x2_3d_rect
!         PUBLIC :: y2_3d_rect
!         PUBLIC :: z2_3d_rect
!         PUBLIC :: px_3d_rect_cd2
!         PUBLIC :: px_3d_rect_cd4
!         PUBLIC :: py_3d_rect_cd2
!         PUBLIC :: py_3d_rect_cd4
!         PUBLIC :: pz_3d_rect_cd2
!         PUBLIC :: pz_3d_rect_cd4
!         PUBLIC :: px2_3d_rect_cd2
!         PUBLIC :: px2_3d_rect_cd4
!         PUBLIC :: py2_3d_rect_cd2
!         PUBLIC :: py2_3d_rect_cd4
!         PUBLIC :: pz2_3d_rect_cd2
!         PUBLIC :: pz2_3d_rect_cd4
!         PUBLIC :: lx_3d_rect_cd2
!         PUBLIC :: lx_3d_rect_cd4
!         PUBLIC :: ly_3d_rect_cd2
!         PUBLIC :: ly_3d_rect_cd4
!         PUBLIC :: lz_3d_rect_cd2
!         PUBLIC :: lz_3d_rect_cd4
!         PUBLIC :: lx2_3d_rect_cd2
!         PUBLIC :: lx2_3d_rect_cd4
!         PUBLIC :: ly2_3d_rect_cd2
!         PUBLIC :: ly2_3d_rect_cd4
!         PUBLIC :: lz2_3d_rect_cd2
!         PUBLIC :: lz2_3d_rect_cd4
!         PUBLIC :: vex_3d_rect
!         PUBLIC :: vmf_3d_rect

         CONTAINS

            REAL FUNCTION l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               l2_norm_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        l2_norm_3d_rect = l2_norm_3d_rect + ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               l2_norm_3d_rect = l2_norm_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION x_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Psi3 )

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

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               x_3d_rect = 0.0 

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb  

                        x_3d_rect = x_3d_rect + X ( j ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               x_3d_rect = x_3d_rect * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION y_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Psi3 )

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

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               y_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        y_3d_rect = y_3d_rect + Y ( k ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               y_3d_rect = y_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION z_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Z , Psi3 )

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

               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               z_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        z_3d_rect = z_3d_rect + Z ( l ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               z_3d_rect = z_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

      END MODULE

! =========================================================================
