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
!     Tuesday, November 25th, 2014
!
! -------------------------------------------------------------------------

      MODULE EVUA

         USE, INTRINSIC :: ISO_FORTRAN_ENV
         USE            :: MATH

         IMPLICIT NONE
         PRIVATE

         REAL, PUBLIC :: l2Norm = 0.0
         REAL, PUBLIC :: avgX = 0.0
         REAL, PUBLIC :: avgX2 = 0.0
         REAL, PUBLIC :: avgX2COM = 0.0
         REAL, PUBLIC :: avgY = 0.0
         REAL, PUBLIC :: avgY2 = 0.0
         REAL, PUBLIC :: avgY2COM = 0.0
         REAL, PUBLIC :: avgZ = 0.0
         REAL, PUBLIC :: avgZ2 = 0.0
         REAL, PUBLIC :: avgZ2COM = 0.0
         REAL, PUBLIC :: avgRxy = 0.0
         REAL, PUBLIC :: avgR2xy = 0.0
         REAL, PUBLIC :: avgRxyz = 0.0
         REAL, PUBLIC :: avgR2xyz = 0.0
         REAL, PUBLIC :: avgPx = 0.0
         REAL, PUBLIC :: avgPx2 = 0.0
         REAL, PUBLIC :: avgPy = 0.0
         REAL, PUBLIC :: avgPy2 = 0.0
         REAL, PUBLIC :: avgPz = 0.0
         REAL, PUBLIC :: avgPz2 = 0.0
         REAL, PUBLIC :: avgLx = 0.0
         REAL, PUBLIC :: avgLxCOM = 0.0
         REAL, PUBLIC :: avgLx2 = 0.0
         REAL, PUBLIC :: avgLx2COM = 0.0
         REAL, PUBLIC :: avgLy = 0.0
         REAL, PUBLIC :: avgLyCOM = 0.0 
         REAL, PUBLIC :: avgLy2 = 0.0
         REAL, PUBLIC :: avgLy2COM = 0.0 
         REAL, PUBLIC :: avgLz = 0.0
         REAL, PUBLIC :: avgLzCOM = 0.0
         REAL, PUBLIC :: avgLz2 = 0.0
         REAL, PUBLIC :: avgLz2COM = 0.0
         REAL, PUBLIC :: avgFx = 0.0
         REAL, PUBLIC :: avgFy = 0.0
         REAL, PUBLIC :: avgFz = 0.0
         REAL, PUBLIC :: avgTauX = 0.0
         REAL, PUBLIC :: avgTauXCOM = 0.0
         REAL, PUBLIC :: avgTauY = 0.0
         REAL, PUBLIC :: avgTauYCOM = 0.0
         REAL, PUBLIC :: avgTauZ = 0.0
         REAL, PUBLIC :: avgTauZCOM = 0.0
         REAL, PUBLIC :: avgVex = 0.0
         REAL, PUBLIC :: avgVmf = 0.0
         REAL, PUBLIC :: avgIxx = 0.0
         REAL, PUBLIC :: avgIyy = 0.0
         REAL, PUBLIC :: avgIzz = 0.0
         REAL, PUBLIC :: avgIxy = 0.0
         REAL, PUBLIC :: avgIyz = 0.0
         REAL, PUBLIC :: avgIxz = 0.0
         REAL, PUBLIC :: avgIxxCOM = 0.0
         REAL, PUBLIC :: avgIyyCOM = 0.0
         REAL, PUBLIC :: avgIzzCOM = 0.0
         REAL, PUBLIC :: avgIxyCOM = 0.0
         REAL, PUBLIC :: avgIyzCOM = 0.0
         REAL, PUBLIC :: avgIxzCOM = 0.0
         REAL, PUBLIC :: avgTx = 0.0
         REAL, PUBLIC :: avgTy = 0.0
         REAL, PUBLIC :: avgTz = 0.0
         REAL, PUBLIC :: avgE = 0.0
         REAL, PUBLIC :: avgL2 = 0.0
         REAL, PUBLIC :: avgMu = 0.0
         REAL, PUBLIC :: sigX = 0.0
         REAL, PUBLIC :: sigY = 0.0
         REAL, PUBLIC :: sigZ = 0.0
         REAL, PUBLIC :: sigRxy = 0.0
         REAL, PUBLIC :: sigPx = 0.0
         REAL, PUBLIC :: sigPy = 0.0
         REAL, PUBLIC :: sigPz = 0.0
         REAL, PUBLIC :: sigLx = 0.0
         REAL, PUBLIC :: sigLy = 0.0
         REAL, PUBLIC :: sigLz = 0.0

         PUBLIC :: l2_norm_3d_rect
         PUBLIC :: x_3d_rect
         PUBLIC :: y_3d_rect
         PUBLIC :: z_3d_rect
         PUBLIC :: r_xy_3d_rect
         PUBLIC :: r_xyz_3d_rect
         PUBLIC :: x2_3d_rect
         PUBLIC :: y2_3d_rect
         PUBLIC :: z2_3d_rect
         PUBLIC :: r2_xy_3d_rect
         PUBLIC :: r2_xyz_3d_rect
         PUBLIC :: px_3d_rect_cd2
         PUBLIC :: px_3d_rect_cd4
         PUBLIC :: py_3d_rect_cd2
         PUBLIC :: py_3d_rect_cd4
         PUBLIC :: pz_3d_rect_cd2
         PUBLIC :: pz_3d_rect_cd4
!         PUBLIC :: pr_3d_rect_cd2
!         PUBLIC :: pr_3d_rect_cd4
         PUBLIC :: px2_3d_rect_cd2
         PUBLIC :: px2_3d_rect_cd4
         PUBLIC :: py2_3d_rect_cd2
         PUBLIC :: py2_3d_rect_cd4
         PUBLIC :: pz2_3d_rect_cd2
         PUBLIC :: pz2_3d_rect_cd4
!         PUBLIC :: pr2_3d_rect_cd2
!         PUBLIC :: pr2_3d_rect_cd4
         PUBLIC :: lx_3d_rect_cd2
         PUBLIC :: lx_3d_rect_cd4
         PUBLIC :: ly_3d_rect_cd2
         PUBLIC :: ly_3d_rect_cd4
         PUBLIC :: lz_3d_rect_cd2
         PUBLIC :: lz_3d_rect_cd4
         PUBLIC :: lx2_3d_rect_cd2
         PUBLIC :: lx2_3d_rect_cd4
         PUBLIC :: ly2_3d_rect_cd2
         PUBLIC :: ly2_3d_rect_cd4
         PUBLIC :: lz2_3d_rect_cd2
         PUBLIC :: lz2_3d_rect_cd4
         PUBLIC :: fx_3d_rect_cd2
         PUBLIC :: fx_3d_rect_cd4
         PUBLIC :: fy_3d_rect_cd2
         PUBLIC :: fy_3d_rect_cd4
         PUBLIC :: fz_3d_rect_cd2
         PUBLIC :: fz_3d_rect_cd4
!         PUBLIC :: fr_3d_rect_cd2
!         PUBLIC :: fr_3d_rect_cd4
         PUBLIC :: taux_3d_rect_cd2
         PUBLIC :: taux_3d_rect_cd4
         PUBLIC :: tauy_3d_rect_cd2
         PUBLIC :: tauy_3d_rect_cd4
         PUBLIC :: tauz_3d_rect_cd2
         PUBLIC :: tauz_3d_rect_cd4
         PUBLIC :: ixx_3d_rect
         PUBLIC :: iyy_3d_rect
         PUBLIC :: izz_3d_rect
         PUBLIC :: ixy_3d_rect
         PUBLIC :: iyz_3d_rect
         PUBLIC :: ixz_3d_rect
         PUBLIC :: vex_3d_rect
         PUBLIC :: vmf_3d_rect

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

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : l2_norm_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : l2_norm_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        l2_norm_3d_rect = l2_norm_3d_rect + ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               l2_norm_3d_rect = l2_norm_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION x_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               x_3d_rect = 0.0 

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : x_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : x_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb  

                        x_3d_rect = x_3d_rect + ( X ( j ) - xO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               x_3d_rect = x_3d_rect * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION y_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               y_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : y_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : y_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        y_3d_rect = y_3d_rect + ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               y_3d_rect = y_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION z_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               z_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : z_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : z_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        z_3d_rect = z_3d_rect + ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               z_3d_rect = z_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               r_xy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r_xy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r_xy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r_xy_3d_rect = r_xy_3d_rect + SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               r_xy_3d_rect = r_xy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r_xyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , X , Y , Z , Psi3 )

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

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               r_xyz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r_xyz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r_xyz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r_xyz_3d_rect = r_xyz_3d_rect + SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               r_xyz_3d_rect = r_xyz_3d_rect * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l 

               x2_3d_rect = 0.0 

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : x2_3d_rect )
               DO l = nZa , nZb 

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : x2_3d_rect )
                  DO k = nYa , nYb 

                     DO j = nXa , nXb  

                        x2_3d_rect = x2_3d_rect + ( X ( j ) - xO )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               x2_3d_rect = x2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               y2_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : y2_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : y2_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        y2_3d_rect = y2_3d_rect + ( Y ( k ) - yO )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               y2_3d_rect = y2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               z2_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : z2_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : z2_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        z2_3d_rect = z2_3d_rect + ( Z ( l ) - zO )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               z2_3d_rect = z2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r2_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO ,  dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               r2_xy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r2_xy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r2_xy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r2_xy_3d_rect = r2_xy_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               r2_xy_3d_rect = r2_xy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r2_xyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , X , Y , Z , Psi3 )

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

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               r2_xyz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r2_xyz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : r2_xyz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r2_xyz_3d_rect = r2_xyz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               r2_xyz_3d_rect = r2_xyz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION px_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )

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

               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               px_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px_3d_rect_cd2 = px_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               px_3d_rect_cd2 = -0.5 * px_3d_rect_cd2 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION px_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )

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

               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               px_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px_3d_rect_cd4 ) 
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px_3d_rect_cd4 ) 
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px_3d_rect_cd4 = px_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               px_3d_rect_cd4 = -0.08333333333333333 * px_3d_rect_cd4 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION py_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )

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
               REAL, INTENT ( IN ) :: dZ

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               py_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py_3d_rect_cd2 = py_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               py_3d_rect_cd2 = -0.5 * py_3d_rect_cd2 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION py_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )

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
               REAL, INTENT ( IN ) :: dZ

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               py_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py_3d_rect_cd4 = py_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k -1 , l ) ) + Psi3 ( j , k - 2 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               py_3d_rect_cd4 = -0.08333333333333333 * py_3d_rect_cd4 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION pz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )

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

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               pz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz_3d_rect_cd2 ) 
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz_3d_rect_cd2 = pz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               pz_3d_rect_cd2 = -0.5 * pz_3d_rect_cd2 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION pz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )

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

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l 

               pz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz_3d_rect_cd4 = pz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               pz_3d_rect_cd4 = -0.08333333333333333 * pz_3d_rect_cd4 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION px2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               px2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px2_3d_rect_cd2 ) 
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px2_3d_rect_cd2 = px2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               px2_3d_rect_cd2 = -px2_3d_rect_cd2 * dY * dZ / dX

               RETURN

            END FUNCTION

            REAL FUNCTION px2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               px2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : px2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px2_3d_rect_cd4 = px2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - Psi3 ( j - 2 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               px2_3d_rect_cd4 = -0.08333333333333333 * px2_3d_rect_cd4 * dY * dZ / dX

               RETURN

            END FUNCTION

            REAL FUNCTION py2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               py2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py2_3d_rect_cd2 = py2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               py2_3d_rect_cd2 = -py2_3d_rect_cd2 * dX * dZ / dY

               RETURN

            END FUNCTION

            REAL FUNCTION py2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               py2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : py2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py2_3d_rect_cd4 = py2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - Psi3 ( j , k - 2 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               py2_3d_rect_cd4 = -0.08333333333333333 * py2_3d_rect_cd4 * dX * dZ / dY

               RETURN
       
            END FUNCTION

            REAL FUNCTION pz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               pz2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz2_3d_rect_cd2 = pz2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) ) 

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               pz2_3d_rect_cd2 = -pz2_3d_rect_cd2 * dX * dY / dZ

               RETURN

            END FUNCTION

            REAL FUNCTION pz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               pz2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : pz2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz2_3d_rect_cd4 = pz2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - Psi3 ( j , k , l - 2 ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               pz2_3d_rect_cd4 = -0.08333333333333333 * pz2_3d_rect_cd4 * dX * dY / dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               lx_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx_3d_rect_cd2 = lx_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Y ( k ) - yO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - CMPLX ( ( Z ( l ) - zO ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lx_3d_rect_cd2 = -0.5 * lx_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               lx_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx_3d_rect_cd4 = lx_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Y ( k ) - yO ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) - CMPLX ( ( Z ( l ) - zO ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lx_3d_rect_cd4 = -0.08333333333333333 * lx_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly_3d_rect_cd2 = ly_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Z ( l ) - zO ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - CMPLX ( ( X ( j ) - xO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ly_3d_rect_cd2 = -0.5 * ly_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly_3d_rect_cd4 = ly_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Z ( l ) - zO ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) - CMPLX ( ( X ( j ) - xO ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ly_3d_rect_cd4 = -0.08333333333333333 * ly_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz_3d_rect_cd2 = lz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( X ( j ) - xO ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - CMPLX ( ( Y ( k ) - yO ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lz_3d_rect_cd2 = -0.5 * lz_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz_3d_rect_cd4 = lz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( X ( j ) - xO ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) - CMPLX ( ( Y ( k ) - yO ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lz_3d_rect_cd4 = -0.08333333333333333 * lz_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               lx2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx2_3d_rect_cd2 = lx2_3d_rect_cd2 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( ( Y ( k ) - yO ) / dZ )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * ( Y ( k ) - yO ) / dY , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                           &    CMPLX ( 0.5 * ( Y ( k ) - yO ) * ( Z ( l ) - zO ) / ( dY * dZ ) , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l + 1 ) - Psi3 ( j , k + 1 , l - 1 ) - &
                           &       Psi3 ( j , k - 1 , l + 1 ) + Psi3 ( j , k - 1 , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * ( Z ( l ) - zO ) / dZ , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + &
                           &    CMPLX ( ( ( Z ( l ) - zO ) / dY )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lx2_3d_rect_cd2 = -lx2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lx2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lx2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx2_3d_rect_cd4 = lx2_3d_rect_cd4 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( ( Y ( k ) - yO ) / dZ )**2 , 0.0 ) * ( -Psi3 ( j , k , l - 2 ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l + 2 ) ) + &
                           &    CMPLX ( ( Y ( k ) - yO ) / dY , 0.0 ) * ( Psi3 ( j , k + 2 , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                           &       Psi3 ( j , k - 2 , l ) ) - &
                           &    CMPLX ( ( Y ( k ) - yO ) * ( Z ( l ) - zO ) / ( 6.0 * dZ * dY ) , 0.0 ) * ( Psi3 ( j , k - 2 , l - 2 ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l - 2 ) - Psi3 ( j , k - 1 , l - 2 ) ) - &
                           &       Psi3 ( j , k + 2 , l - 2 ) - CMPLX ( 8.0 , 0.0 ) * Psi3 ( j , k - 2 , l - 1 ) - &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l - 1 ) - Psi3 ( j , k - 1 , l - 1 ) ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 2 , l - 1 ) + Psi3 ( j , k - 2 , l + 1 ) ) + &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l + 1 ) - Psi3 ( j , k - 1 , l + 1 ) ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * Psi3 ( j , k + 2 , l + 1 ) - Psi3 ( j , k - 2 , l + 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l + 2 ) - Psi3 ( j , k - 1 , l + 2 ) ) + & 
                           &       Psi3 ( j , k + 2 , l + 2 ) ) + &
                           &    CMPLX ( ( Z ( l ) - zO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                           &       Psi3 ( j , k , l - 2 ) ) + &
                           &    CMPLX ( ( ( Z ( l ) - zO ) / dY )**2 , 0.0 ) * ( -Psi3 ( j , k - 2 , l ) + &
                           &    CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &    CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - Psi3 ( j , k + 2 , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lx2_3d_rect_cd4 = -0.08333333333 * lx2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly2_3d_rect_cd2 = ly2_3d_rect_cd2 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( ( Z ( l ) - zO ) / dX )**2 , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) - &
                           &    CMPLX ( 0.5 * ( Z ( l ) - zO ) / dZ , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * ( X ( j ) - xO ) * ( Z ( l ) - zO ) / ( dX * dZ ) , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l + 1 ) - Psi3 ( j - 1 , k , l + 1 ) - &
                           &       Psi3 ( j + 1 , k , l - 1 ) + Psi3 ( j - 1 , k , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * ( X ( j ) - xO ) / dX , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + &
                           &    CMPLX ( ( ( X ( j ) - xO ) / dZ )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) ) )

                     END DO
       
                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ly2_3d_rect_cd2 = -ly2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ly2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly2_3d_rect_cd4 = ly2_3d_rect_cd4 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( ( Z ( l ) - zO ) / dX )**2 , 0.0 ) * ( -Psi3 ( j - 2 , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - CMPLX ( 30.0, 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - Psi3 ( j + 2 , k , l ) ) + &
                           &    CMPLX ( ( Z ( l ) - zO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                           &       Psi3 ( j , k , l - 2 ) ) - &
                           &    CMPLX ( ( Z ( l ) - zO ) * ( X ( j ) - xO ) / ( 6.0 * dX * dZ ) , 0.0 ) * ( Psi3 ( j - 2 , k , l - 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l - 2 ) - Psi3 ( j + 1 , k , l - 2 ) ) - &
                           &       Psi3 ( j + 2 , k , l - 2 ) - CMPLX ( 8.0 , 0.0 ) * Psi3 ( j - 2 , k , l - 1 ) + &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l - 1 ) - Psi3 ( j + 1 , k , l - 1 ) ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 2 , k , l - 1 ) + Psi3 ( j - 2 , k , l + 1 ) ) - &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l + 1 ) - Psi3 ( j + 1 , k , l + 1 ) ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * Psi3 ( j + 2 , k , l + 1 ) - Psi3 ( j - 2 , k , l + 2 ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l + 2 ) - Psi3 ( j + 1 , k , l + 2 ) ) + &
                           &       Psi3 ( j + 2 , k , l + 2 ) ) + &
                           &    CMPLX ( ( X ( j ) - xO ) / dX , 0.0 ) * ( Psi3 ( j + 2 , k , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                           &       Psi3 ( j - 2 , k , l ) ) + &
                           &    CMPLX ( ( ( X ( j ) - xO ) / dZ )**2 , 0.0 ) * ( -Psi3 ( j , k , l - 2 ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l + 2 ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ly2_3d_rect_cd4 = -0.08333333333 * ly2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz2_3d_rect_cd2 = lz2_3d_rect_cd2 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( ( X ( j ) - xO ) / dY )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) - &
                           &    CMPLX ( 0.5 *  ( X ( j ) - xO ) / dX , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                           &    CMPLX ( 0.5 * ( X ( j ) - xO ) * ( Y ( k ) - yO ) / ( dX * dY ) , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k + 1 , l ) - Psi3 ( j + 1 , k - 1 , l ) - &
                           &       Psi3 ( j - 1 , k + 1 , l ) + Psi3 ( j - 1 , k - 1 , l ) ) - &
                           &    CMPLX ( 0.5 * ( Y ( k ) - yO ) / dY , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + &
                           &    CMPLX ( ( ( Y ( k ) - yO ) / dX )**2 , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lz2_3d_rect_cd2 = -lz2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : lz2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz2_3d_rect_cd4 = lz2_3d_rect_cd4 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( ( X ( j ) - xO ) / dY )**2 , 0.0 ) * ( -Psi3 ( j , k - 2 , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - Psi3 ( j , k + 2 , l ) ) + &
                           &    CMPLX ( ( X ( j ) - xO ) / dX , 0.0 ) * ( Psi3 ( j + 2 , k , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                           &       Psi3 ( j - 2 , k , l ) ) - &
                           &    CMPLX ( ( X ( j ) - xO ) * ( Y ( k ) - yO ) / ( 6.0 * dY * dX ) ) * ( Psi3 ( j - 2 , k - 2 , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 2 , k - 1 , l ) - Psi3 ( j - 2 , k + 1 , l ) ) - &
                           &       Psi3 ( j - 2 , k + 2 , l ) - CMPLX ( 8.0 , 0.0 ) * Psi3 ( j - 1 , k - 2 , l ) + &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j - 1 , k - 1 , l ) - Psi3 ( j - 1 , k + 1 , l ) ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 1 , k + 2 , l ) + Psi3 ( j + 1 , k - 2 , l ) ) - &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j + 1 , k - 1 , l ) - Psi3 ( j + 1 , k + 1 , l ) ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * Psi3 ( j + 1 , k + 2 , l ) - Psi3 ( j + 2 , k - 2 , l ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 2 , k - 1 , l ) - Psi3 ( j + 2 , k + 1 , l ) ) + &
                           &       Psi3 ( j + 2 , k + 2 , l ) ) + &
                           &    CMPLX ( ( Y ( k ) - yO ) / dY , 0.0 ) * ( Psi3 ( j , k + 2 , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                           &       Psi3 ( j , k - 2 , l ) ) + &
                           &    CMPLX ( ( ( Y ( k ) - yO ) / dX )**2 , 0.0 ) * ( -Psi3 ( j - 2 , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - Psi3 ( j + 2 , k , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               lz2_3d_rect_cd4 = -0.08333333333 * lz2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )

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

               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fx_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fx_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fx_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fx_3d_rect_cd2 = fx_3d_rect_cd2 + ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               fx_3d_rect_cd2 = 0.5 * fx_3d_rect_cd2 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )

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

               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fx_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fx_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fx_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fx_3d_rect_cd4 = fx_3d_rect_cd4 + ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               fx_3d_rect_cd4 = 0.08333333333333333 * fx_3d_rect_cd4 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )

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
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fy_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fy_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fy_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fy_3d_rect_cd2 = fy_3d_rect_cd2 + ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               fy_3d_rect_cd2 = 0.5 * fy_3d_rect_cd2 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )

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
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fy_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fy_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fy_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fy_3d_rect_cd4 = fy_3d_rect_cd4 + ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               fy_3d_rect_cd4 = 0.08333333333333333 * fy_3d_rect_cd4 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )

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

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fz_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fz_3d_rect_cd2 = fz_3d_rect_cd2 + ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               fz_3d_rect_cd2 = 0.5 * fz_3d_rect_cd2 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION fz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )

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

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : fz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fz_3d_rect_cd4 = fz_3d_rect_cd4 + ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               fz_3d_rect_cd4 = 0.08333333333333333 * fz_3d_rect_cd4 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION taux_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Vex3 , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               taux_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : taux_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : taux_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        taux_3d_rect_cd2 = taux_3d_rect_cd2 + ( ( Y ( k ) - yO ) * ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) / dZ - ( Z ( l ) - zO ) * ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               taux_3d_rect_cd2 = 0.5 * taux_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION taux_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Vex3 , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               taux_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : taux_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : taux_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        taux_3d_rect_cd4 = taux_3d_rect_cd4 + ( ( Y ( k ) - yO ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ - ( Z ( l ) - zO ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               taux_3d_rect_cd4 = 0.08333333333333333 * taux_3d_rect_cd4 * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION tauy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Vex3 , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauy_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauy_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauy_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauy_3d_rect_cd2 = tauy_3d_rect_cd2 + ( ( Z ( l ) - zO ) * ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) / dX - ( X ( j ) - xO ) * ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               tauy_3d_rect_cd2 = 0.5 * tauy_3d_rect_cd2 * dX * dY * dZ

            END FUNCTION

            REAL FUNCTION tauy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Vex3 , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauy_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauy_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauy_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauy_3d_rect_cd4 = tauy_3d_rect_cd4 + ( ( Z ( l ) - zO ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX - ( X ( j ) - xO ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               tauy_3d_rect_cd4 = 0.08333333333333333 * tauy_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION tauz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Vex3 , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauz_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauz_3d_rect_cd2 = tauz_3d_rect_cd2 + ( ( X ( j ) - xO ) * ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) / dY - ( Y ( k ) - yO ) * ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               tauz_3d_rect_cd2 = 0.5 * tauz_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION tauz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Vex3 , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : tauz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauz_3d_rect_cd4 = tauz_3d_rect_cd4 + ( ( X ( j ) - xO ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2, l ) ) / dY - ( Y ( k ) - yO ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               tauz_3d_rect_cd4 = 0.08333333333333333 * tauz_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ixx_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ixx_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ixx_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ixx_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ixx_3d_rect = ixx_3d_rect + ( ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ixx_3d_rect = ixx_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION iyy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               iyy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : iyy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : iyy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        iyy_3d_rect = iyy_3d_rect + ( ( X ( j ) - xO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               iyy_3d_rect = iyy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION izz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               izz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : izz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : izz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        izz_3d_rect = izz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               izz_3d_rect = izz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ixy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ixy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ixy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ixy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ixy_3d_rect = ixy_3d_rect + ( X ( j ) - xO ) * ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ixy_3d_rect = -ixy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION iyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               iyz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : iyz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : iyz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        iyz_3d_rect = iyz_3d_rect + ( Y ( k ) - yO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               iyz_3d_rect = -iyz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ixz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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
               REAL, INTENT ( IN ) :: zO
               REAL, INTENT ( IN ) :: dX
               REAL, INTENT ( IN ) :: dY
               REAL, INTENT ( IN ) :: dZ

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ixz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ixz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : ixz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ixz_3d_rect = ixz_3d_rect + ( X ( j ) - xO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               ixz_3d_rect = -ixz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION vex_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Vex3 , Psi3 )

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

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               vex_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : vex_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : vex_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        vex_3d_rect = vex_3d_rect + Vex3 ( j , k , l ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               vex_3d_rect = vex_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION vmf_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , gS , Psi3 )

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
               REAL, INTENT ( IN ) :: gS

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               vmf_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : vmf_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : vmf_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        vmf_3d_rect = vmf_3d_rect + ABS ( Psi3 ( j , k , l ) )**4

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               vmf_3d_rect = 0.5 * gS * vmf_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

      END MODULE

! =========================================================================
