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
!     Thursday, October 2nd, 2014
!
! -------------------------------------------------------------------------

      MODULE EVUA

         USE, INTRINSIC :: ISO_FORTRAN_ENV
         USE            :: GRID, ONLY: nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ
         USE            :: MATH

         IMPLICIT NONE
         PRIVATE

         REAL, PUBLIC :: l2Norm = 0.0
         REAL, PUBLIC :: avgX = 0.0
         REAL, PUBLIC :: avgX2 = 0.0
         REAL, PUBLIC :: avgY = 0.0
         REAL, PUBLIC :: avgY2 = 0.0
         REAL, PUBLIC :: avgZ = 0.0
         REAL, PUBLIC :: avgZ2 = 0.0
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
         REAL, PUBLIC :: avgLx2 = 0.0
         REAL, PUBLIC :: avgLy = 0.0
         REAL, PUBLIC :: avgLy2 = 0.0
         REAL, PUBLIC :: avgLz = 0.0
         REAL, PUBLIC :: avgLz2 = 0.0
         REAL, PUBLIC :: avgFx = 0.0
         REAL, PUBLIC :: avgFy = 0.0
         REAL, PUBLIC :: avgFz = 0.0
         REAL, PUBLIC :: avgTauX = 0.0
         REAL, PUBLIC :: avgTauY = 0.0
         REAL, PUBLIC :: avgTauZ = 0.0
         REAL, PUBLIC :: avgVex = 0.0
         REAL, PUBLIC :: avgVmf = 0.0
         REAL, PUBLIC :: avgIxxCM = 0.0
         REAL, PUBLIC :: avgIyyCM = 0.0
         REAL, PUBLIC :: avgIzzCM = 0.0
         REAL, PUBLIC :: avgIxyCM = 0.0
         REAL, PUBLIC :: avgIyzCM = 0.0
         REAL, PUBLIC :: avgIxzCM = 0.0
         REAL, PUBLIC :: avgIxxO = 0.0
         REAL, PUBLIC :: avgIyyO = 0.0
         REAL, PUBLIC :: avgIzzO = 0.0
         REAL, PUBLIC :: avgIxyO = 0.0
         REAL, PUBLIC :: avgIyzO = 0.0
         REAL, PUBLIC :: avgIxzO = 0.0
         REAL, PUBLIC :: avgTx = 0.0
         REAL, PUBLIC :: avgTy = 0.0
         REAL, PUBLIC :: avgTz = 0.0
         REAL, PUBLIC :: avgE = 0.0
         REAL, PUBLIC :: avgL2 = 0.0
         REAL, PUBLIC :: mu = 0.0
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

            REAL FUNCTION l2_norm_3d_rect ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               l2_norm_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : l2_norm_3d_rect )
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

            REAL FUNCTION x_3d_rect ( X , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               x_3d_rect = 0.0 

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : x_3d_rect )
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

            REAL FUNCTION y_3d_rect ( Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               y_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : y_3d_rect )
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

            REAL FUNCTION z_3d_rect ( Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               z_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : z_3d_rect )
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

            REAL FUNCTION r_xy_3d_rect ( X , Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               r_xy_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : r_xy_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r_xy_3d_rect = r_xy_3d_rect + SQRT ( X ( j )**2 + Y ( k )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               r_xy_3d_rect = r_xy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r_xyz_3d_rect ( X , Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               r_xyz_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : r_xyz_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r_xyz_3d_rect = r_xyz_3d_rect + SQRT ( X ( j )**2 + Y ( k )**2 + Z ( l )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               r_xyz_3d_rect = r_xyz_3d_rect * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION x2_3d_rect ( X , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l 

               x2_3d_rect = 0.0 

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : x2_3d_rect )
               DO l = nZa , nZb 

                  DO k = nYa , nYb 

                     DO j = nXa , nXb  

                        x2_3d_rect = x2_3d_rect + X ( j )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               x2_3d_rect = x2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION y2_3d_rect ( Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               y2_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : y2_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        y2_3d_rect = y2_3d_rect + Y ( k )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               y2_3d_rect = y2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION z2_3d_rect ( Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               z2_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : z2_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        z2_3d_rect = z2_3d_rect + Z ( l )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               z2_3d_rect = z2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r2_xy_3d_rect ( X , Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               r2_xy_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : r2_xy_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r2_xy_3d_rect = r2_xy_3d_rect + ( X ( j )**2 + Y ( k )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               r2_xy_3d_rect = r2_xy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION r2_xyz_3d_rect ( X , Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               r2_xyz_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : r2_xyz_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        r2_xyz_3d_rect = r2_xyz_3d_rect + ( X ( j )**2 + Y ( k )**2 + Z ( l )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               r2_xyz_3d_rect = r2_xyz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION


            REAL FUNCTION px_3d_rect_cd2 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               px_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : px_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px_3d_rect_cd2 = px_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               px_3d_rect_cd2 = -0.5 * px_3d_rect_cd2 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION px_3d_rect_cd4 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               px_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : px_3d_rect_cd4 ) 
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px_3d_rect_cd4 = px_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               px_3d_rect_cd4 = -0.08333333333333333 * px_3d_rect_cd4 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION py_3d_rect_cd2 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               py_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : py_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py_3d_rect_cd2 = py_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               py_3d_rect_cd2 = -0.5 * py_3d_rect_cd2 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION py_3d_rect_cd4 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               py_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : py_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py_3d_rect_cd4 = py_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k -1 , l ) ) + Psi3 ( j , k - 2 , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               py_3d_rect_cd4 = -0.08333333333333333 * py_3d_rect_cd4 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION pz_3d_rect_cd2 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               pz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : pz_3d_rect_cd2 ) 
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz_3d_rect_cd2 = pz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               pz_3d_rect_cd2 = -0.5 * pz_3d_rect_cd2 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION pz_3d_rect_cd4 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l 

               pz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : pz_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz_3d_rect_cd4 = pz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               pz_3d_rect_cd4 = -0.08333333333333333 * pz_3d_rect_cd4 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION px2_3d_rect_cd2 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               px2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : px2_3d_rect_cd2 ) 
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px2_3d_rect_cd2 = px2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               px2_3d_rect_cd2 = -px2_3d_rect_cd2 * dY * dZ / dX

               RETURN

            END FUNCTION

            REAL FUNCTION px2_3d_rect_cd4 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               px2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : px2_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        px2_3d_rect_cd4 = px2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - Psi3 ( j - 2 , k , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               px2_3d_rect_cd4 = -0.08333333333333333 * px2_3d_rect_cd4 * dY * dZ / dX

               RETURN

            END FUNCTION

            REAL FUNCTION py2_3d_rect_cd2 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               py2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : py2_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py2_3d_rect_cd2 = py2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               py2_3d_rect_cd2 = -py2_3d_rect_cd2 * dX * dZ / dY

               RETURN

            END FUNCTION

            REAL FUNCTION py2_3d_rect_cd4 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               py2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : py2_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        py2_3d_rect_cd4 = py2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - Psi3 ( j , k - 2 , l ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               py2_3d_rect_cd4 = -0.08333333333333333 * py2_3d_rect_cd4 * dX * dZ / dY

               RETURN
       
            END FUNCTION

            REAL FUNCTION pz2_3d_rect_cd2 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               pz2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : pz2_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz2_3d_rect_cd2 = pz2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) ) 

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL 

               pz2_3d_rect_cd2 = -pz2_3d_rect_cd2 * dX * dY / dZ

               RETURN

            END FUNCTION

            REAL FUNCTION pz2_3d_rect_cd4 ( Psi3 )

               IMPLICIT NONE

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               pz2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : pz2_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        pz2_3d_rect_cd4 = pz2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - Psi3 ( j , k , l - 2 ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               pz2_3d_rect_cd4 = -0.08333333333333333 * pz2_3d_rect_cd4 * dX * dY / dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx_3d_rect_cd2 ( Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               lx_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lx_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx_3d_rect_cd2 = lx_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( Y ( k ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - CMPLX ( Z ( l ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lx_3d_rect_cd2 = -0.5 * lx_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx_3d_rect_cd4 ( Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               lx_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lx_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx_3d_rect_cd4 = lx_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( Y ( k ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) - CMPLX ( Z ( l ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lx_3d_rect_cd4 = -0.08333333333333333 * lx_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly_3d_rect_cd2 ( X , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ly_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly_3d_rect_cd2 = ly_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( Z ( l ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - CMPLX ( X ( j ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ly_3d_rect_cd2 = -0.5 * ly_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly_3d_rect_cd4 ( X , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ly_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly_3d_rect_cd4 = ly_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( Z ( l ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) - CMPLX ( X ( j ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ly_3d_rect_cd4 = -0.08333333333333333 * ly_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz_3d_rect_cd2 ( X , Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lz_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz_3d_rect_cd2 = lz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( X ( j ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - CMPLX ( Y ( k ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lz_3d_rect_cd2 = -0.5 * lz_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz_3d_rect_cd4 ( X , Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lz_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz_3d_rect_cd4 = lz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( X ( j ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) - CMPLX ( Y ( k ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lz_3d_rect_cd4 = -0.08333333333333333 * lz_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx2_3d_rect_cd2 ( Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               lx2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lx2_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx2_3d_rect_cd2 = lx2_3d_rect_cd2 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( Y ( k ) / dZ )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * Y ( k ) / dY , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                           &    CMPLX ( 0.5 * Y ( k ) * Z ( l ) / ( dY * dZ ) , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l + 1 ) - Psi3 ( j , k + 1 , l - 1 ) - &
                           &       Psi3 ( j , k - 1 , l + 1 ) + Psi3 ( j , k - 1 , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * Z ( l ) / dZ , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + &
                           &    CMPLX ( ( Z ( l ) / dY )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lx2_3d_rect_cd2 = -lx2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lx2_3d_rect_cd4 ( Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lx2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lx2_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lx2_3d_rect_cd4 = lx2_3d_rect_cd4 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( Y ( k ) / dZ )**2 , 0.0 ) * ( -Psi3 ( j , k , l - 2 ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l + 2 ) ) + &
                           &    CMPLX ( Y ( k ) / dY , 0.0 ) * ( Psi3 ( j , k + 2 , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                           &       Psi3 ( j , k - 2 , l ) ) - &
                           &    CMPLX ( Y ( k ) * Z ( l ) / ( 6.0 * dZ * dY ) , 0.0 ) * ( Psi3 ( j , k - 2 , l - 2 ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l - 2 ) - Psi3 ( j , k - 1 , l - 2 ) ) - &
                           &       Psi3 ( j , k + 2 , l - 2 ) - CMPLX ( 8.0 , 0.0 ) * Psi3 ( j , k - 2 , l - 1 ) - &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l - 1 ) - Psi3 ( j , k - 1 , l - 1 ) ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 2 , l - 1 ) + Psi3 ( j , k - 2 , l + 1 ) ) + &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l + 1 ) - Psi3 ( j , k - 1 , l + 1 ) ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * Psi3 ( j , k + 2 , l + 1 ) - Psi3 ( j , k - 2 , l + 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l + 2 ) - Psi3 ( j , k - 1 , l + 2 ) ) + & 
                           &       Psi3 ( j , k + 2 , l + 2 ) ) + &
                           &    CMPLX ( Z ( l ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                           &       Psi3 ( j , k , l - 2 ) ) + &
                           &    CMPLX ( ( Z ( l ) / dY )**2 , 0.0 ) * ( -Psi3 ( j , k - 2 , l ) + &
                           &    CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &    CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - Psi3 ( j , k + 2 , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lx2_3d_rect_cd4 = -0.08333333333 * lx2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly2_3d_rect_cd2 ( X , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ly2_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly2_3d_rect_cd2 = ly2_3d_rect_cd2 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( Z ( l ) / dX )**2 , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) - &
                           &    CMPLX ( 0.5 * Z ( l ) / dZ , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * X ( j ) * Z ( l ) / ( dX * dZ ) , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l + 1 ) - Psi3 ( j - 1 , k , l + 1 ) - &
                           &       Psi3 ( j + 1 , k , l - 1 ) + Psi3 ( j - 1 , k , l - 1 ) ) - &
                           &    CMPLX ( 0.5 * X ( j ) / dX , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + &
                           &    CMPLX ( ( X ( j ) / dZ )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) ) )

                     END DO
       
                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ly2_3d_rect_cd2 = -ly2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ly2_3d_rect_cd4 ( X , Z , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ly2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ly2_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ly2_3d_rect_cd4 = ly2_3d_rect_cd4 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( Z ( l ) / dX )**2 , 0.0 ) * ( -Psi3 ( j - 2 , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - CMPLX ( 30.0, 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - Psi3 ( j + 2 , k , l ) ) + &
                           &    CMPLX ( Z ( l ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                           &       Psi3 ( j , k , l - 2 ) ) - &
                           &    CMPLX ( Z ( l ) * X ( j ) / ( 6.0 * dX * dZ ) , 0.0 ) * ( Psi3 ( j - 2 , k , l - 2 ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l - 2 ) - Psi3 ( j + 1 , k , l - 2 ) ) - &
                           &       Psi3 ( j + 2 , k , l - 2 ) - CMPLX ( 8.0 , 0.0 ) * Psi3 ( j - 2 , k , l - 1 ) + &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l - 1 ) - Psi3 ( j + 1 , k , l - 1 ) ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 2 , k , l - 1 ) + Psi3 ( j - 2 , k , l + 1 ) ) - &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l + 1 ) - Psi3 ( j + 1 , k , l + 1 ) ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * Psi3 ( j + 2 , k , l + 1 ) - Psi3 ( j - 2 , k , l + 2 ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 1 , k , l + 2 ) - Psi3 ( j + 1 , k , l + 2 ) ) + &
                           &       Psi3 ( j + 2 , k , l + 2 ) ) + &
                           &    CMPLX ( X ( j ) / dX , 0.0 ) * ( Psi3 ( j + 2 , k , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                           &       Psi3 ( j - 2 , k , l ) ) + &
                           &    CMPLX ( ( X ( j ) / dZ )**2 , 0.0 ) * ( -Psi3 ( j , k , l - 2 ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l + 2 ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ly2_3d_rect_cd4 = -0.08333333333 * ly2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz2_3d_rect_cd2 ( X , Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lz2_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz2_3d_rect_cd2 = lz2_3d_rect_cd2 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( X ( j ) / dY )**2 , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) - &
                           &    CMPLX ( 0.5 * X ( j ) / dX , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                           &    CMPLX ( 0.5 * X ( j ) * Y ( k ) / ( dX * dY ) , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k + 1 , l ) - Psi3 ( j + 1 , k - 1 , l ) - &
                           &       Psi3 ( j - 1 , k + 1 , l ) + Psi3 ( j - 1 , k - 1 , l ) ) - &
                           &    CMPLX ( 0.5 * Y ( k ) / dY , 0.0 ) * ( &
                           &       Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + &
                           &    CMPLX ( ( Y ( k ) / dX )**2 , 0.0 ) * ( &
                           &       Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lz2_3d_rect_cd2 = -lz2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION lz2_3d_rect_cd4 ( X , Y , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X 
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               lz2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : lz2_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        lz2_3d_rect_cd4 = lz2_3d_rect_cd4 + &
                           & REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( &
                           &    CMPLX ( ( X ( j ) / dY )**2 , 0.0 ) * ( -Psi3 ( j , k - 2 , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - Psi3 ( j , k + 2 , l ) ) + &
                           &    CMPLX ( X ( j ) / dX , 0.0 ) * ( Psi3 ( j + 2 , k , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                           &       Psi3 ( j - 2 , k , l ) ) - &
                           &    CMPLX ( X ( j ) * Y ( k ) / ( 6.0 * dY * dX ) ) * ( Psi3 ( j - 2 , k - 2 , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 2 , k - 1 , l ) - Psi3 ( j - 2 , k + 1 , l ) ) - &
                           &       Psi3 ( j - 2 , k + 2 , l ) - CMPLX ( 8.0 , 0.0 ) * Psi3 ( j - 1 , k - 2 , l ) + &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j - 1 , k - 1 , l ) - Psi3 ( j - 1 , k + 1 , l ) ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j - 1 , k + 2 , l ) + Psi3 ( j + 1 , k - 2 , l ) ) - &
                           &       CMPLX ( 64.0 , 0.0 ) * ( Psi3 ( j + 1 , k - 1 , l ) - Psi3 ( j + 1 , k + 1 , l ) ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * Psi3 ( j + 1 , k + 2 , l ) - Psi3 ( j + 2 , k - 2 , l ) + &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 2 , k - 1 , l ) - Psi3 ( j + 2 , k + 1 , l ) ) + &
                           &       Psi3 ( j + 2 , k + 2 , l ) ) + &
                           &    CMPLX ( Y ( k ) / dY , 0.0 ) * ( Psi3 ( j , k + 2 , l ) - &
                           &       CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                           &       Psi3 ( j , k - 2 , l ) ) + &
                           &    CMPLX ( ( Y ( k ) / dX )**2 , 0.0 ) * ( -Psi3 ( j - 2 , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                           &       CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - Psi3 ( j + 2 , k , l ) ) ) )

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               lz2_3d_rect_cd4 = -0.08333333333 * lz2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fx_3d_rect_cd2 ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fx_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : fx_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fx_3d_rect_cd2 = fx_3d_rect_cd2 + ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               fx_3d_rect_cd2 = 0.5 * fx_3d_rect_cd2 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fx_3d_rect_cd4 ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fx_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : fx_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fx_3d_rect_cd4 = fx_3d_rect_cd4 + ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               fx_3d_rect_cd4 = 0.08333333333333333 * fx_3d_rect_cd4 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fy_3d_rect_cd2 ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fy_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : fy_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fy_3d_rect_cd2 = fy_3d_rect_cd2 + ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               fy_3d_rect_cd2 = 0.5 * fy_3d_rect_cd2 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fy_3d_rect_cd4 ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fy_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : fy_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fy_3d_rect_cd4 = fy_3d_rect_cd4 + ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               fy_3d_rect_cd4 = 0.08333333333333333 * fy_3d_rect_cd4 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION fz_3d_rect_cd2 ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : fz_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fz_3d_rect_cd2 = fz_3d_rect_cd2 + ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               fz_3d_rect_cd2 = 0.5 * fz_3d_rect_cd2 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION fz_3d_rect_cd4 ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               fz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : fz_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        fz_3d_rect_cd4 = fz_3d_rect_cd4 + ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               fz_3d_rect_cd4 = 0.08333333333333333 * fz_3d_rect_cd4 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION taux_3d_rect_cd2 ( Y , Z , Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               taux_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : taux_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        taux_3d_rect_cd2 = taux_3d_rect_cd2 + ( Y ( k ) * ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) / dZ - Z ( l ) * ( Vex3 ( j , k - 1 , l ) -Vex3 ( j , k + 1 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               taux_3d_rect_cd2 = 0.5 * taux_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION taux_3d_rect_cd4 ( Y , Z , Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y 
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z 
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3 

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 

               INTEGER :: j , k , l

               taux_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : taux_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        taux_3d_rect_cd4 = taux_3d_rect_cd4 + ( Y ( k ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ - Z ( l ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               taux_3d_rect_cd4 = 0.08333333333333333 * taux_3d_rect_cd4 * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION tauy_3d_rect_cd2 ( X , Z , Vex3 , Psi3 ) 

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauy_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : tauy_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauy_3d_rect_cd2 = tauy_3d_rect_cd2 + ( Z ( l ) * ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) / dX - X ( j ) * ( Vex3 ( j , k , l - 1 ) -Vex3 ( j , k , l + 1 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               tauy_3d_rect_cd2 = 0.5 * tauy_3d_rect_cd2 * dX * dY * dZ

            END FUNCTION

            REAL FUNCTION tauy_3d_rect_cd4 ( X , Z , Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauy_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : tauy_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauy_3d_rect_cd4 = tauy_3d_rect_cd4 + ( Z ( l ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX - X ( j ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               tauy_3d_rect_cd4 = 0.08333333333333333 * tauy_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION tauz_3d_rect_cd2 ( X , Y , Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : tauz_3d_rect_cd2 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauz_3d_rect_cd2 = tauz_3d_rect_cd2 + ( X ( j ) * ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) / dY - Y ( k ) * ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               tauz_3d_rect_cd2 = 0.5 * tauz_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION tauz_3d_rect_cd4 ( X , Y , Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               tauz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : tauz_3d_rect_cd4 )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        tauz_3d_rect_cd4 = tauz_3d_rect_cd4 + ( X ( j ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2, l ) ) / dY - Y ( k ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               tauz_3d_rect_cd4 = 0.08333333333333333 * tauz_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ixx_3d_rect ( yO , zO , Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ixx_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ixx_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ixx_3d_rect = ixx_3d_rect + ( ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ixx_3d_rect = ixx_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION iyy_3d_rect ( xO , zO , X , Z , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: zO

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               iyy_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : iyy_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        iyy_3d_rect = iyy_3d_rect + ( ( X ( j ) - xO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               iyy_3d_rect = iyy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION izz_3d_rect ( xO , yO , X , Y , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               izz_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : izz_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        izz_3d_rect = izz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               izz_3d_rect = izz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ixy_3d_rect ( xO , yO , X , Y , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: yO

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ixy_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ixy_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ixy_3d_rect = ixy_3d_rect + ( X ( j ) - xO ) * ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ixy_3d_rect = -ixy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION iyz_3d_rect ( yO , zO , Y , Z , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: yO
               REAL, INTENT ( IN ) :: zO

               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               iyz_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : iyz_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        iyz_3d_rect = iyz_3d_rect + ( Y ( k ) - yO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               iyz_3d_rect = -iyz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION ixz_3d_rect ( xO , zO , X , Z , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: xO
               REAL, INTENT ( IN ) :: zO

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               ixz_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : ixz_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        ixz_3d_rect = ixz_3d_rect + ( X ( j ) - xO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               ixz_3d_rect = -ixz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION vex_3d_rect ( Vex3 , Psi3 )

               IMPLICIT NONE

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               vex_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : vex_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        vex_3d_rect = vex_3d_rect + Vex3 ( j , k , l ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               vex_3d_rect = vex_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION vmf_3d_rect ( gS , Psi3 )

               IMPLICIT NONE

               REAL, INTENT ( IN ) :: gS

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               INTEGER :: j , k , l

               vmf_3d_rect = 0.0

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC ) REDUCTION ( + : vmf_3d_rect )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        vmf_3d_rect = vmf_3d_rect + ABS ( Psi3 ( j , k , l ) )**4

                     END DO

                  END DO

               END DO
!$OMP          END DO
!$OMP          END PARALLEL

               vmf_3d_rect = 0.5 * gS * vmf_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

      END MODULE

! =========================================================================
