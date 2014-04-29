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
!     Tuesday, April 29th, 2014
!
! -------------------------------------------------------------------------

      MODULE EVUA

         USE, INTRINSIC :: ISO_FORTRAN_ENV
         USE            :: MATH

         IMPLICIT NONE
         PRIVATE

         PUBLIC :: l2_norm_3d_rect
         PUBLIC :: x_3d_rect
         PUBLIC :: y_3d_rect
         PUBLIC :: z_3d_rect
         PUBLIC :: x2_3d_rect
         PUBLIC :: y2_3d_rect
         PUBLIC :: z2_3d_rect
         PUBLIC :: px_3d_rect_cd2
         PUBLIC :: px_3d_rect_cd4
         PUBLIC :: py_3d_rect_cd2
         PUBLIC :: py_3d_rect_cd4
         PUBLIC :: pz_3d_rect_cd2
         PUBLIC :: pz_3d_rect_cd4
         PUBLIC :: px2_3d_rect_cd2
         PUBLIC :: px2_3d_rect_cd4
         PUBLIC :: py2_3d_rect_cd2
         PUBLIC :: py2_3d_rect_cd4
         PUBLIC :: pz2_3d_rect_cd2
         PUBLIC :: pz2_3d_rect_cd4
         PUBLIC :: lx_3d_rect_cd2
         PUBLIC :: lx_3d_rect_cd4
         PUBLIC :: ly_3d_rect_cd2
         PUBLIC :: ly_3d_rect_cd4
         PUBLIC :: lz_3d_rect_cd2
         PUBLIC :: lz_3d_rect_cd4
         PUBLIC :: lx2_3d_rect_cd2
!         PUBLIC :: lx2_3d_rect_cd4
         PUBLIC :: ly2_3d_rect_cd2
!         PUBLIC :: ly2_3d_rect_cd4
         PUBLIC :: lz2_3d_rect_cd2
!         PUBLIC :: lz2_3d_rect_cd4
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

            REAL FUNCTION x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Psi3 )

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

            REAL FUNCTION y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Psi3 )

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

            REAL FUNCTION z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Z , Psi3 )

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

               px2_3d_rect_cd4 = -0.0833333333333333333 * px2_3d_rect_cd4 * dY * dZ / dX

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

            REAL FUNCTION lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Z , Psi3 )

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

            REAL FUNCTION lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Z , Psi3 )

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

            REAL FUNCTION ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Z , Psi3 )

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

            REAL FUNCTION ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Z , Psi3 )

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

            REAL FUNCTION lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Y , Psi3 )

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

            REAL FUNCTION lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Y , Psi3 )

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

            REAL FUNCTION lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Z , Psi3 )

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

            REAL FUNCTION ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Z , Psi3 )

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

            REAL FUNCTION lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Y , Psi3 )

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
