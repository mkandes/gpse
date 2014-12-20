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
!     Friday, December 19th, 2014
!
! -------------------------------------------------------------------------

      MODULE EVUA

         USE, INTRINSIC :: ISO_FORTRAN_ENV
         USE            :: MPI
         USE            :: MATH

         IMPLICIT NONE
         PRIVATE

         REAL, PRIVATE :: evuaL2Norma = 0.0
         REAL, PRIVATE :: evuaL2Normb = 0.0
         REAL, PRIVATE :: evuaXa = 0.0
         REAL, PRIVATE :: evuaXb = 0.0
         REAL, PRIVATE :: evuaYa = 0.0
         REAL, PRIVATE :: evuaYb = 0.0
         REAL, PRIVATE :: evuaZa = 0.0
         REAL, PRIVATE :: evuaZb = 0.0
         REAL, PRIVATE :: evuaRa = 0.0
         REAL, PRIVATE :: evuaRb = 0.0
         REAL, PRIVATE :: evuaX2a = 0.0
         REAL, PRIVATE :: evuaX2b = 0.0
         REAL, PRIVATE :: evuaY2a = 0.0
         REAL, PRIVATE :: evuaY2b = 0.0
         REAL, PRIVATE :: evuaZ2a = 0.0
         REAL, PRIVATE :: evuaZ2b = 0.0
         REAL, PRIVATE :: evuaR2a = 0.0
         REAL, PRIVATE :: evuaR2b = 0.0
         REAL, PRIVATE :: evuaIxxa = 0.0
         REAL, PRIVATE :: evuaIxxb = 0.0
         REAL, PRIVATE :: evuaIxya = 0.0
         REAL, PRIVATE :: evuaIxyb = 0.0
         REAL, PRIVATE :: evuaIxza = 0.0
         REAL, PRIVATE :: evuaIxzb = 0.0
         REAL, PRIVATE :: evuaIyya = 0.0
         REAL, PRIVATE :: evuaIyyb = 0.0
         REAL, PRIVATE :: evuaIyza = 0.0
         REAL, PRIVATE :: evuaIyzb = 0.0
         REAL, PRIVATE :: evuaIzza = 0.0
         REAL, PRIVATE :: evuaIzzb = 0.0
         REAL, PRIVATE :: evuaVexa = 0.0
         REAL, PRIVATE :: evuaVexb = 0.0
         REAL, PRIVATE :: evuaVmfa = 0.0
         REAL, PRIVATE :: evuaVmfb = 0.0
         REAL, PRIVATE :: evuaPxa = 0.0
         REAL, PRIVATE :: evuaPxb = 0.0
         REAL, PRIVATE :: evuaPya = 0.0
         REAL, PRIVATE :: evuaPyb = 0.0
         REAL, PRIVATE :: evuaPza = 0.0
         REAL, PRIVATE :: evuaPzb = 0.0
         REAL, PRIVATE :: evuaPx2a = 0.0
         REAL, PRIVATE :: evuaPx2b = 0.0
         REAL, PRIVATE :: evuaPy2a = 0.0
         REAL, PRIVATE :: evuaPy2b = 0.0
         REAL, PRIVATE :: evuaPz2a = 0.0
         REAL, PRIVATE :: evuaPz2b = 0.0
         REAL, PRIVATE :: evuaLxa = 0.0
         REAL, PRIVATE :: evuaLxb = 0.0
         REAL, PRIVATE :: evuaLya = 0.0
         REAL, PRIVATE :: evuaLyb = 0.0
         REAL, PRIVATE :: evuaLza = 0.0
         REAL, PRIVATE :: evuaLzb = 0.0
         REAL, PRIVATE :: evuaLx2a = 0.0
         REAL, PRIVATE :: evuaLx2b = 0.0
         REAL, PRIVATE :: evuaLy2a = 0.0
         REAL, PRIVATE :: evuaLy2b = 0.0
         REAL, PRIVATE :: evuaLz2a = 0.0
         REAL, PRIVATE :: evuaLz2b = 0.0
         REAL, PRIVATE :: evuaFxa = 0.0
         REAL, PRIVATE :: evuaFxb = 0.0
         REAL, PRIVATE :: evuaFya = 0.0
         REAL, PRIVATE :: evuaFyb = 0.0
         REAL, PRIVATE :: evuaFza = 0.0
         REAL, PRIVATE :: evuaFzb = 0.0
         REAL, PRIVATE :: evuaTauXa = 0.0
         REAL, PRIVATE :: evuaTauXb = 0.0
         REAL, PRIVATE :: evuaTauYa = 0.0
         REAL, PRIVATE :: evuaTauYb = 0.0
         REAL, PRIVATE :: evuaTauZa = 0.0
         REAL, PRIVATE :: evuaTauZb = 0.0
         REAL, PRIVATE :: evuaE = 0.0
         REAL, PRIVATE :: evuaMu = 0.0
         REAL, PRIVATE :: evuaL2 = 0.0
         REAL, PRIVATE :: evuaTx = 0.0
         REAL, PRIVATE :: evuaTy = 0.0
         REAL, PRIVATE :: evuaTz = 0.0
         REAL, PRIVATE :: evuaSigX = 0.0
         REAL, PRIVATE :: evuaSigY = 0.0
         REAL, PRIVATE :: evuaSigZ = 0.0
         REAL, PRIVATE :: evuaSigPx = 0.0
         REAL, PRIVATE :: evuaSigPy = 0.0
         REAL, PRIVATE :: evuaSigPz = 0.0
         REAL, PRIVATE :: evuaSigLx = 0.0
         REAL, PRIVATE :: evuaSigLy = 0.0
         REAL, PRIVATE :: evuaSigLz = 0.0
         
         PUBLIC :: evua_compute_base
         PUBLIC :: evua_reduce_base
         PUBLIC :: evua_compute_derived
         PUBLIC :: evua_write_all
         PUBLIC :: evua_normalize

         PRIVATE :: evua_l2_norm_3d_rect
         PRIVATE :: evua_x_3d_rect
         PRIVATE :: evua_y_3d_rect
         PRIVATE :: evua_z_3d_rect
         PRIVATE :: evua_r_xy_3d_rect
         PRIVATE :: evua_r_xyz_3d_rect
         PRIVATE :: evua_x2_3d_rect
         PRIVATE :: evua_y2_3d_rect
         PRIVATE :: evua_z2_3d_rect
         PRIVATE :: evua_r2_xy_3d_rect
         PRIVATE :: evua_r2_xyz_3d_rect
         PRIVATE :: evua_px_3d_rect_cd2
         PRIVATE :: evua_px_3d_rect_cd4
         PRIVATE :: evua_py_3d_rect_cd2
         PRIVATE :: evua_py_3d_rect_cd4
         PRIVATE :: evua_pz_3d_rect_cd2
         PRIVATE :: evua_pz_3d_rect_cd4
         PRIVATE :: evua_px2_3d_rect_cd2
         PRIVATE :: evua_px2_3d_rect_cd4
         PRIVATE :: evua_py2_3d_rect_cd2
         PRIVATE :: evua_py2_3d_rect_cd4
         PRIVATE :: evua_pz2_3d_rect_cd2
         PRIVATE :: evua_pz2_3d_rect_cd4
         PRIVATE :: evua_lx_3d_rect_cd2
         PRIVATE :: evua_lx_3d_rect_cd4
         PRIVATE :: evua_ly_3d_rect_cd2
         PRIVATE :: evua_ly_3d_rect_cd4
         PRIVATE :: evua_lz_3d_rect_cd2
         PRIVATE :: evua_lz_3d_rect_cd4
         PRIVATE :: evua_lx2_3d_rect_cd2
         PRIVATE :: evua_lx2_3d_rect_cd4
         PRIVATE :: evua_ly2_3d_rect_cd2
         PRIVATE :: evua_ly2_3d_rect_cd4
         PRIVATE :: evua_lz2_3d_rect_cd2
         PRIVATE :: evua_lz2_3d_rect_cd4
         PRIVATE :: evua_fx_3d_rect_cd2
         PRIVATE :: evua_fx_3d_rect_cd4
         PRIVATE :: evua_fy_3d_rect_cd2
         PRIVATE :: evua_fy_3d_rect_cd4
         PRIVATE :: evua_fz_3d_rect_cd2
         PRIVATE :: evua_fz_3d_rect_cd4
         PRIVATE :: evua_taux_3d_rect_cd2
         PRIVATE :: evua_taux_3d_rect_cd4
         PRIVATE :: evua_tauy_3d_rect_cd2
         PRIVATE :: evua_tauy_3d_rect_cd4
         PRIVATE :: evua_tauz_3d_rect_cd2
         PRIVATE :: evua_tauz_3d_rect_cd4
         PRIVATE :: evua_ixx_3d_rect
         PRIVATE :: evua_iyy_3d_rect
         PRIVATE :: evua_izz_3d_rect
         PRIVATE :: evua_ixy_3d_rect
         PRIVATE :: evua_iyz_3d_rect
         PRIVATE :: evua_ixz_3d_rect
         PRIVATE :: evua_vex_3d_rect
         PRIVATE :: evua_vmf_3d_rect

         CONTAINS

            SUBROUTINE evua_compute_base ( evuaQuadRule , evuaFdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , gS , X , Y , Z , Vex3 , Psi3 )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: evuaQuadRule
               INTEGER, INTENT ( IN ) :: evuaFdOrder
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
               REAL, INTENT ( IN ) :: gS

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
               REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
               REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

               REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3
               
               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

               IF ( evuaQuadRule == 1 ) THEN

                  evuaL2norma = evua_l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                  evuaXa = evua_x_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )
                  evuaYa = evua_y_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )
                  evuaZa = evua_z_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )
                  evuaRa = evua_r_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                  evuaX2a = evua_x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )
                  evuaY2a = evua_y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )
                  evuaZ2a = evua_z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )
                  evuaR2a = evua_r2_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO ,  dX , dY , dZ , X , Y , Psi3 )
                  evuaIxxa = evua_ixx_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )
                  evuaIxya = evua_ixy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                  evuaIxza = evua_ixz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
                  evuaIyya = evua_iyy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
                  evuaIyza = evua_iyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )
                  evuaIzza = evua_izz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                  evuaVexa = evua_vex_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Vex3 , Psi3 )
                  evuaVmfa = evua_vmf_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , gS , Psi3 )

                  IF ( evuaFdOrder == 2 ) THEN

                     evuaPxa = evua_px_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )
                     evuaPya = evua_py_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )
                     evuaPza = evua_pz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )
                     evuaPx2a = evua_px2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                     evuaPy2a = evua_py2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                     evuaPz2a = evua_pz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                     evuaLxa = evua_lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , Z , Psi3 )
                     evuaLya = evua_ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
                     evuaLza = evua_lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                     evuaLx2a = evua_lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )
                     evuaLy2a = evua_ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
                     evuaLz2a = evua_lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                     evuaFxa = evua_fx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )
                     evuaFya = evua_fy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )
                     evuaFza = evua_fz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )
                     evuaTauXa = evua_taux_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Vex3 , Psi3 )
                     evuaTauYa = evua_tauy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Vex3 , Psi3 )
                     evuaTauZa = evua_tauz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Vex3 , Psi3 )

                  ELSE IF ( evuaFdOrder == 4 ) THEN

                     evuaPxa = evua_px_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )
                     evuaPya = evua_py_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )
                     evuaPza = evua_pz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )
                     evuaPx2a = evua_px2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                     evuaPy2a = evua_py2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                     evuaPz2a = evua_pz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
                     evuaLxa = evua_lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , Z , Psi3 )
                     evuaLya = evua_ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
                     evuaLza = evua_lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                     evuaLx2a = evua_lx2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )
                     evuaLy2a = evua_ly2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
                     evuaLz2a = evua_lz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
                     evuaFxa = evua_fx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )
                     evuaFya = evua_fy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )
                     evuaFza = evua_fz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )
                     evuaTauXa = evua_taux_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Vex3 , Psi3 )
                     evuaTauYa = evua_tauy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Vex3 , Psi3 )
                     evuaTauZa = evua_tauz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Vex3 , Psi3 )

                  ELSE

                     WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : evua_compute_base : ERROR - evuaFdOrder is not supported.'
                     STOP

                  END IF

               ELSE

                  WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : evua_compute_base : ERROR - evuaQuadRule is not supported.'
                  STOP

               END IF

               RETURN

            END SUBROUTINE

            SUBROUTINE evua_reduce_base ( mpiMaster , mpiReal , mpiError )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: mpiMaster
               INTEGER, INTENT ( IN ) :: mpiReal
               INTEGER, INTENT ( INOUT ) :: mpiError
 
               CALL MPI_REDUCE ( evuaL2norma , evuaL2normb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaXa , evuaXb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaYa , evuaYb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaZa , evuaZb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaRa , evuaRb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaX2a , evuaX2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaY2a , evuaY2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaZ2a , evuaZ2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaR2a , evuaR2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaIxxa , evuaIxxb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaIxya , evuaIxyb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaIxza , evuaIxzb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaIyya , evuaIyyb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaIyza , evuaIyzb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaIzza , evuaIzzb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaVexa , evuaVexb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaVmfa , evuaVmfb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaPxa , evuaPxb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaPya , evuaPyb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaPza , evuaPzb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaPx2a , evuaPx2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaPy2a , evuaPy2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaPz2a , evuaPz2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaLxa , evuaLxb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaLya , evuaLyb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaLza , evuaLzb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaLx2a , evuaLx2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaLy2a , evuaLy2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaLz2a , evuaLz2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaFxa , evuaFxb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaFya , evuaFyb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaFza , evuaFzb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaTauXa , evuaTauXb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaTauYa , evuaTauYb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_REDUCE ( evuaTauZa , evuaTauZb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

               RETURN

            END SUBROUTINE

            SUBROUTINE evua_compute_derived ( mpiRank , mpiMaster , wX , wY , wZ )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: mpiRank
               INTEGER, INTENT ( IN ) :: mpiMaster
 
               REAL, INTENT ( IN ) :: wX
               REAL, INTENT ( IN ) :: wY
               REAL, INTENT ( IN ) :: wZ

               IF ( mpiRank == mpiMaster ) THEN

!                 Kinetic energy expectation values

                  evuaTx = 0.5 * evuaPx2b
                  evuaTy = 0.5 * evuaPy2b
                  evuaTz = 0.5 * evuaPz2b

!                 Energy expectation value

                  evuaE = evuaTx + evuaTy + evuaTz + evuaVexb + evuaVmfb - wX * evuaLxb - wY * evuaLyb - wZ * evuaLzb

!                 Chemical potential

                  evuaMu =  evuaTx + evuaTy + evuaTz + evuaVexb + 2.0 * evuaVmfb - wX * evuaLxb - wY * evuaLyb - wZ * evuaLzb

!                 Squared angular momentum expectation value

                  evuaL2 = evuaLx2b + evuaLy2b + evuaLz2b

!                 Position, momentum and angular momentum uncertainty 
               
                  evuaSigX  = SQRT ( evuaX2b  - evuaXb**2  )
                  evuaSigY  = SQRT ( evuaY2b  - evuaYb**2  )
                  evuaSigZ  = SQRT ( evuaZ2b  - evuaZb**2  )
                  evuaSigPx = SQRT ( evuaPx2b - evuaPxb**2 )
                  evuaSigPy = SQRT ( evuaPy2b - evuaPyb**2 )
                  evuaSigPz = SQRT ( evuaPz2b - evuaPzb**2 )
                  evuaSigLx = SQRT ( evuaLx2b - evuaLxb**2 )
                  evuaSigLy = SQRT ( evuaLy2b - evuaLyb**2 )
                  evuaSigLz = SQRT ( evuaLz2b - evuaLzb**2 )

               END IF

               RETURN

            END SUBROUTINE

            SUBROUTINE evua_write_all ( mpiRank , mpiMaster , tN , wX , wY , wZ )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: mpiRank
               INTEGER, INTENT ( IN ) :: mpiMaster

               REAL, INTENT ( IN ) :: tN
               REAL, INTENT ( IN ) :: wX
               REAL, INTENT ( IN ) :: wY
               REAL, INTENT ( IN ) :: wZ

!              Write expectation values, uncertainties and uncertainty relations to file from MPI_MASTER

               IF ( mpiRank == mpiMaster ) THEN

                  WRITE ( UNIT = OUTPUT_UNIT , FMT = '(59(F23.15))' ) tN , wX , wY , wZ , evuaL2normb , evuaE , evuaMu , evuaL2 , evuaTx , evuaTy , evuaTz , evuaVexb , evuaVmfb , evuaXb , evuaYb , evuaZb , evuaRb , evuaPxb , evuaPyb , evuaPzb , evuaLxb , evuaLyb , evuaLzb , evuaFxb , evuaFyb , evuaFzb , evuaTauXb , evuaTauYb , evuaTauZb , evuaIxxb , evuaIxyb , evuaIxzb , evuaIyyb , evuaIyzb , evuaIzzb , evuaX2b , evuaY2b , evuaZ2b , evuaSigX , evuaSigY , evuaSigZ , evuaPx2b , evuaPy2b , evuaPz2b , evuaSigPx , evuaSigPy , evuaSigPz , evuaLx2b , evuaLy2b , evuaLz2b , evuaSigLx , evuaSigLy , evuaSigLz , evuaSigX * evuaSigPx , evuaSigY * evuaSigPy , evuaSigZ * evuaSigPz , evuaSigLx * evuaSigLy , evuaSigLy * evuaSigLz , evuaSigLz * evuaSigLx

               END IF

               RETURN

            END SUBROUTINE

            SUBROUTINE evua_normalize ( mpiMaster , mpiReal , mpiError , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

               IMPLICIT NONE

               INTEGER, INTENT ( IN ) :: mpiMaster
               INTEGER, INTENT ( IN ) :: mpiReal
               INTEGER, INTENT ( IN ) :: mpiError
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

               COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3 

               evuaL2norma = evua_l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ ,Psi3 )
               CALL MPI_REDUCE ( evuaL2Norma , evuaL2Normb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
               CALL MPI_BCAST ( evuaL2Normb , 1 , mpiReal , mpiMaster , MPI_COMM_WORLD , mpiError )
               Psi3 = Psi3 / SQRT ( evuaL2normb )

               RETURN

            END SUBROUTINE

            REAL FUNCTION evua_l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_l2_norm_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_l2_norm_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_l2_norm_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_l2_norm_3d_rect = evua_l2_norm_3d_rect + ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_l2_norm_3d_rect = evua_l2_norm_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_x_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )

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

               evua_x_3d_rect = 0.0 

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb  

                        evua_x_3d_rect = evua_x_3d_rect + ( X ( j ) - xO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_x_3d_rect = evua_x_3d_rect * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION evua_y_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )

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

               evua_y_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_y_3d_rect = evua_y_3d_rect + ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_y_3d_rect = evua_y_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_z_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )

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

               evua_z_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_z_3d_rect = evua_z_3d_rect + ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_z_3d_rect = evua_z_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_r_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_r_xy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_r_xy_3d_rect = evua_r_xy_3d_rect + SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_r_xy_3d_rect = evua_r_xy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_r_xyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , X , Y , Z , Psi3 )

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

               evua_r_xyz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xyz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xyz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_r_xyz_3d_rect = evua_r_xyz_3d_rect + SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_r_xyz_3d_rect = evua_r_xyz_3d_rect * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION evua_x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )

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

               evua_x2_3d_rect = 0.0 

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x2_3d_rect )
               DO l = nZa , nZb 

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x2_3d_rect )
                  DO k = nYa , nYb 

                     DO j = nXa , nXb  

                        evua_x2_3d_rect = evua_x2_3d_rect + ( X ( j ) - xO )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_x2_3d_rect = evua_x2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )

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

               evua_y2_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y2_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y2_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_y2_3d_rect = evua_y2_3d_rect + ( Y ( k ) - yO )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_y2_3d_rect = evua_y2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )

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

               evua_z2_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z2_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z2_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_z2_3d_rect = evua_z2_3d_rect + ( Z ( l ) - zO )**2 * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_z2_3d_rect = evua_z2_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_r2_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO ,  dX , dY , dZ , X , Y , Psi3 )

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

               evua_r2_xy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_r2_xy_3d_rect = evua_r2_xy_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_r2_xy_3d_rect = evua_r2_xy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_r2_xyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , X , Y , Z , Psi3 )

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

               evua_r2_xyz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xyz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xyz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_r2_xyz_3d_rect = evua_r2_xyz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_r2_xyz_3d_rect = evua_r2_xyz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_px_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )

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

               evua_px_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_px_3d_rect_cd2 = evua_px_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_px_3d_rect_cd2 = -0.5 * evua_px_3d_rect_cd2 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_px_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )

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

               evua_px_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd4 ) 
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd4 ) 
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_px_3d_rect_cd4 = evua_px_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_px_3d_rect_cd4 = -0.08333333333333333 * evua_px_3d_rect_cd4 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_py_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )

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

               evua_py_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_py_3d_rect_cd2 = evua_py_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_py_3d_rect_cd2 = -0.5 * evua_py_3d_rect_cd2 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_py_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )

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

               evua_py_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_py_3d_rect_cd4 = evua_py_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k -1 , l ) ) + Psi3 ( j , k - 2 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_py_3d_rect_cd4 = -0.08333333333333333 * evua_py_3d_rect_cd4 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_pz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )

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

               evua_pz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd2 ) 
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_pz_3d_rect_cd2 = evua_pz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_pz_3d_rect_cd2 = -0.5 * evua_pz_3d_rect_cd2 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION evua_pz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )

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

               evua_pz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_pz_3d_rect_cd4 = evua_pz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_pz_3d_rect_cd4 = -0.08333333333333333 * evua_pz_3d_rect_cd4 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION evua_px2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_px2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd2 ) 
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_px2_3d_rect_cd2 = evua_px2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_px2_3d_rect_cd2 = -evua_px2_3d_rect_cd2 * dY * dZ / dX

               RETURN

            END FUNCTION

            REAL FUNCTION evua_px2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_px2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_px2_3d_rect_cd4 = evua_px2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - Psi3 ( j - 2 , k , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_px2_3d_rect_cd4 = -0.08333333333333333 * evua_px2_3d_rect_cd4 * dY * dZ / dX

               RETURN

            END FUNCTION

            REAL FUNCTION evua_py2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_py2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_py2_3d_rect_cd2 = evua_py2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_py2_3d_rect_cd2 = -evua_py2_3d_rect_cd2 * dX * dZ / dY

               RETURN

            END FUNCTION

            REAL FUNCTION evua_py2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_py2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_py2_3d_rect_cd4 = evua_py2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - Psi3 ( j , k - 2 , l ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_py2_3d_rect_cd4 = -0.08333333333333333 * evua_py2_3d_rect_cd4 * dX * dZ / dY

               RETURN
       
            END FUNCTION

            REAL FUNCTION evua_pz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_pz2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_pz2_3d_rect_cd2 = evua_pz2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) ) 

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_pz2_3d_rect_cd2 = -evua_pz2_3d_rect_cd2 * dX * dY / dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_pz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )

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

               evua_pz2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_pz2_3d_rect_cd4 = evua_pz2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - Psi3 ( j , k , l - 2 ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_pz2_3d_rect_cd4 = -0.08333333333333333 * evua_pz2_3d_rect_cd4 * dX * dY / dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , Z , Psi3 )

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

               evua_lx_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lx_3d_rect_cd2 = evua_lx_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Y ( k ) - yO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - CMPLX ( ( Z ( l ) - zO ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_lx_3d_rect_cd2 = -0.5 * evua_lx_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               evua_lx_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lx_3d_rect_cd4 = evua_lx_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Y ( k ) - yO ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) - CMPLX ( ( Z ( l ) - zO ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_lx_3d_rect_cd4 = -0.08333333333333333 * evua_lx_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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

               evua_ly_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ly_3d_rect_cd2 = evua_ly_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Z ( l ) - zO ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - CMPLX ( ( X ( j ) - xO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_ly_3d_rect_cd2 = -0.5 * evua_ly_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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

               evua_ly_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ly_3d_rect_cd4 = evua_ly_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( Z ( l ) - zO ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) - CMPLX ( ( X ( j ) - xO ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_ly_3d_rect_cd4 = -0.08333333333333333 * evua_ly_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_lz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lz_3d_rect_cd2 = evua_lz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( X ( j ) - xO ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - CMPLX ( ( Y ( k ) - yO ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_lz_3d_rect_cd2 = -0.5 * evua_lz_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_lz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lz_3d_rect_cd4 = evua_lz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( CMPLX ( ( X ( j ) - xO ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) - CMPLX ( ( Y ( k ) - yO ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) ) )

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_lz_3d_rect_cd4 = -0.08333333333333333 * evua_lz_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               evua_lx2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lx2_3d_rect_cd2 = evua_lx2_3d_rect_cd2 + &
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

               evua_lx2_3d_rect_cd2 = -evua_lx2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lx2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               evua_lx2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lx2_3d_rect_cd4 = evua_lx2_3d_rect_cd4 + &
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

               evua_lx2_3d_rect_cd4 = -0.08333333333 * evua_lx2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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

               evua_ly2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ly2_3d_rect_cd2 = evua_ly2_3d_rect_cd2 + &
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

               evua_ly2_3d_rect_cd2 = -evua_ly2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ly2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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

               evua_ly2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ly2_3d_rect_cd4 = evua_ly2_3d_rect_cd4 + &
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

               evua_ly2_3d_rect_cd4 = -0.08333333333 * evua_ly2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_lz2_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lz2_3d_rect_cd2 = evua_lz2_3d_rect_cd2 + &
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

               evua_lz2_3d_rect_cd2 = -evua_lz2_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_lz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_lz2_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_lz2_3d_rect_cd4 = evua_lz2_3d_rect_cd4 + &
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

               evua_lz2_3d_rect_cd4 = -0.08333333333 * evua_lz2_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_fx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )

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

               evua_fx_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_fx_3d_rect_cd2 = evua_fx_3d_rect_cd2 + ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_fx_3d_rect_cd2 = 0.5 * evua_fx_3d_rect_cd2 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_fx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )

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

               evua_fx_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_fx_3d_rect_cd4 = evua_fx_3d_rect_cd4 + ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_fx_3d_rect_cd4 = 0.08333333333333333 * evua_fx_3d_rect_cd4 * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_fy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )

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

               evua_fy_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_fy_3d_rect_cd2 = evua_fy_3d_rect_cd2 + ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_fy_3d_rect_cd2 = 0.5 * evua_fy_3d_rect_cd2 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_fy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )

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

               evua_fy_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_fy_3d_rect_cd4 = evua_fy_3d_rect_cd4 + ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_fy_3d_rect_cd4 = 0.08333333333333333 * evua_fy_3d_rect_cd4 * dX * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_fz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )

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

               evua_fz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_fz_3d_rect_cd2 = evua_fz_3d_rect_cd2 + ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_fz_3d_rect_cd2 = 0.5 * evua_fz_3d_rect_cd2 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION evua_fz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )

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

               evua_fz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_fz_3d_rect_cd4 = evua_fz_3d_rect_cd4 + ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_fz_3d_rect_cd4 = 0.08333333333333333 * evua_fz_3d_rect_cd4 * dX * dY

               RETURN

            END FUNCTION

            REAL FUNCTION evua_taux_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Vex3 , Psi3 )

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

               evua_taux_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_taux_3d_rect_cd2 = evua_taux_3d_rect_cd2 + ( ( Y ( k ) - yO ) * ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) / dZ - ( Z ( l ) - zO ) * ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_taux_3d_rect_cd2 = 0.5 * evua_taux_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_taux_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Vex3 , Psi3 )

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

               evua_taux_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_taux_3d_rect_cd4 = evua_taux_3d_rect_cd4 + ( ( Y ( k ) - yO ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ - ( Z ( l ) - zO ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_taux_3d_rect_cd4 = 0.08333333333333333 * evua_taux_3d_rect_cd4 * dX * dY * dZ 

               RETURN

            END FUNCTION

            REAL FUNCTION evua_tauy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Vex3 , Psi3 )

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

               evua_tauy_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_tauy_3d_rect_cd2 = evua_tauy_3d_rect_cd2 + ( ( Z ( l ) - zO ) * ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) / dX - ( X ( j ) - xO ) * ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_tauy_3d_rect_cd2 = 0.5 * evua_tauy_3d_rect_cd2 * dX * dY * dZ

            END FUNCTION

            REAL FUNCTION evua_tauy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Vex3 , Psi3 )

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

               evua_tauy_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_tauy_3d_rect_cd4 = evua_tauy_3d_rect_cd4 + ( ( Z ( l ) - zO ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX - ( X ( j ) - xO ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_tauy_3d_rect_cd4 = 0.08333333333333333 * evua_tauy_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_tauz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Vex3 , Psi3 )

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

               evua_tauz_3d_rect_cd2 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd2 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd2 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_tauz_3d_rect_cd2 = evua_tauz_3d_rect_cd2 + ( ( X ( j ) - xO ) * ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) / dY - ( Y ( k ) - yO ) * ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_tauz_3d_rect_cd2 = 0.5 * evua_tauz_3d_rect_cd2 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_tauz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Vex3 , Psi3 )

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

               evua_tauz_3d_rect_cd4 = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd4 )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd4 )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_tauz_3d_rect_cd4 = evua_tauz_3d_rect_cd4 + ( ( X ( j ) - xO ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2, l ) ) / dY - ( Y ( k ) - yO ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_tauz_3d_rect_cd4 = 0.08333333333333333 * evua_tauz_3d_rect_cd4 * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ixx_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               evua_ixx_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixx_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixx_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ixx_3d_rect = evua_ixx_3d_rect + ( ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_ixx_3d_rect = evua_ixx_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_iyy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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

               evua_iyy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_iyy_3d_rect = evua_iyy_3d_rect + ( ( X ( j ) - xO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_iyy_3d_rect = evua_iyy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_izz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_izz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_izz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_izz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_izz_3d_rect = evua_izz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_izz_3d_rect = evua_izz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ixy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )

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

               evua_ixy_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixy_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixy_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ixy_3d_rect = evua_ixy_3d_rect + ( X ( j ) - xO ) * ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_ixy_3d_rect = -evua_ixy_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_iyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )

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

               evua_iyz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_iyz_3d_rect = evua_iyz_3d_rect + ( Y ( k ) - yO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_iyz_3d_rect = -evua_iyz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_ixz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )

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

               evua_ixz_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixz_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixz_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_ixz_3d_rect = evua_ixz_3d_rect + ( X ( j ) - xO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_ixz_3d_rect = -evua_ixz_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_vex_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Vex3 , Psi3 )

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

               evua_vex_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vex_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vex_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_vex_3d_rect = evua_vex_3d_rect + Vex3 ( j , k , l ) * ABS ( Psi3 ( j , k , l ) )**2

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_vex_3d_rect = evua_vex_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

            REAL FUNCTION evua_vmf_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , gS , Psi3 )

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

               evua_vmf_3d_rect = 0.0

!$OMP          PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vmf_3d_rect )
               DO l = nZa , nZb

!$OMP             PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vmf_3d_rect )
                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        evua_vmf_3d_rect = evua_vmf_3d_rect + ABS ( Psi3 ( j , k , l ) )**4

                     END DO

                  END DO
!$OMP             END PARALLEL DO

               END DO
!$OMP          END PARALLEL DO

               evua_vmf_3d_rect = 0.5 * gS * evua_vmf_3d_rect * dX * dY * dZ

               RETURN

            END FUNCTION

      END MODULE

! =========================================================================
