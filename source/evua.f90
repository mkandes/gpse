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
!     Monday, February 10th, 2014
!
! -------------------------------------------------------------------------

      MODULE EVUA

         USE, INTRINSIC :: ISO_FORTRAN_ENV

         USE MATH

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
         PUBLIC :: lx2_3d_rect_cd4
         PUBLIC :: ly2_3d_rect_cd2
         PUBLIC :: ly2_3d_rect_cd4
         PUBLIC :: lz2_3d_rect_cd2
         PUBLIC :: lz2_3d_rect_cd4
         PUBLIC :: vex_3d_rect
         PUBLIC :: vmf_3d_rect

         CONTAINS

      END MODULE

! =========================================================================
