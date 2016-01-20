! ==================================================================================================================================
! NAME
!
!     evua [ evua ] - Expectation Value and Uncertainty Analysis Module
!
! SYNOPSIS
!
!     USE :: EVUA
!
! DESCRIPTION  
!
!     EVUA is a custom Fortran module written to compute expectation values and other useful quantities, such as uncertainty 
!        relations, from the wave functions generated during program execution.
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
!     Copyright (c) 2014, 2015, 2016 Martin Charles Kandes
!
! LAST UPDATED
!
!     Saturday, January 16th, 2016
!
! ----------------------------------------------------------------------------------------------------------------------------------

      MODULE EVUA

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MPI
      USE            :: MATH

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
!     MPI is the standard Message Passing Interface module.
!
!     MATH is a custom Fortran module written to define well-know mathematical constants and compute specialized functions. 
!
! ----------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE

!      INCLUDE 'mpif.h' ! if MPI module is not supported, use header file.

      PRIVATE

! --- VARIABLE DECLARATIONS --------------------------------------------------------------------------------------------------------

      REAL, PRIVATE :: evuaL2Norm = 0.0 
      REAL, PRIVATE :: evuaX      = 0.0
      REAL, PRIVATE :: evuaY      = 0.0
      REAL, PRIVATE :: evuaZ      = 0.0
      REAL, PRIVATE :: evuaR      = 0.0
      REAL, PRIVATE :: evuaX2a    = 0.0
      REAL, PRIVATE :: evuaX2b    = 0.0 
      REAL, PRIVATE :: evuaY2a    = 0.0 
      REAL, PRIVATE :: evuaY2b    = 0.0 
      REAL, PRIVATE :: evuaZ2a    = 0.0 
      REAL, PRIVATE :: evuaZ2b    = 0.0 
      REAL, PRIVATE :: evuaR2a    = 0.0 
      REAL, PRIVATE :: evuaR2b    = 0.0 
      REAL, PRIVATE :: evuaIxxA   = 0.0 
      REAL, PRIVATE :: evuaIxxB   = 0.0 
      REAL, PRIVATE :: evuaIxyA   = 0.0 
      REAL, PRIVATE :: evuaIxyB   = 0.0 
      REAL, PRIVATE :: evuaIxzA   = 0.0 
      REAL, PRIVATE :: evuaIxzB   = 0.0 
      REAL, PRIVATE :: evuaIyyA   = 0.0 
      REAL, PRIVATE :: evuaIyyB   = 0.0 
      REAL, PRIVATE :: evuaIyzA   = 0.0 
      REAL, PRIVATE :: evuaIyzB   = 0.0 
      REAL, PRIVATE :: evuaIzzA   = 0.0 
      REAL, PRIVATE :: evuaIzzB   = 0.0 
      REAL, PRIVATE :: evuaVex    = 0.0 
      REAL, PRIVATE :: evuaVmf    = 0.0 
      REAL, PRIVATE :: evuaPx     = 0.0 
      REAL, PRIVATE :: evuaPy     = 0.0 
      REAL, PRIVATE :: evuaPz     = 0.0 
      REAL, PRIVATE :: evuaPx2    = 0.0 
      REAL, PRIVATE :: evuaPy2    = 0.0 
      REAL, PRIVATE :: evuaPz2    = 0.0 
      REAL, PRIVATE :: evuaLxA    = 0.0 
      REAL, PRIVATE :: evuaLxB    = 0.0 
      REAL, PRIVATE :: evuaLyA    = 0.0 
      REAL, PRIVATE :: evuaLyB    = 0.0 
      REAL, PRIVATE :: evuaLzA    = 0.0 
      REAL, PRIVATE :: evuaLzB    = 0.0 
      REAL, PRIVATE :: evuaLx2a   = 0.0 
      REAL, PRIVATE :: evuaLx2b   = 0.0 
      REAL, PRIVATE :: evuaLy2a   = 0.0 
      REAL, PRIVATE :: evuaLy2b   = 0.0
      REAL, PRIVATE :: evuaLz2a   = 0.0 
      REAL, PRIVATE :: evuaLz2b   = 0.0 
      REAL, PRIVATE :: evuaFx     = 0.0 
      REAL, PRIVATE :: evuaFy     = 0.0 
      REAL, PRIVATE :: evuaFz     = 0.0 
      REAL, PRIVATE :: evuaTauXa  = 0.0 
      REAL, PRIVATE :: evuaTauXb  = 0.0 
      REAL, PRIVATE :: evuaTauYa  = 0.0 
      REAL, PRIVATE :: evuaTauYb  = 0.0 
      REAL, PRIVATE :: evuaTauZa  = 0.0 
      REAL, PRIVATE :: evuaTauZb  = 0.0 
      REAL, PRIVATE :: evuaEa     = 0.0 
      REAL, PRIVATE :: evuaEb     = 0.0 
      REAL, PRIVATE :: evuaMuA    = 0.0 
      REAL, PRIVATE :: evuaMuB    = 0.0 
      REAL, PRIVATE :: evuaL2a    = 0.0 
      REAL, PRIVATE :: evuaL2b    = 0.0 
      REAL, PRIVATE :: evuaTx     = 0.0 
      REAL, PRIVATE :: evuaTy     = 0.0 
      REAL, PRIVATE :: evuaTz     = 0.0 
      REAL, PRIVATE :: evuaSigXa  = 0.0 
      REAL, PRIVATE :: evuaSigXb  = 0.0 
      REAL, PRIVATE :: evuaSigYa  = 0.0 
      REAL, PRIVATE :: evuaSigYb  = 0.0
      REAL, PRIVATE :: evuaSigZa  = 0.0 
      REAL, PRIVATE :: evuaSigZb  = 0.0 
      REAL, PRIVATE :: evuaSigPx  = 0.0 
      REAL, PRIVATE :: evuaSigPy  = 0.0 
      REAL, PRIVATE :: evuaSigPz  = 0.0 
      REAL, PRIVATE :: evuaSigLxA = 0.0 
      REAL, PRIVATE :: evuaSigLxB = 0.0 
      REAL, PRIVATE :: evuaSigLyA = 0.0 
      REAL, PRIVATE :: evuaSigLyB = 0.0 
      REAL, PRIVATE :: evuaSigLzA = 0.0 
      REAL, PRIVATE :: evuaSigLzB = 0.0 

! --- VARIABLE DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     evuaL2norm is a REAL, PRIVATE variable used to store the L^2-norm of the condensate wave function.
!
!     evuaX is a REAL, PRIVATE variable used to store the condensate's average position along the X-axis of the system. It may also
!        be equivalently interpreted as the X-component of the condensate's center of mass.
!
!     evuaY is a REAL, PRIVATE variable used to store the condensate's average position along the Y-axis of the system. It may also
!        be equivalently interpreted as the Y-component of the condensate's center of mass.
!
!     evuaZ is a REAL, PRIVATE variable used to store the condensate's average position along the Z-axis of the system. It may also
!        be equivalently interpreted as the Z-component of the condensate's center of mass.
!
!     evuaR is a REAL, PRIVATE variable used to store the condensate's average radial position in the XY-plane with respect to the
!        origin of the computational domain.
!
!     evuaX2a is a REAL, PRIVATE variable used to store the condensate's average squared position along the X-axis of the system 
!        with respect to the origin of the computation domain.
!
!     evuaX2b is a REAL, PRIVATE variable used to store the condensate's average squared position along the X-axis of the system 
!        with respect to its center of mass.
!
!     evuaY2a is a REAL, PRIVATE variable used to store the condensate's average squared position along the Y-axis of the system 
!        with respect to the origin of the computation domain.
!
!     evuaY2b is a REAL, PRIVATE variable used to store the condensate's average squared position along the Y-axis of the system 
!        with respect to its center of mass.
!
!     evuaZ2a is a REAL, PRIVATE variable used to store the condensate's average squared position along the Z-axis of the system 
!        with respect to the origin of the computation domain.
!
!     evuaZ2b is a REAL, PRIVATE variable used to store the condensate's average squared position along the Z-axis of the system 
!        with respect to its center of mass.
!
!     evuaR2a is a REAL, PRIVATE variable used to store the condensate's average squared radial position in the XY-plane with 
!        respect to the origin of the computational domain.
!
!     evuaR2b is a REAL, PRIVATE variable used to store the condensate's average squared radial position in the XY-plane with 
!        respect to its center of mass.
!
!     evuaIxxA is a REAL, PRIVATE variable used to store the average XX-component of the condensate's moment of inertia tensor with
!        respect to the origin of the computational domian.
!
!     evuaIxxB is a REAL, PRIVATE variable used to store the average XX-component of the condensate's moment of inertia tensor with
!        respect to its center of mass.
!
!     evuaIxyA is a REAL, PRIVATE variable used to store the average XY-component of the condensate's moment of inertia tensor with
!        respect to the origin of the computational domian.
!
!     evuaIxyB is a REAL, PRIVATE variable used to store the average XY-component of the condensate's moment of inertia tensor with
!        respect to its center of mass.
!
!     evuaIxzA is a REAL, PRIVATE variable used to store the average  XZ-component of the condensate's moment of inertia tensor with
!        respect to the origin of the computational domian.
!
!     evuaIxzB is a REAL, PRIVATE variable used to store the average XZ-component of the condensate's moment of inertia tensor with
!        respect to its center of mass.
!
!     evuaIyyA is a REAL, PRIVATE variable used to store the average YY-component of the condensate's moment of inertia tensor with
!        respect to the origin of the computational domian.
!
!     evuaIyyB is a REAL, PRIVATE variable used to store the average YY-component of the condensate's moment of inertia tensor with
!        respect to its center of mass.
!
!     evuaIyzA is a REAL, PRIVATE variable used to store the average YZ-component of the condensate's moment of inertia tensor with
!        respect to the origin of the computational domian.
!
!     evuaIyzB is a REAL, PRIVATE variable used to store the average YZ-component of the condensate's moment of inertia tensor with
!        respect to its center of mass.
!
!     evuaIzzA is a REAL, PRIVATE variable used to store the average ZZ-component of the condensate's moment of inertia tensor with
!        respect to the origin of the computational domian.
!
!     evuaIzzB is a REAL, PRIVATE variable used to store the average ZZ-component of the condensate's moment of inertia tensor with
!        respect to its center of mass.
!
!     evuaVex is a REAL, PRIVATE variable used to store the average potential energy of the condensate due to an external potential
!        acting on it.
!
!     evuaVmf is a REAL, PRIVATE variable used to store the average potential energy of the condensate due to the mean-field 
!        interaction between its atoms.
!
!     evuaPx is a REAL, PRIVATE variable used to store the condensate's average momentum along the X-axis of the system.
!
!     evuaPy is a REAL, PRIVATE variable used to store the condensate's average momentum along the Y-axis of the system.
!
!     evuaPz is a REAL, PRIVATE variable used to store the condensate's average momentum along the Z-axis of the system.
!
!     evuaPx2 is a REAL, PRIVATE variable used to store the condensate's average squared moemntum along the X-axis of the system.
!
!     evuaPy2 is a REAL, PRIVATE variable used to store the condensate's average squared moemntum along the Y-axis of the system.
!
!     evuaPz2 is a REAL, PRIVATE variable used to store the condensate's average squared moemntum along the Z-axis of the system.
!
!     evuaLxA is a REAL, PRIVATE variable used to store the condensate's average angular momentum about the X-axis of the system 
!       with respect to the origin of the computational domain.
!
!     evuaLxB is a REAL, PRIVATE variable used to store the condensate's average angular momentum about the X-axis of the system 
!        with respect to the condensate's center of mass.
!
!     evuaLyA is a REAL, PRIVATE variable used to store the condensate's average angular momentum about the Y-axis of the system 
!        with respect to the origin of the computational domain.
!
!     evuaLyB is a REAL, PRIVATE variable used to store the condensate's average angular momentum about the Y-axis of the system 
!        with respect to the condensate's center of mass.
!
!     evuaLzA is a REAL, PRIVATE variable used to store the condensate's average angular momentum about the Z-axis of the system 
!        with respect to the origin of the computational domain.
!
!     evuaLzB is a REAL, PRIVATE variable used to store the condensate's average angular momentum about the Z-axis of the system 
!        with respect to the condensate's center of mass.
!
!     evuaLx2a is a REAL, PRIVATE variable used to store the condensate's average squared angular momentum about the X-axis of the 
!        system with respect to the origin of the computational domain.
!
!     evuaLx2b is a REAL, PRIVATE variable used to store the condensate's average squared angular momentum about the X-axis of the 
!        system with respect to the condensate's center of mass.
!
!     evuaLy2a is a REAL, PRIVATE variable used to store the condensate's average squared angular momentum about the Y-axis of the 
!        system with respect to the origin of the computational domain.
!
!     evuaLy2b is a REAL, PRIVATE variable used to store the condensate's average squared angular momentum about the Y-axis of the 
!        system with respect to the condensate's center of mass.
!
!     evuaLz2a is a REAL, PRIVATE variable used to store the condensate's average squared angular momentum about the Z-axis of the 
!        system with respect to the origin of the computational domain.
!
!     evuaLz2b is a REAL, PRIVATE variable used to store the condensate's average squared angular momentum about the Z-axis of the 
!        system with respect to the condensate's center of mass.
!
!     evuaFx is a REAL, PRIVATE variable used to store the average force on the condensate along the X-axis of the system due to an
!        external potential acting on the condenate.
!
!     evuaFy is a REAL, PRIVATE variable used to store the average force on the condensate along the Y-axis of the system due to an
!        external potential acting on the condenate.
!
!     evuaFz is a REAL, PRIVATE variable used to store the average force on the condensate along the Z-axis of the system due to an
!        external potential acting on the condenate.
!
!     evuaTauXa is a REAL, PRIVATE variable used to store the average torque on the condensate about the X-axis of the system due to
!        and external potential acting on the condensate and measured with respect to the origin of the computational domain.
!
!     evuaTauXb is a REAL, PRIVATE variable used to store the average torque on the condensate about the X-axis of the system due to
!        and external potential acting on the condensate and measured with respect to its center of mass.
!
!     evuaTauYa is a REAL, PRIVATE variable used to store the average torque on the condensate about the Y-axis of the system due to
!        and external potential acting on the condensate and measured with respect to the origin of the computational domain.
!
!     evuaTauYb is a REAL, PRIVATE variable used to store the average torque on the condensate about the Y-axis of the system due to
!        and external potential acting on the condensate and measured with respect to its center of mass.
!
!     evuaTauZa is a REAL, PRIVATE variable used to store the average torque on the condensate about the Z-axis of the system due to
!        and external potential acting on the condensate and measured with respect to the origin of the computational domain.
!
!     evuaTauZb is a REAL, PRIVATE variable used to store the average torque on the condensate about the Z-axis of the system due to
!        and external potential acting on the condensate and measured with respect to its center of mass.
!
!     evuaEa is a REAL, PRIVATE variable used to store the condensate's average total energy measured with respect to the origin of
!        the computational domain.
!
!     evuaEb is a REAL, PRIVATE variable used to store the condensate's average total energy measured with respect to its center of
!        mass.
!
!     evuaMuA is a REAL, PRIVATE variable used to store the chemical potential of the condensate as measured with respect to the 
!        origin of the computational domain.
!      
!     evuaMuB is a REAL, PRIVATE variable used to store the chemical potential of the condensate as measured with respect to its
!        center of mass. 
!
!     evuaL2a is a REAL, PRIVATE variable used to store the average squared total angular momentum of the condensate as measured
!        with respect to the origin of the computational domain.
!
!     evuaL2b is a REAL, PRIVATE variable used to store the average squared total angular momentum of the condensate as measured
!        with respect to its center of mass.
!
!     evuaTx is a REAL, PRIVATE variable used to store the average kinetic energy of the condensate along the X-axis of the system.
!
!     evuaTy is a REAL, PRIVATE variable used to store the average kinetic energy of the condensate along the Y-axis of the system.
!
!     evuaTz is a REAL, PRIVATE variable used to store the average kinetic energy of the condensate along the Z-axis of the system.
!
!     evuaSigXa is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the X-axis of the system and 
!        measured with respect to the origin of the computational domain.
!
!     evuaSigXa is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the X-axis of the system and 
!        measured with respect to the origin of the computational domain.
!
!     evuaSigXb is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the X-axis of the system and 
!        measured with respect to the condensate's center of mass.
!
!     evuaSigYa is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the Y-axis of the system and 
!        measured with respect to the origin of the computational domain.
!
!     evuaSigYb is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the Y-axis of the system and 
!        measured with respect to the condensate's center of mass.
!
!     evuaSigZa is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the Z-axis of the system and 
!        measured with respect to the origin of the computational domain.
!
!     evuaSigZb is a REAL, PRIVATE variable used to store the condensate's position uncertainty along the Z-axis of the system and 
!        measured with respect to the condensate's center of mass.
!
!     evuaSigPx is a REAL, PRIVATE variable used to store the condensate's momentum uncertainty along the X-axis of the system.
!
!     evuaSigPy is a REAL, PRIVATE variable used to store the condensate's momentum uncertainty along the Y-axis of the system.
!
!     evuaSigPz is a REAL, PRIVATE variable used to store the condensate's momentum uncertainty along the Z-axis of the system.
!
!     evuaSigLxA is a REAL, PRIVATE variable used to store the uncertainty in the angular momentum of the condensate about the 
!        X-axis of the system and measured with respect to the origin of the computational domain.
!
!     evuaSigLxB is a REAL, PRIVATE variable used to store the uncertainty in the angular momentum of the condensate about the 
!        X-axis of the system and measured with respect to the condensate's center of mass.
!
!     evuaSigLyA is a REAL, PRIVATE variable used to store the uncertainty in the angular momentum of the condensate about the 
!        Y-axis of the system and measured with respect to the origin of the computational domain.
!
!     evuaSigLyB is a REAL, PRIVATE variable used to store the uncertainty in the angular momentum of the condensate about the 
!        Y-axis of the system and measured with respect to the condensate's center of mass.
!
!     evuaSigLzA is a REAL, PRIVATE variable used to store the uncertainty in the angular momentum of the condensate about the 
!        Z-axis of the system and measured with respect to the origin of the computational domain.
!
!     evuaSigLzB is a REAL, PRIVATE variable used to store the uncertainty in the angular momentum of the condensate about the 
!        Z-axis of the system and measured with respect to the condensate's center of mass.
!
! --- SUBROUTINE DECLARATIONS ------------------------------------------------------------------------------------------------------
       
      PUBLIC :: evua_compute_base
      PUBLIC :: evua_compute_derived
      PUBLIC :: evua_write_all
      PUBLIC :: evua_normalize

! --- SUBROUTINE DEFINITIONS -------------------------------------------------------------------------------------------------------
!
!     evua_compute_base is a PUBLIC SUBROUTINE that calls functions which computes base quantities associated with a wave function.
!
!     evua_compute_derived is a PUBLIC SUBROUTINE that computes quantities which are derived from base quantities associated with a
!        wave function.
!
!     evua_write_all is a PUBLIC SUBROUTINE that formats and writes all quantities, both base and derived, to standard output.
!
!     evua_normalize is a PULIC SUBROUTINE that may be used to normalize a wave function. It is currently used only to re-normalize
!        the wave function after every time step when performing imaginary time propagation.
!
! --- FUNCTION DECLARATIONS --------------------------------------------------------------------------------------------------------

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

! --- FUNCTION DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     evua_l2_norm_3d_rect is a PRIVATE, REAL scalar function that computes the L^2-norm of a three-dimensional wave function 
!        defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_x_3d_rect is a PRIVATE, REAL scalar function that computes the x-component of the position expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_y_3d_rect is a PRIVATE, REAL scalar function that computes the y-component of the position expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_z_3d_rect is a PRIVATE, REAL scalar function that computes the z-component of the position expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_r_xy_3d_rect is a PRIVATE, REAL scalar function that computes the planar (xy-plane) radial position expectation value 
!        for a three-dimensinal wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_r_xyz_3d_rect is a PRIVATE, REAL scalar function that computes the spherical radial position expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_x2_3d_rect is a PRIVATE, REAL scalar function that computes the x-component of the squared position expectation value 
!        for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_y2_3d_rect is a PRIVATE, REAL scalar function that computes the y-component of the squared position expectation value 
!        for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_z2_3d_rect is a PRIVATE, REAL scalar function that computes the z-component of the squared position expectation value 
!        for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_r2_xy_3d_rect is a PRIVATE, REAL scalar function that computes the squared planar (xy-plane) radial position expectation
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_r2_xyz_3d_rect is a PRIVATE, REAL scalar function that computes the squared spherical radial position expectation value 
!        for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_px_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the x-component of the momentum expectation value for a
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order central 
!        differences.
!
!     evua_px_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the x-component of the momentum expectation value for a
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order central 
!        differences.
!
!     evua_py_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the y-component of the momentum expectation value for a
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order central 
!        differences.
!
!     evua_py_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the y-component of the momentum expectation value for a
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order central 
!        differences.
!
!     evua_pz_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the z-component of the momentum expectation value for a
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order central 
!        differences.
!
!     evua_pz_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the z-component of the momentum expectation value for a
!        three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order central 
!        differences.
!
!     evua_px2_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the x-component of the squared momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order 
!        central differences.
!
!     evua_px2_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the x-component of the squared momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order 
!        central differences.
!
!     evua_py2_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the y-component of the squared momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order 
!        central differences.
!
!     evua_py2_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the y-component of the squared momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order
!        central differences.
!
!     evua_pz2_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the z-component of the squared momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order 
!        central differences.
!
!     evua_pz2_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the z-component of the squared momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order 
!        central differences.
!
!     evua_lx_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the x-component of the angular momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order
!        central differences.
!
!     evua_lx_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the x-component of the angular momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order
!        central differences.
!
!     evua_ly_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the y-component of the angular momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order 
!        central differences.
!
!     evua_ly_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the y-component of the angular momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order 
!        central differences.
!
!     evua_lz_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the z-component of the angular momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 2nd-order
!        central differences.
!
!     evua_lz_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the z-component of the angular momentum expectation 
!        value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 4th-order
!        central differences.
!
!     evua_lx2_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the x-component of the squared angular momentum 
!        expectation value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 
!        2nd-order central differences.
!
!     evua_lx2_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the x-component of the squared angular momentum 
!        expectation value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 
!        4th-order central differences.
!
!     evua_ly2_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the y-component of the squared angular momentum 
!        expectation value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and
!        2nd-order central differences.
!
!     evua_ly2_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the y-component of the squared angular momentum 
!        expectation value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 
!        4th-order central differences.
!
!     evua_lz2_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the z-component of the squared angular momentum 
!        expectation value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 
!        2nd-order central differences.
!
!     evua_lz2_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the z-component of the squared angular momentum 
!        expectation value for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule and 
!        4th-order central differences.
!
!     evua_fx_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the x-component of the force expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid and acted upon by an applied external potential. It 
!        utilizes the rectangle rule and 2nd-order central differences.
!
!     evua_fx_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the x-component of the force expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid and acted upon by an applied external potential. It 
!        utilizes the rectangle rule and 4th-order central differences.
!
!     evua_fy_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the y-component of the force expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid and acted upon by an applied external potential. It 
!        utilizes the rectangle rule and 2nd-order central differences.
!
!     evua_fy_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the y-component of the force expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid and acted upon by an applied external potential. It 
!        utilizes the rectangle rule and 4th-order central differences. 
!
!     evua_fz_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the z-component of the force expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid and acted upon by an applied external potential. It 
!        utilizes the rectangle rule and 2nd-order central differences.
!
!     evua_fz_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the z-component of the force expectation value for a 
!        three-dimensional wave function defined on a regular Cartesian grid and acted upon by an applied external potential. It
!        utilizes the rectangle rule and 4th-order central differences. 
!
!     evua_taux_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the x-component of the torque expectation value due to
!        an applied external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the 
!        rectangle rule and 2nd-order central differences.
!
!     evua_taux_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the x-component of the torque expectation value due to
!        an applied external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the 
!        rectangle rule and 4th-order central differences.
!
!     evua_tauy_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the y-component of the torque expectation value due to
!        an applied external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the 
!        rectangle rule and 2nd-order central differences.
!
!     evua_tauy_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the y-component of the torque expectation value due to
!        an applied external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the 
!        rectangle rule and 4th-order central differences.
!
!     evua_tauz_3d_rect_cd2 is a PRIVATE, REAL scalar function that computes the z-component of the torque expectation value due to
!        an applied external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the 
!        rectangle rule and 2nd-order central differences.
!
!     evua_tauz_3d_rect_cd4 is a PRIVATE, REAL scalar function that computes the z-component of the torque expectation value due to
!        an applied external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the 
!        rectangle rule and 4th-order central differences.
!
!     evua_ixx_3d_rect is a PRIVATE, REAL scalar function that computes the expectation value of the xx-component of the moment of 
!        inertia tensor for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule.
!
!     evua_iyy_3d_rect is a PRIVATE, REAL scalar function that computes the expectation value of the yy-component of the moment of 
!        inertia tensor for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_izz_3d_rect is a PRIVATE, REAL scalar function that computes the expectation value of the zz-component of the moment of 
!        inertia tensor for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_ixy_3d_rect is a PRIVATE, REAL scalar function that computes the expectation value of the xy-component of the moment of 
!        inertia tensor for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_iyz_3d_rect is a PRIVATE, REAL scalar function that computes the expectation value of the yz-component of the moment of 
!        inertia tensor for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
!
!     evua_ixz_3d_rect is a PRIVATE, REAL scalar function that computes the expectation value of the xz-component of the moment of 
!        inertia tensor for a three-dimensional wave function defined on a regular Cartesian grid using the rectangle rule. 
! 
!     evua_vex_3d_rect is a PRIVATE, REAL scalar function that computes the potential energy expectation value due to an applied 
!        external potential acting on a three-dimensional wave function defined on a regular Cartesian grid using the rectangle 
!        rule.
!
!     evua_vmf_3d_rect is a PRIVATE, REAL scalar function that computes the potential energy expectation value due to the mean-field
!        interactions between atoms in a Bose-Einstein condensate for a three-dimensional condensate wave function defined on a 
!        regular Cartesian grid using the rectangle rule.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE evua_compute_base ( mpiMaster , mpiReal , mpiError , evuaQuadRule , evuaFdOrder , nXa , nXb , nXbc , nYa , nYb , &
         & nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , gS , X , Y , Z , Vex3 , Psi3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: mpiMaster
      INTEGER, INTENT ( IN ) :: mpiReal
      INTEGER, INTENT ( INOUT ) :: mpiError
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

      REAL :: temp

      IF ( evuaQuadRule == 1 ) THEN

         temp = evua_l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
         CALL MPI_REDUCE ( temp , evuaL2Norm , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_x_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )
         CALL MPI_REDUCE ( temp , evuaX , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
         CALL MPI_BCAST ( evuaX , 1 , mpiReal , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , dX , dY , dZ , X , Psi3 )
         CALL MPI_REDUCE ( temp , evuaX2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , dX , dY , dZ , X , Psi3 )
         CALL MPI_REDUCE ( temp , evuaX2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_y_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaY , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
         CALL MPI_BCAST ( evuaY , 1 , mpiReal , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , dX , dY , dZ , Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaY2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , dX , dY , dZ , Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaY2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_z_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )
         CALL MPI_REDUCE ( temp , evuaZ , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
         CALL MPI_BCAST ( evuaZ , 1 , mpiReal , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , zO , dX , dY , dZ , Z , Psi3 )
         CALL MPI_REDUCE ( temp , evuaZ2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaZ , dX , dY , dZ , Z , Psi3 )
         CALL MPI_REDUCE ( temp , evuaZ2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_r_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaR , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
         CALL MPI_BCAST ( evuaR , 1 , mpiReal , mpiMaster , MPI_COMM_WORLD , mpiError ) ! is this BCAST needed? 

         temp = evua_r2_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaR2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_r2_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , X , &
            & Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaR2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_ixx_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )
         CALL MPI_REDUCE ( temp , evuaIxxA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_ixx_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ , dX , dY , dZ , Y , Z , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaIxxB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_ixy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaIxyA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_ixy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , X , Y , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaIxyB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_ixz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 )
         CALL MPI_REDUCE ( temp , evuaIxzA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_ixz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , X , Z , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaIxzB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_iyy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , Psi3 ) 
         CALL MPI_REDUCE ( temp , evuaIyyA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_iyy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , X , Z , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaIyyB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_iyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , Psi3 )
         CALL MPI_REDUCE ( temp , evuaIyzA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_iyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ , dX , dY , dZ , Y , Z , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaIyzB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_izz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , Psi3 )
         CALL MPI_REDUCE ( temp , evuaIzzA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_izz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , X , Y , &
            & Psi3 )
         CALL MPI_REDUCE ( temp , evuaIzzB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_vex_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Vex3 , Psi3 )
         CALL MPI_REDUCE ( temp , evuaVex , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         temp = evua_vmf_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , gS , Psi3 )
         CALL MPI_REDUCE ( temp , evuaVmf , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         IF ( evuaFdOrder == 2 ) THEN

            temp = evua_px_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPx , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_py_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPy , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_pz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPz , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_px2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPx2 , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_py2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPy2 , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_pz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPz2 , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , Z ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLxA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ ,  dX , dY , dZ , &
               & Y , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLxB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , &
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLyA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , &
               & X , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLyB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , &
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLzA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , &
               & X , Y , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLzB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLx2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ , dX , dY , dZ , &
               & Y , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLx2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLy2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , &
               & X , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLy2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLz2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , &
               & X , Y , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLz2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_fx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaFx , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_fy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaFy , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_fz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaFz , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_taux_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , &
               & Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauXa , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_taux_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ , dX , dY , dZ , &
               & Y , Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauXb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , &
               & Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauYa , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , &
               & X , Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauYb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , &
               & Y , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauZa , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , &
               & X , Y , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauZb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

         ELSE IF ( evuaFdOrder == 4 ) THEN

            temp = evua_px_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPx , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_py_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPy , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_pz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPz , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_px2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPx2 , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_py2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPy2 , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_pz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3 )
            CALL MPI_REDUCE ( temp , evuaPz2 , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , Z ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLxA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ ,  dX , dY , dZ , &
               & Y , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLxB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , &
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLyA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , &
               & X , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLyB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , &
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLzA , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , &
               & X , Y , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLzB , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLx2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lx2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ , dX , dY , dZ , &
               & Y , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLx2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLy2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_ly2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , &
               & X , Z , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLy2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y ,&
               & Psi3 )
            CALL MPI_REDUCE ( temp , evuaLz2a , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_lz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , &
               & X , Y , Psi3 )
            CALL MPI_REDUCE ( temp , evuaLz2b , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_fx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaFx , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_fy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaFy , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_fz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaFz , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_taux_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , &
               & Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauXa , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_taux_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaY , evuaZ , dX , dY , dZ , &
               & Y , Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauXb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , &
               & Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauYa , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaZ , dX , dY , dZ , &
               & X , Z , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauYb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , &
               & Y , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauZa , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

            temp = evua_tauz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , evuaX , evuaY , dX , dY , dZ , &
               & X , Y , Vex3 , Psi3 )
            CALL MPI_REDUCE ( temp , evuaTauZb , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )

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

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE evua_compute_derived ( mpiRank , mpiMaster , wX , wY , wZ )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: mpiRank
      INTEGER, INTENT ( IN ) :: mpiMaster
 
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ

      IF ( mpiRank == mpiMaster ) THEN

!        Kinetic energy expectation values

         evuaTx = 0.5 * evuaPx2
         evuaTy = 0.5 * evuaPy2
         evuaTz = 0.5 * evuaPz2

!        Energy expectation value

         evuaEa = evuaTx + evuaTy + evuaTz + evuaVex + evuaVmf - wX * evuaLxA - wY * evuaLyA - wZ * evuaLzA
         evuaEb = evuaTx + evuaTy + evuaTz + evuaVex + evuaVmf - wX * evuaLxB - wY * evuaLyB - wZ * evuaLzB

!        Chemical potential

         evuaMuA =  evuaTx + evuaTy + evuaTz + evuaVex + 2.0 * evuaVmf - wX * evuaLxA - wY * evuaLyA - wZ * evuaLzA
         evuaMuB =  evuaTx + evuaTy + evuaTz + evuaVex + 2.0 * evuaVmf - wX * evuaLxB - wY * evuaLyB - wZ * evuaLzB

!        Squared angular momentum expectation value

         evuaL2a = evuaLx2a + evuaLy2a + evuaLz2a
         evuaL2b = evuaLx2b + evuaLy2b + evuaLz2b

!        Position, momentum and angular momentum uncertainty 
           
         evuaSigXa = SQRT ( evuaX2a  - evuaX**2 )
         evuaSigXb = SQRT ( evuaX2b  - evuaX**2 )

         evuaSigYa = SQRT ( evuaY2a  - evuaY**2 )
         evuaSigYb = SQRT ( evuaY2b  - evuaY**2 )

         evuaSigZa = SQRT ( evuaZ2a  - evuaZ**2 )
         evuaSigZb = SQRT ( evuaZ2b  - evuaZ**2 )

         evuaSigPx = SQRT ( evuaPx2 - evuaPx**2 )
         evuaSigPy = SQRT ( evuaPy2 - evuaPy**2 )
         evuaSigPz = SQRT ( evuaPz2 - evuaPz**2 )

         evuaSigLxA = SQRT ( evuaLx2a - evuaLxA**2 )
         evuaSigLxB = SQRT ( evuaLx2b - evuaLxB**2 )

         evuaSigLyA = SQRT ( evuaLy2a - evuaLyA**2 )
         evuaSigLyB = SQRT ( evuaLy2b - evuaLyB**2 )

         evuaSigLzA = SQRT ( evuaLz2a - evuaLzA**2 )
         evuaSigLzB = SQRT ( evuaLz2b - evuaLzB**2 )

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE evua_write_all ( mpiRank , mpiMaster , tN , wX , wY , wZ )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: mpiRank
      INTEGER, INTENT ( IN ) :: mpiMaster

      REAL, INTENT ( IN ) :: tN
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ

!     Write expectation values, uncertainties and uncertainty relations to file from mpiMaster

      IF ( mpiRank == mpiMaster ) THEN
     
         WRITE ( UNIT = OUTPUT_UNIT , FMT = '(82(F23.15))' ) tN , wX , wY , wZ , evuaL2Norm , evuaEa , evuaEb , evuaMuA , evuaMuB ,&
            & evuaL2a , evuaL2b , evuaTx , evuaTy , evuaTz , evuaVex , evuaVmf , evuaX , evuaX2a , evuaSigXa , evuaX2b , &
            & evuaSigXb , evuaPx , evuaPx2 , evuaSigPx , evuaLxA , evuaLx2a , evuaSigLxA , evuaLxB , evuaLx2b , evuaSigLxB , &
            & evuaFx , evuaTauXa , evuaTauXb , evuaY , evuaY2a , evuaSigYa , evuaY2b , evuaSigYb , evuaPy , evuaPy2 , evuaSigPy , &
            & evuaLyA , evuaLy2a , evuaSigLyA , evuaLyB , evuaLy2b , evuaSigLyB , evuaFy , evuaTauYa , evuaTauYb , evuaZ , &
            & evuaZ2a , evuaSigZa , evuaZ2b , evuaSigZb , evuaPz , evuaPz2 , evuaSigPz , evuaLzA , evuaLz2a , evuaSigLzA , &
            & evuaLzB , evuaLz2b , evuaSigLzB , evuaFz , evuaTauZa , evuaTauZb , evuaR , evuaR2a , evuaR2b , evuaIxxA , evuaIxyA , &
            & evuaIxzA , evuaIyyA , evuaIyzA , evuaIzzA , evuaIxxB , evuaIxyB , evuaIxzB , evuaIyyB , evuaIyzB , evuaIzzB

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE evua_normalize ( mpiMaster , mpiReal , mpiError , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , &
         & dY , dZ , Psi3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: mpiMaster
      INTEGER, INTENT ( IN ) :: mpiReal
      INTEGER, INTENT ( INOUT ) :: mpiError
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

      REAL :: temp

      temp = evua_l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ ,Psi3 )
      CALL MPI_REDUCE ( temp , evuaL2Norm , 1 , mpiReal , MPI_SUM , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( evuaL2Norm , 1 , mpiReal , mpiMaster , MPI_COMM_WORLD , mpiError )
      Psi3 = Psi3 / SQRT ( evuaL2Norm )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_l2_norm_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_l2_norm_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_l2_norm_3d_rect = evua_l2_norm_3d_rect + ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_l2_norm_3d_rect = evua_l2_norm_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb  

               evua_x_3d_rect = evua_x_3d_rect + ( X ( j ) - xO ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_x_3d_rect = evua_x_3d_rect * dX * dY * dZ 

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_y_3d_rect = evua_y_3d_rect + ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_y_3d_rect = evua_y_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_z_3d_rect = evua_z_3d_rect + ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_z_3d_rect = evua_z_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_r_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xy_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xy_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_r_xy_3d_rect = evua_r_xy_3d_rect + SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * &
                  & ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_r_xy_3d_rect = evua_r_xy_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_r_xyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , X ,&
         & Y , Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xyz_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r_xyz_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_r_xyz_3d_rect = evua_r_xyz_3d_rect + SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) *&
                  & ABS ( Psi3 ( j , k , l ) )**2 

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_r_xyz_3d_rect = evua_r_xyz_3d_rect * dX * dY * dZ 

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x2_3d_rect )
      DO l = nZa , nZb 

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_x2_3d_rect )
         DO k = nYa , nYb 

            DO j = nXa , nXb  

               evua_x2_3d_rect = evua_x2_3d_rect + ( X ( j ) - xO )**2 * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_x2_3d_rect = evua_x2_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y2_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_y2_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_y2_3d_rect = evua_y2_3d_rect + ( Y ( k ) - yO )**2 * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_y2_3d_rect = evua_y2_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z2_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_z2_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_z2_3d_rect = evua_z2_3d_rect + ( Z ( l ) - zO )**2 * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_z2_3d_rect = evua_z2_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_r2_xy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO ,  dX , dY , dZ , X , Y ,&
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xy_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xy_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_r2_xy_3d_rect = evua_r2_xy_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * &
                  & ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_r2_xy_3d_rect = evua_r2_xy_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_r2_xyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , dX , dY , dZ , &
         & X , Y , Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xyz_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_r2_xyz_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_r2_xyz_3d_rect = evua_r2_xyz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * &
                  & ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_r2_xyz_3d_rect = evua_r2_xyz_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_px_3d_rect_cd2 = evua_px_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - &
                  & Psi3 ( j - 1 , k , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_px_3d_rect_cd2 = -0.5 * evua_px_3d_rect_cd2 * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd4 ) 
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px_3d_rect_cd4 ) 
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_px_3d_rect_cd4 = evua_px_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + &
                  & CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_px_3d_rect_cd4 = -0.08333333333333333 * evua_px_3d_rect_cd4 * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_py_3d_rect_cd2 = evua_py_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - &
                  & Psi3 ( j , k - 1 , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_py_3d_rect_cd2 = -0.5 * evua_py_3d_rect_cd2 * dX * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_py_3d_rect_cd4 = evua_py_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + &
                  & CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k -1 , l ) ) + Psi3 ( j , k - 2 , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_py_3d_rect_cd4 = -0.08333333333333333 * evua_py_3d_rect_cd4 * dX * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd2 ) 
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_pz_3d_rect_cd2 = evua_pz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - &
                  & Psi3 ( j , k , l - 1 ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_pz_3d_rect_cd2 = -0.5 * evua_pz_3d_rect_cd2 * dX * dY

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_pz_3d_rect_cd4 = evua_pz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + &
                  & CMPLX ( 8.0 , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_pz_3d_rect_cd4 = -0.08333333333333333 * evua_pz_3d_rect_cd4 * dX * dY

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd2 ) 
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_px2_3d_rect_cd2 = evua_px2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j + 1 , k , l ) - &
                  & CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j - 1 , k , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_px2_3d_rect_cd2 = -evua_px2_3d_rect_cd2 * dY * dZ / dX

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_px2_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_px2_3d_rect_cd4 = evua_px2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j + 2 , k , l ) + &
                  & CMPLX ( 16.0 , 0.0 ) * Psi3 ( j + 1 , k , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                  & CMPLX ( 16.0 , 0.0 ) * Psi3 ( j - 1 , k , l ) - Psi3 ( j - 2 , k , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_px2_3d_rect_cd4 = -0.08333333333333333 * evua_px2_3d_rect_cd4 * dY * dZ / dX

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_py2_3d_rect_cd2 = evua_py2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k + 1 , l ) - &
                  & CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k - 1 , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_py2_3d_rect_cd2 = -evua_py2_3d_rect_cd2 * dX * dZ / dY

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_py2_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_py2_3d_rect_cd4 = evua_py2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k + 2 , l ) + &
                  & CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k + 1 , l ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                  & CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k - 1 , l ) - Psi3 ( j , k - 2 , l ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_py2_3d_rect_cd4 = -0.08333333333333333 * evua_py2_3d_rect_cd4 * dX * dZ / dY

      RETURN
    
      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_pz2_3d_rect_cd2 = evua_pz2_3d_rect_cd2 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( Psi3 ( j , k , l + 1 ) - &
                  & CMPLX ( 2.0 , 0.0 ) * Psi3 ( j , k , l ) + Psi3 ( j , k , l - 1 ) ) ) 

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_pz2_3d_rect_cd2 = -evua_pz2_3d_rect_cd2 * dX * dY / dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_pz2_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_pz2_3d_rect_cd4 = evua_pz2_3d_rect_cd4 + REAL ( CONJG ( Psi3 ( j , k , l ) ) * ( -Psi3 ( j , k , l + 2 ) + &
                  & CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l + 1 ) - CMPLX ( 30.0 , 0.0 ) * Psi3 ( j , k , l ) + &
                  & CMPLX ( 16.0 , 0.0 ) * Psi3 ( j , k , l - 1 ) - Psi3 ( j , k , l - 2 ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_pz2_3d_rect_cd4 = -0.08333333333333333 * evua_pz2_3d_rect_cd4 * dX * dY / dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO ,  dX , dY , dZ , Y , &
         & Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_lx_3d_rect_cd2 = evua_lx_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( &
                  & CMPLX ( ( Y ( k ) - yO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) - &
                  & CMPLX ( ( Z ( l ) - zO ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lx_3d_rect_cd2 = -0.5 * evua_lx_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lx_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z ,&
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_lx_3d_rect_cd4 = evua_lx_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( &
                  & CMPLX ( ( Y ( k ) - yO ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( &
                  & Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) - &
                  & CMPLX ( ( Z ( l ) - zO ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( &
                  & Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lx_3d_rect_cd4 = -0.08333333333333333 * evua_lx_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z ,&
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_ly_3d_rect_cd2 = evua_ly_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( &
                  & CMPLX ( ( Z ( l ) - zO ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) - &
                  & CMPLX ( ( X ( j ) - xO ) / dZ , 0.0 ) * ( Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ly_3d_rect_cd2 = -0.5 * evua_ly_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ly_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z ,&
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_ly_3d_rect_cd4 = evua_ly_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( &
                  & CMPLX ( ( Z ( l ) - zO ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( &
                  & Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) - &
                  & CMPLX ( ( X ( j ) - xO ) / dZ , 0.0 ) * ( -Psi3 ( j , k , l + 2 ) + CMPLX ( 8.0 , 0.0 ) * ( &
                  & Psi3 ( j , k , l + 1 ) - Psi3 ( j , k , l - 1 ) ) + Psi3 ( j , k , l - 2 ) ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ly_3d_rect_cd4 = -0.08333333333333333 * evua_ly_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y ,&
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_lz_3d_rect_cd2 = evua_lz_3d_rect_cd2 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( &
                  & CMPLX ( ( X ( j ) - xO ) / dY , 0.0 ) * ( Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) - &
                  & CMPLX ( ( Y ( k ) - yO ) / dX , 0.0 ) * ( Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lz_3d_rect_cd2 = -0.5 * evua_lz_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y ,&
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_lz_3d_rect_cd4 = evua_lz_3d_rect_cd4 + REAL ( I * CONJG ( Psi3 ( j , k , l ) ) * ( &
                  & CMPLX ( ( X ( j ) - xO ) / dY , 0.0 ) * ( -Psi3 ( j , k + 2 , l ) + CMPLX ( 8.0 , 0.0 ) * ( &
                  & Psi3 ( j , k + 1 , l ) - Psi3 ( j , k - 1 , l ) ) + Psi3 ( j , k - 2 , l ) ) - &
                  & CMPLX ( ( Y ( k ) - yO ) / dX , 0.0 ) * ( -Psi3 ( j + 2 , k , l ) + CMPLX ( 8.0 , 0.0 ) * ( &
                  & Psi3 ( j + 1 , k , l ) - Psi3 ( j - 1 , k , l ) ) + Psi3 ( j - 2 , k , l ) ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lz_3d_rect_cd4 = -0.08333333333333333 * evua_lz_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , &
         & Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd2 )
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
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lx2_3d_rect_cd2 = -evua_lx2_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lx2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , &
         & Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lx2_3d_rect_cd4 )
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
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lx2_3d_rect_cd4 = -0.08333333333 * evua_lx2_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , &
         & Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd2 )
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
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ly2_3d_rect_cd2 = -evua_ly2_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ly2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , &
         & Z , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ly2_3d_rect_cd4 )
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
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ly2_3d_rect_cd4 = -0.08333333333 * evua_ly2_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , &
         & Y , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd2 )
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
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lz2_3d_rect_cd2 = -evua_lz2_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_lz2_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , &
         & Y , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_lz2_3d_rect_cd4 )
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
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_lz2_3d_rect_cd4 = -0.08333333333 * evua_lz2_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_fx_3d_rect_cd2 = evua_fx_3d_rect_cd2 + ( Vex3 ( j - 1 , k , l ) - Vex3 ( j + 1 , k , l ) ) * &
                  & ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_fx_3d_rect_cd2 = 0.5 * evua_fx_3d_rect_cd2 * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fx_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_fx_3d_rect_cd4 = evua_fx_3d_rect_cd4 + ( Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - &
                  & Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_fx_3d_rect_cd4 = 0.08333333333333333 * evua_fx_3d_rect_cd4 * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_fy_3d_rect_cd2 = evua_fy_3d_rect_cd2 + ( Vex3 ( j , k - 1 , l ) - Vex3 ( j , k + 1 , l ) ) * &
                  & ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_fy_3d_rect_cd2 = 0.5 * evua_fy_3d_rect_cd2 * dX * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fy_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_fy_3d_rect_cd4 = evua_fy_3d_rect_cd4 + ( Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - &
                  & Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2 , l ) ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_fy_3d_rect_cd4 = 0.08333333333333333 * evua_fy_3d_rect_cd4 * dX * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_fz_3d_rect_cd2 = evua_fz_3d_rect_cd2 + ( Vex3 ( j , k , l - 1 ) - Vex3 ( j , k , l + 1 ) ) * &
                  & ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_fz_3d_rect_cd2 = 0.5 * evua_fz_3d_rect_cd2 * dX * dY

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_fz_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_fz_3d_rect_cd4 = evua_fz_3d_rect_cd4 + ( Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - &
                  & Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_fz_3d_rect_cd4 = 0.08333333333333333 * evua_fz_3d_rect_cd4 * dX * dY

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_taux_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , &
         & Z , Vex3 , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_taux_3d_rect_cd2 = evua_taux_3d_rect_cd2 + ( ( Y ( k ) - yO ) * ( Vex3 ( j , k , l - 1 ) - &
                  & Vex3 ( j , k , l + 1 ) ) / dZ - ( Z ( l ) - zO ) * ( Vex3 ( j , k - 1 , l ) - &
                  & Vex3 ( j , k + 1 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_taux_3d_rect_cd2 = 0.5 * evua_taux_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_taux_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , &
         & Z , Vex3 , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_taux_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_taux_3d_rect_cd4 = evua_taux_3d_rect_cd4 + ( ( Y ( k ) - yO ) * ( Vex3 ( j , k , l + 2 ) - 8.0 * ( &
                  & Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - Vex3 ( j , k , l - 2 ) ) / dZ - ( Z ( l ) - zO ) * ( &
                  & Vex3 ( j , k + 2 , l ) - 8.0 * ( Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - &
                  & Vex3 ( j , k - 2 , l ) ) / dY ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_taux_3d_rect_cd4 = 0.08333333333333333 * evua_taux_3d_rect_cd4 * dX * dY * dZ 

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_tauy_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , &
         & Z , Vex3 , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_tauy_3d_rect_cd2 = evua_tauy_3d_rect_cd2 + ( ( Z ( l ) - zO ) * ( Vex3 ( j - 1 , k , l ) - &
                  & Vex3 ( j + 1 , k , l ) ) / dX - ( X ( j ) - xO ) * ( Vex3 ( j , k , l - 1 ) - &
                  & Vex3 ( j , k , l + 1 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_tauy_3d_rect_cd2 = 0.5 * evua_tauy_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_tauy_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , &
         & Z , Vex3 , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauy_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_tauy_3d_rect_cd4 = evua_tauy_3d_rect_cd4 + ( ( Z ( l ) - zO ) * ( Vex3 ( j + 2 , k , l ) - 8.0 * ( &
                  & Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - Vex3 ( j - 2 , k , l ) ) / dX - ( X ( j ) - xO ) * ( &
                  & Vex3 ( j , k , l + 2 ) - 8.0 * ( Vex3 ( j , k , l + 1 ) - Vex3 ( j , k , l - 1 ) ) - &
                  & Vex3 ( j , k , l - 2 ) ) / dZ ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_tauy_3d_rect_cd4 = 0.08333333333333333 * evua_tauy_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_tauz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , &
         & Y , Vex3 , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd2 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd2 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_tauz_3d_rect_cd2 = evua_tauz_3d_rect_cd2 + ( ( X ( j ) - xO ) * ( Vex3 ( j , k - 1 , l ) - &
                  & Vex3 ( j , k + 1 , l ) ) / dY - ( Y ( k ) - yO ) * ( Vex3 ( j - 1 , k , l ) - &
                  & Vex3 ( j + 1 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_tauz_3d_rect_cd2 = 0.5 * evua_tauz_3d_rect_cd2 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_tauz_3d_rect_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , &
         & Y , Vex3 , Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd4 )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_tauz_3d_rect_cd4 )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_tauz_3d_rect_cd4 = evua_tauz_3d_rect_cd4 + ( ( X ( j ) - xO ) * ( Vex3 ( j , k + 2 , l ) - 8.0 * ( &
                  & Vex3 ( j , k + 1 , l ) - Vex3 ( j , k - 1 , l ) ) - Vex3 ( j , k - 2, l ) ) / dY - ( Y ( k ) - yO ) * ( &
                  & Vex3 ( j + 2 , k , l ) - 8.0 * ( Vex3 ( j + 1 , k , l ) - Vex3 ( j - 1 , k , l ) ) - &
                  & Vex3 ( j - 2 , k , l ) ) / dX ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_tauz_3d_rect_cd4 = 0.08333333333333333 * evua_tauz_3d_rect_cd4 * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ixx_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixx_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixx_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_ixx_3d_rect = evua_ixx_3d_rect + ( ( Y ( k ) - yO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2 

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ixx_3d_rect = evua_ixx_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_iyy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyy_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyy_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_iyy_3d_rect = evua_iyy_3d_rect + ( ( X ( j ) - xO )**2 + ( Z ( l ) - zO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_iyy_3d_rect = evua_iyy_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_izz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_izz_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_izz_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_izz_3d_rect = evua_izz_3d_rect + ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_izz_3d_rect = evua_izz_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ixy_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , dX , dY , dZ , X , Y , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixy_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixy_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_ixy_3d_rect = evua_ixy_3d_rect + ( X ( j ) - xO ) * ( Y ( k ) - yO ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ixy_3d_rect = -evua_ixy_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_iyz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , yO , zO , dX , dY , dZ , Y , Z , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyz_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_iyz_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_iyz_3d_rect = evua_iyz_3d_rect + ( Y ( k ) - yO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_iyz_3d_rect = -evua_iyz_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      REAL FUNCTION evua_ixz_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , zO , dX , dY , dZ , X , Z , &
         & Psi3 )

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixz_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_ixz_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_ixz_3d_rect = evua_ixz_3d_rect + ( X ( j ) - xO ) * ( Z ( l ) - zO ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

        END DO
!$OMP   END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_ixz_3d_rect = -evua_ixz_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vex_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vex_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_vex_3d_rect = evua_vex_3d_rect + Vex3 ( j , k , l ) * ABS ( Psi3 ( j , k , l ) )**2

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_vex_3d_rect = evua_vex_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

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

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vmf_3d_rect )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC ) REDUCTION ( + : evua_vmf_3d_rect )
         DO k = nYa , nYb

            DO j = nXa , nXb

               evua_vmf_3d_rect = evua_vmf_3d_rect + ABS ( Psi3 ( j , k , l ) )**4

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      evua_vmf_3d_rect = 0.5 * gS * evua_vmf_3d_rect * dX * dY * dZ

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------------------------------------------------------------------

      END MODULE

! ==================================================================================================================================
