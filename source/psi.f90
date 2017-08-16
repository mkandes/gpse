! ==================================================================================================================================
! NAME
!
!     psi [ (p)sī ] - Psi Module
!
! SYNOPSIS
!
!     USE :: PSI
!
! DESCRIPTION  
!
!     PSI is a custom Fortran module written to compute analytic solutions of the Schrodinger and Gross-Pitaevskii equations, which
!        may be used as initial conditions for simulations.
!
! OPTIONS
!
! SEE ALSO
!
! BUGS
!
! HISTORY
!
! AUTHOR
!
!     Marty Kandes, Ph.D.
!     Computational & Data Science Research Specialist
!     User Services Group
!     San Diego Supercomputer Center
!     University of California, San Diego
!
! COPYRIGHT
!     
!     Copyright (c) 2014, 2015, 2016, 2017 Martin Charles Kandes
!
! LAST UPDATED
!
!     Wednesday, August 16th, 2017
!
! ----------------------------------------------------------------------------------------------------------------------------------

      MODULE PSI

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MATH

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
!     MATH is a custom Fortran module written to define well-know mathematical constants and compute specialized functions. The
!        module only has dependency on the ISO_FORTRAN_ENV module. 
!
! ----------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      PRIVATE

! --- VARIABLE DECLARATIONS --------------------------------------------------------------------------------------------------------

      INTEGER, PUBLIC :: psiInput  = -1
      INTEGER, PUBLIC :: psiOutput = -1
      INTEGER, PUBLIC :: psiFileNo = -1
      INTEGER ( KIND = 8 ), PUBLIC :: psiFilePos = -1
      INTEGER, PUBLIC :: psiInit   = -1

! --- VARIABLE DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     psiInput is a PUBLIC, INTEGER-valued variable that sets the file format of the initial wave function to be read in at the 
!        start of program execution. 0 = No input wave function; 1 = Read wave function from .bin file; 2 = Read wave function 
!        from .vtk file. NOTE: VTK READER NOT AVAILABLE YET.
!
!     psiOutput is a PUBLIC, INTEGER-valued variable that sets the file format of the wave functions to be written out to disk 
!        during program execution. 0 = No output wave function; 1 = Write wave function to .bin file; 2 = Write wave function to 
!        .vtk file.
!   
!     psiFileNo is a PUBLIC, INTEGER-valued variable that sets the unit number of the file to be read as the initial wave function 
!        file.
!
!     psiFilePos is a PUBLIC, INTEGER-valued variable that tracks the file position when reading and writing the input and output 
!        wave function files, respectively. NOTE: THIS POSITION VARIABLE IS CURRENTLY HARD-CODED WITH KIND = 8 TO AVOID INTEGER 
!        OVERFLOWS WHEN WRITING OUT LARGE FILES.
!
!     psiInit is a PUBLIC, INTEGER-valued variable that determines which analytical wave function to use as an initial condition to
!        the simulation when no input wave function file is used. 0 = Isotropic 3D SHO; 1 = Anisotropic 3D SHO; 2 = Axisymmetric 3D
!        SHO; 3 = Approx 3D SHOR.
! 
! --- SUBROUTINE DECLARATIONS ------------------------------------------------------------------------------------------------------

      PUBLIC :: psi_init
      PUBLIC :: psi_normalize
      PUBLIC :: psi_boost

      PRIVATE :: psi_3d_se_sho_ani
      PRIVATE :: psi_3d_se_sho_axi
      PRIVATE :: psi_3d_se_shor_axi
      PRIVATE :: psi_3d_se_sho_iso

! --- SUBROUTINE DEFINITIONS -------------------------------------------------------------------------------------------------------
!
!     psi_init is a PUBLIC SUBROUTINE that selects which analytic wave wave function will be used as an initial condition to the
!        simulation if no input wave function file is provided.
!
!     psi_normalize is NOT IMPLEMENTED YET.
!
!     psi_boost is a PUBLIC SUBROUTINE that applies a linear momentum boost to a three-dimensional wave function defined on a 
!        regular grid. It is currently used to provide an initial linear momentum to the initial wave function of the simulation.
!
!     psi_3d_se_sho_ani is a PRIVATE SUBROUTINE that computes the analytic solution of the time-independent Schrodinger equation 
!        for a three-dimensional, anisotropic simple harmonic oscillator potential.
!
!     psi_3d_se_sho_axi is a PRIVATE SUBROUTINE that computes the analytic solution of the time-independent Schrodinger equation 
!        for a three-dimensional, axisymmetric simple harmonic oscillator potential. 
!        
!     psi_3d_se_shor_axi is a PRIVATE SUBROUTINE that computes an analytic approximation to the solution of the time-dependent 
!        Schrodinger equation for a three-dimensional, simple harmonic oscillator ring potential.
!
!     psi_3d_se_sho_iso is NOT IMPLEMENTED YET.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_init ( psiInit , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nX , nY , nZ , nR , mL , xO , yO , &
         & zO , rO , wX , wY , wZ , wR , X , Y , Z , Psi3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: psiInit
      INTEGER, INTENT ( IN ) :: nXa 
      INTEGER, INTENT ( IN ) :: nXb 
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa 
      INTEGER, INTENT ( IN ) :: nYb 
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: nZa 
      INTEGER, INTENT ( IN ) :: nZb 
      INTEGER, INTENT ( IN ) :: nZbc
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nZ
      INTEGER, INTENT ( IN ) :: nR
      INTEGER, INTENT ( IN ) :: mL

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: rO
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: wR 

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

      IF ( psiInit == 0 ) THEN

         ! ?

      ELSE IF ( psiInit == 1 ) THEN 

         CALL psi_3d_se_sho_ani ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nX , nY , nZ , xO , yO , zO , wX , wY , &
            & wZ , X , Y , Z , Psi3 )

      ELSE IF ( psiInit == 2 ) THEN 

         CALL psi_3d_se_sho_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nR , mL , nZ , xO , yO , zO , wR , wZ , &
            & X , Y , Z , Psi3 )

      ELSE IF ( psiInit == 3 ) THEN

         CALL psi_3d_se_shor_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nR , mL , nZ , xO , yO , zO , rO , wR ,&
            & wZ , X , Y , Z , Psi3 )

      ELSE 

         WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : psi : psi_init : ERROR - psiInit not supported.'
         STOP

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_3d_se_sho_ani ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nX , nY , nZ , xO , yO , zO , wX , &
         & wY , wZ , X , Y , Z , Psi3 )
 
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
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nZ

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Psi3 ( j , k , l )  = CMPLX ( & 
                  & ( 1.0 / SQRT ( REAL ( 2**nX * factorial ( nX ) ) ) ) * SQRT ( SQRT ( wX / PI ) ) * & 
                  & ( 1.0 / SQRT ( REAL ( 2**nY * factorial ( nY ) ) ) ) * SQRT ( SQRT ( wY / PI ) ) * & 
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * & 
                  & hermite ( nX , SQRT ( wX ) * ( X ( j ) - xO ) ) * EXP ( -0.5 * wX * ( X ( j ) - xO )**2 ) * & 
                  & hermite ( nY , SQRT ( wY ) * ( Y ( k ) - yO ) ) * EXP ( -0.5 * wY * ( Y ( k ) - yO )**2 ) * & 
                  & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_3d_se_sho_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nR , mL , nZ , xO , yO , zO , wR , &
         & wZ , X , Y , Z , Psi3 )

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
      INTEGER, INTENT ( IN ) :: nR
      INTEGER, INTENT ( IN ) :: mL
      INTEGER, INTENT ( IN ) :: nZ

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: wR
      REAL, INTENT ( IN ) :: wZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb  

               Psi3 ( j , k , l ) = CMPLX ( &
                  & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / & 
                  & ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * & 
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * & 
                  & SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 )**ABS ( mL ) * & 
                  & alaguerre ( nR , ABS ( mL ) , wR * ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) ) * & 
                  & EXP ( -0.5 * wR * ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) ) * & 
                  & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 ) * & 
                  & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( Y ( k ) - yO , X ( j ) - xO ) ) ) 

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_3d_se_shor_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nR , mL , nZ , xO , yO , zO , rO , &
         & wR , wZ , X , Y , Z , Psi3 )

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
      INTEGER, INTENT ( IN ) :: nR
      INTEGER, INTENT ( IN ) :: mL
      INTEGER, INTENT ( IN ) :: nZ

      REAL, INTENT ( IN ) :: xO
      REAL, INTENT ( IN ) :: yO
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: rO
      REAL, INTENT ( IN ) :: wR
      REAL, INTENT ( IN ) :: wZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

      INTEGER :: j , k , l 

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Psi3 ( j , k , l ) = CMPLX ( &
                  & SQRT ( ( wR**( ABS ( mL ) + 1 ) * REAL ( factorial ( nR ) ) ) / &
                  & ( PI * REAL ( factorial ( nR + ABS ( mL ) ) ) ) ) * &
                  & ( 1.0 / SQRT ( REAL ( 2**nZ * factorial ( nZ ) ) ) ) * SQRT ( SQRT ( wZ / PI ) ) * &
                  & SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 )**ABS ( mL ) * &
                  & alaguerre ( nR , ABS ( mL ) , wR * ( SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) - rO )**2 ) * &
                  & EXP ( -0.5 * wR * ( SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) - rO )**2 ) * &
                  & hermite ( nZ , SQRT ( wZ ) * ( Z ( l ) - zO ) ) * EXP ( -0.5 * wZ * ( Z ( l ) - zO )**2 ) , 0.0 ) * &
                  & EXP ( CMPLX ( 0.0 , REAL ( mL ) * ATAN2 ( Y ( k ) - yO , X ( j ) - xO ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_3d_se_sho_iso ( )

      IMPLICIT NONE

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_normalize ( ) 

      IMPLICIT NONE

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE psi_boost ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , pX , pY , pZ , X , Y , Z , &
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
      REAL, INTENT ( IN ) :: zO
      REAL, INTENT ( IN ) :: pX
      REAL, INTENT ( IN ) :: pY
      REAL, INTENT ( IN ) :: pZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa - nZbc , nZb + nZbc ! if you don't include boundary points, they will not be properly boosted with interior points

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa - nYbc , nYb + nYbc

            DO j = nXa - nXbc , nXb + nXbc

               Psi3 ( j , k , l ) = Psi3 ( j , k , l ) * EXP ( CMPLX ( 0.0 , pX * ( X ( j ) - xO ) + pY * ( Y ( k ) - yO ) + &
                  & pZ * ( Z ( l ) - zO ) ) )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      END MODULE

! ==================================================================================================================================
