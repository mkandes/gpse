! ==================================================================================================================================
! NAME
!
!     vex [ veks ] - External Potential Module
!
! SYNOPSIS
!
!     USE :: VEX
!
! DESCRIPTION  
!
!     VEX is a custom Fortran module written to compute analytic external potentials.
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
!     High-Performance Computing User Services Group
!     San Diego Supercomputer Center
!     University of California, San Diego
!
! COPYRIGHT
!     
!     Copyright (c) 2014, 2015, 2016, 2017, 2018, 2019, 2020 Martin Charles Kandes
!
! LAST UPDATED
!
!     Tuesday, April 2nd, 2019
!
! ----------------------------------------------------------------------------------------------------------------------------------

      MODULE VEX

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      PRIVATE

! --- VARIABLE DECLARATIONS --------------------------------------------------------------------------------------------------------

      INTEGER, PUBLIC :: vexInput  = -1
      INTEGER, PUBLIC :: vexOutput = -1
      INTEGER, PUBLIC :: vexFileNo = -1
      INTEGER, PUBLIC :: vexInit = -1

! --- VARIABLE DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     vexInput is a PUBLIC, INTEGER-valued input variable intended to designate the file format of the initial external potential to
!       be read in at the start of program execution. THIS HAS NOT BEEN IMPLEMENTED YET. 0 = No input external potential.
!
!     vexOutput is a PUBLIC, INTEGER-valued input variable intented to designate the file format of the external potentials to be 
!        written out to disk during program execution. THE HAS NOT BEEN IMPLEMENTED YET. 0 = No output external potential.
!       
!     vexFileNo is a PUBLIC, INTEGER-valued input variable intented to designate the unit number of the file to be read in as the 
!        initial external potential file. THIS HAS NOT BEEN IMPLEMENTED YET. 
!
!     vexInit is a PUBLIC, INTEGER-valued input variable that determines which time-independent external potential is applied 
!        throughout the course of a simulation when no input external potential is provided by the user. 0 =  Linear ; 1 = SHO ; 
!        2 = SHOR.
!
! --- SUBROUTINE DECLARATIONS ------------------------------------------------------------------------------------------------------

      PUBLIC :: vex_init

      PRIVATE :: vex_3d_lin
      PRIVATE :: vex_3d_sho
      PRIVATE :: vex_3d_shor

! --- SUBROUTINE DEFINITIONS -------------------------------------------------------------------------------------------------------
!
!     vex_init is a PUBLIC SUBROUTINE that selects which time-independent, external potential will be applied throughout the course
!        of a simulation if no input external poential is provided by the user at the beginning of program execution.
!
!     vex_3d_lin is a PRIVATE SUBROUTINE that computes the analytic form of a three-dimensional, linear external potential.
!
!     vex_3d_sho is a PRIVATE SUBROUTINE that computes the analytic form of a three-dimensional, simple harmonic oscillator 
!        potential.
!
!     vex_3d_shor is a PRIVATE SUBROUTINE that computes the analytic form of a three-dimensional, simple harmonic oscillator ring 
!        potential.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE vex_init ( vexInit , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , rO , fX , fY , fZ , &
         & wX , wY , wZ , wR , X , Y , Z , Vex3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: vexInit
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
      REAL, INTENT ( IN ) :: rO
      REAL, INTENT ( IN ) :: fX
      REAL, INTENT ( IN ) :: fY
      REAL, INTENT ( IN ) :: fZ
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ
      REAL, INTENT ( IN ) :: wR 

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

      IF ( vexInit == 0 ) THEN 

         CALL vex_3d_lin ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , fX , fY , fZ , X , Y , Z , Vex3 )

      ELSE IF ( vexInit == 1 ) THEN 

         CALL vex_3d_sho ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , wX , wY , wZ , X , Y , Z , Vex3 )

      ELSE IF ( vexInit == 2 ) THEN 

         CALL vex_3d_shor ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , rO , wR , wZ , X , Y , Z , &
            & Vex3 )

      ELSE 

         WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : vex : vex_init : ERROR - vexInit not supported.'

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE vex_3d_lin ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , fX , fY , fZ , X , Y , Z , &
         & Vex3 )

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
      REAL, INTENT ( IN ) :: fX
      REAL, INTENT ( IN ) :: fY
      REAL, INTENT ( IN ) :: fZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + fX * ( X ( j ) - xO ) + fY * ( Y ( k ) - yO ) + fZ * ( Z ( l ) - zO )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE vex_3d_sho ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , wX , wY , wZ , X , Y , Z , &
         & Vex3 )

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
      REAL, INTENT ( IN ) :: wX
      REAL, INTENT ( IN ) :: wY
      REAL, INTENT ( IN ) :: wZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

      INTEGER :: j , k , l

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb

            DO j = nXa , nXb

               Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( ( wX * ( X ( j ) - xO ) )**2 + ( wY * ( Y ( k ) - yO ) )**2 + &
                  & ( wZ * ( Z ( l ) - zO ) )**2 )

            END DO

         END DO
!$OMP    END PARALLEL DO

      END DO
!$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE vex_3d_shor ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xO , yO , zO , rO , wR , wZ , X , Y , Z , &
         & Vex3 )

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
      REAL, INTENT ( IN ) :: rO
      REAL, INTENT ( IN ) :: wR
      REAL, INTENT ( IN ) :: wZ

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

      INTEGER :: j , k , l 

!$OMP PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
      DO l = nZa , nZb 

!$OMP    PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
         DO k = nYa , nYb 

            DO j = nXa , nXb

               Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( wR * ( SQRT ( ( X ( j ) - xO )**2 + &
                  & ( Y ( k ) - yO )**2 ) - rO )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

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
