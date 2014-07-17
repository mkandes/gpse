! ==========================================================================
! NAME
!
!     psi [ (p)sī ] - Psi Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     psi is a Fortran module ... 
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
!     Wednesday, July 9th, 2014
!
! -------------------------------------------------------------------------

      MODULE PSI

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: GRID, ONLY: nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc
      USE            :: MATH

      IMPLICIT NONE
      PRIVATE

      INTEGER, PARAMETER, PRIVATE :: unitPsiIn   = 502
      INTEGER, PARAMETER, PRIVATE :: unitPsiInit = 999 

      LOGICAL, PUBLIC :: readPsi = .FALSE. ! Read initial wave function from file ( init.psi ) ? .TRUE. = Yes ; .FALSE. = No
      LOGICAL, PUBLIC :: writePsi    = .FALSE. ! Write wave function to file? .TRUE. = Yes ; .FALSE. = No

      INTEGER, PUBLIC :: psiInit = -1 ! 0 = Isotropic 3D SHO ; 1 = Anisotropic 3D SHO ; 2 = Axially-Symmetric 3D SHO
      INTEGER, PUBLIC :: psiNx   = 0  ! Degree of Hermite polynomial used to define anisotropic SHO wave function along x-axis
      INTEGER, PUBLIC :: psiNy   = 0  ! Degree of Hermite polynomial used to define anisotropic SHO wave function along y-axis
      INTEGER, PUBLIC :: psiNz   = 0  ! Degree of Hermite polynomial used to define both anisotropic and axially-symmetric SHO wave functions along z-axis
      INTEGER, PUBLIC :: psiNr   = 0  ! Degree of (associated) Laguerre polynomials used to define radial components of isotropic and axi-symmetric SHO 
      INTEGER, PUBLIC :: psiMl   = 0  ! Projection of orbital angular momentum along z-axis for axially-symmetric SHO wave function
      INTEGER, PUBLIC :: fmtWritePsi = -1 ! File format for output wave functions? 0 = Binary ; 1 = GPI ; 2 = VTK ; 3 = VTK_XML

      REAL, PUBLIC :: psiXo = 0.0 ! X-coordinate of origin used to define initial wave function
      REAL, PUBLIC :: psiYo = 0.0 ! Y-coordinate of origin used to define initial wave function
      REAL, PUBLIC :: psiZo = 0.0 ! Z-coordinate of origin used to define initial wave function
      REAL, PUBLIC :: psiWx = 0.0 ! Angular frequency of SHO potential along x-axis used to define anisotropic SHO wave function
      REAL, PUBLIC :: psiWy = 0.0 ! Angular frequency of SHO potential along y-axis used to define anisotropic SHO wave function
      REAL, PUBLIC :: psiWz = 0.0 ! Angular frequency of SHO potential along z-axis used to define both anisotropic and axially-symmetric SHO wave functions
      REAL, PUBLIC :: psiWr = 0.0 ! Angular frequency of isotropic (radially-symmetric) SHO potential used to define ...    

      PUBLIC :: psi_read_inputs
      PUBLIC :: psi_read_init
      PUBLIC :: psi_compute_init
      PUBLIC :: psi_normalize

      PRIVATE :: psi_3d_se_sho_ani
      PRIVATE :: psi_3d_se_sho_axi
      PRIVATE :: psi_3d_se_sho_iso

      NAMELIST /nmlPsiIn/ readPsi , psiInit , psiNx , psiNy , psiNz , psiNr , psiMl , psiXo , psiYo , psiZo , psiWx , psiWy , psiWz , psiWr

      CONTAINS

         SUBROUTINE psi_read_inputs ( )

            IMPLICIT NONE

            OPEN ( UNIT = unitPsiIn, FILE = 'psi.in' , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
               READ ( UNIT = unitPsiIn , NML = nmlPsiIn )
            CLOSE ( UNIT = unitPsiIn , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_read_init ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_compute_init ( X , Y , Z , Psi3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            IF ( psiInit == 0 ) THEN 

               ! Error: Isotropic SHO not supported yet.

            ELSE IF ( psiInit == 1 ) THEN 

               CALL psi_3d_se_sho_ani ( X , Y , Z , Psi3 )

            ELSE IF ( psiInit == 2 ) THEN 

               CALL psi_3d_se_sho_axi ( X , Y , Z , Psi3 )

            ELSE 

               ! Error: psiInit not defined.

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_3d_se_sho_ani ( X , Y , Z , Psi3 )
 
            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Psi3 ( j , k , l )  = CMPLX ( & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNx * factorial ( psiNx ) ) ) ) * SQRT ( SQRT ( psiWx / PI ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNy * factorial ( psiNy ) ) ) ) * SQRT ( SQRT ( psiWy / PI ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNz * factorial ( psiNz ) ) ) ) * SQRT ( SQRT ( psiWz / PI ) ) * & 
                        & hermite ( psiNx , SQRT ( psiWx ) * ( X ( j ) - psiXo ) ) * EXP ( -0.5 * psiWx * ( X ( j ) - psiXo )**2 ) * & 
                        & hermite ( psiNy , SQRT ( psiWy ) * ( Y ( k ) - psiYo ) ) * EXP ( -0.5 * psiWy * ( Y ( k ) - psiYo )**2 ) * & 
                        & hermite ( psiNz , SQRT ( psiWz ) * ( Z ( l ) - psiZo ) ) * EXP ( -0.5 * psiWz * ( Z ( l ) - psiZo )**2 ) , 0.0 )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_3d_se_sho_axi ( X , Y , Z , Psi3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb  

                     Psi3 ( j , k , l ) = CMPLX ( &
                        & SQRT ( ( psiWr**( ABS ( psiMl ) + 1 ) * REAL ( factorial ( psiNr ) ) ) / & 
                        & ( PI * REAL ( factorial ( psiNr + ABS ( psiMl ) ) ) ) ) * & 
                        & ( 1.0 / SQRT ( REAL ( 2**psiNz * factorial ( psiNz ) ) ) ) * SQRT ( SQRT ( psiWz / PI ) ) * & 
                        & SQRT ( ( X ( j ) - psiXo )**2 + ( Y ( k ) - psiYo )**2 )**ABS ( psiMl ) * & 
                        & alaguerre ( psiNr , ABS ( psiMl ) , psiWr * ( ( X ( j ) - psiXo )**2 + ( Y ( k ) - psiYo )**2 ) ) * & 
                        & EXP ( -0.5 * psiWr * ( ( X ( j ) - psiXo )**2 + ( Y ( k ) - psiYo )**2 ) ) * & 
                        & hermite ( psiNz , SQRT ( psiWz ) * ( Z ( l ) - psiZo ) ) * EXP ( -0.5 * psiWz * ( Z ( l ) - psiZo )**2 ) , 0.0 ) * & 
                        & EXP ( CMPLX ( 0.0 , REAL ( psiMl ) * ATAN2 ( Y ( k ) - psiYo , X ( j ) - psiXo ) ) ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_3d_se_sho_iso ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE psi_normalize ( ) 

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
