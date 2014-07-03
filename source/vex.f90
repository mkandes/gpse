! ==========================================================================
! NAME
!
!     vex [ veks ] - VEX Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     VEX is a Fortran module ... 
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
!     Wednesday, July 2nd, 2014
!
! -------------------------------------------------------------------------

      MODULE VEX

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MPI
      USE            :: GRID, ONLY: nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc

      IMPLICIT NONE
      PRIVATE

      INTEGER, PARAMETER, PRIVATE :: vexUnitIn = 503
      INTEGER, PARAMETER, PRIVATE :: vexUnitInit = 998

      LOGICAL, PUBLIC :: vexRead   = .FALSE.
      LOGICAL, PUBLIC :: vexWrite  = .FALSE.

      INTEGER, PRIVATE :: vexFmtIn  = -1      ! 0 = Unformatted ( Binary ) ; 1 = Formatted ( GPI ) 
      INTEGER, PRIVATE :: vexFmtOut = -1      ! 0 = Unformatted ( Binary ) ; 1 = Formatted ( GPI ) ; 2 = VTK ; 4 = VTK_XML
      INTEGER, PRIVATE :: vexInit = -1 ! 0 = Linear ; 1 = Simple Harmonic Oscillator ; 2 = Simple Harmonic Oscillator Ring

      REAL, PRIVATE :: xO = 0.0
      REAL, PRIVATE :: yO = 0.0
      REAL, PRIVATE :: zO = 0.0
      REAL, PRIVATE :: rO = 0.0
      REAL, PRIVATE :: fX = 0.0  
      REAL, PRIVATE :: fY = 0.0  
      REAL, PRIVATE :: fZ = 0.0  
      REAL, PRIVATE :: wX = 0.0  
      REAL, PRIVATE :: wY = 0.0  
      REAL, PRIVATE :: wZ = 0.0
      REAL, PRIVATE :: wR = 0.0

      PUBLIC :: vex_read_inputs
      PUBLIC :: vex_write_inputs
      PUBLIC :: vex_mpi_bcast_inputs
      PUBLIC :: vex_read_init
      PUBLIC :: vex_compute_init

      PUBLIC :: vex_3d_lin
      PUBLIC :: vex_3d_sho
      PUBLIC :: vex_3d_shor

      NAMELIST /nmlVexIn/ vexRead , vexWrite , vexFmtIn , vexFmtOut , vexInit , xO , yO , zO , rO , fX , fY , fZ , wX , wY , wZ , wR 

      CONTAINS

         SUBROUTINE vex_read_inputs ( )

            IMPLICIT NONE

            OPEN ( UNIT = vexUnitIn, FILE = 'vex.in' , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
               READ ( UNIT = vexUnitIn , NML = nmlVexIn )
            CLOSE ( UNIT = vexUnitIn , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_write_inputs ( )

            IMPLICIT NONE

            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# vexRead   = ', vexRead
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# vexWrite  = ', vexWrite
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# vexFmtIn  = ', vexFmtIn
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# vexFmtOut = ', vexFmtOut
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# vexInit   = ', vexInit
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xO = ', xO
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yO = ', yO
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zO = ', zO
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        rO = ', rO
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fX = ', fX
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fY = ', fY
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fZ = ', fZ
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wX = ', wX
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wY = ', wY
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZ = ', wZ
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wR = ', wR

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_mpi_bcast_inputs ( mpiMaster , mpiInt , mpiReal , mpiCmplx , mpiError )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: mpiMaster
            INTEGER, INTENT ( IN    ) :: mpiInt
            INTEGER, INTENT ( IN    ) :: mpiReal
            INTEGER, INTENT ( IN    ) :: mpiCmplx
            INTEGER, INTENT ( INOUT ) :: mpiError

            CALL MPI_BCAST ( vexRead   , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexWrite  , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFmtIn  , 1 , mpiInt  , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFmtOut , 1 , mpiInt  , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexInit   , 1 , mpiInt  , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( xO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( yO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( zO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( rO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fX        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fY        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fZ        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wX        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wY        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wZ        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wR        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_read_init ( )

            IMPLICIT NONE

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_compute_init ( X , Y , Z , Vex3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            IF ( vexInit == 0 ) THEN 

               CALL vex_3d_lin ( X , Y , Z , Vex3 )

            ELSE IF ( vexInit == 1 ) THEN 

               CALL vex_3d_sho ( X , Y , Z , Vex3 )

            ELSE IF ( vexInit == 2 ) THEN 

               CALL vex_3d_shor ( X , Y , Z , Vex3 )

            ELSE 

               ! Error: vexInit not defined.

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_3d_lin ( X , Y , Z , Vex3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + fX * ( X ( j ) - xO ) + fY * ( Y ( k ) - yO ) + fZ * ( Z ( l ) - zO )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_3d_sho ( X , Y , Z , Vex3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( ( wX * ( X ( j ) - xO ) )**2 + ( wY * ( Y ( k ) - yO ) )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE vex_3d_shor ( X , Y , Z , Vex3 )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l 

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb 

               DO k = nYa , nYb 

                  DO j = nXa , nXb

                     Vex3 ( j , k , l ) = Vex3 ( j , k , l ) + 0.5 * ( wR * ( SQRT ( ( X ( j ) - xO )**2 + ( Y ( k ) - yO )**2 ) - rO )**2 + ( wZ * ( Z ( l ) - zO ) )**2 )

                  END DO
                 
               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

      END MODULE

! =========================================================================
