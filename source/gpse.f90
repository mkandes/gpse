! ==========================================================================
! NAME
!
!     gpse [ jip-see ] - 
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     gpse is a Fortran program ... 
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
!     Thursday, July 3rd, 2014
!
! -------------------------------------------------------------------------

      PROGRAM GPSE

! --- MODULE DECLARATIONS -------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MPI
      USE            :: EVUA
      USE            :: GRID
      USE            :: IO
      USE            :: MATH
      USE            :: PSI
      USE            :: VEX

! --- MODULE DEFINITIONS --------------------------------------------------
! -------------------------------------------------------------------------
      
      IMPLICIT NONE

! --- PARAMETER DECLARATIONS  ---------------------------------------------

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.1.8'
      CHARACTER ( LEN = * ), PARAMETER :: GPSE_LAST_UPDATED = 'Thursday, July 3rd, 2014'

      INTEGER, PARAMETER :: INT_DEFAULT_KIND   = KIND ( 0 ) 
      INTEGER, PARAMETER :: REAL_DEFAULT_KIND  = KIND ( 0.0 )
      INTEGER, PARAMETER :: CMPLX_DEFAULT_KIND = KIND ( CMPLX ( 0.0 , 0.0 ) )
      INTEGER, PARAMETER :: MAX_LEN_FILENAME   = 255
      INTEGER, PARAMETER :: MPI_MASTER         = 0

! --- PARAMETER DEFINITIONS -----------------------------------------------
! --- VARIABLE DECLARATIONS -----------------------------------------------

      LOGICAL :: itpOn = .FALSE.

      CHARACTER ( LEN = 8  ) :: startDate = 'NONE'
      CHARACTER ( LEN = 10 ) :: startTime = 'NONE'
      CHARACTER ( LEN = 5  ) :: startZone = 'NONE'
      CHARACTER ( LEN = 8  ) :: stopDate  = 'NONE'
      CHARACTER ( LEN = 10 ) :: stopTime  = 'NONE'
      CHARACTER ( LEN = 5  ) :: stopZone  = 'NONE'

      INTEGER :: rk4Lambda     = 0
      INTEGER :: fdOrder       = 0
      INTEGER :: quadRule      = 0
      INTEGER :: mpiCmplx      = 0 
      INTEGER :: mpiError      = 0 
      INTEGER :: mpiErrorCode  = 0 
      INTEGER :: mpiErrorClass = 0 
      INTEGER :: mpiInt        = 0 
      INTEGER :: mpiProcesses  = 0 
      INTEGER :: mpiProvided   = 0 
      INTEGER :: mpiRank       = 0
      INTEGER :: mpiReal       = 0 
      INTEGER :: ompThreads    = 0 
      INTEGER :: ompThreadID   = 0
      INTEGER :: nTsteps       = 0
      INTEGER :: nTwrite       = 0
      INTEGER :: j , k , l , n , m ! Reserved loop counters. 

      REAL :: tN = 0.0 ! Simulation time at nth time step
      REAL :: t0 = 0.0 ! Time at the beginning of the simulation
      REAL :: dT = 0.0 ! Interval of a time step
      REAL :: wX = 0.0
      REAL :: wY = 0.0
      REAL :: wZ = 0.0
      REAL :: gS = 0.0 

      COMPLEX :: dTz = CMPLX ( 0.0 , 0.0 )

! --- VARIABLE DEFINITIONS ------------------------------------------------
! --- ARRAY DECLARATIONS --------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIStatus
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartVals
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopVals

      REAL, ALLOCATABLE, DIMENSION ( : , : , : ) :: Vex3
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: X
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Y
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Z

      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K4
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3a
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3b
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3

! --- ARRAY DEFINITIONS ---------------------------------------------------
! --- FUNCTION AND SUBROUTINE DECLARATIONS --------------------------------

      INTEGER :: OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM

! --- NAMELIST DECLARATIONS -----------------------------------------------

      NAMELIST /gpseIn/ itpOn , rk4Lambda , fdOrder , quadRule , nTsteps , nTwrite , t0 , dT , wX , wY , wZ , gS

! --- NAMELIST DEFINITIONS ------------------------------------------------

! --- FUNCTION AND SUBROUTINE DEFINITIONS ---------------------------------
! --- MAIN PROGRAM --------------------------------------------------------

      ALLOCATE ( MPIStatus ( MPI_STATUS_SIZE ) ) 
      CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiProvided , mpiError )
      IF ( mpiError /= MPI_SUCCESS ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'gpse: ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
         CALL MPI_ABORT ( MPI_COMM_WORLD , mpiErrorCode , mpiError )

      END IF
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiProcesses , mpiError )
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRank , mpiError )
      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!$OMP PARALLEL DEFAULT ( SHARED )
      ompThreads = OMP_GET_NUM_THREADS ( )
!$OMP END PARALLEL

      ALLOCATE ( StartVals ( 8 ) )
      CALL DATE_AND_TIME ( startDate , startTime , startZone , StartVals )
      
      CALL select_mpi_default_kinds ( )

      IF ( mpiRank == MPI_MASTER ) THEN

         CALL read_inputs ( )
         CALL write_header ( )

         CALL grid_read_inputs ( )
         CALL grid_write_inputs ( )
         CALL grid_bound_cond_size ( fdOrder )

         CALL vex_read_inputs ( )
         CALL vex_write_inputs ( )

         CALL psi_read_inputs ( )
         CALL psi_write_inputs ( )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RANGE CHECKING INPUT PARAMETERS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     BROADCASTING INPUT PARAMETERS TO ALL MPI PROCESSES ... '

      END IF

      CALL gpse_mpi_bcast_inputs ( MPI_MASTER , mpiInt , mpiReal , mpiCmplx , mpiError )
      CALL grid_mpi_bcast_inputs ( MPI_MASTER , mpiInt , mpiReal , mpiCmplx , mpiError )
      CALL vex_mpi_bcast_inputs  ( MPI_MASTER , mpiInt , mpiReal , mpiCmplx , mpiError )
      CALL psi_mpi_bcast_inputs  ( MPI_MASTER , mpiInt , mpiReal , mpiCmplx , mpiError )

      IF ( itpOn .EQV. .TRUE. ) THEN 

         dTz = CMPLX ( 0.0 , -dT )

      ELSE

         dTz = CMPLX ( dT , 0.0 )

      END IF

      nXa = 1
      nXb = nX
      nYa = 1
      nYb = nY
      nZa = 1 + mpiRank * FLOOR ( REAL ( nZ / mpiProcesses ) )

      IF ( ( mpiRank + 1 ) == mpiProcesses ) THEN ! include any remainder z-points on last MPI process ...

         nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcesses ) ) + MODULO ( nZ , mpiProcesses )

      ELSE ! all other MPI processes have same number of z-points ... 

         nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcesses ) )

      END IF

      ALLOCATE ( X    ( nXa - nXbc : nXb + nXbc                                                     ) )
      ALLOCATE ( Y    ( nYa - nYbc : nYb + nYbc                                                     ) )
      ALLOCATE ( Z    ( nZa - nZbc : nZb + nZbc                                                     ) )
      ALLOCATE ( Vex3 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

      ALLOCATE ( K1    ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K2    ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K3    ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K4    ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Psi3a ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Psi3b ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

!      IF ( mpiRank == MPI_MASTER ) THEN
!
!         ALLOCATE ( Psi3 ( -nXbc : nX + nXbc , -nYbc : nY + nYbc , -nZbc : nZ + nZbc ) ) 
!
!      END IF

      Vex3 = 0.0
      X   = 0.0
      Y   = 0.0
      Z   = 0.0

      K1   = CMPLX ( 0.0 , 0.0 )
      K2   = CMPLX ( 0.0 , 0.0 )
      K3   = CMPLX ( 0.0 , 0.0 )
      K4   = CMPLX ( 0.0 , 0.0 )
      Psi3a = CMPLX ( 0.0 , 0.0 )
      Psi3b = CMPLX ( 0.0 , 0.0 )
      
!      IF ( mpiRank == MPI_MASTER ) THEN 
! 
!         Psi3 = CMPLX ( 0.0 , 0.0 )
!
!      END IF

      CALL grid_regular ( X , Y , Z )

      IF ( vexRead .EQV. .TRUE. ) THEN ! initialize external potential from file
     
         CALL vex_read_init ( )

      ELSE ! initialize external potential using built-in module functions and/or subroutines

         CALL vex_compute_init ( X , Y , Z , Vex3 )

      END IF

      IF ( readPsiInit .EQV. .TRUE. ) THEN ! read initial wave function from file

         CALL psi_read_init ( ) 

      ELSE ! initialize wave function using built-in module functions and/or subroutines

         CALL psi_compute_init ( X , Y , Z , Psi3a )

      END IF

      CALL mpi_exchange_ghosts ( Psi3a )

! --- BEGIN MAIN TIME PROPAGATION LOOP ------------------------------------

      tN = t0 ! initialize simulation time

      DO n = 0 , nTsteps

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

         IF ( MODULO ( n , nTwrite ) == 0 ) THEN

!           Compute partial expectation values locally on each MPI process

            IF ( quadRule == 1 ) THEN ! use rectangle rule
            
               normL2Local = l2_norm_3d_rect ( Psi3a )
               avgXLocal   = x_3d_rect ( X , Psi3a )
               avgX2Local  = x2_3d_rect ( X , Psi3a ) 
               avgYLocal   = y_3d_rect ( Y , Psi3a )
               avgY2Local  = y2_3d_rect ( Y , Psi3a )
               avgZLocal   = z_3d_rect (  Z , Psi3a )
               avgZ2Local  = z2_3d_rect ( Z , Psi3a )
               avgVexLocal = vex_3d_rect ( Vex3 , Psi3a )
               avgVmfLocal = vmf_3d_rect ( gS , Psi3a )

               IF ( fdOrder == 2 ) THEN ! use 2nd-order central differences

                  avgPxLocal  = px_3d_rect_cd2 ( Psi3a )
                  avgPx2Local = px2_3d_rect_cd2 ( Psi3a )
                  avgPyLocal  = py_3d_rect_cd2 ( Psi3a )
                  avgPy2Local = py2_3d_rect_cd2 ( Psi3a )
                  avgPzLocal  = pz_3d_rect_cd2 ( Psi3a )
                  avgPz2Local = pz2_3d_rect_cd2 ( Psi3a )
                  avgLxLocal  = lx_3d_rect_cd2 ( Y , Z , Psi3a )
                  avgLx2Local = lx2_3d_rect_cd2 ( Y , Z , Psi3a )
                  avgLyLocal  = ly_3d_rect_cd2 ( X , Z , Psi3a )
                  avgLy2Local = ly2_3d_rect_cd2 ( X , Z , Psi3a )
                  avgLzLocal  = lz_3d_rect_cd2 ( X , Y , Psi3a )
                  avgLz2Local = lz2_3d_rect_cd2 ( X , Y , Psi3a )

               ELSE IF ( fdOrder == 4 ) THEN ! use 4th-order central differences

                  avgPxLocal  = px_3d_rect_cd4 ( Psi3a )
                  avgPx2Local = px2_3d_rect_cd4 ( Psi3a )
                  avgPyLocal  = py_3d_rect_cd4 ( Psi3a )
                  avgPy2Local = py2_3d_rect_cd4 ( Psi3a )
                  avgPzLocal  = pz_3d_rect_cd4 ( Psi3a )
                  avgPz2Local = pz2_3d_rect_cd4 ( Psi3a )
                  avgLxLocal  = lx_3d_rect_cd4 ( Y , Z , Psi3a )
                  !avgLx2Local = lx2_3d_rect_cd4 ( Y , Z , Psi3a )
                  avgLyLocal  = ly_3d_rect_cd4 ( X , Z , Psi3a )
                  !avgLy2Local = ly2_3d_rect_cd4 ( X , Z , Psi3a )
                  avgLzLocal  = lz_3d_rect_cd4 ( X , Y , Psi3a )
                  !avgLz2Local = lz2_3d_rect_cd4 ( X , Y , Psi3a )

               ELSE ! fdOrder not supported

               END IF

            ELSE ! quadRule not supported

            END IF

!           Reduce sum of partial expectation values from all MPI processes on MPI_MASTER

            CALL MPI_REDUCE ( normL2Local , normL2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgXLocal   , avgX   , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgX2Local  , avgX2  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgYLocal   , avgY   , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgY2Local  , avgY2  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgZLocal   , avgZ   , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgZ2Local  , avgZ2  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgVexLocal , avgVex , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgVmfLocal , avgVmf , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPxLocal  , avgPx  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPx2Local , avgPx2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPyLocal  , avgPy  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPy2Local , avgPy2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPzLocal  , avgPz  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPz2Local , avgPz2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLxLocal  , avgLx  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLx2Local , avgLx2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLyLocal  , avgLy  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLy2Local , avgLy2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLzLocal  , avgLz  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLz2Local , avgLz2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )

            IF ( mpiRank == MPI_MASTER ) THEN

!              Compute uncertainty in position, momentum and angular momentum

               sigX  = SQRT ( avgX2 - avgX**2 )
               sigY  = SQRT ( avgY2 - avgY**2 )
               sigZ  = SQRT ( avgZ2 - avgZ**2 )
               sigPx = SQRT ( avgPx2 - avgPx**2 )
               sigPy = SQRT ( avgPy2 - avgPy**2 )
               sigPz = SQRT ( avgPz2 - avgPz**2 )
               sigLx = SQRT ( avgLx2 - avgLx**2 )
               sigLy = SQRT ( avgLy2 - avgLy**2 )
               sigLz = SQRT ( avgLz2 - avgLz**2 )

!              Compute expectation value of the square of the total angular momentum

               avgL2 = avgLx2 + avgLy2 + avgLz2

!              Compute expectation values of kinetic energy

               avgTx = 0.5 * avgPx2
               avgTy = 0.5 * avgPy2
               avgTz = 0.5 * avgPz2

!              Compute expectation value of the total energy

               avgE = avgTx + avgTy + avgTz + avgVex + avgVmf

!              Write expectation values, uncertainties and uncertainty relations to file from MPI_MASTER

               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) &
                  & tN , normL2 , avgX  , avgX2  , sigX  ,                   &
                  &               avgY  , avgY2  , sigY  ,                   &
                  &               avgZ  , avgZ2  , sigZ  ,                   & 
                  &               avgPx , avgPx2 , sigPx ,                   & 
                  &               avgPy , avgPy2 , sigPy ,                   & 
                  &               avgPz , avgPz2 , sigPz ,                   & 
                  &               avgLx , avgLx2 , sigLx ,                   &
                  &               avgLy , avgLy2 , sigLy ,                   & 
                  &               avgLz , avgLz2 , sigLz ,                   &
                  & avgL2 ,                                                  &
                  &               avgTx , avgTy  , avgTz , avgVex , avgVmf , &
                  & avgE ,                                                   &
                  &         sigX  * sigPx , sigY  * sigPy , sigZ  * sigPz ,  &
                  &         sigLx * sigLy , sigLy * sigLz , sigLz * sigLx

            END IF

            ! Write wave function to file from MPI_MASTER

         END IF

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!        Compute 1st stage of generalized 4th-order Runge-Kutta (GRK4): k_1 = f ( t_n , y_n )

         CALL compute_f ( X , Y , Z , Vex3 , Psi3a , K1 )

!        Compute the intermediate wave function for the 2nd stage of the GRK4 method: y_n + 0.5 * dT * k_1
!        Note that the intermediate wave function is only computed on the interior grid points assigned to each MPI process.

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + CMPLX ( 0.5 , 0.0 ) * dTz * K1 ( j , k , l )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

!        Exchange boundary condition information among nearest neighbor processes

         CALL mpi_exchange_ghosts ( Psi3b )

!        Compute V ( x , t_n + dT / 2 ), W ( t_n + dT / 2)

         tN = t0 + ( REAL ( n ) + 0.5 ) * dT

!        Compute 2nd stage of GRK4: k_2 = f ( t_n + 0.5 * dT , y_n + 0.5 * dT * k_1 )

         CALL compute_f ( X , Y , Z , Vex3 , Psi3b , K2 )

!        Compute intermediate wave function for 3rd stage of GRK4: y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( lambda / 2 ) * dT * k_2

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + dTz * ( CMPLX ( 0.5 - 1.0 / REAL( rk4Lambda ) , 0.0 ) * K1 ( j , k , l ) + CMPLX ( 1.0 / REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

         CALL mpi_exchange_ghosts ( Psi3b )

!        Compute stage three ... k_3 = f ( t_n + 0.5 * dT , y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2

         CALL compute_f ( X , Y , Z , Vex3 , Psi3b , K3 )

!        Compute intermediate wave function for 4th stage of GRK4: y_n + ( 1 - lambda / 2 ) * dT * k_2 + ( lambda / 2 ) * dT * k_3

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + dTz * ( CMPLX ( 1.0 - 0.5 * REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( 0.5 * REAL ( rk4Lambda ) , 0.0 ) * K3 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

         CALL mpi_exchange_ghosts ( Psi3b )

!        Calculate V ( x , t_n + dT ), W ( t_n + dT )

         tN = t0 + REAL ( n + 1 ) * dT

!        Compute the fourth stage ... k_4 = f ( t_n + dT , y_n + ( 1 - lamda / 2 ) * dT * k_2 + ( lamda / 2) * dT * k_3 

         CALL compute_f ( X , Y , Z , Vex3 , Psi3b , K4 )

!        Compute wave function at nth+1 time step ... y_{ n + 1 } = y_n + ( dT / 6 ) * [ k_1 + ( 4 - lambda ) * k_2 + lambda * k_3 + k_4 ]

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  Psi3b ( j , k , l ) = Psi3a ( j , k , l ) + CMPLX ( 1.0 / 6.0 , 0.0 ) * dTz * ( K1 ( j , k , l ) + CMPLX ( 4.0 - REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( REAL ( rk4Lambda ) , 0.0 ) * K3 ( j , k , l ) + K4 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

         CALL mpi_exchange_ghosts ( Psi3b )

!        y_{ n + 1 } ---> y_n        
         Psi3a = Psi3b

      END DO

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

! --- END MAIN TIME PROPAGATION LOOP / BEGIN CLEAN UP TO STOP -------------

      DEALLOCATE ( Vex3 )
      DEALLOCATE ( Z )
      DEALLOCATE ( Y )
      DEALLOCATE ( X ) 

      IF ( mpiRank == MPI_MASTER ) THEN

         ALLOCATE ( StopVals ( 8 ) ) 
         CALL DATE_AND_TIME ( stopDate , stopTime , stopZone , StopVals )
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RUN STOPPED @ ', stopTime, ' ON ', stopDate, ' ... '
         DEALLOCATE ( StopVals )
         DEALLOCATE ( StartVals )

      END IF
      
      CALL MPI_FINALIZE ( mpiError )

      DEALLOCATE ( MPIStatus )

! --- FORMAT STATEMENTS ---------------------------------------------------

      STOP

      CONTAINS

         SUBROUTINE select_mpi_default_kinds ( )

            IMPLICIT NONE

            SELECT CASE ( INT_DEFAULT_KIND )

               CASE ( INT8 )

                  mpiInt = MPI_INTEGER1

               CASE ( INT16 )

                  mpiInt = MPI_INTEGER2

               CASE ( INT32 )

                  mpiInt = MPI_INTEGER

               CASE ( INT64 )

                  mpiInt = MPI_INTEGER8

               CASE DEFAULT

                  mpiInt = -1

            END SELECT

            SELECT CASE ( REAL_DEFAULT_KIND )

               CASE ( REAL32 )

                  mpiReal = MPI_REAL

               CASE ( REAL64 )

                  mpiReal= MPI_DOUBLE_PRECISION

               CASE ( REAL128 )

                  mpiReal = MPI_REAL16

               CASE DEFAULT

                  mpiReal = -1

            END SELECT

            SELECT CASE ( CMPLX_DEFAULT_KIND )

               CASE ( REAL32 )

                  mpiCmplx = MPI_COMPLEX

               CASE ( REAL64 )

                  mpiCmplx = MPI_DOUBLE_COMPLEX

               CASE ( REAL128 )

                  mpiCmplx = MPI_COMPLEX32

               CASE DEFAULT

                  mpiCmplx = -1

            END SELECT

            RETURN

         END SUBROUTINE

         SUBROUTINE read_inputs ( )

            IMPLICIT NONE

            OPEN ( UNIT = 500 , FILE = 'gpse.in' , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
               READ ( UNIT = 500 , NML = gpseIn )
            CLOSE ( UNIT = 500 , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE write_header ( )

            IMPLICIT NONE

            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# =========================================================================='
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     GPSE VERSION ', GPSE_VERSION_NUMBER
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        Compiled by ', COMPILER_VERSION ( ) , ' using the options ', COMPILER_OPTIONS ( )
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     AUTHOR(S)'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         Marty Kandes'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         Computational Science Research Center'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         San Diego State University'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         5500 Campanile Drive'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         San Diego, California 92182'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     COPYRIGHT'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'     
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         Copyright (c) 2014 Martin Charles Kandes'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     LAST UPDATED'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         ', GPSE_LAST_UPDATED
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# -------------------------------------------------------------------------'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     STARTING GPSE ... '
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RUN STARTED @ ', startTime, ' ON ', startDate, ' ... '
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RUNNING ON ', mpiProcesses, ' MPI PROCESSES WITH ', ompThreads , ' OPENMP THREADS PER PROCESS ... '
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     CHECKING MACHINE/COMPILER-SPECIFIC DATA TYPE SUPPORT AND CONFIGURATION ... '
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        INTEGER_KINDS SUPPORTED ... ', INTEGER_KINDS
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        REAL_KINDS SUPPORTED ...    ', REAL_KINDS
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        DEFAULT INTEGER KIND ...    ', INT_DEFAULT_KIND
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        DEFAULT REAL KIND ...       ', REAL_DEFAULT_KIND
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        DEFAULT COMPLEX KIND ...    ', CMPLX_DEFAULT_KIND
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     INPUT PARAMETERS ... '
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        itpOn   = ', itpOn
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        rk4Lambda  = ', rk4Lambda
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fdOrder    = ', fdOrder
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        quadRule   = ', quadRule
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nTsteps    = ', nTsteps
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nTwrite    = ', nTwrite
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        t0         = ', t0
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dT         = ', dT
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wX         = ', wX
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wY         = ', wY
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZ         = ', wZ
            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        gS         = ', gS

            RETURN

         END SUBROUTINE

         SUBROUTINE gpse_mpi_bcast_inputs ( mpiMaster , mpiInt , mpiReal , mpiCmplx , mpiError )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: mpiMaster
            INTEGER, INTENT ( IN ) :: mpiInt
            INTEGER, INTENT ( IN ) :: mpiReal
            INTEGER, INTENT ( IN ) :: mpiCmplx
            INTEGER, INTENT ( IN ) :: mpiError

            CALL MPI_BCAST ( itpOn      , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( rk4Lambda  , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fdOrder    , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( quadRule   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nTsteps    , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nTwrite    , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( t0         , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dT         , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wX         , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wY         , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wZ         , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( gS         , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            RETURN

         END SUBROUTINE

         SUBROUTINE compute_f ( X , Y , Z , Vex3 , Psi3 , F )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            IF ( fdOrder == 2 ) THEN 

               CALL f_gp_3d_rrf_cd2 ( X , Y , Z , Vex3 , Psi3 , F )

            ELSE IF ( fdOrder == 4 ) THEN

               CALL f_gp_3d_rrf_cd4 ( X , Y , Z , Vex3 , Psi3 , F )

            ELSE IF ( fdOrder == 6 ) THEN

               CALL f_gp_3d_rrf_cd6 ( X , Y , Z , Vex3 , Psi3 , F )

            ELSE IF ( fdOrder == 8 ) THEN

               CALL f_gp_3d_rrf_cd8 ( X , Y , Z , Vex3 , Psi3 , F )

            ELSE

               ! fdOrder not supported ...

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd2 ( X , Y , Z , Vex3 , Psi3 , F ) 

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb 

               DO k = nYa , nYb 

                  DO j = nXa , nXb 

                     F ( j , k , l ) = &
                        & CMPLX ( 0.5 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 0.5 / dZ**2 ) * Psi3 ( j , k , l - 1 ) + & 
                        & CMPLX ( 0.5 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 0.5 / dY**2 ) * Psi3 ( j , k - 1 , l ) + & 
                        & CMPLX ( 0.5 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 0.5 / dX**2 ) * Psi3 ( j - 1 , k , l ) - & 
                        & CMPLX ( 0.0 , 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 + Vex3 ( j , k , l ) + & 
                        &    gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + & 
                        & CMPLX ( 0.5 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 0.5 / dX**2 ) * Psi3 ( j + 1 , k , l ) + & 
                        & CMPLX ( 0.5 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 0.5 / dY**2 ) * Psi3 ( j , k + 1 , l ) + & 
                        & CMPLX ( 0.5 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 0.5 / dZ**2 ) * Psi3 ( j , k , l + 1 ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd4 ( X , Y , Z , Vex3 , Psi3 , F ) 

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb 

               DO k = nYa , nYb 

                  DO j = nXa , nXb 

                     F ( j , k , l ) = & 
                        & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi3 ( j , k , l - 2 ) + & 
                        & CMPLX ( 0.75 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi3 ( j , k , l - 1 ) + & 
                        & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi3 ( j , k - 2 , l ) + & 
                        & CMPLX ( 0.75 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi3 ( j , k - 1 , l ) + & 
                        & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi3 ( j - 2 , k , l ) + & 
                        & CMPLX ( 0.75 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi3 ( j - 1 , k , l ) - & 
                        & CMPLX ( 0.0 , 1.25 * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex3 ( j , k , l ) + & 
                        &    gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + & 
                        & CMPLX ( 0.75 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi3 ( j + 1 , k , l ) + & 
                        & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi3 ( j + 2 , k , l ) + & 
                        & CMPLX ( 0.75 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi3 ( j , k + 1 , l ) + & 
                        & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi3 ( j , k + 2 , l ) + & 
                        & CMPLX ( 0.75 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi3 ( j , k , l + 1 ) + & 
                        & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi3 ( j , k , l + 2 ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd6 ( X , Y , Z , Vex3 , Psi3 , F )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     F ( j , k , l ) = &
                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi3 ( j , k , l - 3 ) + &
                        & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * Psi3 ( j , k , l - 2 ) + &
                        & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * Psi3 ( j , k , l - 1 ) + &
                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi3 ( j , k - 3 , l ) + &
                        & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * Psi3 ( j , k - 2 , l ) + &
                        & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * Psi3 ( j , k - 1 , l ) + &
                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi3 ( j - 3 , k , l ) + &
                        & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * Psi3 ( j - 2 , k , l ) + &
                        & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * Psi3 ( j - 1 , k , l ) - &
                        & CMPLX ( 0.0 , ( 49.0 / 36.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex3 ( j , k , l ) + gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
                        & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * Psi3 ( j + 1 , k , l ) + &
                        & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * Psi3 ( j + 2 , k , l ) + &
                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi3 ( j + 3 , k , l ) + &
                        & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * Psi3 ( j , k + 1 , l ) + &
                        & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * Psi3 ( j , k + 2 , l ) + &
                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi3 ( j , k + 3 , l ) + &
                        & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * Psi3 ( j , k , l + 1 ) + &
                        & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * Psi3 ( j , k , l + 2 ) + &
                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi3 ( j , k , l + 3 )   

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd8 ( X , Y , Z , Vex3 , Psi3 , F )

            IMPLICIT NONE

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     F ( j , k , l ) = &
                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi3 ( j , k , l - 4 ) + &
                        & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 105.0 * dZ ) ,  4.0 / ( 315.0 * dZ**2 ) ) * Psi3 ( j , k , l - 3 ) + &
                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dZ ) , -1.0 / ( 10.0 * dZ**2 ) ) * Psi3 ( j , k , l - 2 ) + &
                        & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dZ ) ,  4.0 / ( 5.0 * dZ**2 ) ) * Psi3 ( j , k , l - 1 ) + &
                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 280.0 * dY ) , -1.0 / ( 1120 * dY**2 ) ) * Psi3 ( j , k - 4 , l ) + &
                        & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi3 ( j , k - 3 , l ) + &
                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi3 ( j , k - 2 , l ) + &
                        & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * Psi3 ( j , k - 1 , l ) + &
                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * Psi3 ( j - 4 , k , l ) + &
                        & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * Psi3 ( j - 3 , k , l ) + &
                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi3 ( j - 2 , k , l ) + &
                        & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi3 ( j - 1 , k , l ) - &
                        & CMPLX ( 0.0 , ( 205.0 / 144.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex3 ( j , k , l ) + gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
                        & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi3 ( j + 1 , k , l ) + &
                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi3 ( j + 2 , k , l ) + &
                        & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * Psi3 ( j + 3 , k , l ) + &
                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * Psi3 ( j + 4 , k , l ) + &
                        & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * Psi3 ( j , k + 1 , l ) + &
                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi3 ( j , k + 2 , l ) + & 
                        & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi3 ( j , k + 3 , l ) + &
                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 280.0 * dY ) , -1.0 / ( 1120.0 * dY**2 ) ) * Psi3 ( j , k + 4 , l ) + &
                        & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi3 ( j , k , l + 1 ) + &
                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi3 ( j , k , l + 2 ) + &
                        & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi3 ( j , k , l + 3 ) + &
                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi3 ( j , k , l + 4 ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_exchange_ghosts ( Psi3 ) 

            IMPLICIT NONE

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN 

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 0 , MPI_COMM_WORLD, mpiError )

               END IF

               IF ( mpiRank /= MPI_MASTER ) THEN 

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 1 , MPI_COMM_WORLD , MPIStatus , mpiError )

               END IF

            ELSE 

               CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 0 , MPI_COMM_WORLD , MPIStatus , mpiError )

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 1 , MPI_COMM_WORLD , mpiError )

               END IF

            END IF

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN

               IF ( mpiRank /= MPI_MASTER ) THEN

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 2 , MPI_COMM_WORLD , mpiError )

               END IF

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 3 , MPI_COMM_WORLD , MPIStatus , mpiError )

               END IF

            ELSE

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 2 , MPI_COMM_WORLD , MPIStatus , mpiError )

               END IF
               CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 3 , MPI_COMM_WORLD , mpiError )

            END IF

            RETURN

         END SUBROUTINE

      END PROGRAM

! =========================================================================
