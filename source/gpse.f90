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
!     Friday, December 19th, 2014
!
! -------------------------------------------------------------------------

      PROGRAM GPSE

! --- MODULE DECLARATIONS -------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MPI
      USE            :: EVUA
      USE            :: GRID
      USE            :: GRK4
      USE            :: IO
      USE            :: MATH
      USE            :: PSI
      USE            :: VEX

! --- MODULE DEFINITIONS --------------------------------------------------
! -------------------------------------------------------------------------
      
      IMPLICIT NONE

! --- PARAMETER DECLARATIONS  ---------------------------------------------

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.4.2'
      CHARACTER ( LEN = * ), PARAMETER :: GPSE_LAST_UPDATED = 'Friday, December 19th, 2014'

      INTEGER, PARAMETER :: MPI_MASTER = 0

! --- PARAMETER DEFINITIONS -----------------------------------------------
! --- VARIABLE DECLARATIONS -----------------------------------------------

      LOGICAL :: itpOn   = .FALSE. ! Perform imaginary time propagation? .TRUE. = Yes ; .FALSE. = No

      CHARACTER ( LEN = 8  ) :: startDate = 'NONE'
      CHARACTER ( LEN = 10 ) :: startTime = 'NONE'
      CHARACTER ( LEN = 5  ) :: startZone = 'NONE'
      CHARACTER ( LEN = 8  ) :: stopDate  = 'NONE'
      CHARACTER ( LEN = 10 ) :: stopTime  = 'NONE'
      CHARACTER ( LEN = 5  ) :: stopZone  = 'NONE'

      INTEGER :: rk4Lambda     = -1 ! 1 = Tan-Chen Lambda-1 ; 2 = Classical 4th-Order Runge-Kutta ; 3 = Tan-Chen Lambda-3 ; 4 = England ; 5 = Tan-Chen Lambda-5
      INTEGER :: fdOrder       = -1 ! 2 = 2nd-Order Central Differences ( CD ); 4 = 4th-Order CD ; 6 = 6th-Order CD; 8 = 8th-Order CD
      INTEGER :: quadRule      = 1  ! 1 = Rectangle Rule ( only quadrature rule currently implemented ) 
      INTEGER :: nTsteps       = -1 ! Number of time steps in simulation
      INTEGER :: nTwrite       = -1 ! Period of writes to disk; i.e., number of time steps between writes to disk
      INTEGER :: nX            = -1 ! Number of grid points along the x-axis
      INTEGER :: nXa           = -1 !
      INTEGER :: nXb           = -1 !
      INTEGER :: nXbc          = -1 !
      INTEGER :: nY            = -1 ! Number of grid points along the y-axis
      INTEGER :: nYa           = -1 !
      INTEGER :: nYb           = -1 !
      INTEGER :: nYbc          = -1 !
      INTEGER :: nZ            = -1 ! Number of grid points along the z-axis
      INTEGER :: nZa           = -1 !
      INTEGER :: nZb           = -1 !
      INTEGER :: nZbc          = -1 !
      INTEGER :: nXpsi         = -1 ! Degree of Hermite polynomial used to define anisotropic SHO wave function along x-axis
      INTEGER :: nYpsi         = -1 ! Degree of Hermite polynomial used to define anisotropic SHO wave function along y-axis
      INTEGER :: nZpsi         = -1 ! Degree of Hermite polynomial used to define both anisotropic and axially-symmetric SHO wave functions along z-axis
      INTEGER :: nRpsi         = -1 ! Degree of Laguerre polynomials used to define radial components of isotropic and axially-symmetric SHO wave functions
      INTEGER :: mLpsi         = 0  ! Projection of orbital angular momentum along z-axis for axially-symmetric SHO wave function
      INTEGER :: mpiCmplx      = -1 
      INTEGER :: mpiError      = -1 
      INTEGER :: mpiErrorCode  = -1 
      INTEGER :: mpiErrorClass = -1
      INTEGER :: mpiInt        = -1 
      INTEGER :: mpiProcesses  = -1 
      INTEGER :: mpiProvided   = -1 
      INTEGER :: mpiRank       = -1
      INTEGER :: mpiSource     = -1
      INTEGER :: mpiDestination = -1
      INTEGER :: mpiReal       = -1 
      INTEGER :: mpiFile       = -1
      INTEGER :: ompThreads    = -1 ! Number of threads in OpenMP PARALLEL regions 
      INTEGER :: ompThreadID   = -1 
      INTEGER :: fileNumber    = -1
      INTEGER :: filePosX      = -1
      INTEGER :: filePosY      = -1
      INTEGER :: filePosZ      = -1
      INTEGER :: filePosRePsi  = -1
      INTEGER :: filePosImPsi  = -1
      INTEGER :: psiUnit = -1
      INTEGER :: vexUnit = -1
      INTEGER :: j , k , l , m , n  ! Reserved loop counters
      
      REAL :: tN    = 0.0 ! Time at nth time step
      REAL :: t0    = 0.0 ! Time at the beginning of the simulation
      REAL :: xO    = 0.0 !
      REAL :: yO    = 0.0 !
      REAL :: zO    = 0.0 !
      REAL :: dT    = 0.0 ! Interval of a time step
      REAL :: dX    = 0.0 !
      REAL :: dY    = 0.0 !
      REAL :: dZ    = 0.0 !
      REAL :: wX    = 0.0 ! X-component of the rotating reference frame's angular velocity vector
      REAL :: wY    = 0.0 ! Y-component of the rotating reference frame's angular velocity vector
      REAL :: wZ    = 0.0 ! Z-component of the rotating reference frame's angular velocity vector
      REAL :: gS    = 0.0 ! Nonlinear atom-atom interaction coupling constant
      REAL :: xOpsi = 0.0 !
      REAL :: yOpsi = 0.0 !
      REAL :: zOpsi = 0.0 !
      REAL :: wXpsi = 0.0 !
      REAL :: wYpsi = 0.0 !
      REAL :: wZpsi = 0.0 !
      REAL :: wRpsi = 0.0 !
      REAL :: pXpsi = 0.0 !
      REAL :: pYpsi = 0.0 !
      REAL :: pZpsi = 0.0 !
      REAL :: xOvex = 0.0 !
      REAL :: yOvex = 0.0 !
      REAL :: zOvex = 0.0 !
      REAL :: rOvex = 0.0 !
      REAL :: fXvex = 0.0 !
      REAL :: fYvex = 0.0 !
      REAL :: fZvex = 0.0 !
      REAL :: wXvex = 0.0 !
      REAL :: wYvex = 0.0 !
      REAL :: wZvex = 0.0 !
      REAL :: wRvex = 0.0 !

      COMPLEX :: dTz = CMPLX ( 0.0 , 0.0 ) ! Stores simulation time step in complex form; Used for imaginary time propagation

! --- VARIABLE DEFINITIONS ------------------------------------------------
! --- ARRAY DECLARATIONS --------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MpiStatus

      REAL, ALLOCATABLE, DIMENSION ( : ) :: Xa
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Xb
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Ya
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Yb
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Za
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Zb
      REAL, ALLOCATABLE, DIMENSION ( : , : , : ) :: Vex3a
      REAL, ALLOCATABLE, DIMENSION ( : , : , : ) :: Vex3b

      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K4
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3a
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3b

! --- ARRAY DEFINITIONS ---------------------------------------------------

! --- FUNCTION AND SUBROUTINE DECLARATIONS --------------------------------

      INTEGER :: OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM

! --- NAMELIST DECLARATIONS -----------------------------------------------

      NAMELIST /gpseIn/ itpOn , rk4Lambda , fdOrder , nTsteps , nTwrite , nX , nY , nZ , t0 , xO , yO , zO , dT , dX , dY , dZ , wX , wY , wZ , gS , psiInput , psiOutput , psiFileNo , psiInit , nXpsi , nYpsi , nZpsi , nRpsi , mLpsi , xOpsi , yOpsi , zOpsi , wXpsi , wYpsi , wZpsi , wRpsi , pXpsi , pYpsi , pZpsi , vexInput , vexOutput , vexFileNo , vexInit , xOvex , yOvex , zOvex , rOvex , fXvex , fYvex , fZvex , wXvex , wYvex , wZvex , wRvex

! --- NAMELIST DEFINITIONS ------------------------------------------------

! --- FUNCTION AND SUBROUTINE DEFINITIONS ---------------------------------

! --- MAIN PROGRAM --------------------------------------------------------

      ALLOCATE ( StartValues ( 8 ) )
      ALLOCATE ( StopValues  ( 8 ) )
      CALL DATE_AND_TIME ( startDate , startTime , startZone , StartValues )

      ALLOCATE ( MpiStatus ( MPI_STATUS_SIZE ) )
      CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiProvided , mpiError )
      IF ( mpiError /= MPI_SUCCESS ) THEN

         WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
         CALL MPI_ABORT ( MPI_COMM_WORLD , mpiErrorCode , mpiError )
         STOP

      END IF
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiProcesses , mpiError )
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRank , mpiError )

!$OMP PARALLEL DEFAULT ( SHARED )

      ompThreads = OMP_GET_NUM_THREADS ( )

!$OMP END PARALLEL

      IF ( mpiRank == MPI_MASTER ) THEN

         CALL gpse_read_stdin ( 'gpse.input' , 500 )
         CALL grid_boundary_condition_size ( fdOrder , nXbc , nYbc , nZbc )
         CALL gpse_write_stdout_header ( )

      END IF
      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      CALL mpi_select_default_kinds ( )
      CALL mpi_bcast_inputs ( MPI_MASTER , mpiInt , mpiReal , mpiCmplx , mpiError )

      IF ( itpOn .EQV. .TRUE. ) THEN ! run simulation using imaginary time propagation

         dTz = CMPLX ( 0.0 , -dT )

      ELSE ! run simulation normally

         dTz = CMPLX ( dT , 0.0 )

      END IF

      nXa = 1
      nXb = nX
      nYa = 1
      nYb = nY
      nZa = 1 + mpiRank * FLOOR ( REAL ( nZ / mpiProcesses ) )

      IF ( ( mpiRank + 1 ) == mpiProcesses ) THEN ! include any remaining z-axis grid points on last MPI process

         nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcesses ) ) + MODULO ( nZ , mpiProcesses )

      ELSE ! all other MPI processes have same number of z-axis grid points

         nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcesses ) )

      END IF

      ALLOCATE ( Xa ( nXa - nXbc : nXb + nXbc ) )
      ALLOCATE ( Ya ( nYa - nYbc : nYb + nYbc ) )
      ALLOCATE ( Za ( nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Xb ( nXa - nXbc : nXb + nXbc ) )
      ALLOCATE ( Yb ( nYa - nYbc : nYb + nYbc ) )
      ALLOCATE ( Zb ( nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Vex3a ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Vex3b ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K1 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K2 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K3 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K4 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Psi3a ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Psi3b ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

      Xa = 0.0
      Xb = 0.0
      Ya = 0.0
      Yb = 0.0
      Za = 0.0
      Zb = 0.0
      Vex3a = 0.0
      Vex3b = 0.0
      K1 = CMPLX ( 0.0 , 0.0 )
      K2 = CMPLX ( 0.0 , 0.0 )
      K3 = CMPLX ( 0.0 , 0.0 )
      K4 = CMPLX ( 0.0 , 0.0 )
      Psi3a = CMPLX ( 0.0 , 0.0 )
      Psi3b = CMPLX ( 0.0 , 0.0 )

      CALL grid_regular_axis ( nX , nXa , nXb , nXbc , xO , dX , Xa ) 
      CALL grid_regular_axis ( nY , nYa , nYb , nYbc , yO , dY , Ya ) 
      CALL grid_regular_axis ( nZ , nZa , nZb , nZbc , zO , dZ , Za ) 

      IF ( psiInput == 0 ) THEN ! compute initial wave function from available analytic expression

         CALL psi_init ( psiInit , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nXpsi , nYpsi , nZpsi , nRpsi , mLpsi , xOpsi , yOpsi , zOpsi , wXpsi , wYpsi , wZpsi , wRpsi , Xa , Ya , Za , Psi3a )

      ELSE IF ( psiInput == 1 ) THEN ! read initial wave function from binary file on MPI_MASTER

         psiFilePos = 1 ! initialize file position

         IF ( mpiRank == MPI_MASTER ) THEN ! read and initialize first block of wave function values from binary file on MPI_MASTER

            CALL io_read_psi3 ( psiFileNo , psiFilePos , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3a )

         END IF

         DO mpiDestination = 1 , mpiProcesses - 1 ! loop over mpiDesination ranks that are not equal to MPI_MASTER

            IF ( mpiRank == MPI_MASTER ) THEN ! read in next block of wave function values from binary file on MPI_MASTER

               CALL io_read_psi3 ( psiFileNo , psiFilePos , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )

            END IF
!           send read block from MPI_MASTER to mpiRank equal to mpiDestination
            CALL mpi_copy_psi3 ( mpiRank , MPI_MASTER , mpiDestination , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )

         END DO

         IF ( mpiRank /= MPI_MASTER ) THEN ! after all blocks have been sent, copy wave function values into Psi3a; MPI_MASTER block was read into Psi3a directly
         
            Psi3a = Psi3b

         END IF

      ELSE

         IF ( mpiRank == MPI_MASTER ) THEN

            WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : ERROR - psiInput not reognized.'
            STOP

         END IF

      END IF
      CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3a )
      CALL psi_boost ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOpsi , yOpsi , zOpsi , pXpsi , pYpsi , pZpsi , Xa , Ya , Za , Psi3a )

      IF ( vexInput == 0 ) THEN ! compute initial external potential from known analytic expression

         CALL vex_init ( vexInit , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOvex , yOvex , zOvex , rOvex , fXvex , fYvex , fZvex , wXvex , wYvex , wZvex , wRvex , Xa , Ya , Za , Vex3a )

      ELSE

         IF ( mpiRank == MPI_MASTER ) THEN

            WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : ERROR - vexInput not recognized.'
            STOP

         END IF 

      END IF

      psiFileNo = 1000
      vexFileNo = 1000

! --- BEGIN MAIN TIME PROPAGATION LOOP ------------------------------------

      tN = t0 ! initialize simulation time

      DO n = 0 , nTsteps

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

         IF ( MODULO ( n , nTwrite ) == 0 ) THEN ! compute partial 
!        relations to file from MPI_MASTER; write wave function and 
!        external potential to file from MPI_MASTER

            CALL evua_compute_base ( 1 , fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , 0.0 , 0.0 , 0.0 , dX , dY , dZ , gS , Xa , Ya , Za , Vex3a , Psi3a ) ! Compute partial base expectation values locally on each MPI process
            CALL evua_reduce_base ( MPI_MASTER , mpiReal , mpiError ) ! Reduce partial base expectation values from all MPI processes to MPI_MASTER to get full base expectation values
            CALL evua_compute_derived ( mpiRank , MPI_MASTER , wX , wY , wZ ) ! Compute derived expectation values, uncertainties from base expectation values 
            CALL evua_write_all ( mpiRank , MPI_MASTER , tN , wX , wY , wZ )  ! Write expectation values, undertainties and uncertainty relations to file from MPI_MASTER

            IF ( psiOutput == 1 ) THEN ! Write wave function to file from MPI_MASTER using streaming I/O binary with partial reduce

               Psi3b = Psi3a
               psiFilePos = 1
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_psi3 ( mpiRank , mpiSource , MPI_MASTER , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )
                  IF ( mpiRank == MPI_MASTER ) THEN

                     CALL io_write_psi3 ( psiFileNo , psiFilePos , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )

                  END IF

               END DO
               psiFileNo = psiFileNo + 1

            ELSE IF ( psiOutput == 2 ) THEN ! Write wave function to file from MPI_MASTER using streaming I/O vtk format with partial reduce

               Zb = Za
               Psi3b = Psi3a
               IF ( mpiRank == MPI_MASTER ) THEN

                  CALL io_write_vtk_header ( 'psi-', psiFileNo , psiFilePos , nX , nY , nZ )
                  CALL io_write_vtk_xcoordinates ( 'psi-' , psiFileNo , psiFilePos , nX , nXa , nXb , nXbc , Xa )
                  CALL io_write_vtk_ycoordinates ( 'psi-' , psiFileNo , psiFilePos , nY , nYa , nYb , nYbc , Ya )

               END IF
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_q ( mpiRank , mpiSource , MPI_MASTER , nZ , nZa , nZb , nZbc , Zb )
                  IF ( mpiRank == MPI_MASTER ) THEN

                     CALL io_write_vtk_zcoordinates ( 'psi-' , psiFileNo , psiFilePos , mpiSource , nZ , nZa , nZb , nZbc , Zb )

                  END IF

               END DO
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_psi3 ( mpiRank , mpiSource , MPI_MASTER , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )
                  IF ( mpiRank == MPI_MASTER ) THEN

                     CALL io_write_vtk_repsi ( 'psi-' , psiFileNo , psiFilePos , mpiSource , nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )

                  END IF

               END DO
               Psi3b = Psi3a
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_psi3 ( mpiRank , mpiSource , MPI_MASTER , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )
                  IF ( mpiRank == MPI_MASTER ) THEN 

                     CALL io_write_vtk_impsi ( 'psi-' , psiFileNo , psiFilePos , mpiSource , nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )

                  END IF

               END DO
               psiFileNo = psiFileNo + 1

            END IF

         END IF

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!        Compute 1st stage of generalized 4th-order Runge-Kutta (GRK4L): k_1 = f ( t_n , y_n )
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3a , K1 )
!        Compute the intermediate wave function for the 2nd stage of the GRK4L method: y_n + 0.5 * dT * k_1
!        Note that the intermediate wave function is only computed on the interior grid points assigned to each MPI process.
         CALL grk4_y_3d_stgx ( 1 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Compute V ( x , t_n + dT / 2 ), W ( t_n + dT / 2)
         tN = t0 + ( REAL ( n ) + 0.5 ) * dT
!        Compute 2nd stage of GRK4L: k_2 = f ( t_n + 0.5 * dT , y_n + 0.5 * dT * k_1 )
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3b , K2 )
!        Compute intermediate wave function for 3rd stage of GRK4L: y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2
         CALL grk4_y_3d_stgx ( 2 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Compute stage three ... k_3 = f ( t_n + 0.5 * dT , y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3b , K3 )
!        Compute intermediate wave function for 4th stage of GRK4L: y_n + ( 1 - lambda / 2 ) * dT * k_2 + ( lambda / 2 ) * dT * k_3
         CALL grk4_y_3d_stgx ( 3 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Calculate V ( x , t_n + dT ), W ( t_n + dT )
         tN = t0 + REAL ( n + 1 ) * dT
!        Compute the fourth stage ... k_4 = f ( t_n + dT , y_n + ( 1 - lamda / 2 ) * dT * k_2 + ( lamda / 2) * dT * k_3 
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3b , K4 )
!        Compute wave function at nth+1 time step ... y_{ n + 1 } = y_n + ( dT / 6 ) * [ k_1 + ( 4 - lambda ) * k_2 + lambda * k_3 + k_4 ]
         CALL grk4_y_3d_stgx ( 4 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Update wave function each time step: y_{ n + 1 } ---> y_n        
         Psi3a = Psi3b
!        If running in ITP mode, then also renormalize condensate wave function each time step
         IF ( itpOn .EQV. .TRUE. ) THEN

             CALL evua_normalize ( MPI_MASTER , mpiReal , mpiError , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Psi3a )

         END IF

      END DO

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

! --- END MAIN TIME PROPAGATION LOOP / BEGIN CLEAN UP TO STOP -------------

      !DEALLOCATE ( Vex3L )
      !DEALLOCATE ( ZL )
      !DEALLOCATE ( YL )
      !DEALLOCATE ( XL ) 

      DEALLOCATE ( MpiStatus )

      CALL DATE_AND_TIME ( stopDate , stopTime , stopZone , StopValues )
      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RUN STOPPED @ ', stopTime, ' ON ', stopDate, ' ... '

      END IF
      DEALLOCATE ( StopValues )
      DEALLOCATE ( StartValues )
      
      CALL MPI_FINALIZE ( mpiError )

! --- FORMAT STATEMENTS ---------------------------------------------------

900   FORMAT(1X,74(F23.15))

      STOP

      CONTAINS

         SUBROUTINE gpse_read_stdin ( fileName , fileUnit )

            IMPLICIT NONE

            CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
            
            INTEGER, INTENT ( IN ) :: fileUnit

            OPEN  ( UNIT = fileUnit , FILE = fileName , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )

               READ  ( UNIT = fileUnit , NML = gpseIn )

            CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

            RETURN

         END SUBROUTINE

         SUBROUTINE gpse_write_stdout_header ( )

            IMPLICIT NONE

            IF ( mpiRank == MPI_MASTER ) THEN

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
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        itpOn     = ', itpOn
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        rk4Lambda = ', rk4Lambda
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fdOrder   = ', fdOrder
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nTsteps   = ', nTsteps
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nTwrite   = ', nTwrite
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nX        = ', nX
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nY        = ', nY
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nZ        = ', nZ
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        t0        = ', t0
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xO        = ', xO
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yO        = ', yO
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zO        = ', zO
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dT        = ', dT
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dX        = ', dX
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dY        = ', dY
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dZ        = ', dZ
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wX        = ', wX
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wY        = ', wY
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZ        = ', wZ
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        gS        = ', gS
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nXpsi     = ', nXpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nYpsi     = ', nYpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nZpsi     = ', nZpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nRpsi     = ', nRpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        mLpsi     = ', mLpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xOpsi     = ', xOpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yOpsi     = ', yOpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zOpsi     = ', zOpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wXpsi     = ', wXpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wYpsi     = ', wYpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZpsi     = ', wZpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wRpsi     = ', wRpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        pXpsi     = ', pXpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        pYpsi     = ', pYpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        pZpsi     = ', pZpsi
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xOvex     = ', xOvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yOvex     = ', yOvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zOvex     = ', zOvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        rOvex     = ', rOvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fXvex     = ', fXvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fYvex     = ', fYvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fZvex     = ', fZvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wXvex     = ', wXvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wYvex     = ', wYvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZvex     = ', wZvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wRvex     = ', wRvex
               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     BROADCASTING INPUT PARAMETERS TO ALL MPI PROCESSES ... '

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_select_default_kinds ( )

            IMPLICIT NONE

            SELECT CASE ( INT_DEFAULT_KIND )

               CASE ( INT_8 )

                  mpiInt = MPI_INTEGER1

               CASE ( INT_16 )

                  mpiInt = MPI_INTEGER2

               CASE ( INT_32 )

                  mpiInt = MPI_INTEGER

               CASE ( INT_64 )

                  mpiInt = MPI_INTEGER8

               CASE DEFAULT

                  mpiInt = -1

            END SELECT

            SELECT CASE ( REAL_DEFAULT_KIND )

               CASE ( REAL_32 )

                  mpiReal = MPI_REAL

               CASE ( REAL_64 )

                  mpiReal= MPI_DOUBLE_PRECISION

               CASE ( REAL_128 )

                  mpiReal = MPI_REAL16

               CASE DEFAULT

                  mpiReal = -1

            END SELECT

            SELECT CASE ( CMPLX_DEFAULT_KIND )

               CASE ( REAL_32 )

                  mpiCmplx = MPI_COMPLEX

               CASE ( REAL_64 )

                  mpiCmplx = MPI_DOUBLE_COMPLEX

               CASE ( REAL_128 )

                  mpiCmplx = MPI_COMPLEX32

               CASE DEFAULT

                  mpiCmplx = -1

            END SELECT

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_bcast_inputs ( mpiMaster , mpiInt , mpiReal , mpiCmplx , mpiError )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: mpiMaster
            INTEGER, INTENT ( IN    ) :: mpiInt
            INTEGER, INTENT ( IN    ) :: mpiReal
            INTEGER, INTENT ( IN    ) :: mpiCmplx
            INTEGER, INTENT ( INOUT ) :: mpiError

            CALL MPI_BCAST ( itpOn     , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( rk4Lambda , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fdOrder   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( quadRule  , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nTsteps   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nTwrite   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nX        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nXbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nY        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nYbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nZ        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nZbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( t0        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( xO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( yO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( zO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dT        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dX        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dY        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dZ        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wX        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wY        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wZ        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( gS        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            CALL MPI_BCAST ( psiInput  , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiOutput , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiFileNo , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiInit   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nXpsi     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nYpsi     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nZpsi     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nRpsi     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( mLpsi     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( xOpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( yOpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( zOpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wXpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wYpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wZpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wRpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( pXpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( pYpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( pZpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            CALL MPI_BCAST ( vexInput  , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexOutput , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFileNo , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexInit   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( xOvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( yOvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( zOvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( rOvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fXvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fYvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fZvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wXvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wYvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wZvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wRvex     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_copy_q ( mpiRank , mpiSource , mpiDestination , nQ , nQa , nQb , nQbc , Q )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: mpiRank
            INTEGER, INTENT ( IN ) :: mpiSource
            INTEGER, INTENT ( IN ) :: mpiDestination
            INTEGER, INTENT ( IN ) :: nQ
            INTEGER, INTENT ( IN ) :: nQa
            INTEGER, INTENT ( IN ) :: nQb
            INTEGER, INTENT ( IN ) :: nQbc 

            REAL, DIMENSION ( nQa - nQbc : nQb + nQbc ), INTENT ( INOUT ) :: Q

            INTEGER :: j

            IF ( ( mpiRank == mpiDestination ) .AND. ( mpiRank /= mpiSource ) ) THEN

               DO j = nQa , nQb

                  CALL MPI_RECV ( Q ( j ) , 1 , mpiReal , mpiSource , mpiSource , MPI_COMM_WORLD , MpiStatus , mpiError )

               END DO

            ELSE IF ( ( mpiRank == mpiSource ) .AND. ( mpiRank /= mpiDestination ) ) THEN

               DO j = nQa , nQb

                  CALL MPI_SSEND ( Q ( j ) , 1 , mpiReal , mpiDestination , mpiRank , MPI_COMM_WORLD , mpiError )

               END DO

            END IF

            CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

            RETURN

         END SUBROUTINE 

         SUBROUTINE mpi_copy_vex3 ( mpiRank , mpiSource , mpiDestination , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 ) 

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: mpiRank
            INTEGER, INTENT ( IN ) :: mpiSource
            INTEGER, INTENT ( IN ) :: mpiDestination
            INTEGER, INTENT ( IN ) :: nXa
            INTEGER, INTENT ( IN ) :: nXb
            INTEGER, INTENT ( IN ) :: nXbc 
            INTEGER, INTENT ( IN ) :: nYa
            INTEGER, INTENT ( IN ) :: nYb
            INTEGER, INTENT ( IN ) :: nYbc 
            INTEGER, INTENT ( IN ) :: nZa
            INTEGER, INTENT ( IN ) :: nZb
            INTEGER, INTENT ( IN ) :: nZbc 

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Vex3

            INTEGER :: j , k , l

            IF ( ( mpiRank == mpiDestination ) .AND. ( mpiRank /= mpiSource ) ) THEN

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        CALL MPI_RECV ( Vex3 ( j , k , l ) , 1 , mpiReal , mpiSource , mpiSource , MPI_COMM_WORLD , MpiStatus , mpiError )

                     END DO ! possible alternative: no loop over j; call MPI_RECV ( Real3b ( nXa , nYa , l ) , nXb - nXa + 1 , mpiReal ... )

                  END DO

               END DO

            ELSE IF ( ( mpiRank == mpiSource ) .AND. ( mpiRank /= mpiDestination ) ) THEN

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        CALL MPI_SSEND ( Vex3 ( j , k , l ) , 1 , mpiReal , mpiDestination , mpiRank , MPI_COMM_WORLD , mpiError )

                     END DO ! possible alternative: no loop over j; call MPI_SSEND ( Real3a ( nXa , nYa , l ) , nXb - nXa + 1 , mpiReal ... )

                  END DO

               END DO

            END IF

            CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

            RETURN
 
         END SUBROUTINE

         SUBROUTINE mpi_copy_psi3 ( mpiRank , mpiSource , mpiDestination , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: mpiRank
            INTEGER, INTENT ( IN ) :: mpiSource
            INTEGER, INTENT ( IN ) :: mpiDestination
            INTEGER, INTENT ( IN ) :: nXa
            INTEGER, INTENT ( IN ) :: nXb
            INTEGER, INTENT ( IN ) :: nXbc 
            INTEGER, INTENT ( IN ) :: nYa
            INTEGER, INTENT ( IN ) :: nYb
            INTEGER, INTENT ( IN ) :: nYbc 
            INTEGER, INTENT ( IN ) :: nZa
            INTEGER, INTENT ( IN ) :: nZb
            INTEGER, INTENT ( IN ) :: nZbc 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

            INTEGER :: j , k , l

            IF ( ( mpiRank == mpiDestination ) .AND. ( mpiRank /= mpiSource ) ) THEN

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        CALL MPI_RECV ( Psi3 ( j , k , l ) , 1 , mpiCmplx , mpiSource , mpiSource , MPI_COMM_WORLD , MpiStatus , mpiError )

                     END DO ! possible alternative: no loop over j; call MPI_RECV ( Real3b ( nXa , nYa , l ) , nXb - nXa + 1 , mpiReal ... )

                  END DO

               END DO

            ELSE IF ( ( mpiRank == mpiSource ) .AND. ( mpiRank /= mpiDestination ) ) THEN

               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        CALL MPI_SSEND ( Psi3 ( j , k , l ) , 1 , mpiCmplx , mpiDestination , mpiRank , MPI_COMM_WORLD , mpiError )

                     END DO ! possible alternative: no loop over j; call MPI_SSEND ( Real3a ( nXa , nYa , l ) , nXb - nXa + 1 , mpiReal ... )

                  END DO

               END DO

            END IF

            CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_gather_real3 ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Real3p , Real3f )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nXa
         INTEGER, INTENT ( IN ) :: nXb
         INTEGER, INTENT ( IN ) :: nXbc 
         INTEGER, INTENT ( IN ) :: nY 
         INTEGER, INTENT ( IN ) :: nYa
         INTEGER, INTENT ( IN ) :: nYb
         INTEGER, INTENT ( IN ) :: nYbc 
         INTEGER, INTENT ( IN ) :: nZ 
         INTEGER, INTENT ( IN ) :: nZa
         INTEGER, INTENT ( IN ) :: nZb
         INTEGER, INTENT ( IN ) :: nZbc

         REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Real3p
         REAL, DIMENSION ( 1 : nX , 1 : nY , 1 : nZ ), INTENT ( INOUT ) :: Real3f

         INTEGER :: mpiSource
         INTEGER :: nZaSource
         INTEGER :: nZbSource
         INTEGER :: j , k , l

         IF ( mpiRank == MPI_MASTER ) THEN ! receive external potential and wave function data from all MPI processes

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Real3f ( j , k , l ) = Real3p ( j , k , l )

                  END DO

               END DO
!$OMP          END PARALLEL DO

            END DO
!$OMP       END PARALLEL DO

            DO mpiSource = 1 , mpiProcesses - 1

               CALL MPI_RECV ( nZaSource , 1 , mpiInt , mpiSource , 0 , MPI_COMM_WORLD , MpiStatus , mpiError )
               CALL MPI_RECV ( nZbSource , 1 , mpiInt , mpiSource , 1 , MPI_COMM_WORLD , MpiStatus , mpiError )
               DO l = nZaSource , nZbSource

                  DO k = 1 , nY

                     DO j = 1 , nX

                        CALL MPI_RECV ( Real3f ( j , k , l ) , 1 , mpiReal , mpiSource , 2 , MPI_COMM_WORLD , MpiStatus , mpiError )

                     END DO

                  END DO

               END DO

            END DO

         ELSE ! send real-valued, three-dimensional array data to MPI_MASTER

            CALL MPI_SSEND ( nZa , 1 , mpiInt , MPI_MASTER , 0 , MPI_COMM_WORLD , mpiError )
            CALL MPI_SSEND ( nZb , 1 , mpiInt , MPI_MASTER , 1 , MPI_COMM_WORLD , mpiError )
            DO l = nZa , nZb

               DO k = 1 , nY

                  DO j = 1 , nX

                     CALL MPI_SSEND ( Real3p ( j , k , l ) , 1 , mpiReal , MPI_MASTER , 2 , MPI_COMM_WORLD , mpiError )

                  END DO

               END DO

            END DO

         END IF

         RETURN

         END SUBROUTINE

         SUBROUTINE mpi_gather_cmplx3 ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Cmplx3p , Cmplx3f )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nXa
         INTEGER, INTENT ( IN ) :: nXb
         INTEGER, INTENT ( IN ) :: nXbc
         INTEGER, INTENT ( IN ) :: nY
         INTEGER, INTENT ( IN ) :: nYa
         INTEGER, INTENT ( IN ) :: nYb
         INTEGER, INTENT ( IN ) :: nYbc
         INTEGER, INTENT ( IN ) :: nZ
         INTEGER, INTENT ( IN ) :: nZa
         INTEGER, INTENT ( IN ) :: nZb
         INTEGER, INTENT ( IN ) :: nZbc

         COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Cmplx3p
         COMPLEX, DIMENSION ( 1 : nX , 1 : nY , 1 : nZ ), INTENT ( INOUT ) :: Cmplx3f

         INTEGER :: mpiSource
         INTEGER :: nZaSource
         INTEGER :: nZbSource
         INTEGER :: j , k , l

         IF ( mpiRank == MPI_MASTER ) THEN ! receive complex-valued, three-dimensional array data from all MPI processes

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Cmplx3f ( j , k , l ) = Cmplx3p ( j , k , l )

                  END DO

               END DO
!$OMP          END PARALLEL DO

            END DO
!$OMP       END PARALLEL DO
            DO mpiSource = 1 , mpiProcesses - 1

               CALL MPI_RECV ( nZaSource , 1 , mpiInt , mpiSource , 0 , MPI_COMM_WORLD , MpiStatus , mpiError )
               CALL MPI_RECV ( nZbSource , 1 , mpiInt , mpiSource , 1 , MPI_COMM_WORLD , MpiStatus , mpiError )

               DO l = nZaSource , nZbSource

                  DO k = 1 , nY

                     DO j = 1 , nX

                        CALL MPI_RECV ( Cmplx3f ( j , k , l ) , 1 , mpiCmplx , mpiSource , 2 , MPI_COMM_WORLD , MpiStatus , mpiError )

                     END DO

                  END DO

               END DO

            END DO

         ELSE ! send complex-valued, three-dimensional array data to MPI_MASTER

            CALL MPI_SSEND ( nZa , 1 , mpiInt , MPI_MASTER , 0 , MPI_COMM_WORLD , mpiError )
            CALL MPI_SSEND ( nZb , 1 , mpiInt , MPI_MASTER , 1 , MPI_COMM_WORLD , mpiError )

            DO l = nZa , nZb

               DO k = 1 , nY

                  DO j = 1 , nX

                     CALL MPI_SSEND ( Cmplx3p ( j , k , l ) , 1 , mpiCmplx , MPI_MASTER , 2 , MPI_COMM_WORLD , mpiError )

                  END DO

               END DO

            END DO

         END IF

         RETURN

         END SUBROUTINE

         SUBROUTINE mpi_scatter_real3 ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Real3p , Real3f )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nXa
         INTEGER, INTENT ( IN ) :: nXb
         INTEGER, INTENT ( IN ) :: nXbc
         INTEGER, INTENT ( IN ) :: nY
         INTEGER, INTENT ( IN ) :: nYa
         INTEGER, INTENT ( IN ) :: nYb
         INTEGER, INTENT ( IN ) :: nYbc
         INTEGER, INTENT ( IN ) :: nZ
         INTEGER, INTENT ( IN ) :: nZa
         INTEGER, INTENT ( IN ) :: nZb
         INTEGER, INTENT ( IN ) :: nZbc

         REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Real3p
         REAL, DIMENSION ( 1 : nX , 1 : nY , 1 : nZ ), INTENT ( IN ) :: Real3f

         INTEGER :: mpiDest
         INTEGER :: nZaDest
         INTEGER :: nZbDest
         INTEGER :: j , k , l

         IF ( mpiRank == MPI_MASTER ) THEN ! send real-valued, three-dimensional array data to all MPI processes

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Real3p ( j , k , l ) = Real3f ( j , k , l )

                  END DO

               END DO
!$OMP          END PARALLEL DO 

            END DO
!$OMP       END PARALLEL DO
            DO mpiDest = 1 , mpiProcesses - 1

               CALL MPI_RECV ( nZaDest , 1 , mpiInt , mpiDest , 0 , MPI_COMM_WORLD , MpiStatus , mpiError )
               CALL MPI_RECV ( nZbDest , 1 , mpiInt , mpiDest , 1 , MPI_COMM_WORLD , MpiStatus , mpiError )

               DO l = nZaDest , nZbDest

                  DO k = 1 , nY

                     DO j = 1 , nX

                        CALL MPI_SSEND ( Real3f ( j , k , l ) , 1 , mpiReal , mpiDest , 2 , MPI_COMM_WORLD , mpiError )

                     END DO

                  END DO

               END DO

            END DO

         ELSE ! receive real-valued, three-dimensional array data from MPI_MASTER

            CALL MPI_SSEND ( nZa , 1 , mpiInt , MPI_MASTER , 0 , MPI_COMM_WORLD , mpiError )
            CALL MPI_SSEND ( nZb , 1 , mpiInt , MPI_MASTER , 1 , MPI_COMM_WORLD , mpiError )
            DO l = nZa , nZb

               DO k = 1 , nY

                  DO j = 1 , nY

                     CALL MPI_RECV ( Real3p ( j , k , l ) , 1 , mpiReal , MPI_MASTER , 2 , MPI_COMM_WORLD , MpiStatus , mpiError )

                  END DO

               END DO

            END DO

         END IF

         RETURN

         END SUBROUTINE

         SUBROUTINE mpi_scatter_cmplx3 ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Cmplx3p , Cmplx3f )

         IMPLICIT NONE

         INTEGER, INTENT ( IN ) :: nX
         INTEGER, INTENT ( IN ) :: nXa
         INTEGER, INTENT ( IN ) :: nXb
         INTEGER, INTENT ( IN ) :: nXbc
         INTEGER, INTENT ( IN ) :: nY
         INTEGER, INTENT ( IN ) :: nYa
         INTEGER, INTENT ( IN ) :: nYb
         INTEGER, INTENT ( IN ) :: nYbc
         INTEGER, INTENT ( IN ) :: nZ
         INTEGER, INTENT ( IN ) :: nZa
         INTEGER, INTENT ( IN ) :: nZb
         INTEGER, INTENT ( IN ) :: nZbc

         COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Cmplx3p
         COMPLEX, DIMENSION ( 1 : nX , 1 : nY , 1 : nZ ), INTENT ( IN ) :: Cmplx3f

         INTEGER :: mpiDest
         INTEGER :: nZaDest
         INTEGER :: nZbDest
         INTEGER :: j , k , l

         IF ( mpiRank == MPI_MASTER ) THEN ! send complex-valued, three-dimensional array data to all MPI processes

!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
            DO l = nZa , nZb

!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
               DO k = nYa , nYb

                  DO j = nXa , nXb

                     Cmplx3p ( j , k , l ) = Cmplx3f ( j , k , l )

                  END DO

               END DO
!$OMP          END PARALLEL DO 

            END DO
!$OMP       END PARALLEL DO
            DO mpiDest = 1 , mpiProcesses - 1

               CALL MPI_RECV ( nZaDest , 1 , mpiInt , mpiDest , 0 , MPI_COMM_WORLD , MpiStatus , mpiError )
               CALL MPI_RECV ( nZbDest , 1 , mpiInt , mpiDest , 1 , MPI_COMM_WORLD , MpiStatus , mpiError )

               DO l = nZaDest , nZbDest

                  DO k = 1 , nY

                     DO j = 1 , nX

                        CALL MPI_SSEND ( Cmplx3f ( j , k , l ) , 1 , mpiCmplx , mpiDest , 2 , MPI_COMM_WORLD , mpiError )

                     END DO

                  END DO

               END DO

            END DO

         ELSE ! receive complex-valued, three-dimensional array data from MPI_MASTER

            CALL MPI_SSEND ( nZa , 1 , mpiInt , MPI_MASTER , 0 , MPI_COMM_WORLD , mpiError )
            CALL MPI_SSEND ( nZb , 1 , mpiInt , MPI_MASTER , 1 , MPI_COMM_WORLD , mpiError )
            DO l = nZa , nZb

               DO k = 1 , nY

                  DO j = 1 , nY

                     CALL MPI_RECV ( Cmplx3p ( j , k , l ) , 1 , mpiCmplx , MPI_MASTER , 2 , MPI_COMM_WORLD , MpiStatus , mpiError )

                  END DO

               END DO

            END DO

         END IF

         RETURN

         END SUBROUTINE

!         SUBROUTINE compute_f ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            IMPLICIT NONE
!
!            INTEGER, INTENT ( IN ) :: nXa
!            INTEGER, INTENT ( IN ) :: nXb
!            INTEGER, INTENT ( IN ) :: nXbc
!            INTEGER, INTENT ( IN ) :: nYa
!            INTEGER, INTENT ( IN ) :: nYb
!            INTEGER, INTENT ( IN ) :: nYbc
!            INTEGER, INTENT ( IN ) :: nZa
!            INTEGER, INTENT ( IN ) :: nZb
!            INTEGER, INTENT ( IN ) :: nZbc
!
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3
!
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F3
!
!            IF ( fdOrder == 2 ) THEN 
!
!               CALL f_gp_3d_rrf_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            ELSE IF ( fdOrder == 4 ) THEN
!
!               CALL f_gp_3d_rrf_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            ELSE IF ( fdOrder == 6 ) THEN
!
!               CALL f_gp_3d_rrf_cd6 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            ELSE IF ( fdOrder == 8 ) THEN
!
!               CALL f_gp_3d_rrf_cd8 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            ELSE
!
!               WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : compute_f : ERROR - fdOrder not supported.'
!
!            END IF
!
!            RETURN
!
!         END SUBROUTINE
!
!         SUBROUTINE f_gp_3d_rrf_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 ) 
!
!            IMPLICIT NONE
!
!            INTEGER, INTENT ( IN ) :: nXa
!            INTEGER, INTENT ( IN ) :: nXb
!            INTEGER, INTENT ( IN ) :: nXbc
!            INTEGER, INTENT ( IN ) :: nYa
!            INTEGER, INTENT ( IN ) :: nYb
!            INTEGER, INTENT ( IN ) :: nYbc
!            INTEGER, INTENT ( IN ) :: nZa
!            INTEGER, INTENT ( IN ) :: nZb
!            INTEGER, INTENT ( IN ) :: nZbc
!
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3
!
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F3
!
!            INTEGER :: j , k , l
!
!!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!            DO l = nZa , nZb  
!
!!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!               DO k = nYa , nYb  
!
!                  DO j = nXa , nXb  
!
!                     F3 ( j , k , l ) = &
!                        & CMPLX ( 0.5 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 0.5 / dZ**2 ) * Psi3 ( j , k , l - 1 ) + &  
!                        & CMPLX ( 0.5 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 0.5 / dY**2 ) * Psi3 ( j , k - 1 , l ) + &  
!                        & CMPLX ( 0.5 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 0.5 / dX**2 ) * Psi3 ( j - 1 , k , l ) - &  
!                        & CMPLX ( 0.0 , 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 + Vex3 ( j , k , l ) + gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
!                        & CMPLX ( 0.5 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 0.5 / dX**2 ) * Psi3 ( j + 1 , k , l ) + &  
!                        & CMPLX ( 0.5 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 0.5 / dY**2 ) * Psi3 ( j , k + 1 , l ) + &  
!                        & CMPLX ( 0.5 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 0.5 / dZ**2 ) * Psi3 ( j , k , l + 1 )  
!
!                  END DO
!
!               END DO
!!$OMP          END PARALLEL DO
!
!            END DO
!!$OMP       END PARALLEL DO 
!
!            RETURN
!
!         END SUBROUTINE
!
!         SUBROUTINE f_gp_3d_rrf_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 ) 
!
!            IMPLICIT NONE
!
!            INTEGER, INTENT ( IN ) :: nXa
!            INTEGER, INTENT ( IN ) :: nXb
!            INTEGER, INTENT ( IN ) :: nXbc
!            INTEGER, INTENT ( IN ) :: nYa
!            INTEGER, INTENT ( IN ) :: nYb
!            INTEGER, INTENT ( IN ) :: nYbc
!            INTEGER, INTENT ( IN ) :: nZa
!            INTEGER, INTENT ( IN ) :: nZb
!            INTEGER, INTENT ( IN ) :: nZbc
!
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3 
!
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3 
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F3
!
!            INTEGER :: j , k , l
!
!!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!            DO l = nZa , nZb
!
!!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!               DO k = nYa , nYb
!
!                  DO j = nXa , nXb
!
!                     F3 ( j , k , l ) = &
!                        & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi3 ( j , k , l - 2 ) + &
!                        & CMPLX ( 0.75 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi3 ( j , k , l - 1 ) + &
!                        & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi3 ( j , k - 2 , l ) + &
!                        & CMPLX ( 0.75 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi3 ( j , k - 1 , l ) + &
!                        & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi3 ( j - 2 , k , l ) + &
!                        & CMPLX ( 0.75 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi3 ( j - 1 , k , l ) - &
!                        & CMPLX ( 0.0 , 1.25 * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex3 ( j , k , l ) + &
!                        &    gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
!                        & CMPLX ( 0.75 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi3 ( j + 1 , k , l ) + &
!                        & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi3 ( j + 2 , k , l ) + &
!                        & CMPLX ( 0.75 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi3 ( j , k + 1 , l ) + &
!                        & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi3 ( j , k + 2 , l ) + &
!                        & CMPLX ( 0.75 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi3 ( j , k , l + 1 ) + &
!                        & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi3 ( j , k , l + 2 )
!
!                  END DO
!
!               END DO
!!$OMP          END PARALLEL DO
!
!            END DO
!!$OMP       END PARALLEL DO
!
!            RETURN
!
!         END SUBROUTINE
!
!         SUBROUTINE f_gp_3d_rrf_cd6 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            IMPLICIT NONE
!
!            INTEGER, INTENT ( IN ) :: nXa
!            INTEGER, INTENT ( IN ) :: nXb
!            INTEGER, INTENT ( IN ) :: nXbc
!            INTEGER, INTENT ( IN ) :: nYa
!            INTEGER, INTENT ( IN ) :: nYb
!            INTEGER, INTENT ( IN ) :: nYbc
!            INTEGER, INTENT ( IN ) :: nZa
!            INTEGER, INTENT ( IN ) :: nZb
!            INTEGER, INTENT ( IN ) :: nZbc
!
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3
!
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: F3
!
!            INTEGER :: j , k , l
!
!!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!            DO l = nZa , nZb
!
!!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!               DO k = nYa , nYb
!
!                  DO j = nXa , nXb
!
!                     F3 ( j , k , l ) = &
!                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi3 ( j , k , l - 3 ) + &
!                        & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * Psi3 ( j , k , l - 2 ) + &
!                        & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * Psi3 ( j , k , l - 1 ) + &
!                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi3 ( j , k - 3 , l ) + &
!                        & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * Psi3 ( j , k - 2 , l ) + &
!                        & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * Psi3 ( j , k - 1 , l ) + &
!                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi3 ( j - 3 , k , l ) + &
!                        & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * Psi3 ( j - 2 , k , l ) + &
!                        & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * Psi3 ( j - 1 , k , l ) - &
!                        & CMPLX ( 0.0 , ( 49.0 / 36.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex3 ( j , k , l ) + gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
!                        & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * Psi3 ( j + 1 , k , l ) + &
!                        & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * Psi3 ( j + 2 , k , l ) + &
!                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi3 ( j + 3 , k , l ) + &
!                        & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * Psi3 ( j , k + 1 , l ) + &
!                        & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * Psi3 ( j , k + 2 , l ) + &
!                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi3 ( j , k + 3 , l ) + &
!                        & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * Psi3 ( j , k , l + 1 ) + &
!                        & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * Psi3 ( j , k , l + 2 ) + &
!                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi3 ( j , k , l + 3 )   
!
!                  END DO
!
!               END DO
!!$OMP          END PARALLEL DO
!
!            END DO
!!$OMP       END PARALLEL DO
!
!            RETURN
!
!         END SUBROUTINE
!
!         SUBROUTINE f_gp_3d_rrf_cd8 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , X , Y , Z , Vex3 , Psi3 , F3 )
!
!            IMPLICIT NONE
!
!            INTEGER, INTENT ( IN ) :: nXa
!            INTEGER, INTENT ( IN ) :: nXb
!            INTEGER, INTENT ( IN ) :: nXbc
!            INTEGER, INTENT ( IN ) :: nYa
!            INTEGER, INTENT ( IN ) :: nYb
!            INTEGER, INTENT ( IN ) :: nYbc
!            INTEGER, INTENT ( IN ) :: nZa
!            INTEGER, INTENT ( IN ) :: nZb
!            INTEGER, INTENT ( IN ) :: nZbc
!
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
!            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
!            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
!            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3
!
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3
!            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F3
!
!            INTEGER :: j , k , l
!
!!$OMP       PARALLEL DO IF ( nZa /= nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!            DO l = nZa , nZb
!
!!$OMP          PARALLEL DO IF ( nZa == nZb ) DEFAULT ( SHARED ) SCHEDULE ( STATIC )
!               DO k = nYa , nYb
!
!                  DO j = nXa , nXb
!
!                     F3 ( j , k , l ) = &
!                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi3 ( j , k , l - 4 ) + &
!                        & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 105.0 * dZ ) ,  4.0 / ( 315.0 * dZ**2 ) ) * Psi3 ( j , k , l - 3 ) + &
!                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dZ ) , -1.0 / ( 10.0 * dZ**2 ) ) * Psi3 ( j , k , l - 2 ) + &
!                        & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dZ ) ,  4.0 / ( 5.0 * dZ**2 ) ) * Psi3 ( j , k , l - 1 ) + &
!                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 280.0 * dY ) , -1.0 / ( 1120 * dY**2 ) ) * Psi3 ( j , k - 4 , l ) + &
!                        & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi3 ( j , k - 3 , l ) + &
!                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi3 ( j , k - 2 , l ) + &
!                        & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * Psi3 ( j , k - 1 , l ) + &
!                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * Psi3 ( j - 4 , k , l ) + &
!                        & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * Psi3 ( j - 3 , k , l ) + &
!                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi3 ( j - 2 , k , l ) + &
!                        & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi3 ( j - 1 , k , l ) - &
!                        & CMPLX ( 0.0 , ( 205.0 / 144.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex3 ( j , k , l ) + gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
!                        & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi3 ( j + 1 , k , l ) + &
!                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi3 ( j + 2 , k , l ) + &
!                        & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * Psi3 ( j + 3 , k , l ) + &
!                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * Psi3 ( j + 4 , k , l ) + &
!                        & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * Psi3 ( j , k + 1 , l ) + &
!                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi3 ( j , k + 2 , l ) + & 
!                        & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi3 ( j , k + 3 , l ) + &
!                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 280.0 * dY ) , -1.0 / ( 1120.0 * dY**2 ) ) * Psi3 ( j , k + 4 , l ) + &
!                        & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi3 ( j , k , l + 1 ) + &
!                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi3 ( j , k , l + 2 ) + &
!                        & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi3 ( j , k , l + 3 ) + &
!                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi3 ( j , k , l + 4 ) 
!
!                  END DO
!
!               END DO
!!$OMP          END PARALLEL DO
!
!            END DO
!!$OMP       END PARALLEL DO
!
!            RETURN
!
!         END SUBROUTINE
!
         SUBROUTINE mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3 ) 

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: nX
            INTEGER, INTENT ( IN ) :: nXa
            INTEGER, INTENT ( IN ) :: nXb
            INTEGER, INTENT ( IN ) :: nXbc
            INTEGER, INTENT ( IN ) :: nY
            INTEGER, INTENT ( IN ) :: nYa
            INTEGER, INTENT ( IN ) :: nYb
            INTEGER, INTENT ( IN ) :: nYbc
            INTEGER, INTENT ( IN ) :: nZ
            INTEGER, INTENT ( IN ) :: nZa
            INTEGER, INTENT ( IN ) :: nZb
            INTEGER, INTENT ( IN ) :: nZbc

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN 

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 0 , MPI_COMM_WORLD, mpiError )

               END IF

               IF ( mpiRank /= MPI_MASTER ) THEN 

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 1 , MPI_COMM_WORLD , MpiStatus , mpiError )

               END IF

            ELSE 

               CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 0 , MPI_COMM_WORLD , MpiStatus , mpiError )

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 1 , MPI_COMM_WORLD , mpiError )

               END IF

            END IF

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN

               IF ( mpiRank /= MPI_MASTER ) THEN

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 2 , MPI_COMM_WORLD , mpiError )

               END IF

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 3 , MPI_COMM_WORLD , MpiStatus , mpiError )

               END IF

            ELSE

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 2 , MPI_COMM_WORLD , MpiStatus , mpiError )

               END IF
               CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 3 , MPI_COMM_WORLD , mpiError )

            END IF

            RETURN

         END SUBROUTINE

      END PROGRAM

! =========================================================================
