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
!     Thursday, July 17th, 2014
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

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.2.0'
      CHARACTER ( LEN = * ), PARAMETER :: GPSE_LAST_UPDATED = 'Thursday, July 17th, 2014'

      INTEGER, PARAMETER :: INT_DEFAULT_KIND   = KIND ( 0 ) 
      INTEGER, PARAMETER :: REAL_DEFAULT_KIND  = KIND ( 0.0 )
      INTEGER, PARAMETER :: CMPLX_DEFAULT_KIND = KIND ( CMPLX ( 0.0 , 0.0 ) )
      INTEGER, PARAMETER :: MAX_LEN_FILENAME   = 255
      INTEGER, PARAMETER :: MPI_MASTER         = 0

! --- PARAMETER DEFINITIONS -----------------------------------------------
! --- VARIABLE DECLARATIONS -----------------------------------------------

      LOGICAL :: itpOn  = .FALSE. ! Perform imaginary time propagation? .TRUE. = Yes ; .FALSE. = No
      LOGICAL :: fullIO = .FALSE. 

      CHARACTER ( LEN = 8  ) :: startDate = 'NONE'
      CHARACTER ( LEN = 10 ) :: startTime = 'NONE'
      CHARACTER ( LEN = 5  ) :: startZone = 'NONE'
      CHARACTER ( LEN = 8  ) :: stopDate  = 'NONE'
      CHARACTER ( LEN = 10 ) :: stopTime  = 'NONE'
      CHARACTER ( LEN = 5  ) :: stopZone  = 'NONE'

      INTEGER :: fmtIO         = -1
      INTEGER :: rk4Lambda     = -1 ! 1 = Tan-Chen Lambda-1 ; 2 = Classical 4th-Order Runge-Kutta ; 3 = Tan-Chen Lambda-3 ; 4 = England ; 5 = Tan-Chen Lambda-5
      INTEGER :: fdOrder       = -1 ! 2 = 2nd-Order Central Differences ( CD ); 4 = 4th-Order CD ; 6 = 6th-Order CD; 8 = 8th-Order CD
      INTEGER :: quadRule      = -1 ! 1 = Rectangle Rule ; 2 = Trapezoidal Rule
      INTEGER :: mpiCmplx      = -1 
      INTEGER :: mpiError      = -1 
      INTEGER :: mpiErrorCode  = -1 
      INTEGER :: mpiErrorClass = -1
      INTEGER :: mpiInt        = -1 
      INTEGER :: mpiProcesses  = -1 
      INTEGER :: mpiProvided   = -1 
      INTEGER :: mpiRank       = -1
      INTEGER :: mpiReal       = -1 
      INTEGER :: ompThreads    = -1 ! Number of threads in OpenMP parallel regions 
      INTEGER :: ompThreadID   = -1 
      INTEGER :: nTsteps       = -1 ! Number of time steps in simulation
      INTEGER :: nTwrite       = -1 ! Period of writes to disk; i.e., number of time steps between writes to disk
      INTEGER :: fileNumber    = -1
      INTEGER :: j , k , l , m , n  ! Reserved loop counters

      REAL :: tN = 0.0 ! Simulation time at nth time step
      REAL :: t0 = 0.0 ! Simulation time at the beginning of the simulation
      REAL :: dT = 0.0 ! Interval of a time step
      REAL :: wX = 0.0 ! X-component of the rotating reference frame's angular velocity vector
      REAL :: wY = 0.0 ! Y-component of the rotating reference frame's angular velocity vector
      REAL :: wZ = 0.0 ! Z-component of the rotating reference frame's angular velocity vector
      REAL :: gS = 0.0 ! Nonlinear atom-atom interaction coupling constant

      COMPLEX :: zDt = CMPLX ( 0.0 , 0.0 ) ! Stores simulation time step in complex form; Used for imaginary time propagation

! --- VARIABLE DEFINITIONS ------------------------------------------------
! --- ARRAY DECLARATIONS --------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIstatus
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopValues

      REAL, ALLOCATABLE, DIMENSION ( :         ) :: X
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Y
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Z
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Xm
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Ym
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Zm
      REAL, ALLOCATABLE, DIMENSION ( : , : , : ) :: Vex
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: VexM

      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K4
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: PsiA
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: PsiB
      COMPLEX, ALLOCATABLE, DIMENSION ( :         ) :: PsiM

! --- ARRAY DEFINITIONS ---------------------------------------------------

! --- FUNCTION AND SUBROUTINE DECLARATIONS --------------------------------

      INTEGER :: OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM

! --- NAMELIST DECLARATIONS -----------------------------------------------

      NAMELIST /nmlGPSEin/ itpOn , fullIO , fmtIO , rk4Lambda , fdOrder , quadRule , nTsteps , nTwrite , t0 , dT , wX , wY , wZ , gS

! --- NAMELIST DEFINITIONS ------------------------------------------------

! --- FUNCTION AND SUBROUTINE DEFINITIONS ---------------------------------

! --- MAIN PROGRAM --------------------------------------------------------

      ALLOCATE ( MPIstatus ( MPI_STATUS_SIZE ) ) 
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

      ALLOCATE ( StartValues ( 8 ) )
      ALLOCATE ( StopValues  ( 8 ) )
      CALL DATE_AND_TIME ( startDate , startTime , startZone , StartValues )
      
      IF ( mpiRank == MPI_MASTER ) THEN

         OPEN  ( UNIT = 500 , FILE = 'gpse.in' , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
         READ  ( UNIT = 500 , NML = nmlGPSEin )
         CLOSE ( UNIT = 500 , STATUS = 'KEEP' )

         CALL grid_read_inputs ( )
         CALL grid_bound_cond_size ( fdOrder )
         CALL vex_read_inputs ( )
         CALL psi_read_inputs ( )

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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fullIO    = ', fullIO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fmtIO     = ', fmtIO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        rk4Lambda = ', rK4Lambda
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        fdOrder   = ', fdOrder
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        quadRule  = ', quadRule
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nTsteps   = ', nTsteps
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nTwrite   = ', nTwrite
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        t0        = ', t0
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dT        = ', dT
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wX        = ', wX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wY        = ', wY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZ        = ', wZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        gS        = ', gS
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nX        = ', nX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nY        = ', nY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nZ        = ', nZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xO        = ', xO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yO        = ', yO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zO        = ', zO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dX        = ', dX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dY        = ', dY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dZ        = ', dZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexRead   = ', vexRead
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexFmtIn  = ', vexFmtIn
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexInit   = ', vexInit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexXo     = ', vexXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexYo     = ', vexYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexZo     = ', vexZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexRo     = ', vexRo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexFx     = ', vexFx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexFy     = ', vexFy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexFz     = ', vexFz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexWx     = ', vexWx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexWy     = ', vexWy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexWz     = ', vexWz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexWr     = ', vexWr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiInit   = ', psiInit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiNx     = ', psiNx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiNy     = ', psiNy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiNz     = ', psiNz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiNr     = ', psiNr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiMl     = ', psiMl
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiXo     = ', psiXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiYo     = ', psiYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiZo     = ', psiZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiWx     = ', psiWx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiWy     = ', psiWy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiWz     = ', psiWz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiWr     = ', psiWr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RANGE CHECKING INPUT PARAMETERS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     BROADCASTING INPUT PARAMETERS TO ALL MPI PROCESSES ... '

      END IF

      CALL mpi_select_default_kinds ( )
      CALL mpi_bcast_inputs ( MPI_MASTER , mpiInt , mpiReal , mpiCmplx , mpiError )

      IF ( itpOn .EQV. .TRUE. ) THEN ! run simulation using imaginary time propagation

         zDt = CMPLX ( 0.0 , -dT )

      ELSE ! run simulation normally

         zDt = CMPLX ( dT , 0.0 )

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

      ALLOCATE ( X    ( nXa - nXbc : nXb + nXbc ) )
      ALLOCATE ( Y    ( nYa - nYbc : nYb + nYbc ) )
      ALLOCATE ( Z    ( nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Vex  ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K1   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K2   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K3   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K4   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( PsiA ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( PsiB ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

      IF ( ( mpiRank == MPI_MASTER ) .AND. ( fullIO .EQV. .TRUE. ) ) THEN

         ALLOCATE ( Xm   ( nX           ) )
         ALLOCATE ( Ym   ( nY           ) )
         ALLOCATE ( Zm   ( nZ           ) )
         ALLOCATE ( VexM ( nX * nY * nZ ) )
         ALLOCATE ( PsiM ( nX * nY * nZ ) ) 

      END IF

      X    = 0.0
      Y    = 0.0
      Z    = 0.0
      Vex  = 0.0
      K1   = CMPLX ( 0.0 , 0.0 )
      K2   = CMPLX ( 0.0 , 0.0 )
      K3   = CMPLX ( 0.0 , 0.0 )
      K4   = CMPLX ( 0.0 , 0.0 )
      PsiA = CMPLX ( 0.0 , 0.0 )
      PsiB = CMPLX ( 0.0 , 0.0 )
      
      IF ( ( mpiRank == MPI_MASTER ) .AND. ( fullIO .EQV. .TRUE. ) ) THEN 

         Xm   = 0.0
         Ym   = 0.0
         Zm   = 0.0
         VexM = 0.0
         PsiM = CMPLX ( 0.0 , 0.0 )

      END IF

      CALL grid_regular_axis ( nX , nXa , nXb , nXbc , xO , dX , X ) 
      CALL grid_regular_axis ( nY , nYa , nYb , nYbc , yO , dY , Y ) 
      CALL grid_regular_axis ( nZ , nZa , nZb , nZbc , zO , dZ , Z ) 

      IF ( ( mpiRank == MPI_MASTER ) .AND. ( fullIO .EQV. .TRUE. ) ) THEN

         CALL grid_regular_axis ( nX , 1 , nX , 0 , xO , dX , Xm )
         CALL grid_regular_axis ( nY , 1 , nY , 0 , yO , dY , Ym )
         CALL grid_regular_axis ( nZ , 1 , nZ , 0 , zO , dZ , Zm )

      END IF

      IF ( vexRead .EQV. .TRUE. ) THEN ! read initial potential from file
     
         CALL vex_read_init ( )

      ELSE ! compute initial potential 

         CALL vex_compute_init ( X , Y , Z , Vex )

      END IF

      IF ( readPsi .EQV. .TRUE. ) THEN ! read initial wave function from file

         CALL psi_read_init ( ) 

      ELSE ! compute initial wave function

         CALL psi_compute_init ( X , Y , Z , PsiA )

      END IF

      CALL mpi_exchange_ghosts ( PsiA )

! --- BEGIN MAIN TIME PROPAGATION LOOP ------------------------------------

      tN = t0 ! initialize simulation time

      DO n = 0 , nTsteps

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

         IF ( MODULO ( n , nTwrite ) == 0 ) THEN ! compute partial 
!        expectation values locally on each MPI process; reduce sum of 
!        partial expectation values from all MPI processes on MPI_MASTER;
!        write expectation values, uncertainties and uncertainty 
!        relations to file from MPI_MASTER; write wave function and 
!        external potential to file from MPI_MASTER

            IF ( quadRule == 1 ) THEN ! use rectangle rule
            
               normL2L = l2_norm_3d_rect ( PsiA        )
               avgXL   = x_3d_rect       ( X    , PsiA )
               avgX2L  = x2_3d_rect      ( X    , PsiA ) 
               avgYL   = y_3d_rect       ( Y    , PsiA )
               avgY2L  = y2_3d_rect      ( Y    , PsiA )
               avgZL   = z_3d_rect       ( Z    , PsiA )
               avgZ2L  = z2_3d_rect      ( Z    , PsiA )
               avgVexL = vex_3d_rect     ( Vex  , PsiA )
               avgVmfL = vmf_3d_rect     ( gS   , PsiA )

               IF ( fdOrder == 2 ) THEN ! use 2nd-order CD 

                  avgPxL  = px_3d_rect_cd2  ( PsiA            )
                  avgPx2L = px2_3d_rect_cd2 ( PsiA            )
                  avgPyL  = py_3d_rect_cd2  ( PsiA            )
                  avgPy2L = py2_3d_rect_cd2 ( PsiA            )
                  avgPzL  = pz_3d_rect_cd2  ( PsiA            )
                  avgPz2L = pz2_3d_rect_cd2 ( PsiA            )
                  avgLxL  = lx_3d_rect_cd2  ( Y    , Z , PsiA )
                  avgLx2L = lx2_3d_rect_cd2 ( Y    , Z , PsiA )
                  avgLyL  = ly_3d_rect_cd2  ( X    , Z , PsiA )
                  avgLy2L = ly2_3d_rect_cd2 ( X    , Z , PsiA )
                  avgLzL  = lz_3d_rect_cd2  ( X    , Y , PsiA )
                  avgLz2L = lz2_3d_rect_cd2 ( X    , Y , PsiA )

               ELSE IF ( fdOrder == 4 ) THEN ! use 4th-order CD 

                  avgPxL  = px_3d_rect_cd4  ( PsiA            )
                  avgPx2L = px2_3d_rect_cd4 ( PsiA            )
                  avgPyL  = py_3d_rect_cd4  ( PsiA            )
                  avgPy2L = py2_3d_rect_cd4 ( PsiA            )
                  avgPzL  = pz_3d_rect_cd4  ( PsiA            )
                  avgPz2L = pz2_3d_rect_cd4 ( PsiA            )
                  avgLxL  = lx_3d_rect_cd4  ( Y    , Z , PsiA )
                  !avgLx2L = lx2_3d_rect_cd4 ( Y , Z , Psi3La )
                  avgLyL  = ly_3d_rect_cd4  ( X    , Z , PsiA )
                  !avgLy2L = ly2_3d_rect_cd4 ( X , Z , Psi3La )
                  avgLzL  = lz_3d_rect_cd4  ( X    , Y , PsiA )
                  !avgLz2L = lz2_3d_rect_cd4 ( X , Y , Psi3La )

               ELSE ! fdOrder not supported

               END IF

            ELSE ! quadRule not supported

            END IF

!           Reduce sum of partial expectation values from all MPI 
!           processes on MPI_MASTER

            CALL MPI_REDUCE ( normL2L , normL2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgXL   , avgX   , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgX2L  , avgX2  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgYL   , avgY   , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgY2L  , avgY2  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgZL   , avgZ   , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgZ2L  , avgZ2  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgVexL , avgVex , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgVmfL , avgVmf , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPxL  , avgPx  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPx2L , avgPx2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPyL  , avgPy  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPy2L , avgPy2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPzL  , avgPz  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgPz2L , avgPz2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLxL  , avgLx  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLx2L , avgLx2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLyL  , avgLy  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLy2L , avgLy2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLzL  , avgLz  , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgLz2L , avgLz2 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )

!           Gather partial potentials and wave functions from all MPI processes; will fail for nZ odd - fix in future ... need same count for 

            IF ( mpiRank == MPI_MASTER ) THEN

!              Position, momentum and angular momentum uncertainty

               sigX  = SQRT ( avgX2  - avgX**2  )
               sigY  = SQRT ( avgY2  - avgY**2  )
               sigZ  = SQRT ( avgZ2  - avgZ**2  )
               sigPx = SQRT ( avgPx2 - avgPx**2 )
               sigPy = SQRT ( avgPy2 - avgPy**2 )
               sigPz = SQRT ( avgPz2 - avgPz**2 )
               sigLx = SQRT ( avgLx2 - avgLx**2 )
               sigLy = SQRT ( avgLy2 - avgLy**2 )
               sigLz = SQRT ( avgLz2 - avgLz**2 )

!              Squared angular momentum expectation value

               avgL2 = avgLx2 + avgLy2 + avgLz2

!              Kinetic energy expectation values

               avgTx = 0.5 * avgPx2
               avgTy = 0.5 * avgPy2
               avgTz = 0.5 * avgPz2

!              Energy expectation value

               avgE = avgTx + avgTy + avgTz + avgVex + avgVmf

!              Write expectation values, uncertainties and uncertainty 
!              relations to file from MPI_MASTER

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

!              Write wave function to file from MPI_MASTER

            END IF

            IF ( fullIO .EQV. .TRUE. ) THEN ! gather external potential and wave function from all MPI processes to MPI_MASTER

               CALL mpi_gather_custom ( Vex , PsiA )

               IF ( mpiRank == MPI_MASTER ) THEN ! write external potential and wave function to file from MPI_MASTER

                  fileNumber = fileNumber + 1
                  CALL write_vtk ( 'gpse-vex-psi' , fileNumber , Xm , Ym , Zm , VexM , PsiM )

               END IF

            END IF

         END IF

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!        Compute 1st stage of generalized 4th-order Runge-Kutta (GRK4L): k_1 = f ( t_n , y_n )

         CALL compute_f ( X , Y , Z , Vex , PsiA , K1 )

!        Compute the intermediate wave function for the 2nd stage of the GRK4L method: y_n + 0.5 * dT * k_1
!        Note that the intermediate wave function is only computed on the interior grid points assigned to each MPI process.

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + CMPLX ( 0.5 , 0.0 ) * zDt * K1 ( j , k , l )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

!        Exchange boundary condition information among nearest neighbor processes

         CALL mpi_exchange_ghosts ( PsiB )

!        Compute V ( x , t_n + dT / 2 ), W ( t_n + dT / 2)

         tN = t0 + ( REAL ( n ) + 0.5 ) * dT

!        Compute 2nd stage of GRK4L: k_2 = f ( t_n + 0.5 * dT , y_n + 0.5 * dT * k_1 )

         CALL compute_f ( X , Y , Z , Vex , PsiB , K2 )

!        Compute intermediate wave function for 3rd stage of GRK4L: y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + zDt * ( CMPLX ( 0.5 - 1.0 / REAL( rk4Lambda ) , 0.0 ) * K1 ( j , k , l ) + CMPLX ( 1.0 / REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

         CALL mpi_exchange_ghosts ( PsiB )

!        Compute stage three ... k_3 = f ( t_n + 0.5 * dT , y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2

         CALL compute_f ( X , Y , Z , Vex , PsiB , K3 )

!        Compute intermediate wave function for 4th stage of GRK4L: y_n + ( 1 - lambda / 2 ) * dT * k_2 + ( lambda / 2 ) * dT * k_3

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + zDt * ( CMPLX ( 1.0 - 0.5 * REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( 0.5 * REAL ( rk4Lambda ) , 0.0 ) * K3 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

         CALL mpi_exchange_ghosts ( PsiB )

!        Calculate V ( x , t_n + dT ), W ( t_n + dT )

         tN = t0 + REAL ( n + 1 ) * dT

!        Compute the fourth stage ... k_4 = f ( t_n + dT , y_n + ( 1 - lamda / 2 ) * dT * k_2 + ( lamda / 2) * dT * k_3 

         CALL compute_f ( X , Y , Z , Vex , PsiB , K4 )

!        Compute wave function at nth+1 time step ... y_{ n + 1 } = y_n + ( dT / 6 ) * [ k_1 + ( 4 - lambda ) * k_2 + lambda * k_3 + k_4 ]

!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + CMPLX ( 1.0 / 6.0 , 0.0 ) * zDt * ( K1 ( j , k , l ) + CMPLX ( 4.0 - REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( REAL ( rk4Lambda ) , 0.0 ) * K3 ( j , k , l ) + K4 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL

         CALL mpi_exchange_ghosts ( PsiB )

!        y_{ n + 1 } ---> y_n        
         PsiA = PsiB

      END DO

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

! --- END MAIN TIME PROPAGATION LOOP / BEGIN CLEAN UP TO STOP -------------

      !DEALLOCATE ( Vex3L )
      !DEALLOCATE ( ZL )
      !DEALLOCATE ( YL )
      !DEALLOCATE ( XL ) 

      IF ( mpiRank == MPI_MASTER ) THEN

         CALL DATE_AND_TIME ( stopDate , stopTime , stopZone , StopValues )
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RUN STOPPED @ ', stopTime, ' ON ', stopDate, ' ... '
         DEALLOCATE ( StopValues )
         DEALLOCATE ( StartValues )

      END IF
      
      CALL MPI_FINALIZE ( mpiError )

      DEALLOCATE ( MPIstatus )

! --- FORMAT STATEMENTS ---------------------------------------------------

      STOP

      CONTAINS

         SUBROUTINE mpi_select_default_kinds ( )

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

         SUBROUTINE mpi_bcast_inputs ( mpiMaster , mpiInt , mpiReal , mpiCmplx , mpiError )

            IMPLICIT NONE

            INTEGER, INTENT ( IN    ) :: mpiMaster
            INTEGER, INTENT ( IN    ) :: mpiInt
            INTEGER, INTENT ( IN    ) :: mpiReal
            INTEGER, INTENT ( IN    ) :: mpiCmplx
            INTEGER, INTENT ( INOUT ) :: mpiError

            CALL MPI_BCAST ( itpOn     , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fullIO    , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fmtIO     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( rk4Lambda , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( fdOrder   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( quadRule  , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nTsteps   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nTwrite   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( t0        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dT        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wX        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wY        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( wZ        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( gS        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            CALL MPI_BCAST ( nX        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nXbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nY        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nYbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )    
            CALL MPI_BCAST ( nZ        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( nZbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( xO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( yO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( zO        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dX        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dY        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( dZ        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            CALL MPI_BCAST ( vexRead   , 1 , MPI_LOGICAL , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFmtIn  , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexInit   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexXo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexYo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexZo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexRo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFx     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFy     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexFz     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexWx     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexWy     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexWz     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( vexWr     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            CALL MPI_BCAST ( psiInit   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiNx     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiNy     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiNz     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiNr     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiMl     , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiXo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiYo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiZo     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiWx     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiWy     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiWz     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
            CALL MPI_BCAST ( psiWr     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_gather_custom ( Vex3 , Psi3 )

            IMPLICIT NONE

            INTEGER :: mpiSource
            INTEGER :: nZaSource
            INTEGER :: nZbSource

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

            IF ( mpiRank == MPI_MASTER ) THEN ! receive external potential and wave function data from all MPI processes

!$OMP          PARALLEL DEFAULT ( SHARED )
!$OMP          DO SCHEDULE ( STATIC )
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        VexM ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) = Vex3 ( j , k , l )
                        PsiM ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) = Psi3 ( j , k , l ) 

                     END DO 

                  END DO 

               END DO 
!$OMP          END DO
!$OMP          END PARALLEL

               DO mpiSource = 1 , mpiProcesses - 1

                  CALL MPI_RECV ( nZaSource , 1 , mpiInt , mpiSource , 0 , MPI_COMM_WORLD , MPIStatus , mpiError )
                  CALL MPI_RECV ( nZbSource , 1 , mpiInt , mpiSource , 1 , MPI_COMM_WORLD , MPIStatus , mpiError )
                  
                  DO l = nZaSource , nZbSource

                     DO k = 1 , nY

                        DO j = 1 , nX

                           CALL MPI_RECV ( VexM ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) , 1 , mpiCmplx , mpiSource , 2 , MPI_COMM_WORLD , MPIStatus , mpiError )
                           CALL MPI_RECV ( PsiM ( j + nX * ( ( k - 1 ) + nY * ( l - 1 ) ) ) , 1 , mpiCmplx , mpiSource , 3 , MPI_COMM_WORLD , MPIStatus , mpiError )

                        END DO

                     END DO

                  END DO

               END DO

            ELSE ! send external potential and wave function data to MPI_MASTER

               CALL MPI_SSEND ( nZa , 1 , mpiInt , MPI_MASTER , 0 , MPI_COMM_WORLD , mpiError )
               CALL MPI_SSEND ( nZb , 1 , mpiInt , MPI_MASTER , 1 , MPI_COMM_WORLD , mpiError )

               DO l = nZa , nZb

                  DO k = 1 , nY

                     DO j = 1 , nY
   
                        CALL MPI_SSEND ( Vex3 ( j , k , l ) , 1 , mpiCmplx , MPI_MASTER , 2 , MPI_COMM_WORLD , mpiError )
                        CALL MPI_SSEND ( Psi3 ( j , k , l ) , 1 , mpiCmplx , MPI_MASTER , 3 , MPI_COMM_WORLD , mpiError )

                     END DO

                  END DO

               END DO

            END IF 

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
                        & CMPLX ( 0.0 , 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 + Vex3 ( j , k , l ) + gS * ABS ( Psi3 ( j , k , l ) )**2 ) * Psi3 ( j , k , l ) + &
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

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 1 , MPI_COMM_WORLD , MPIstatus , mpiError )

               END IF

            ELSE 

               CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 0 , MPI_COMM_WORLD , MPIstatus , mpiError )

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 1 , MPI_COMM_WORLD , mpiError )

               END IF

            END IF

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN

               IF ( mpiRank /= MPI_MASTER ) THEN

                  CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 2 , MPI_COMM_WORLD , mpiError )

               END IF

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 3 , MPI_COMM_WORLD , MPIstatus , mpiError )

               END IF

            ELSE

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 2 , MPI_COMM_WORLD , MPIstatus , mpiError )

               END IF
               CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 3 , MPI_COMM_WORLD , mpiError )

            END IF

            RETURN

         END SUBROUTINE

      END PROGRAM

! =========================================================================
