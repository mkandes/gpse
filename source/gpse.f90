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
!     Tuesday, June 17th, 2014
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

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.1.7'

      INTEGER, PARAMETER :: INT_DEFAULT_KIND   = KIND ( 0 ) 
      INTEGER, PARAMETER :: REAL_DEFAULT_KIND  = KIND ( 0.0 )
      INTEGER, PARAMETER :: CMPLX_DEFAULT_KIND = KIND ( CMPLX ( 0.0 , 0.0 ) )
      INTEGER, PARAMETER :: MAX_LEN_FILENAME   = 255
      INTEGER, PARAMETER :: MPI_MASTER         = 0

! --- PARAMETER DEFINITIONS -----------------------------------------------
! --- VARIABLE DECLARATIONS -----------------------------------------------

      LOGICAL :: itpOn = .FALSE.
      LOGICAL :: psiRead  = .FALSE.
      LOGICAL :: psiWrite = .FALSE.
      LOGICAL :: vexRead  = .FALSE.
      LOGICAL :: vexWrite = .FALSE.
      LOGICAL :: vexLin   = .FALSE.
      LOGICAL :: vexSho   = .FALSE.
      LOGICAL :: vexShor  = .FALSE.

      CHARACTER ( LEN = 7  ) :: mainInFile   = 'gpse.in'
      CHARACTER ( LEN = 8  ) :: mainOutFile  = 'gpse.out'
      CHARACTER ( LEN = 6  ) :: psiInFile    = 'psi.in'
      CHARACTER ( LEN = 21 ) :: fIntName     = 'NONE'
      CHARACTER ( LEN = 27 ) :: fRealName    = 'NONE'
      CHARACTER ( LEN = 29 ) :: fCmplxName   = 'NONE'
      CHARACTER ( LEN = 28 ) :: mpiIntName   = 'NONE'
      CHARACTER ( LEN = 34 ) :: mpiRealName  = 'NONE'
      CHARACTER ( LEN = 36 ) :: mpiCmplxName = 'NONE'
      CHARACTER ( LEN = 8  ) :: startDate    = 'NONE'
      CHARACTER ( LEN = 10 ) :: startTime    = 'NONE'
      CHARACTER ( LEN = 5  ) :: startZone    = 'NONE'
      CHARACTER ( LEN = 8  ) :: stopDate     = 'NONE'
      CHARACTER ( LEN = 10 ) :: stopTime     = 'NONE'
      CHARACTER ( LEN = 5  ) :: stopZone     = 'NONE'

      INTEGER :: rk4Lambda      = 0
      INTEGER :: fdOrder        = 0
      INTEGER :: quadRule       = 0
      INTEGER :: fCmplx     = 0 
      INTEGER :: fInt       = 0
      INTEGER :: fReal      = 0
      INTEGER :: mpiCmplx   = 0 
      INTEGER :: mpiDestination = 0 
      INTEGER :: mpiError       = 0 
      INTEGER :: mpiErrorCode   = 0 
      INTEGER :: mpiErrorClass  = 0 
      INTEGER :: mpiInteger     = 0 
      INTEGER :: mpiProcesses   = 0 
      INTEGER :: mpiProvided    = 0 
      INTEGER :: mpiRank        = 0
      INTEGER :: mpiReal    = 0 
      INTEGER :: mpiSource      = 0
      INTEGER :: ompThreads     = 0 
      INTEGER :: ompThreadID    = 0
      INTEGER :: mainInUnit     = 500
      INTEGER :: mainOutUnit    = 600
      INTEGER :: nTsteps        = 0
      INTEGER :: nTwrite        = 0
      INTEGER :: nTgS           = 0 
      INTEGER :: nTomega        = 0 
      INTEGER :: nTpsi          = 0 
      INTEGER :: nTvex          = 0  
      INTEGER :: nX             = 0 
      INTEGER :: nXa            = 0 
      INTEGER :: nXb            = 0 
      INTEGER :: nXbc           = 0 
      INTEGER :: nY             = 0 
      INTEGER :: nYa            = 0 
      INTEGER :: nYb            = 0 
      INTEGER :: nYbc           = 0 
      INTEGER :: nZ             = 0
      INTEGER :: nZa            = 0 
      INTEGER :: nZb            = 0 
      INTEGER :: nZbc           = 0
      INTEGER :: psiInFmt       = 0
      INTEGER :: psiInUnit      = 0
      INTEGER :: psiOutFmt      = 0
      INTEGER :: psiOutUnit     = 0
      INTEGER :: psiInit        = 0
      INTEGER :: psiNx          = 0
      INTEGER :: psiNy          = 0
      INTEGER :: psiNz          = 0
      INTEGER :: psiNr          = 0
      INTEGER :: psiMl          = 0
      INTEGER :: vexInFmt       = 0
      INTEGER :: vexInUnit      = 0
      INTEGER :: vexOutFmt      = 0
      INTEGER :: vexOutUnit     = 0
      INTEGER :: j , k , l , n , m ! Reserved loop counters. 

      REAL :: tN    = 0.0 ! Simulation time at nth time step
      REAL :: t0    = 0.0 ! Time at the beginning of the simulation
      REAL :: xO    = 0.0 ! X-coordinate of origin used to define computational grid
      REAL :: yO    = 0.0 ! Y-coordinate of origin used to define computational grid
      REAL :: zO    = 0.0 ! Z-coordinate of origin used to define computational grid
      REAL :: dT    = 0.0 ! Interval of a time step
      REAL :: dX    = 0.0 !
      REAL :: dY    = 0.0
      REAL :: dZ    = 0.0
      REAL :: wX    = 0.0
      REAL :: wY    = 0.0
      REAL :: wZ    = 0.0
      REAL :: gS    = 0.0 
      REAL :: psiXo = 0.0
      REAL :: psiYo = 0.0
      REAL :: psiZo = 0.0
      REAL :: psiWx = 0.0
      REAL :: psiWy = 0.0
      REAL :: psiWz = 0.0
      REAL :: psiWr = 0.0
      REAL :: vexLinXo = 0.0
      REAL :: vexLinYo = 0.0
      REAL :: vexLinZo = 0.0
      REAL :: vexLinFx = 0.0 
      REAL :: vexLinFy = 0.0 
      REAL :: vexLinFz = 0.0 
      REAL :: vexShoXo = 0.0 
      REAL :: vexShoYo = 0.0 
      REAL :: vexShoZo = 0.0
      REAL :: vexShoWx = 0.0 
      REAL :: vexShoWy = 0.0 
      REAL :: vexShoWz = 0.0
      REAL :: vexShorXo = 0.0 
      REAL :: vexShorYo = 0.0 
      REAL :: vexShorZo = 0.0
      REAL :: vexShorRo = 0.0
      REAL :: vexShorWr = 0.0
      REAL :: vexShorWz = 0.0
      REAL :: normL2 = 0.0
      REAL :: normL20 = 0.0
      REAL :: avgX = 0.0
      REAL :: avgX0 = 0.0
      REAL :: avgY = 0.0
      REAL :: avgY0 = 0.0
      REAL :: avgZ = 0.0
      REAL :: avgZ0 = 0.0
      REAL :: avgPx = 0.0
      REAL :: avgPx0 = 0.0
      REAL :: avgPy = 0.0
      REAL :: avgPy0 = 0.0
      REAL :: avgPz = 0.0
      REAL :: avgPz0 = 0.0
      REAL :: avgLx = 0.0
      REAL :: avgLxSum = 0.0
      REAL :: avgLy = 0.0
      REAL :: avgLySum = 0.0
      REAL :: avgLz = 0.0
      REAL :: avgLzSum = 0.0
      REAL :: avgTx = 0.0
      REAL :: avgTy = 0.0
      REAL :: avgTz = 0.0
      REAL :: avgVex = 0.0
      REAL :: avgVexSum = 0.0
      REAL :: avgVmf = 0.0
      REAL :: avgVmfSum = 0.0
      REAL :: avgE = 0.0
      REAL :: avgX2 = 0.0
      REAL :: avgX20 = 0.0
      REAL :: avgY2 = 0.0
      REAL :: avgY20 = 0.0
      REAL :: avgZ2 = 0.0
      REAL :: avgZ20 = 0.0
      REAL :: avgPx2 = 0.0
      REAL :: avgPx2sum = 0.0
      REAL :: avgPy2 = 0.0
      REAL :: avgPy2sum = 0.0
      REAL :: avgPz2 = 0.0
      REAL :: avgPz20 = 0.0
      REAL :: avgLx2 = 0.0
      REAL :: avgLx20 = 0.0
      REAL :: avgLy2 = 0.0
      REAL :: avgLy20 = 0.0
      REAL :: avgLz2 = 0.0
      REAL :: avgLz20 = 0.0
      REAL :: avgL2 = 0.0
      REAL :: sigX = 0.0
      REAL :: sigY = 0.0
      REAL :: sigZ = 0.0
      REAL :: sigPx = 0.0
      REAL :: sigPy = 0.0
      REAL :: sigPz = 0.0
      REAL :: sigLx = 0.0
      REAL :: sigLy = 0.0
      REAL :: sigLz = 0.0

      COMPLEX :: dTz = CMPLX ( 0.0 , 0.0 )

! --- VARIABLE DEFINITIONS ------------------------------------------------
! --- ARRAY DECLARATIONS --------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIStatus
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartVals
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopVals

      REAL, ALLOCATABLE, DIMENSION ( : , : , : ) :: Vex
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: X
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Y
      REAL, ALLOCATABLE, DIMENSION ( :         ) :: Z

      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: F
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K4
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: PsiA
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: PsiB
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi

! --- ARRAY DEFINITIONS ---------------------------------------------------
! --- FUNCTION AND SUBROUTINE DECLARATIONS --------------------------------

      INTEGER :: OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM

! --- NAMELIST DECLARATIONS -----------------------------------------------

      NAMELIST /gpseIn/ itpOn , rk4Lambda , fdOrder , quadRule , nTsteps , nTwrite , nX , nXbc , nY , nYbc , nZ , nZbc , t0 , &
         & xO , yO , zO , dT , dX , dY , dZ , gS , psiRead , psiInFile , psiInFmt , psiInUnit , psiWrite , &
         & psiOutFmt , psiOutUnit , psiInit , psiNx , psiNy , psiNz , psiNr , psiMl , psiXo , psiYo , psiZo , psiWx , psiWy , & 
         & psiWz , psiWr , vexRead , vexWrite , vexInFmt , vexOutFmt , vexInUnit , vexOutUnit , vexLin , vexLinXo , vexLinYo , & 
         & vexLinZo , vexLinFx , vexLinFy , vexLinFz , vexSho , vexShoXo , vexShoYo , vexShoZo , vexShoWx , vexShoWy , vexShoWz , &
         & vexShor , vexShorXo , vexShorYo , vexShorZo , vexShorRo , vexShorWr , vexShorWz

! --- NAMELIST DEFINITIONS ------------------------------------------------

! --- FUNCTION AND SUBROUTINE DEFINITIONS ---------------------------------
! --- MAIN PROGRAM --------------------------------------------------------

      ALLOCATE ( MPIStatus ( MPI_STATUS_SIZE ) ) 

      CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiProvided , mpiError )
!     CALL mpi_init_thread_errchk
      IF ( mpiError /= MPI_SUCCESS ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'gpse: ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
         CALL MPI_ABORT ( MPI_COMM_WORLD , mpiErrorCode , mpiError )

      END IF
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiProcesses , mpiError )
!     CALL mpi_comm_size_errchk
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRank , mpiError )
!     CALL mpi_comm_rank_errchk
      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!     CALL mpi_barrier_errchk
!$OMP PARALLEL DEFAULT ( SHARED )
!      ompThreads = OMP_GET_NUM_THREADS ( )
!$OMP END PARALLEL

      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '! =========================================================================='
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     GPSE VERSION ', GPSE_VERSION_NUMBER
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        Compiled by ', COMPILER_VERSION ( ) , ' using the options ', COMPILER_OPTIONS ( )
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     AUTHOR(S)'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         Marty Kandes'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         Computational Science Research Center'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         San Diego State University'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         5500 Campanile Drive'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         San Diego, California 92182'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     COPYRIGHT'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'     
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         Copyright (c) 2014 Martin Charles Kandes'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     LAST UPDATED'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         Tuesday, June 17th, 2014'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '! -------------------------------------------------------------------------'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     STARTING GPSE ... '

         ALLOCATE ( StartVals ( 8 ) )
         CALL DATE_AND_TIME ( startDate , startTime , startZone , StartVals )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RUN STARTED @ ', startTime, ' ON ', startDate, ' ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RUNNING ON ', mpiProcesses, ' MPI PROCESSES WITH ', ompThreads , ' OPENMP THREADS PER PROCESS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     CHECKING MACHINE/COMPILER-SPECIFIC DATA TYPE SUPPORT AND CONFIGURATION ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        INTEGER_KINDS SUPPORTED ... ', INTEGER_KINDS
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        REAL_KINDS SUPPORTED ...    ', REAL_KINDS

      END IF

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      SELECT CASE ( INT_DEFAULT_KIND )

         CASE ( INT8 )

            fIntName = 'INTEGER*1'
            mpiIntName = 'MPI_INTEGER1'
            mpiInteger = MPI_INTEGER1

         CASE ( INT16 )

            fIntName = 'INTEGER*2'
            mpiIntName = 'MPI_INTEGER2'
            mpiInteger = MPI_INTEGER2

         CASE ( INT32 )

            fIntName = 'INTEGER ( INTEGER*4 )'
            mpiIntName = 'MPI_INTEGER ( MPI_INTEGER4 )'
            mpiInteger = MPI_INTEGER

         CASE ( INT64 )

            fIntName = 'INTEGER*8'
            mpiIntName = 'MPI_INTEGER8'
            mpiInteger = MPI_INTEGER8

         CASE DEFAULT

            mpiInteger = -1

      END SELECT

      SELECT CASE ( REAL_DEFAULT_KIND )

         CASE ( REAL32 )

            fRealName = 'REAL ( REAL*4 )'
            mpiRealName = 'MPI_REAL ( MPI_REAL4 )'
            mpiReal = MPI_REAL

         CASE ( REAL64 )

            fRealName = 'DOUBLE PRECISION ( REAL*8 )'
            mpiRealName = 'MPI_DOUBLE_PRECISION ( MPI_REAL8 )'
            mpiReal= MPI_DOUBLE_PRECISION

         CASE ( REAL128 )

            fRealName = 'REAL*16'
            mpiRealName = 'MPI_REAL16'
            mpiReal = MPI_REAL16

         CASE DEFAULT

            mpiReal = -1

      END SELECT

      SELECT CASE ( CMPLX_DEFAULT_KIND )

         CASE ( REAL32 )

            fCmplxName = 'COMPLEX ( COMPLEX*8 )'
            mpiCmplxName = 'MPI_COMPLEX ( MPI_COMPLEX8 )'
            mpiCmplx = MPI_COMPLEX

         CASE ( REAL64 )

            fCmplxName = 'DOUBLE COMPLEX ( COMPLEX*16 )'
            mpiCmplxName = 'MPI_DOUBLE_COMPLEX ( MPI_COMPLEX16 )'
            mpiCmplx = MPI_DOUBLE_COMPLEX

         CASE ( REAL128 )

            fCmplxName = 'COMPLEX*32'
            mpiCmplxName = 'MPI_COMPLEX32'
            mpiCmplx = MPI_COMPLEX32

         CASE DEFAULT

            mpiCmplx = -1            

      END SELECT

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        DEFAULT INTEGER KIND ...    ', fIntName
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        DEFAULT REAL KIND ...       ', fRealName
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        DEFAULT COMPLEX KIND ...    ', fCmplxName
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        MPI_INTEGER TYPE ...        ', mpiIntName
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        MPI_REAL TYPE ...           ', mpiRealName
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        MPI_COMPLEX TYPE ...        ', mpiCmplxName
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     READING MAIN INPUT FILE ... ', mainInFile
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     MAIN INPUT FILE UNIT NUMBER ... ', mainInUnit

         OPEN ( UNIT = mainInUnit , FILE = mainInFile , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )
         READ ( UNIT = mainInUnit , NML = gpseIn )
         CLOSE ( UNIT = mainInUnit , STATUS = 'KEEP' )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     INPUT PARAMETERS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        itpOn   = ', itpOn
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        rk4Lambda  = ', rk4Lambda
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        fdOrder    = ', fdOrder
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        quadRule   = ', quadRule
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        nTsteps    = ', nTsteps
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        nTwrite    = ', nTwrite
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        nX         = ', nX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        nY         = ', nY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        nZ         = ', nZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        t0         = ', t0
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        xO         = ', xO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        yO         = ', yO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        zO         = ', zO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dT         = ', dT
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dX         = ', dX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dY         = ', dY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dZ         = ', dZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        gS         = ', gS
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiRead    = ', psiRead
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiInFmt   = ', psiInFmt
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiInUnit  = ', psiInUnit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiWrite   = ', psiWrite
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiOutFmt  = ', psiOutFmt
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiOutUnit = ', psiOutUnit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiInit    = ', psiInit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiNx      = ', psiNx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiNy      = ', psiNy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiNz      = ', psiNz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiNr      = ', psiNr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiMl      = ', psiMl
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiXo      = ', psiXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiYo      = ', psiYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiZo      = ', psiZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiWx      = ', psiWx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiWy      = ', psiWy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiWz      = ', psiWz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        psiWr      = ', psiWr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexRead    = ', vexRead
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexInFmt   = ', vexInFmt
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexInUnit  = ', vexInUnit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexWrite   = ', vexWrite
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexOutFmt  = ', vexOutFmt
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexOutUnit = ', vexOutUnit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLin     = ', vexLin
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLinXo   = ', vexLinXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLinYo   = ', vexLinYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLinZo   = ', vexLinZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLinFx   = ', vexLinFx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLinFy   = ', vexLinFy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexLinFz   = ', vexLinFz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexSho     = ', vexSho
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShoXo   = ', vexShoXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShoYo   = ', vexShoYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShoZo   = ', vexShoZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShoWx   = ', vexShoWx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShoWy   = ', vexShoWy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShoWz   = ', vexShoWz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShor    = ', vexShor
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShorXo  = ', vexShorXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShorYo  = ', vexShorYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShorZo  = ', vexShorZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShorRo  = ', vexShorRo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShorWr  = ', vexShorWr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexShorWz  = ', vexShorWz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RANGE CHECKING INPUT PARAMETERS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     BROADCASTING INPUT PARAMETERS TO ALL MPI PROCESSES ... '

      END IF

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      CALL MPI_BCAST ( itpOn   , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( rk4Lambda  , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( fdOrder    , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( quadRule   , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nTsteps    , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nTwrite    , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nX         , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nXbc       , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nY         , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nYbc       , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )    
      CALL MPI_BCAST ( nZ         , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nZbc       , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( t0         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( xO         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( yO         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( zO         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dT         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dX         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dY         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dZ         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( gS         , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiRead    , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiInFmt   , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiInUnit  , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWrite   , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiOutFmt  , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiOutUnit , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiInit    , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNx      , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNy      , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNz      , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNr      , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiMl      , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiXo      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiYo      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiZo      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWx      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWy      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWz      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWr      , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexRead    , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexInFmt   , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexInUnit  , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexWrite   , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexOutFmt  , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexOutUnit , 1 , mpiInteger  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLin     , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLinXo   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLinYo   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLinZo   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLinFx   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLinFy   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLinFz   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexSho     , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShoXo   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShoYo   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShoZo   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShoWx   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShoWy   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShoWz   , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShor    , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShorXo  , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShorYo  , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShorZo  , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShorWr  , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexShorWz  , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiError )

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

! --- INITIALIZING ... 

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

      ALLOCATE ( Vex  ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( X    ( nXa - nXbc : nXb + nXbc                                                     ) )
      ALLOCATE ( Y    ( nYa - nYbc : nYb + nYbc                                                     ) )
      ALLOCATE ( Z    ( nZa - nZbc : nZb + nZbc                                                     ) )

      ALLOCATE ( F    ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K1   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K2   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K3   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K4   ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( PsiA ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( PsiB ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

!      IF ( mpiRank == MPI_MASTER ) THEN
!
!         ALLOCATE ( Psi ( -nXbc : nX + nXbc , -nYbc : nY + nYbc , -nZbc : nZ + nZbc ) ) 
!
!      END IF

      Vex = 0.0
      X   = 0.0
      Y   = 0.0
      Z   = 0.0

      F    = CMPLX ( 0.0 , 0.0 )
      K1   = CMPLX ( 0.0 , 0.0 )
      K2   = CMPLX ( 0.0 , 0.0 )
      K3   = CMPLX ( 0.0 , 0.0 )
      K4   = CMPLX ( 0.0 , 0.0 )
      PsiA = CMPLX ( 0.0 , 0.0 )
      PsiB = CMPLX ( 0.0 , 0.0 )
      
!      IF ( mpiRank == MPI_MASTER ) THEN 
! 
!         Psi = CMPLX ( 0.0 , 0.0 )
!
!      END IF

      CALL regular_grid_axis ( nX , nXa , nXb , nXbc , xO , dX , X )
      CALL regular_grid_axis ( nY , nYa , nYb , nYbc , yO , dY , Y )
      CALL regular_grid_axis ( nZ , nZa , nZb , nZbc , zO , dZ , Z )

      IF ( psiRead .EQV. .TRUE. ) THEN

         ! Read initial wave function from input file.

      ELSE

         IF ( psiInit == 1 ) THEN

            ! Isotropic SHO not supported yet.

         ELSE IF ( psiInit == 2 ) THEN

            CALL psi_3d_se_sho_ani ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , psiNx , psiNy , psiNz , psiXo , &
               & psiYo , psiZo , psiWx , psiWy , psiWz , X , Y , Z , PsiA )

         ELSE IF ( psiInit == 3 ) THEN

            CALL psi_3d_se_sho_axi ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , psiNr , psiMl , psiNz , psiXo , &
               & psiYo , psiZo , psiWr , psiWz , X , Y , Z , PsiA )

         ELSE

            ! ERROR: psiInit not defined.

         END IF

      END IF

      CALL mpi_ghost_exchange ( PsiA )

      IF ( vexRead .EQV. .TRUE. ) THEN ! initialize external potential from file ...

      END IF

      IF ( vexLin .EQV. .TRUE. ) THEN ! add linear potential to external potential ...

         CALL vex_3d_lin ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , vexLinXo , vexLinYo , vexLinZo , vexLinFx , vexLinFy , vexLinFz , X , Y , Z , Vex )

      END IF

      IF ( vexSho .EQV. .TRUE. ) THEN ! add simple harmonic oscillator potential to external potential ...

         CALL vex_3d_sho ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , vexShoXo , vexShoYo , vexShoZo , vexShoWx , vexShoWy , vexShoWz , X , Y , Z , Vex )

      END IF

      IF ( vexShor .EQV. .TRUE. ) THEN ! add simple harmonic oscillator ring potential to external potential ...

         CALL vex_3d_shor ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , vexShorXo , vexShorYo , vexShorZo , vexShorRo , vexShorWr , vexShorWz , X , Y , Z , Vex )

      END IF

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

! --- BEGIN MAIN TIME PROPAGATION LOOP ------------------------------------

      tN = t0 ! Initialize simulation time. 
!     Read initial wave function from file OR calculate initial wave function at tN = t0 using a subroutine available in psi module.
!     Read initial external potential from file OR calculate inital external potential at tN = t0 using a subroutine available in vex module.
!     Calculate initial values of angular velocity vector at tN = t0. 

      DO n = 0 , nTsteps

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

         IF ( MODULO ( n , nTwrite ) == 0 ) THEN

!           Calculate expectation values.

            IF ( quadRule == 1 ) THEN 
            
               normL2 = l2_norm_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , PsiA )
               avgX = x_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , PsiA )
               avgX2 = x2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , PsiA ) 
               avgY = y_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , PsiA )
               avgY2 = y2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , PsiA )
               avgZ = z_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Z , PsiA )
               avgZ2 = z2_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Z , PsiA )
               avgVex = vex_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Vex , PsiA )
               avgVmf = vmf_3d_rect ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , gS , PsiA )

!               IF ( fdOrder == 2 ) THEN

                  avgPx = px_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dY , dZ , PsiA )
                  avgPx2 = px2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , PsiA )
                  avgPy = py_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dZ , PsiA )
                  avgPy2 = py2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , PsiA )
                  avgPz = pz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , PsiA )
                  avgPz2 = pz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , PsiA )
                  avgLx = lx_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Z , PsiA )
                  avgLx2 = lx2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , Y , Z , PsiA )
                  avgLy = ly_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Z , PsiA )
                  avgLy2 = ly2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Z , PsiA )
                  avgLz = lz_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Y , PsiA )
                  avgLz2 = lz2_3d_rect_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , X , Y , PsiA )

!               ELSE IF ( fdOrder == 4 ) THEN
!               ELSE
!               END IF

            ELSE ! quadRule is not supported. 

            END IF

!           ... 

            CALL MPI_REDUCE ( normL2 , normL20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
            CALL MPI_REDUCE ( avgX , avgX0 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgX2 , avgX20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgY , avgY0 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgY2 , avgY20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgZ , avgZ0 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgZ2 , avgZ20,  1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgVex , avgVexSum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgVmf , avgVmfSum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgPx , avgPx0 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgPx2 , avgPx2sum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgPy , avgPy0 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgPy2 , avgPy2sum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgPz , avgPz0 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgPz2 , avgPz20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgLx , avgLxSum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgLx2 , avgLx20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgLy , avgLySum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgLy2 , avgLy20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgLz , avgLzSum , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )
            CALL MPI_REDUCE ( avgLz2 , avgLz20 , 1 , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD )

            IF ( mpiRank == MPI_MASTER ) THEN

!              Uncertainties in position, momentum and angular momentum.

               sigX  = SQRT ( avgX20 - avgX0**2 )
               sigY  = SQRT ( avgY20 - avgY0**2 )
               sigZ  = SQRT ( avgZ20 - avgZ0**2 )
               sigPx = SQRT ( avgPx2sum - avgPx0**2 )
               sigPy = SQRT ( avgPy2sum - avgPy0**2 )
               sigPz = SQRT ( avgPz20 - avgPz0**2 )
               sigLx = SQRT ( avgLx20 - avgLxSum**2 )
               sigLy = SQRT ( avgLy20 - avgLySum**2 )
               sigLz = SQRT ( avgLz20 - avgLzSum**2 )

!              Expectation value of the square of the total angular momentum. 

               avgL2 = avgLx20 + avgLy20 + avgLz20

!              Expectation values of kinetic energy.

               avgTx = 0.5 * avgPx2sum
               avgTy = 0.5 * avgPy2sum
               avgTz = 0.5 * avgPz20

!              Expectation value of the total energy.

               avgE = avgTx + avgTy + avgTz + avgVexSum + avgVmfSum

!              Write expectation values, uncertainties and uncertainty relations to file from MPI_MASTER.

               WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) tN , normL20 , avgX0 , avgX20 , sigX , avgY0 , avgY20 , sigY , avgZ0 , avgZ20 , sigZ , avgPx0 , avgPx2sum , sigPx , avgPy0 , avgPy2sum , sigPy , avgPz0 , avgPz20 , sigPz , avgLxSum , avgLx20 , sigLx , avgLySum , avgLy20 , sigLy , avgLzSum , avgLz20 , sigLz , avgL2 , avgTx , avgTy , avgTz , avgVexSum , avgVmfSum , avgE , sigX * sigPx , sigY * sigPy , sigZ * sigPz , sigLx * sigLy , sigLy * sigLz , sigLz * sigLx

            END IF
            ! Write wave function to file from MPI_MASTER.

         END IF

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!     -- GRK4 ALGORITHM STEP #1 -------------------------------------------
!
!        k_1 = f ( t_n , y_n )
!
!        Compute first stage of generalized 4th-order Runge-Kutta (GRK4) 
!        method using either 2nd- or 4th-order central differences.
!
         CALL compute_f ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , PsiA , K1 )

!     -- GRK4 ALGORITHM STEP #2 -------------------------------------------
!
!        y_n + 0.5 * dT * k_1
!
!        Compute the intermediate wave function for the second stage of the 
!        GRK4 method. Note that the intermediate wave function is only 
!        computed on the interior grid points assigned to each MPI process.
!
!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + CMPLX ( 0.5 , 0.0 ) * dTz * K1 ( j , k , l )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL
!
!     -- GRK4 ALGORITHM STEP #3 -------------------------------------------
!
!        Send boundary condition information to nearest neighbor processes.
!

         CALL mpi_ghost_exchange ( PsiB )

!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 0 , MPI_COMM_WORLD, mpiError )
!
!            END IF
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 1 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 0 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 1 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!
!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 2 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 3 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 2 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!            CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 3 , MPI_COMM_WORLD , mpiError )
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!     -- GRK4 ALGORITHM STEP #4 -------------------------------------------
!
!        Calculate V ( x , t_n + dT / 2 ), W ( t_n + dT / 2)

         tN = t0 + ( REAL ( n ) + 0.5 ) * dT

!     -- GRK4 ALGORITHM STEP #5 -------------------------------------------
!
!        k_2 = f ( t_n + 0.5 * dT , y_n + 0.5 * dT * k_1 )
!
!        Calculate second stage of GRK4 ...

         CALL compute_f ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , PsiB , K2 )

!     -- GRK4 ALGORITHM STEP #6 -------------------------------------------
!
!        y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( lambda / 2 ) * dT * k_2
!
!        Compute the intermediate wave function for the third stage of the
!        GRK4 method. Again, note that the wave function is computed only
!        on the interior grid points assigned to each MPI processes.
!
!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + dTz * ( CMPLX ( 0.5 - 1.0 / REAL( rk4Lambda ) , 0.0 ) * K1 ( j , k , l ) + CMPLX ( 1.0 / REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL
!
!     -- GRK4 ALGORITHM STEP #7 -------------------------------------------
!
!        Distribute boundary conditions ... again ...
!

         CALL mpi_ghost_exchange ( PsiB )

!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 4 , MPI_COMM_WORLD, mpiError )
!
!            END IF
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 5 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 4 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 5 , MPI_COMM_WORLD , mpiError !)
!
!            END IF
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!
!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= MPI_MASTER ) !THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 6 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!            IF ( mpiRank /= mpiProcesses  - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 7 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 6 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!            CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 7 , MPI_COMM_WORLD , mpiError )
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!!
!     -- GRK4 ALGORITHM STEP #8 -------------------------------------------
!
!        k_3 = f ( t_n + 0.5 * dT , y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2
!
!        Compute stage three ...

         CALL compute_f ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , PsiB , K3 )

!     -- GRK4 ALGORITHM STEP #9 ------------------------------------------
!
!        y_n + ( 1 - lambda / 2 ) * dT * k_2 + ( lambda / 2 ) * dT * k_3
!
!        Compute intermediate wave function for the fourth stage of the
!        GRK4 method ...
!
!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + dTz * ( CMPLX ( 1.0 - 0.5 * REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( 0.5 * REAL ( rk4Lambda ) , 0.0 ) * K3 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL
!
!     -- GRK4 ALGORITHM STEP #10 ------------------------------------------
!
!        Distribute boundary conditions ...
!

         CALL mpi_ghost_exchange ( PsiB )

!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 8 , MPI_COMM_WORLD, mpiError )
!
!            END IF
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 9 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 8 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 9 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!
!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 10 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 11 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 10 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!            CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 11 , MPI_COMM_WORLD , mpiError )
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!!
!     -- GRK4 ALGORITHM STEP #11 ------------------------------------------
!
!        Calculate V ( x , t_n + dT ), W ( t_n + dT )

         tN = t0 + REAL ( n + 1 ) * dT

!     -- GRK4 ALGORITHM STEP #12 ------------------------------------------
!
!        k_4 = f ( t_n + dT , y_n + ( 1 - lamda / 2 ) * dT * k_2 + ( lamda / 2) * dT * k_3 
!
!        Compute the fourth stage ...

         CALL compute_f ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , PsiB , K4 )

!     -- GRK4 ALGORITHM STEP #13 ------------------------------------------
!
!        y_{ n + 1 } = y_n + ( dT / 6 ) * [ k_1 + ( 4 - lambda ) * k_2 + lambda * k_3 + k_4 ]
!
!        Compute wave function at nth+1 time step ...
!
!$OMP    PARALLEL DEFAULT ( SHARED )
!$OMP    DO SCHEDULE ( STATIC )
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  PsiB ( j , k , l ) = PsiA ( j , k , l ) + CMPLX ( 1.0 / 6.0 , 0.0 ) * dTz * ( K1 ( j , k , l ) + CMPLX ( 4.0 - REAL ( rk4Lambda ) , 0.0 ) * K2 ( j , k , l ) + CMPLX ( REAL ( rk4Lambda ) , 0.0 ) * K3 ( j , k , l ) + K4 ( j , k , l ) )

               END DO

            END DO

         END DO
!$OMP    END DO
!$OMP    END PARALLEL
!
!     -- GRK4 ALGORITHM STEP #14 ------------------------------------------
!
!        Distribute boundary conditions ...
!

         CALL mpi_ghost_exchange ( PsiB )

!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 12 , MPI_COMM_WORLD, mpiError )
!!
!            END IF
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 13 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            CALL MPI_RECV ( PsiB ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 12 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 13 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!         END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!
!         IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN
!
!            IF ( mpiRank /= MPI_MASTER ) THEN
!
!               CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 14 , MPI_COMM_WORLD , mpiError )
!
!            END IF
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 15 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!         ELSE
!
!            IF ( mpiRank /= mpiProcesses - 1 ) THEN
!
!               CALL MPI_RECV ( PsiB ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 14 , MPI_COMM_WORLD , MPIStatus , mpiError )
!
!            END IF
!
!            CALL MPI_SEND ( PsiB ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 15 , MPI_COMM_WORLD , mpiError )
!
!        END IF
!
!         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )
!!
!     -- GRK4 ALGORITHM STEP #15 ------------------------------------------
!
!        y_{ n + 1 } ---> y_n        
!
         PsiA = PsiB
!
!     ---------------------------------------------------------------------
         
         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )        

      END DO

! --- END MAIN TIME PROPAGATION LOOP / BEGIN CLEAN UP TO STOP -------------

      DEALLOCATE ( Vex )
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

         SUBROUTINE compute_f ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: fdOrder
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
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: gS

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            IF ( fdOrder == 2 ) THEN 

               CALL f_gp_3d_rrf_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

            ELSE IF ( fdOrder == 4 ) THEN

               CALL f_gp_3d_rrf_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

            ELSE IF ( fdOrder == 6 ) THEN

               CALL f_gp_3d_rrf_cd6 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

            ELSE IF ( fdOrder == 8 ) THEN

               CALL f_gp_3d_rrf_cd8 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

            ELSE

               ! fdOrder not supported ...

            END IF

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd2 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F ) 

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
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: gS

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi 
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb 

               DO k = nYa , nYb 

                  DO j = nXa , nXb 

                     F ( j , k , l ) = &
                        & CMPLX ( 0.5 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 0.5 / dZ**2 ) * Psi ( j , k , l - 1 ) + & 
                        & CMPLX ( 0.5 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 0.5 / dY**2 ) * Psi ( j , k - 1 , l ) + & 
                        & CMPLX ( 0.5 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 0.5 / dX**2 ) * Psi ( j - 1 , k , l ) - & 
                        & CMPLX ( 0.0 , 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 + Vex ( j , k , l ) + & 
                        &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + & 
                        & CMPLX ( 0.5 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 0.5 / dX**2 ) * Psi ( j + 1 , k , l ) + & 
                        & CMPLX ( 0.5 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 0.5 / dY**2 ) * Psi ( j , k + 1 , l ) + & 
                        & CMPLX ( 0.5 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 0.5 / dZ**2 ) * Psi ( j , k , l + 1 ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd4 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F ) 

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
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: gS

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex 

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi 
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb 

               DO k = nYa , nYb 

                  DO j = nXa , nXb 

                     F ( j , k , l ) = & 
                        & CMPLX ( ( wX * Y ( k ) - wY * X ( j ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi ( j , k , l - 2 ) + & 
                        & CMPLX ( 0.75 * ( wY * X ( j ) - wX * Y ( k ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi ( j , k , l - 1 ) + & 
                        & CMPLX ( ( wZ * X ( j ) - wX * Z ( l ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi ( j , k - 2 , l ) + & 
                        & CMPLX ( 0.75 * ( wX * Z ( l ) - wZ * X ( j ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi ( j , k - 1 , l ) + & 
                        & CMPLX ( ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi ( j - 2 , k , l ) + & 
                        & CMPLX ( 0.75 * ( wZ * Y ( k ) - wY * Z ( l ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi ( j - 1 , k , l ) + & 
                        & CMPLX ( 0.0 , -1.25 * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) - Vex ( j , k , l ) - & 
                        &    gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + & 
                        & CMPLX ( 0.75 * ( wY * Z ( l ) - wZ * Y ( k ) ) / dX , 2.0 / ( 3.0 * dX**2 ) ) * Psi ( j + 1 , k , l ) + & 
                        & CMPLX ( ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 12.0 * dX ) , -1.0 / ( 24.0 * dX**2 ) ) * Psi ( j + 2 , k , l ) + & 
                        & CMPLX ( 0.75 * ( wZ * X ( j ) - wX * Z ( l ) ) / dY , 2.0 / ( 3.0 * dY**2 ) ) * Psi ( j , k + 1 , l ) + & 
                        & CMPLX ( ( wX * Z ( l ) - wZ * X ( j ) ) / ( 12.0 * dY ) , -1.0 / ( 24.0 * dY**2 ) ) * Psi ( j , k + 2 , l ) + & 
                        & CMPLX ( 0.75 * ( wX * Y ( k ) - wY * X ( j ) ) / dZ , 2.0 / ( 3.0 * dZ**2 ) ) * Psi ( j , k , l + 1 ) + & 
                        & CMPLX ( ( wY * X ( j ) - wX * Y ( k ) ) / ( 12.0 * dZ ) , -1.0 / ( 24.0 * dZ**2 ) ) * Psi ( j , k , l + 2 ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd6 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

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
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: gS

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     F ( j , k , l ) = &
                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi ( j , k , l - 3 ) + &
                        & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * Psi ( j , k , l - 2 ) + &
                        & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * Psi ( j , k , l - 1 ) + &
                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi ( j , k - 3 , l ) + &
                        & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * Psi ( j , k - 2 , l ) + &
                        & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * Psi ( j , k - 1 , l ) + &
                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi ( j - 3 , k , l ) + &
                        & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * Psi ( j - 2 , k , l ) + &
                        & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * Psi ( j - 1 , k , l ) - &
                        & CMPLX ( 0.0 , ( 49.0 / 36.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex ( j , k , l ) + gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                        & CMPLX ( 3.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 4.0  * dX ) ,  3.0 / ( 4.0 * dX**2 ) ) * Psi ( j + 1 , k , l ) + &
                        & CMPLX ( 3.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 20.0 * dX ) , -3.0 / ( 40.0 * dX**2 ) ) * Psi ( j + 2 , k , l ) + &
                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 60.0 * dX ) ,  1.0 / ( 180.0 * dX**2 ) ) * Psi ( j + 3 , k , l ) + &
                        & CMPLX ( 3.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 4.0  * dY ) ,  3.0 / ( 4.0 * dY**2 ) ) * Psi ( j , k + 1 , l ) + &
                        & CMPLX ( 3.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 20.0 * dY ) , -3.0 / ( 40.0 * dY**2 ) ) * Psi ( j , k + 2 , l ) + &
                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 60.0 * dY ) ,  1.0 / ( 180.0 * dY**2 ) ) * Psi ( j , k + 3 , l ) + &
                        & CMPLX ( 3.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 4.0  * dZ ) ,  3.0 / ( 4.0 * dZ**2 ) ) * Psi ( j , k , l + 1 ) + &
                        & CMPLX ( 3.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 20.0 * dZ ) , -3.0 / ( 40.0 * dZ**2 ) ) * Psi ( j , k , l + 2 ) + &
                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 60.0 * dZ ) ,  1.0 / ( 180.0 * dZ**2 ) ) * Psi ( j , k , l + 3 )   

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE f_gp_3d_rrf_cd8 ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , dY , dZ , wX , wY , wZ , gS , X , Y , Z , Vex , Psi , F )

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
            REAL, INTENT ( IN ) :: wX
            REAL, INTENT ( IN ) :: wY
            REAL, INTENT ( IN ) :: wZ
            REAL, INTENT ( IN ) :: gS

            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
            REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
            REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z
            REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi
            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( OUT ) :: F

            INTEGER :: j , k , l

!$OMP       PARALLEL DEFAULT ( SHARED )
!$OMP       DO SCHEDULE ( STATIC )
            DO l = nZa , nZb

               DO k = nYa , nYb

                  DO j = nXa , nXb

                     F ( j , k , l ) = &
                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi ( j , k , l - 4 ) + &
                        & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 105.0 * dZ ) ,  4.0 / ( 315.0 * dZ**2 ) ) * Psi ( j , k , l - 3 ) + &
                        & CMPLX (       ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dZ ) , -1.0 / ( 10.0 * dZ**2 ) ) * Psi ( j , k , l - 2 ) + &
                        & CMPLX ( 4.0 * ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dZ ) ,  4.0 / ( 5.0 * dZ**2 ) ) * Psi ( j , k , l - 1 ) + &
                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 280.0 * dY ) , -1.0 / ( 1120 * dY**2 ) ) * Psi ( j , k - 4 , l ) + &
                        & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi ( j , k - 3 , l ) + &
                        & CMPLX (       ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi ( j , k - 2 , l ) + &
                        & CMPLX ( 4.0 * ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * Psi ( j , k - 1 , l ) + &
                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * Psi ( j - 4 , k , l ) + &
                        & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * Psi ( j - 3 , k , l ) + &
                        & CMPLX (       ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi ( j - 2 , k , l ) + &
                        & CMPLX ( 4.0 * ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi ( j - 1 , k , l ) - &
                        & CMPLX ( 0.0 , ( 205.0 / 144.0 ) * ( 1.0 / dX**2 + 1.0 / dY**2 + 1.0 / dZ**2 ) + Vex ( j , k , l ) + gS * ABS ( Psi ( j , k , l ) )**2 ) * Psi ( j , k , l ) + &
                        & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi ( j + 1 , k , l ) + &
                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dX**2 ) ) * Psi ( j + 2 , k , l ) + &
                        & CMPLX ( 4.0 * ( wY * Z ( l ) - wZ * Y ( k ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dX**2 ) ) * Psi ( j + 3 , k , l ) + &
                        & CMPLX (       ( wZ * Y ( k ) - wY * Z ( l ) ) / ( 280.0 * dX ) , -1.0 / ( 1120.0 * dX**2 ) ) * Psi ( j + 4 , k , l ) + &
                        & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 5.0   * dY ) ,  4.0 / ( 5.0 * dY**2 ) ) * Psi ( j , k + 1 , l ) + &
                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 5.0   * dY ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi ( j , k + 2 , l ) + & 
                        & CMPLX ( 4.0 * ( wZ * X ( j ) - wX * Z ( l ) ) / ( 105.0 * dY ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi ( j , k + 3 , l ) + &
                        & CMPLX (       ( wX * Z ( l ) - wZ * X ( j ) ) / ( 280.0 * dY ) , -1.0 / ( 1120.0 * dY**2 ) ) * Psi ( j , k + 4 , l ) + &
                        & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 5.0   * dX ) ,  4.0 / ( 5.0 * dX**2 ) ) * Psi ( j , k , l + 1 ) + &
                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 5.0   * dX ) , -1.0 / ( 10.0 * dY**2 ) ) * Psi ( j , k , l + 2 ) + &
                        & CMPLX ( 4.0 * ( wX * Y ( k ) - wY * X ( j ) ) / ( 105.0 * dX ) ,  4.0 / ( 315.0 * dY**2 ) ) * Psi ( j , k , l + 3 ) + &
                        & CMPLX (       ( wY * X ( j ) - wX * Y ( k ) ) / ( 280.0 * dZ ) , -1.0 / ( 1120.0 * dZ**2 ) ) * Psi ( j , k , l + 4 ) 

                  END DO

               END DO

            END DO
!$OMP       END DO
!$OMP       END PARALLEL

            RETURN

         END SUBROUTINE

         SUBROUTINE mpi_ghost_exchange ( Psi ) 

            IMPLICIT NONE

            COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN 

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 0 , MPI_COMM_WORLD, mpiError )

               END IF

               IF ( mpiRank /= MPI_MASTER ) THEN 

                  CALL MPI_RECV ( Psi ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 1 , MPI_COMM_WORLD , MPIStatus , mpiError )

               END IF

            ELSE 

               CALL MPI_RECV ( Psi ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 0 , MPI_COMM_WORLD , MPIStatus , mpiError )

               IF ( mpiRank /= mpiProcesses - 1 ) THEN 

                  CALL MPI_SEND ( Psi ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 1 , MPI_COMM_WORLD , mpiError )

               END IF

            END IF

            IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN

               IF ( mpiRank /= MPI_MASTER ) THEN

                  CALL MPI_SEND ( Psi ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 2 , MPI_COMM_WORLD , mpiError )

               END IF

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 3 , MPI_COMM_WORLD , MPIStatus , mpiError )

               END IF

            ELSE

               IF ( mpiRank /= mpiProcesses - 1 ) THEN

                  CALL MPI_RECV ( Psi ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 2 , MPI_COMM_WORLD , MPIStatus , mpiError )

               END IF
               CALL MPI_SEND ( Psi ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 3 , MPI_COMM_WORLD , mpiError )

            END IF

            RETURN

         END SUBROUTINE

      END PROGRAM

! =========================================================================
