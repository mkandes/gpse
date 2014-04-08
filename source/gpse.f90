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
!     Tuesday, April 8th, 2014
!
! -------------------------------------------------------------------------

      PROGRAM GPSE

! --- MODULE DECLARATIONS -------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: GRID
      USE            :: IO
      USE            :: MATH
      USE            :: PSI
      USE            :: VEX

! --- MODULE DEFINITIONS --------------------------------------------------
! -------------------------------------------------------------------------
      
      IMPLICIT NONE

      INCLUDE 'mpif.h'

! --- PARAMETER DECLARATIONS  ---------------------------------------------

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.1.1'

      INTEGER, PARAMETER :: INT_DEFAULT_KIND   = KIND ( 0 ) 
      INTEGER, PARAMETER :: REAL_DEFAULT_KIND  = KIND ( 0.0 )
      INTEGER, PARAMETER :: CMPLX_DEFAULT_KIND = KIND ( CMPLX ( 0.0 , 0.0 ) )
      INTEGER, PARAMETER :: MAX_LEN_FILENAME   = 255
      INTEGER, PARAMETER :: MPI_MASTER         = 0

! --- PARAMETER DEFINITIONS -----------------------------------------------
! --- VARIABLE DECLARATIONS -----------------------------------------------

      LOGICAL :: psiRead  = .FALSE.
      LOGICAL :: psiWrite = .FALSE.
      LOGICAL :: vexRead  = .FALSE.
      LOGICAL :: vexWrite = .FALSE.
      LOGICAL :: vexLin   = .FALSE.
      LOGICAL :: vexSHO   = .FALSE.
      LOGICAL :: vexSHOR  = .FALSE.

      CHARACTER ( LEN = 7  ) :: mainInFile   = 'gpse.in'
      CHARACTER ( LEN = 8  ) :: mainOutFile  = 'gpse.out'
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

      INTEGER :: runMode        = 0
      INTEGER :: dimEq          = 0
      INTEGER :: odeSolve       = 0
      INTEGER :: rk4Lambda      = 0
      INTEGER :: fdOrder        = 0
      INTEGER :: quadRule       = 0
      INTEGER :: fCmplxKind     = 0 
      INTEGER :: fIntKind       = 0
      INTEGER :: fRealKind      = 0
      INTEGER :: mpiCmplxKind   = 0 
      INTEGER :: mpiDestination = 0 
      INTEGER :: mpiError       = 0 
      INTEGER :: mpiErrorCode   = 0 
      INTEGER :: mpiErrorClass  = 0 
      INTEGER :: mpiIntKind     = 0 
      INTEGER :: mpiProcesses   = 0 
      INTEGER :: mpiProvided    = 0 
      INTEGER :: mpiRank        = 0
      INTEGER :: mpiRealKind    = 0 
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

      REAL :: tN    = 0.0 ! Time of simulation (at nth time step)
      REAL :: t0    = 0.0 ! Time at start of simulation
      REAL :: xO    = 0.0 ! X-coordinate of origin used to define computational grid
      REAL :: yO    = 0.0 ! Y-coordinate of origin used to define computational grid
      REAL :: zO    = 0.0 ! Z-coordinate of origin used to define computational grid
      REAL :: dTRe  = 0.0 ! Interval of a real time step
      REAL :: dTIm  = 0.0 ! Interval of an imaginary time step
      REAL :: dX    = 0.0 !
      REAL :: dY    = 0.0
      REAL :: dZ    = 0.0
      REAL :: gSRe  = 0.0 
      REAL :: gSIm  = 0.0
      REAL :: psiXo = 0.0
      REAL :: psiYo = 0.0
      REAL :: psiZo = 0.0
      REAL :: psiWx = 0.0
      REAL :: psiWy = 0.0
      REAL :: psiWz = 0.0
      REAL :: psiWr = 0.0
      REAL :: vexXo = 0.0
      REAL :: vexYo = 0.0
      REAL :: vexZo = 0.0
      REAL :: vexRo = 0.0
      REAL :: vexFx = 0.0
      REAL :: vexFy = 0.0
      REAL :: vexFz = 0.0
      REAL :: vexWx = 0.0
      REAL :: vexWy = 0.0
      REAL :: vexWz = 0.0
      REAL :: vexWr = 0.0

      COMPLEX :: dT = CMPLX ( 0.0 , 0.0 )
      COMPLEX :: gS = CMPLX ( 0.0 , 0.0 )

! --- VARIABLE DEFINITIONS ------------------------------------------------
! --- ARRAY DECLARATIONS --------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIStatus
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: Pvt

      REAL, ALLOCATABLE, DIMENSION ( :             ) :: ImPsi1
      REAL, ALLOCATABLE, DIMENSION ( : , :         ) :: ImPsi2
      REAL, ALLOCATABLE, DIMENSION ( : , : , :     ) :: ImPsi3
      REAL, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: ImPsi4
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: RePsi1
      REAL, ALLOCATABLE, DIMENSION ( : , :         ) :: RePsi2
      REAL, ALLOCATABLE, DIMENSION ( : , : , :     ) :: RePsi3
      REAL, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: RePsi4
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: Wx
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: Wy
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: Wz
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: X
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: Y
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: Z

      COMPLEX, ALLOCATABLE, DIMENSION ( :             ) :: Gs1   ! Note: Same as below.  
      COMPLEX, ALLOCATABLE, DIMENSION ( : , :         ) :: Gs2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , :     ) :: Gs3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: Gs4
      COMPLEX, ALLOCATABLE, DIMENSION ( :             ) :: Psi1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , :         ) :: Psi2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , :     ) :: Psi3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: Psi4
      COMPLEX, ALLOCATABLE, DIMENSION ( :             ) :: Vex1   ! Note: Arrays reserved for storing external potential functions 
      COMPLEX, ALLOCATABLE, DIMENSION ( : , :         ) :: Vex2   ! are declared as COMPLEX to provide the built-in flexability to  
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , :     ) :: Vex3   ! accomodate complex potentials. 
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: Vex4

! --- ARRAY DEFINITIONS ---------------------------------------------------
! --- FUNCTION AND SUBROUTINE DECLARATIONS --------------------------------

      INTEGER :: OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM

! --- NAMELIST DECLARATIONS -----------------------------------------------

      NAMELIST /gpseIn/ runMode , dimEq , odeSolve , rk4Lambda , fdOrder , quadRule , nTsteps , nTwrite , nX , nXbc , nY , nYbc , &
         & nZ , nZbc , t0 , xO , yO , zO , dTRe , dTIm , dX , dY , dZ , gSRe , gSIm , psiRead , psiInFmt , psiInUnit , psiWrite , &
         & psiOutFmt , psiOutUnit , psiInit , psiNx , psiNy , psiNz , psiNr , psiMl , psiXo , psiYo , psiZo , psiWx , psiWy , &
         & psiWz , psiWr , vexRead , vexWrite , vexInFmt , vexOutFmt , vexInUnit , vexOutUnit , vexLin , vexSHO , vexSHOR , vexXo ,&
         & vexYo , vexZo , vexRo , vexFx , vexFy , vexFz , vexWx , vexWy , vexWz , vexWr

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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         Tuesday, April 8th, 2014'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '! -------------------------------------------------------------------------'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     STARTING GPSE ... '

         ALLOCATE ( StartValues ( 8 ) )
         CALL DATE_AND_TIME ( startDate , startTime , startZone , StartValues )

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
            mpiIntKind = MPI_INTEGER1

         CASE ( INT16 )

            fIntName = 'INTEGER*2'
            mpiIntName = 'MPI_INTEGER2'
            mpiIntKind = MPI_INTEGER2

         CASE ( INT32 )

            fIntName = 'INTEGER ( INTEGER*4 )'
            mpiIntName = 'MPI_INTEGER ( MPI_INTEGER4 )'
            mpiIntKind = MPI_INTEGER

         CASE ( INT64 )

            fIntName = 'INTEGER*8'
            mpiIntName = 'MPI_INTEGER8'
            mpiIntKind = MPI_INTEGER8

         CASE DEFAULT

            mpiIntKind = -1

      END SELECT

      SELECT CASE ( REAL_DEFAULT_KIND )

         CASE ( REAL32 )

            fRealName = 'REAL ( REAL*4 )'
            mpiRealName = 'MPI_REAL ( MPI_REAL4 )'
            mpiRealKind = MPI_REAL

         CASE ( REAL64 )

            fRealName = 'DOUBLE PRECISION ( REAL*8 )'
            mpiRealName = 'MPI_DOUBLE_PRECISION ( MPI_REAL8 )'
            mpiRealKind= MPI_DOUBLE_PRECISION

         CASE ( REAL128 )

            fRealName = 'REAL*16'
            mpiRealName = 'MPI_REAL16'
            mpiRealKind = MPI_REAL16

         CASE DEFAULT

            mpiRealKind = -1

      END SELECT

      SELECT CASE ( CMPLX_DEFAULT_KIND )

         CASE ( REAL32 )

            fCmplxName = 'COMPLEX ( COMPLEX*8 )'
            mpiCmplxName = 'MPI_COMPLEX ( MPI_COMPLEX8 )'
            mpiCmplxKind = MPI_COMPLEX

         CASE ( REAL64 )

            fCmplxName = 'DOUBLE COMPLEX ( COMPLEX*16 )'
            mpiCmplxName = 'MPI_DOUBLE_COMPLEX ( MPI_COMPLEX16 )'
            mpiCmplxKind = MPI_DOUBLE_COMPLEX

         CASE ( REAL128 )

            fCmplxName = 'COMPLEX*32'
            mpiCmplxName = 'MPI_COMPLEX32'
            mpiCmplxKind = MPI_COMPLEX32

         CASE DEFAULT

            mpiCmplxKind = -1            

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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        runMode    = ', runMode
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dimEq      = ', dimEq
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        odeSolve   = ', odeSolve
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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dTRe       = ', dTRe
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dTIm       = ', dTIm
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dX         = ', dX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dY         = ', dY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        dZ         = ', dZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        gSRe       = ', gSRe
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        gSIm       = ', gSIm
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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexSHO     = ', vexSHO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexSHOR    = ', vexSHOR
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexXo      = ', vexXo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexYo      = ', vexYo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexZo      = ', vexZo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexRo      = ', vexRo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexFx      = ', vexFx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexFy      = ', vexFy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexFz      = ', vexFz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexWx      = ', vexWx
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexWy      = ', vexWy
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexWz      = ', vexWz
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        vexWr      = ', vexWr
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RANGE CHECKING INPUT PARAMETERS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     BROADCASTING INPUT PARAMETERS TO ALL MPI PROCESSES ... '

      END IF

      CALL MPI_BCAST ( runMode    , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dimEq      , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( odeSolve   , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( rk4Lambda  , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( fdOrder    , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( quadRule   , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nTsteps    , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nTwrite    , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nX         , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nXbc       , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nY         , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nYbc       , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )    
      CALL MPI_BCAST ( nZ         , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nZbc       , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( t0         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( xO         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( yO         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( zO         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dTRe       , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dTIm       , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dX         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dY         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dZ         , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( gSRe       , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( gSIm       , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiRead    , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiInFmt   , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiInUnit  , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWrite   , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiOutFmt  , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiOutUnit , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiInit    , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNx      , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNy      , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNz      , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiNr      , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiMl      , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiXo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiYo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiZo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWx      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWy      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWz      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( psiWr      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexRead    , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexInFmt   , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexInUnit  , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexWrite   , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexOutFmt  , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexOutUnit , 1 , mpiIntKind  , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexLin     , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexSHO     , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexSHOR    , 1 , MPI_LOGICAL , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexXo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexYo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexZo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexRo      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexFx      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexFy      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexFz      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexWx      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexWy      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexWz      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( vexWr      , 1 , mpiRealKind , MPI_MASTER , MPI_COMM_WORLD , mpiError )

      IF ( runMode == 0 ) THEN ! perform imaginary time propagation calculation ...

         dT = CMPLX ( 0.0 , -dTIm )
 
      ELSE IF ( runMode == 1 ) THEN ! perform normal calculation ... 

         dT = CMPLX ( dTRe , 0.0 )

      ELSE

         ! ERROR: Run mode not supported. 

      END IF

      IF ( dimEq == 3 ) THEN

         nXa = 1
         nXb = nX
         nYa = 1
         nYb = nY
         nZa = 1 + mpiRank * FLOOR ( REAL ( nZ / mpiProcesses ) )

         IF ( ( mpiRank + 1 ) == mpiProcesses ) THEN ! include remaining z-points on last MPI process ...

            nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcesses ) ) + MODULO ( nZ , mpiProcesses )

         ELSE ! all MPI processes have same number of z-points ... 

            nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcesses ) )

         END IF

         ALLOCATE ( X ( nXa - nXbc : nXb + nXbc ) )
         ALLOCATE ( Y ( nYa - nYbc : nYb + nYbc ) )
         ALLOCATE ( Z ( nZa - nZbc : nZb + nZbc ) )
         ALLOCATE ( Psi3 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

         CALL regular_grid_axis ( nX , nXa - nXbc , nXb + nXbc , xO , dX , X )
         CALL regular_grid_axis ( nY , nYa - nYbc , nYb + nYbc , yO , dY , Y )
         CALL regular_grid_axis ( nZ , nZa - nZbc , nZb + nZbc , zO , dZ , Z )

      ELSE

         ! ERROR: Dimensionality of equation not supported yet.

      END IF

      CALL MPI_FINALIZE ( mpiError )

      DEALLOCATE ( MPIStatus )

! --- FORMAT STATEMENTS ---------------------------------------------------

      STOP

      END PROGRAM

! =========================================================================
