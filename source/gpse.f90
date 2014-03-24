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
!     Monday, March 24th, 2014
!
! -------------------------------------------------------------------------

      PROGRAM GPSE

! --- MODULE DECLARATIONS -------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: GRID
      USE            :: IO
      USE            :: PSI
      USE            :: VEX

! --- MODULE DEFINITIONS --------------------------------------------------
! -------------------------------------------------------------------------
      
      IMPLICIT NONE

      INCLUDE 'mpif.h'

! --- PARAMETER DECLARATIONS  ---------------------------------------------

      CHARACTER ( LEN = * ), PARAMETER :: VERSION_NUMBER = '0.1.0'

      INTEGER, PARAMETER :: INT_DEFAULT_KIND   = KIND ( 0 ) 
      INTEGER, PARAMETER :: INT_SINGLE_KIND    = SELECTED_INT_KIND ( 7 ) 
      INTEGER, PARAMETER :: INT_DBLE_KIND      = SELECTED_INT_KIND ( 15 )
      INTEGER, PARAMETER :: INT_QUAD_KIND      = SELECTED_INT_KIND ( 31 )
      INTEGER, PARAMETER :: REAL_DEFAULT_KIND  = KIND ( 0.0 )
      INTEGER, PARAMETER :: REAL_SINGLE_KIND   = SELECTED_REAL_KIND ( 6  , 37   )   
      INTEGER, PARAMETER :: REAL_DBLE_KIND     = SELECTED_REAL_KIND ( 15 , 307  )
      INTEGER, PARAMETER :: REAL_QUAD_KIND     = SELECTED_REAL_KIND ( 33 , 4931 )
      INTEGER, PARAMETER :: CMPLX_DEFAULT_KIND = KIND ( CMPLX ( 0.0 , 0.0 ) )
      INTEGER, PARAMETER :: MAX_LEN            = 80
      INTEGER, PARAMETER :: MPI_MASTER         = 0

! --- PARAMETER DEFINITIONS -----------------------------------------------
! --- VARIABLE DECLARATIONS -----------------------------------------------

      CHARACTER ( LEN = MAX_LEN ) :: dimEq
      CHARACTER ( LEN = MAX_LEN ) :: dimPsi
      CHARACTER ( LEN = MAX_LEN ) :: dimVex
      CHARACTER ( LEN = MAX_LEN ) :: nameBC
      CHARACTER ( LEN = MAX_LEN ) :: nameEq
      CHARACTER ( LEN = MAX_LEN ) :: namePsi
      CHARACTER ( LEN = MAX_LEN ) :: nameVex
      CHARACTER ( LEN = MAX_LEN ) :: refFrame
      CHARACTER ( LEN = MAX_LEN ) :: quadRule
      CHARACTER ( LEN = MAX_LEN ) :: fdScheme
      CHARACTER ( LEN = MAX_LEN ) :: fdOrder
      CHARACTER ( LEN = MAX_LEN ) :: fmtInput
      CHARACTER ( LEN = MAX_LEN ) :: fmtOutput
      CHARACTER ( LEN = MAX_LEN ) :: cmdArg
      CHARACTER ( LEN = MAX_LEN ) :: runMode
      CHARACTER ( LEN = MAX_LEN ) :: fIntType = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: fRealType = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: fCmplxType = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: mpiIntType = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: mpiRealType = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: mpiCmplxType = 'NONE'
      CHARACTER ( LEN = 8 ) :: startDate = 'NONE'
      CHARACTER ( LEN = 10 ) :: startTime = 'NONE'
      CHARACTER ( LEN = 5 ) :: startZone = 'NONE'
      CHARACTER ( LEN = 8 ) :: stopDate = 'NONE'
      CHARACTER ( LEN = 10 ) :: stopTime = 'NONE'
      CHARACTER ( LEN = 5 ) :: stopZone = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: envPwd = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: envHostname = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: envHome = 'NONE'
      CHARACTER ( LEN = MAX_LEN ) :: envNodename = 'NONE'

      INTEGER :: cmdCnt
      INTEGER :: cmdLen
      INTEGER :: cmdNum
      INTEGER :: cmdStat
      INTEGER :: fileOut
      INTEGER :: fileInp
      INTEGER :: lambda = 2
      INTEGER :: lpkStat
      INTEGER :: mpiCmplx
      INTEGER :: mpiDst
      INTEGER :: mpiErr
      INTEGER :: mpiErrCode
      INTEGER :: mpiErrClass
      INTEGER :: mpiInt
      INTEGER :: mpiProcs
      INTEGER :: mpiProv
      INTEGER :: mpiRank
      INTEGER :: mpiReal
      INTEGER :: mpiSrc
      INTEGER :: ompTs
      INTEGER :: ompTID
      INTEGER :: nDim
      INTEGER :: nIterates
      INTEGER :: nParts
      INTEGER :: nTgS
      INTEGER :: nTomega
      INTEGER :: nTpsi
      INTEGER :: nTsteps
      INTEGER :: nTvex
      INTEGER :: nTwrite
      INTEGER :: nX 
      INTEGER :: nXa
      INTEGER :: nXb
      INTEGER :: nXbc
      INTEGER :: nY
      INTEGER :: nYa
      INTEGER :: nYb
      INTEGER :: nYbc
      INTEGER :: nZ
      INTEGER :: nZa
      INTEGER :: nZb
      INTEGER :: nZbc
      INTEGER :: j , k , l , n , m ! Reserved loop counters. 

      REAL :: xO = 0.0
      REAL :: yO = 0.0
      REAL :: zO = 0.0
      REAL :: dT = 0.0
      REAL :: dX = 0.0
      REAL :: dY = 0.0
      REAL :: dZ = 0.0
      REAL :: gRe = 0.0
      REAL :: gIm = 0.0
      REAL :: xMin = 0.0
      REAL :: xMax = 0.0
      REAL :: yMin = 0.0
      REAL :: yMax = 0.0
      REAL :: zMin = 0.0
      REAL :: zMax = 0.0

! --- VARIABLE DEFINITIONS ------------------------------------------------
! --- ARRAY DECLARATIONS --------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIStat
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

! --- FUNCTION AND SUBROUTINE DEFINITIONS ---------------------------------
! --- MAIN PROGRAM --------------------------------------------------------

      ALLOCATE ( MPIStat ( MPI_STATUS_SIZE ) ) 

      CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiProv , mpiErr )
!     CALL MPI_INIT_THREAD_ERRCHK
      IF ( mpiErr /= MPI_SUCCESS ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'gpse: ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
         CALL MPI_ABORT ( MPI_COMM_WORLD , mpiErrCode , mpiErr )

      END IF
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiProcs , mpiErr )
!     CALL MPI_COMM_SIZE_ERRCHK
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRank , mpiErr )
!     CALL MPI_COMM_RANK_ERRCHK
      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiErr )
!     CALL MPI_BARRIER_ERRCHK

      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '! =========================================================================='
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     GPSE VERSION ',VERSION_NUMBER
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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!         Monday, March 24th, 2014'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '! -------------------------------------------------------------------------'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     STARTING GPSE ... '

         ALLOCATE ( StartValues ( 8 ) )
         CALL DATE_AND_TIME ( startDate , startTime , startZone , StartValues )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RUN STARTED @ ', startTime , ' ON ' , startDate , ' ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     CHECKING MACHINE/COMPILER-SPECIFIC DATA TYPE SUPPORT AND CONFIGURATION ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        INTEGER_KINDS SUPPORTED ... ', INTEGER_KINDS
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        REAL_KINDS SUPPORTED ...    ', REAL_KINDS

      END IF

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiErr )

      SELECT CASE ( INT_DEFAULT_KIND )

         CASE ( INT8 )

            fIntType = 'INTEGER*1'
            mpiIntType = 'MPI_INTEGER1'
            mpiInt = MPI_INTEGER1

         CASE ( INT16 )

            fIntType = 'INTEGER*2'
            mpiIntType = 'MPI_INTEGER2'
            mpiInt = MPI_INTEGER2

         CASE ( INT32 )

            fIntType = 'INTEGER ( INTEGER*4 )'
            mpiIntType = 'MPI_INTEGER ( MPI_INTEGER4 )'
            mpiInt = MPI_INTEGER

         CASE ( INT64 )

            fIntType = 'INTEGER*8'
            mpiIntType = 'MPI_INTEGER8'
            mpiInt = MPI_INTEGER8

         CASE DEFAULT

            mpiInt = -1

      END SELECT

      SELECT CASE ( REAL_DEFAULT_KIND )

         CASE ( REAL32 )

            fRealType = 'REAL ( REAL*4 )'
            mpiRealType = 'MPI_REAL ( MPI_REAL4 )'
            mpiReal = MPI_REAL

         CASE ( REAL64 )

            fRealType = 'DOUBLE_PRECISION ( REAL*8 )'
            mpiRealType = 'MPI_DOUBLE_PRECISION ( MPI_REAL8 )'
            mpiReal = MPI_DOUBLE_PRECISION

         CASE ( REAL128 )

            fRealType = 'REAL*16'
            mpiRealType = 'MPI_REAL16'
            mpiReal = MPI_REAL16

         CASE DEFAULT

            mpiReal = -1

      END SELECT

      SELECT CASE ( CMPLX_DEFAULT_KIND )

         CASE ( REAL32 )

            fCmplxType = 'COMPLEX ( COMPLEX*8 )'
            mpiCmplxType = 'MPI_COMPLEX ( MPI_COMPLEX8 )'
            mpiCmplx = MPI_COMPLEX

         CASE ( REAL64 )

            fCmplxType = 'DOUBLE COMPLEX ( COMPLEX*16 )'
            mpiCmplxType = 'MPI_DOUBLE_COMPLEX ( MPI_COMPLEX16 )'
            mpiCmplx = MPI_DOUBLE_COMPLEX

         CASE ( REAL128 )

            fCmplxType = 'COMPLEX*32'
            mpiCmplxType = 'MPI_COMPLEX32'
            mpiCmplx = MPI_COMPLEX32

         CASE DEFAULT

            mpiCmplx = -1            

      END SELECT

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiErr )

      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        DEFAULT INTEGER KIND ...    ' , fIntType
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        DEFAULT REAL KIND ...       ' , fRealType
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        DEFAULT COMPLEX KIND ...    ' , fCmplxType
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        MPI_INTEGER TYPE ...        ' , mpiIntType
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        MPI_REAL TYPE ...           ' , mpiRealType
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!        MPI_COMPLEX TYPE ...        ' , mpiCmplxType
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     READING AND PARSING COMMAND-LINE ARGUMENTS ... '
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '!     RANGE CHECKING INPUT PARAMETER VALUES ... '

         nDim = 3
         nX = 256
         nXbc = 2
         nY = 256
         nYbc = 2
         nZ = 256
         nZbc = 2
         dX = 0.015625
         dY = 0.015625
         dZ = 0.015625

      END IF

      CALL MPI_BCAST ( nDim , 1 , MPI_INTEGER , MPI_MASTER , MPI_COMM_WORLD , mpiErr )

      IF ( nDim == 3 ) THEN

         CALL MPI_BCAST ( nX , 1 , mpiInt , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( nXbc , 1 , mpiInt , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( nY , 1 , mpiInt , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( nYbc , 1 , mpiInt , MPI_MASTER , MPI_COMM_WORLD , mpiErr )    
         CALL MPI_BCAST ( nZ , 1 , mpiInt , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( nZbc , 1 , mpiInt , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( dX , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( dY , 1 , mpiReal , MPI_MASTER , MPI_COMM_WORLD , mpiErr )
         CALL MPI_BCAST ( dZ , 1 , mpiReal, MPI_MASTER , MPI_COMM_WORLD , mpiErr )

         nXa = 1
         nXb = nX
         nYa = 1
         nYb = nY
         nZa = 1 + mpiRank * FLOOR ( REAL ( nZ / mpiProcs ) )

         IF ( ( mpiRank + 1 ) == mpiProcs ) THEN ! include any remainder grid points in z-range on  

            nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcs ) ) + MODULO ( nZ , mpiProcs )

         ELSE

            nZb = ( mpiRank + 1 ) * FLOOR ( REAL ( nZ / mpiProcs ) )

         END IF

         ALLOCATE ( X ( nXa - nXbc : nXb + nXbc ) )
         ALLOCATE ( Y ( nYa - nYbc : nYb + nYbc ) )
         ALLOCATE ( Z ( nZa - nZbc : nZb + nZbc ) )
         ALLOCATE ( Psi3 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

         CALL regular_grid_axis ( nX , nXa - nXbc , nXb + nXbc , xO , dX , X )
         CALL regular_grid_axis ( nY , nYa - nYbc , nYb + nYbc , yO , dY , Y )
         CALL regular_grid_axis ( nZ , nZa - nZbc , nZb + nZbc , zO , dZ , Z )

         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  Psi3 ( j , k , l ) = psi_3d_se_sho_ani ( 0 , 0 , 0 , xO , yO , zO , 1.0 , 1.0 , 2.0 , X ( j ) , Y ( k ) , Z ( l ) )

               END DO

            END DO 

         END DO

      ELSE

      END IF

      CALL MPI_FINALIZE ( mpiErr )

      DEALLOCATE ( MPIStat )

! --- FORMAT STATEMENTS ---------------------------------------------------

      STOP

      END PROGRAM

! =========================================================================
