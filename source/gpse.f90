! =====================================================================
! NAME
!
!     gpse [ jip-see ] - Gross-Pitaevskii and Schrodinger Equation Solver 
!
! SYNOPSIS
!
!     export OMP_NUM_THREADS=X
!     mpirun -np Y ./gpse.x > gpse.output
!
! DESCRIPTION  
!
!     GPSE is a Fortran program created to numerically approximate solutions of the time-dependent Gross-Pitaevskii equation (GPE) 
!        in non-uniformly rotating frames of reference for a single component Bose-Einstein condensate (BEC) in three dimensions. 
!        It implements a method-of-lines approach that combines an explicit, generalized 4th-order Runge-Kutta (GRK4) 
!        time-integration scheme* with 2nd- or 4th-order central differences (CD2/CD4) to discretize the spatial derivatives of the
!        equation.
!
!     Zero Dirichlet boundary conditions are assumed to apply at the boundary of the computational domain at all times in any given 
!        simulation. As a result, gpse is best suited for investigating BECs confined within externally applied trapping potentials
!        where the condensate density decreases to zero at the boundary of the computational domain.
!
!     The GRK4 + CDX algorithm is implemented in parallel using a one-dimensional slab decomposition of the computational domain 
!        embedded within a hybrid MPI + OpenMP framework to enable scalable, high-resolution numerical simulations. A two-round 
!        communication pattern, where even- and odd-numbered MPI processes alternate between send and receive calls, is utilized to
!        reduce communication overhead and improve message passing throughput.
!
!     * See 'On A General Formula of Fourth Order Runge-Kutta Method' by D. Tan and Z. Chen in the Journal of Mathematical Sciences
!        & Mathematics Education Vol. 7 No. 2 (2012).
!
! OPTIONS
!
!     All runtime options available to the user, with the exception of the number of MPI processes and OpenMP threads to run on, 
!        must be specified within the input file gpse.input.
!
! SEE ALSO
!
! BUGS
!
! HISTORY
!
! ADDITIONAL NOTES
!
!     1. Add support for rectilinear grids.
!     2. Add support for spatial- and time-varying mean-field coupling.
!     3. Add support for imaginary potentialis, including mean-field.
!     4. Add support for translating reference frames.
!     5. Add support for general reference frames.
!     6. Add support for dipolar condensates; requires complex comm -- multi-two-round?
!     7. Add support for composite condensates.
!     8. Add support for one- and two-dimensions.
!     9. Add iterative methods.
!     10. Add support for bosonic and fermionic many-body systems.
!     11. Add support for accelerator hardware, i.e., OpenACC
!     12. Remodularize source code. Use submodules?
!     13. Rewrite IO module for parallel I/O with MPI I/O.
!     14. Explore other MPI communication patterns.
!     15. Add command-line interface.
!     16. Complete support for 6th- and 8th-order central differences.
!     17. Explore OpenMP scheduling options.
!     18. Create integrated, automated test infrastructure.
!     19. Test automatic differentiation.
!     20. Segment limited I/O output into multiple files, e.g. one file per expectation value.
!     21. Support Cartesian, regular, rectilinear and curvilinear grids.
!     22. Compare performance of array functions vs. passing by reference to subroutines, e.g. in computing potentials
!
! AUTHOR
!
!     Marty Kandes, Ph.D.
!     Distributed High-Throughput Computing Group
!     San Diego Supercomputer Center
!     University of California, San Diego
!
! COPYRIGHT
!     
!     Copyright (c) 2014, 2015, 2016, 2017 Martin Charles Kandes
!
! LAST UPDATED
!
!     Sunday, January 22nd, 2017
!
! ----------------------------------------------------------------------------------------------------------------------------------

      PROGRAM GPSE

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: MPI
      USE            :: EVUA
      USE            :: GRID
      USE            :: GRK4
      USE            :: IO
      USE            :: MATH
      USE            :: PSI
      USE            :: ROT
      USE            :: VEX

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
!     MPI is the standard Message Passing Interface module.
!
!     EVUA is a custom Fortran module written to compute expectation values and other useful quantities, such as uncertainty 
!        relations, from the wave functions generated during program execution. 
!
!     GRID is a custom Fortran module written to build the regular Cartesian grid of points that defines the spatial computational 
!        domain on which external potentials, wave functions and other physical quantities may be computed. 
!
!     GRK4 is a custom Fortran module written to numerically approximate solutions of the time-dependent, three-dimensional 
!        Gross-Pitaveskii equation for a single component Bose-Einstein condensate in a rotating reference frame using the 
!        GRK4 + CDX algorithm.
!
!     IO is a custom Fortran module written to perform input-output operations for gpse.
!
!     MATH is a custom Fortran module written to define well-know mathematical constants and compute specialized functions. 
!
!     PSI is a custom Fortran module written to compute analytic solutions of the Schrodinger and Gross-Pitaevskii equations, which
!        may be used as initial conditions for simulations.
!
!     ROT is a custom Fortran module written to compute the elemental rotation matrices of a Cartesian coordinate system.
!
!     VEX is a custom Fortran module written to compute analytic external potentials.
!     
! ----------------------------------------------------------------------------------------------------------------------------------
      
      IMPLICIT NONE

!      INCLUDE 'mpif.h' ! if MPI module is not supported, use header file.

! --- PARAMETER DECLARATIONS  ------------------------------------------------------------------------------------------------------

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.5.8'
      CHARACTER ( LEN = * ), PARAMETER :: GPSE_LAST_UPDATED = 'Sunday, January 22nd, 2017'

      INTEGER, PARAMETER :: MPI_MASTER = 0

! --- PARAMETER DEFINITIONS --------------------------------------------------------------------------------------------------------
!
!     GPSE_VERSION_NUMBER is a character-valued parameter that defines the current working version of the gpse's source code. It is
!        to be maintained manually by the author(s) of the program.
!
!     GPSE_LAST_UPDATED is a character-valude parameter that defines the last date any portion of the current working version of the 
!        gpse's source code was modified. It is to be maintained manually by the author(s) of the program
!
!     MPI_MASTER is an integer-valued parameter that defines the default rank of the master MPI process when the program begins 
!        execution. 
!
! --- VARIABLE DECLARATIONS --------------------------------------------------------------------------------------------------------

      LOGICAL :: itpOn = .FALSE.

      CHARACTER ( LEN = 8  ) :: startDate     = 'NONE'
      CHARACTER ( LEN = 10 ) :: startTime     = 'NONE'
      CHARACTER ( LEN = 5  ) :: startZone     = 'NONE'
      CHARACTER ( LEN = 8  ) :: stopDate      = 'NONE'
      CHARACTER ( LEN = 10 ) :: stopTime      = 'NONE'
      CHARACTER ( LEN = 5  ) :: stopZone      = 'NONE'
      CHARACTER ( LEN = 4 )  :: fileUnitChar  = 'None'

      INTEGER :: rk4Lambda      = -1 
      INTEGER :: fdOrder        = -1 
      INTEGER :: quadRule       = 1  
      INTEGER :: nTsteps        = -1 
      INTEGER :: nTwrite        = -1 
      INTEGER :: nX             = -1 
      INTEGER :: nXa            = -1 
      INTEGER :: nXb            = -1 
      INTEGER :: nXbc           = -1
      INTEGER :: dNx            = -1
      INTEGER :: nY             = -1 
      INTEGER :: nYa            = -1 
      INTEGER :: nYb            = -1
      INTEGER :: nYbc           = -1
      INTEGER :: dNy            = -1
      INTEGER :: nZ             = -1
      INTEGER :: nZa            = -1
      INTEGER :: nZb            = -1
      INTEGER :: nZbc           = -1
      INTEGER :: dNz            = -1
      INTEGER :: nXpsi          = -1
      INTEGER :: nYpsi          = -1
      INTEGER :: nZpsi          = -1
      INTEGER :: nRpsi          = -1
      INTEGER :: mLpsi          = 0 
      INTEGER :: mpiCmplx       = -1 
      INTEGER :: mpiError       = -1 
      INTEGER :: mpiErrorCode   = -1 
      INTEGER :: mpiInt         = -1 
      INTEGER :: mpiProcesses   = -1 
      INTEGER :: mpiProvided    = -1 
      INTEGER :: mpiRank        = -1
      INTEGER :: mpiSource      = -1
      INTEGER :: mpiDestination = -1
      INTEGER :: mpiReal        = -1
      INTEGER ( KIND = MPI_OFFSET_KIND ) :: mpiOffset = -1
      INTEGER :: mpiFileHandle  = -1
      INTEGER :: ompThreads     = -1  
      INTEGER :: j , k , l , m , n 
      
      REAL :: tN      = 0.0 
      REAL :: t0      = 0.0
      REAL :: tF      = 0.0 
      REAL :: xO      = 0.0 
      REAL :: yO      = 0.0 
      REAL :: zO      = 0.0 
      REAL :: dT      = 0.0 
      REAL :: dX      = 0.0 
      REAL :: dY      = 0.0 
      REAL :: dZ      = 0.0 
      REAL :: xOrrf   = 0.0 
      REAL :: yOrrf   = 0.0 
      REAL :: zOrrf   = 0.0 
      REAL :: wXo     = 0.0 
      REAL :: wYo     = 0.0 
      REAL :: wZo     = 0.0 
      REAL :: wX      = 0.0
      REAL :: wY      = 0.0
      REAL :: wZ      = 0.0 
      REAL :: gS      = 0.0 
      REAL :: xOpsi   = 0.0 
      REAL :: yOpsi   = 0.0 
      REAL :: zOpsi   = 0.0 
      REAL :: rOpsi   = 0.0 
      REAL :: wXpsi   = 0.0 
      REAL :: wYpsi   = 0.0 
      REAL :: wZpsi   = 0.0 
      REAL :: wRpsi   = 0.0 
      REAL :: pXpsi   = 0.0 
      REAL :: pYpsi   = 0.0 
      REAL :: pZpsi   = 0.0 
      REAL :: xOvex   = 0.0 
      REAL :: yOvex   = 0.0 
      REAL :: zOvex   = 0.0 
      REAL :: rOvex   = 0.0 
      REAL :: fXvex   = 0.0 
      REAL :: fYvex   = 0.0 
      REAL :: fZvex   = 0.0 
      REAL :: wXvex   = 0.0 
      REAL :: wYvex   = 0.0 
      REAL :: wZvex   = 0.0 
      REAL :: wRvex   = 0.0 
      REAL :: thetaXo = 0.0
      REAL :: nu      = 0.0 
      REAL :: sigma   = 0.0 
      REAL :: tSigma  = 0.0 

      COMPLEX :: dTz = CMPLX ( 0.0 , 0.0 ) 

! --- VARIABLE DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     itpOn is a LOGICAL input variable used as a flag to determine if imaginary time propagation will be performed.
!
!     startDate is a CHARACTER variable that records the start date of program execution.
!
!     startTime is a CHARACTER variable that records the start time of program execution.
!
!     startZone is a CHARACTER variable that records the time zone of the system.
!
!     stopDate is a CHARACTER variable that records the date when program execution ends.
!
!     stopTime is a CHARACTER variable that records the time when program execution ends.
!
!     stopZone is a CHARACTER variable that records the time zone of the system.
!
!     rk4Lambda is an INTEGER-valued input variable that selects the specific Runge-Kutta scheme used for the time-integration. 
!        1 = Tan-Chen-1 ; 2 = Classical 4th-Order Runge-Kutta ; 3 = Tan-Chen Lambda-3 ; 4 = England ; 5 = Tan-Chen-5
!
!     fdOrder is an INTEGER-valued input variable that selects the order-of-accuracy of the central differences used to compute
!        approximations to the spatial derivatives in the GPE. 
!
!     quadRule is an INTEGER-valued variable that is intended to select which quadrature rule is used to perform numerical 
!        integration in the calculation of expectation values. This option has not been implemented yet.
!
!     nTsteps is an INTEGER-valued input variable that sets the total number of time steps to be performed in a simulation.
!
!     nTwrite is an INTEGER-valued input variable that sets the period with which data is written to disk during a simulation (in 
!        number of time steps).
!
!     nX is an INTEGER-valued input variable that sets the number of grid points along the x-axis of the system.
!
!     nXa is an INTEGER-valued variable that sets the lower bound index of grid points along the x-axis.
!
!     nXb is an INTEGER-valued variable that sets the upper bound index of grid points along the x-axis.
!
!     nXbc is an INTEGER-valued variable that sets the number of grid points in the boundary condition along the x-axis required by
!        the order-of-accuracy of central differences to be used.
!
!     nY is an INTEGER-valued input variable that sets the number of grid points along the y-axis of the system.
!
!     nYa is an INTEGER-valued variable that sets the lower bound index of grid points along the y-axis.
!
!     nYb is an INTEGER-valued variable that sets the upper bound index of grid points along the y-axis.
!
!     nYbc is an INTEGER-valued variable that sets the number of grid points in the boundary condition along the y-axis required by
!        the order-of-accuracy of central differences to be used.
!
!     nZ is an INTEGER-valued input variable that sets the number of grid points along the z-axis of the system.
!
!     nZa is an INTEGER-valued variable that sets the lower bound index of grid points along the z-axis.
!
!     nZb is an INTEGER-valued variable that sets the upper bound index of grid points along the z-axis.
!
!     nZbc is an INTEGER-valued variable that sets the number of grid points in the boundary condition along the z-axis required by
!        the order-of-accuracy of central differences to be used.
!
!     nXpsi is an INTEGER-valued input variable that sets the degree of Hermite polynomial used to define anisotropic SHO wave 
!        function along x-axis, if psiInit = 1.
!
!     nYpsi is an INTEGER-valued input variable that sets the degree of Hermite polynomial used to define anisotropic SHO wave 
!        function along y-axis, if psiInit = 1.
!
!     nZpsi is an INTEGER-valued input variable that sets the degree of Hermite polynomial used to define both the anisotropic and 
!        axially-symmetric SHO wave functions along z-axis, if psiInit = 1 or psiInit = 2.
!
!     nRpsi is an INTEGER-valued input variable that sets the degree of (associated) Laguerre polynomials used to define radial 
!        components of isotropic and axially-symmetric SHO wave functions, if psiInit = 0 or psiInit = 2.
!
!     mLpsi is an INTEGER-valued input variable that sets the orbital angular momentum quantum number along the z-axis for an 
!        axially-symmetric SHO wave function, i.e., when psiInit = 2.
!
!     mpiCmplx is an INTEGER-valued variable that stores the default KIND parameter for MPI COMPLEXs.
!
!     mpiError is an INTEGER-valued variable that stores MPI error codes returned by MPI subroutines.
!
!     mpiErrorCode is an INTEGER-valued variable that is used to return an MPI error code to invoking environment via MPI_ABORT.
!
!     mpiInt is an INTEGER-valued variable that stores the default KIND parameter for MPI INTEGERs.
!
!     mpiProcesses is an INTEGER-valued variable that stores the total number of MPI processes in use.
!
!     mpiProvided is an INTEGER-valued variable that stores the level of thread support that is provided by MPI.
!
!     mpiRank is an INTEGER-valued variable that stores the rank of an MPI process. 
!
!     mpiSource is an INTEGER-valued variable that stores the rank of a source MPI process, where an MPI_SEND call originates.
!
!     mpiDestination is an INTEGER-valued variable that stores the rank of a destination MPI process, where an MPI_RECV call 
!        originates.
!
!     mpiReal is an INTEGER-valued variable that stores the default KIND parameter for MPI REALs.
!
!     ompThreads is an INTEGER-valued variable that stores the total number of OpenMP threads in a PARALLEL region.
!
!     j, k, l , m, and n are INTEGER-valued variables reserved for use as loop counters.
!
!     tN is a REAL-valued input variable that stores the running simulation time.
!
!     t0 is a REAL-valued input variable that stores the initial simulation time.
!
!     tF is a REAL-valued input variable that stores the final simulation time.
!
!     xO is a REAL-valued input variable that stores the x-coordinate of the centre of the computational domain.
!
!     yO is a REAL-valued input variable that stores the y-coordinate of the centre of the computational domain. 
!
!     zO is a REAL-valued input variable that stores the z-coordinate of the centre of the computational domain.
!
!     dT is a REAL-valued input variable that stores the interval of a simulation time step.
!
!     dX is a REAL-valued input variable that stores the distance between grid points along the x-axis of a regular grid. 
!
!     dY is a REAL-valued input variable that stores the distance between grid points along the y-axis of a regular grid.
!
!     dZ is a REAL-valued input variable that stores the distance between grid points along the z-axis of a regular grid.
!
!     xOrrf is a REAL-valued input variable that stores the x-coordinate of the origin of the rotating frame. 
!
!     yOrrf is a REAL-valued input variable that stores the y-coordinate of the origin of the rotating frame.
!
!     zOrrf is a REAL-valued input variable that stores the z-coordinate of the origin of the rotating frame.
!
!     wXo is a REAL-valued variable that stores the initial rotation rate of the reference frame about the x-axis.
!
!     wYo is a REAL-valued variable that stores the initial rotation rate of the reference frame about the y-axis.
!
!     wZo is a REAL-valued variable that stores the initial rotation rate of the reference frame about the z-axis.
!
!     wX is a REAL-valued input variable that sets the initial rotation rate of the reference frame about the x-axis as well as 
!        stores its time-varying value, if applicable.
!
!     wY is a REAL-valued input variable that sets the initial rotation rate of the reference frame about the y-axis as well as 
!        stores its time-varying value, if applicable. 
!
!     wZ is a REAL-valued input variable that sets the initial rotation rate of the reference frame about the z-axis as well as 
!        stores its time-varying value, if applicable. 
!
!     gS is a REAL-valued input variable that determines the nonlinear, mean-field coupling constant set by the strength of 
!        atom-atom interactions within a Bose-Einstein condensate.
!
!     xOpsi is a REAL-valued input variable that sets the x-coordinate of the origin to define an initial wave function built using
!        the psi_init subroutine.
!
!     yOpsi is a REAL-valued input variable that sets the y-coordinate of the origin to define an initial wave function built using
!        the psi_init subroutine.
!
!     zOpsi is a REAL-valued input variable that sets the z-coordinate of the origin to define an initial wave function built using
!        the psi_init subroutine. 
!
!     rOpsi is a REAL-valued input variable that sets the radius of the approximate analytic wave function for a simple harmonic
!        oscillator external potential. See psi_init for psiInit = 3.
!
!     wXpsi is a REAL-valued input variable that sets the angular frequency of a simple harmonic oscillator potential along x-axis
!        used to define an initial anisotropic simple harmonic oscillator wave function built using the psi_init subroutine.
!
!     wYpsi is a REAL-valued input variable that sets the angular frequency of a simple harmonic oscillator potential along y-axis
!        used to define an initial anisotropic simple harmonic oscillator wave function built using the psi_init subroutine.
!
!     wZpsi is a REAL-valued input variable that sets the angular frequency of a simple harmonic oscillator potential along z-axis
!        used to define an initial anisotropic simple harmonic oscillator wave function built using the psi_init subroutine.
!
!     wRpsi is a REAL-valued input variable that sets the radial angular frequency of used to define an initial axially-symmetric 
!        simple harmonic oscillator wave function or the approximate simple harmonic oscillator ring wave function built using the 
!        psi_init subroutine.
!
!     pXpsi is a REAL-valued input variable that sets the magnitude and direction of an initial linear momentum boost to the wave 
!        function along the x-axis.
!
!     pYpsi is a REAL-valued input variable that sets the magnitude and direction of an initial linear momentum boost to the wave 
!        function along the y-axis.
!
!     pZpsi is a REAL-valued input variable that sets the magnitude and direction of an initial linear momentum boost to the wave 
!        function along the z-axis.
!
!     xOvex is a REAL-valued input variable that sets the x-coordinate of the external potential's origin.
!
!     yOvex is a REAL-valued input variable that sets the y-coordinate of the external potential's origin.
!
!     zOvex is a REAL-valued input variable that sets the z-coordinate of the external potential's origin.
!
!     rOvex is a REAL-valued input variable that sets the radius of a simple harmonic oscillator ring potential.
!
!     fXvex is a REAL-valued input variable that sets the magnitude and direction of a constant force along the x-axis.
!
!     fYvex is a REAL-valued input variable that sets the magnitude and direction of a constant force along the y-axis.
!
!     fZvex is a REAL-valued input variable that sets the magnitude and direction of a constant force along the z-axis.
!
!     wXvex is a REAL-valued input variable that sets the angular frequency of a simple harmonic oscillator potential along the 
!        x-axis.
!
!     wYvex is a REAL-valued input variable that sets the angular frequency of a simple harmonic oscillator potential along the 
!        y-axis.
!
!     wZvex is a REAL-valued input variable that sets the angular frequency of a simple harmonic oscillator potential along the 
!        z-axis.
!
!     wRvex is a REAL-valued input variable that sets the radial angular fequency of a simple harmonic oscillator ring potential.
!
!     thetaXo is a REAL-valued variable that sets the maximum angle the rotating frame's angular velocity vector is nutated about 
!        the x-axis.
!
!     nu is a REAL-valued variable that sets the frequency at which the rotating frame's angular velocity vector is nutated about 
!        the x-axis.
!
!     sigma is a REAL-valued variable that sets the rate at which the rotating frame's angular velocity vector is flipped 180 
!        degrees about the x-axis. 
!
!     tSigma is a REAL-valued variable that sets the time at which the rotating frame's angular velocity vector is flipped 180 
!        degrees about the x-axis.
!
!     dTz is a COMPLEX_valued variable that determines the interval of a simulation time.
!
! --- ARRAY DECLARATIONS -----------------------------------------------------------------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StartValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: StopValues
      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MpiStatus ! replace usage with MPI_STATUS_IGNORE?

      REAL, ALLOCATABLE, DIMENSION ( : ) :: Xa
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Ya
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Za
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Zb
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Zc
      REAL, ALLOCATABLE, DIMENSION ( : ) :: Omega
      REAL, ALLOCATABLE, DIMENSION ( : , : , : ) :: Vex3a
      REAL, ALLOCATABLE, DIMENSION ( : , : ) :: Psi2a
      REAL, ALLOCATABLE, DIMENSION ( : , : ) :: Psi2b

      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: K4
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3a
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : ) :: Psi3b

! --- ARRAY DEFINITIONS ------------------------------------------------------------------------------------------------------------
!
!     StartValues is an INTEGER-valued, rank-one array that stores the values returned by the real-time system clock when the 
!        DATE_AND_TIME subroutine is called at the beginning of program execution.
!
!     StopValues is an INTEGER-valued, rank-one array that stores the values returned by the real-time system clock when the
!        DATE_AND_TIME subroutine is called at the end of program execution.
!
!     MpiStatus is an INTEGER-valued, rank-one array that stores the MPI status object. Replace usage with MPI_STATUS_IGNORE?
!
!     Xa is a REAL-valued, rank-one array that stores the x-coordinates of the grid points that define the computational domain.
!
!     Ya is a REAL-valued, rank-one array that stores the y-coordinates of the grid points that define the computational domain.
!
!     Za is a REAL-valued, rank-one array that stores the z-coordinates of the grid points that definte the computational domain.
!
!     Zb is a REAL-valued, rank-one array that stores the z-coordinates of the grid points that definte the computational domain.
!
!     Omega is a REAL-valued, rank-one array that stores the computed values of the time-varying angular velocity vector of the
!        rotating frame.
!
!     Vex3a is a REAL-valued, rank-three array that stores the external potential.
!
!     K1 is a COMPLEX-valued, rank-three array that stores the first increment of the generalized 4th-order Runge-Kutta 
!        time-integration scheme.
!
!     K2 is a COMPLEX-valued, rank-three array that stores the second increment of the generalized 4th-order Runge-Kutta 
!        time-integration scheme.
!
!     K3 is a COMPLEX-valued, rank-three array that stores the third increment of the generalized 4th-order Runge-Kutta 
!        time-integration scheme.
!
!     K4 is a COMPLEX-valued, rank-three array that stores the fourth increment of the generalized 4th-order Runge-Kutta 
!        time-integration scheme.
! 
!     Psi3a is a COMPLEX-valued, rank-three array that stores the wave function.
!
!     Psi3b is a COMPLEX-valued, rank-three array that stores the wave function.
!
! --- FUNCTION AND SUBROUTINE DECLARATIONS -----------------------------------------------------------------------------------------

      INTEGER :: OMP_GET_NUM_THREADS

! --- FUNCTION AND SUBROUTINE DEFINITIONS ------------------------------------------------------------------------------------------
!
!     OMP_GET_NUM_THREADS is an INTEGER scalar function that returns the total number of OpenMP threads executed in a parallel 
!        region during program execution.
!
! --- NAMELIST DECLARATIONS --------------------------------------------------------------------------------------------------------

      NAMELIST /gpseIn/ itpOn , rk4Lambda , fdOrder , nTsteps , nTwrite , nX , nY , nZ , dNx , dNy , dNz , t0 , tF , xO , yO , zO ,&
         & dT , dX , dY , dZ , xOrrf , yOrrf , zOrrf , wX , wY , wZ , gS , psiInput , psiOutput , psiFileNo , psiInit , nXpsi , &
         & nYpsi , nZpsi , nRpsi , mLpsi , xOpsi , yOpsi , zOpsi , rOpsi , wXpsi , wYpsi , wZpsi , wRpsi , pXpsi , pYpsi , pZpsi ,& 
         & vexInput , vexOutput , vexFileNo , vexInit , xOvex , yOvex , zOvex , rOvex , fXvex , fYvex , fZvex , wXvex , wYvex , & 
         & wZvex , wRvex

! --- NAMELIST DEFINITIONS ---------------------------------------------------------------------------------------------------------
!
!     gpseIn is the input variable namelist, which must be defined in the file gpse.input.
!
! --- BEGIN MAIN PROGRAM -----------------------------------------------------------------------------------------------------------

      CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiProvided , mpiError )
      IF ( mpiError /= MPI_SUCCESS ) THEN

         WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
         CALL MPI_ABORT ( MPI_COMM_WORLD , mpiErrorCode , mpiError )
         STOP

      END IF
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiProcesses , mpiError )
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRank , mpiError )

      ALLOCATE ( StartValues ( 8 ) )
      ALLOCATE ( StopValues  ( 8 ) )
      CALL DATE_AND_TIME ( startDate , startTime , startZone , StartValues )

      ALLOCATE ( MpiStatus ( MPI_STATUS_SIZE ) )

!$OMP PARALLEL DEFAULT ( SHARED )

      !ompThreads = OMP_GET_NUM_THREADS ( )

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

      wXo = wX
      wYo = wY
      wZo = wZ

      ALLOCATE ( Xa ( nXa - nXbc : nXb + nXbc ) )
      ALLOCATE ( Ya ( nYa - nYbc : nYb + nYbc ) )
      ALLOCATE ( Za ( nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Zb ( nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Zc ( 1 - nZbc : nZ + nZbc ) )
      ALLOCATE ( Omega ( 3 ) )
      ALLOCATE ( Vex3a ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K1 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K2 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K3 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( K4 ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Psi3a ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )
      ALLOCATE ( Psi3b ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ) )

      Xa = 0.0
      Ya = 0.0
      Za = 0.0
      Zb = 0.0
      Zc = 0.0
      Omega = 0.0
      Vex3a = 0.0
      K1 = CMPLX ( 0.0 , 0.0 )
      K2 = CMPLX ( 0.0 , 0.0 )
      K3 = CMPLX ( 0.0 , 0.0 )
      K4 = CMPLX ( 0.0 , 0.0 )
      Psi3a = CMPLX ( 0.0 , 0.0 )
      Psi3b = CMPLX ( 0.0 , 0.0 )

      CALL grid_regular_axis ( nX , nXa , nXb , nXbc , xO , dX , Xa ) 
      CALL grid_regular_axis ( nY , nYa , nYb , nYbc , yO , dY , Ya ) 
      CALL grid_regular_axis ( nZ , nZa , nZb , nZbc , zO , dZ , Za )
      CALL grid_regular_axis ( nZ , 1 , nZ , nZbc , zO , dZ , Zc )

      IF ( psiInput == 0 ) THEN ! compute initial wave function from available analytic expression

         CALL psi_init ( psiInit , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , nXpsi , nYpsi , nZpsi , nRpsi , mLpsi ,&
            & xOpsi , yOpsi , zOpsi , rOpsi , wXpsi , wYpsi , wZpsi , wRpsi , Xa , Ya , Za , Psi3a )

      ELSE IF ( psiInput == 1 ) THEN ! read initial wave function from binary file on MPI_MASTER

         psiFilePos = 1 ! initialize file position

         IF ( mpiRank == MPI_MASTER ) THEN ! read and initialize first block of wave function values from binary file on MPI_MASTER

            CALL io_read_bin_psi ( psiFileNo , psiFilePos , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3a )

         END IF

         DO mpiDestination = 1 , mpiProcesses - 1 ! loop over mpiDesination ranks that are not equal to MPI_MASTER

            IF ( mpiRank == MPI_MASTER ) THEN ! read in next block of wave function values from binary file on MPI_MASTER

               CALL io_read_bin_psi ( psiFileNo , psiFilePos , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3b )

            END IF
!           send read block from MPI_MASTER to mpiRank equal to mpiDestination
            CALL mpi_copy_psi ( mpiRank , MPI_MASTER , mpiDestination , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , &
               & Psi3b )

         END DO

         IF ( mpiRank /= MPI_MASTER ) THEN ! after all blocks have been sent, copy wave function values into Psi3a; 
                                           ! MPI_MASTER block was read into Psi3a directly
            Psi3a = Psi3b

         END IF

      ELSE

         IF ( mpiRank == MPI_MASTER ) THEN

            WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : ERROR - psiInput not reognized.'
            STOP

         END IF

      END IF
      CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3a )
      CALL psi_boost ( nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOpsi , yOpsi , zOpsi , pXpsi , pYpsi , pZpsi , &
         & Xa , Ya , Za , Psi3a )

      IF ( vexInput == 0 ) THEN ! compute initial external potential from known analytic expression

         CALL vex_init ( vexInit , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOvex , yOvex , zOvex , rOvex , fXvex ,&
            & fYvex , fZvex , wXvex , wYvex , wZvex , wRvex , Xa , Ya , Za , Vex3a )

      ELSE

         IF ( mpiRank == MPI_MASTER ) THEN

            WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'gpse : ERROR - vexInput not recognized.'
            STOP

         END IF 

      END IF

      psiFileNo = 1000
      vexFileNo = 1000

! --- BEGIN MAIN TIME PROPAGATION LOOP ---------------------------------------------------------------------------------------------

      tN = t0 ! initialize simulation time

      DO n = 0 , nTsteps

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

         IF ( MODULO ( n , nTwrite ) == 0 ) THEN ! compute partial 
!        relations to file from MPI_MASTER; write wave function and 
!        external potential to file from MPI_MASTER

            ! Compute partial base expectation values locally on each MPI process; reduce partial base expectation values from all 
            ! MPI processes to MPI_MASTER to get full base expectation values
            CALL evua_compute_base ( MPI_MASTER , mpiReal , mpiError , 1 , fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , &
               & nZb , nZbc , xO , yO , zO , dX , dY , dZ , gS , Xa , Ya , Za , Vex3a , Psi3a )
            ! Compute derived expectation values, uncertainties from base expectation values
            CALL evua_compute_derived ( mpiRank , MPI_MASTER , wX , wY , wZ )
            ! Write expectation values, undertainties and uncertainty relations to file from MPI_MASTER
            CALL evua_write_all ( mpiRank , MPI_MASTER , tN , wX , wY , wZ )

            IF ( psiOutput == 1 ) THEN ! Write wave function to file from MPI_MASTER using streaming I/O binary with partial reduce

               Psi3b = Psi3a
               psiFilePos = 1 
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_psi ( mpiRank , mpiSource , MPI_MASTER , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , &
                     & Psi3b )
                  IF ( mpiRank == MPI_MASTER ) THEN

                     CALL io_write_bin_psi ( psiFileNo , psiFilePos , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , &
                        & Psi3b )

                  END IF

               END DO
               psiFileNo = psiFileNo + 1

            ELSE IF ( psiOutput == 2 ) THEN ! Write the wave function to a legacy VTK file from MPI_MASTER using streaming I/O
                                            ! with partial reduce
               Zb = Za
               Psi3b = Psi3a
               psiFilePos = 1 ! initialize file position
               IF ( mpiRank == MPI_MASTER ) THEN

                  CALL io_write_vtk_header ( 'psi-', psiFileNo , psiFilePos , nX , nY , nZ , dNx , dNy , dNz )
                  CALL io_write_vtk_xcoordinates ( 'psi-' , psiFileNo , psiFilePos , nX , nXa , nXb , nXbc , dNx , Xa )
                  CALL io_write_vtk_ycoordinates ( 'psi-' , psiFileNo , psiFilePos , nY , nYa , nYb , nYbc , dNy , Ya )

               END IF
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_q ( mpiRank , mpiSource , MPI_MASTER , nZ , nZa , nZb , nZbc , Zb )
                  IF ( mpiRank == MPI_MASTER ) THEN

                     CALL io_write_vtk_zcoordinates ( 'psi-' , psiFileNo , psiFilePos , mpiSource , nZ , nZa , nZb , nZbc , dNz , Zb )

                  END IF

               END DO
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_psi ( mpiRank , mpiSource , MPI_MASTER , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , &
                     & Psi3b )
                  IF ( mpiRank == MPI_MASTER ) THEN

                     CALL io_write_vtk_repsi ( 'psi-' , psiFileNo , psiFilePos , mpiSource , nX , nXa , nXb , nXbc , dNx , nY , &
                        & nYa , nYb , nYbc , dNy , nZ , nZa , nZb , nZbc , dNz , Psi3b )

                  END IF

               END DO
               Psi3b = Psi3a
               DO mpiSource = 0 , mpiProcesses - 1

                  CALL mpi_copy_psi ( mpiRank , mpiSource , MPI_MASTER , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , &
                     & Psi3b )
                  IF ( mpiRank == MPI_MASTER ) THEN 

                     CALL io_write_vtk_impsi ( 'psi-' , psiFileNo , psiFilePos , mpiSource , nX , nXa , nXb , nXbc , dNx , nY , &
                        & nYa , nYb , nYbc , dNy , nZ , nZa , nZb , nZbc , dNz , Psi3b )

                  END IF

               END DO
               psiFileNo = psiFileNo + 1

            ELSE IF ( psiOutput == 3 ) THEN ! Use simple parallel I/O to write out the complete wave function simultaneously to many
                                            ! legacy VTK files that only contain a single, one-dimensional slab of the wave function
                                            ! overseen by an MPI process

               CALL io_write_vtk ( 'psi-', psiFileNo, mpiRank, nX, nXa,&
                  & nXb, nXbc, dNx, nY, nYa, nYb, nYbc, dNy, nZ, nZa,  &
                  & nZb, nZbc, dNz, Xa, Ya, Za, Psi3a )
               psiFileNo = psiFileNo + 1

            ELSE IF ( psiOutput == 4 ) THEN ! Use MPI-I/O to write the wave function out to a single binary file in parallel

               WRITE ( UNIT = fileUnitChar, FMT = '(I4.4)' ) psiFileNo
               CALL MPI_FILE_OPEN ( MPI_COMM_WORLD, &
                  & TRIM ( 'psi-'//fileUnitChar//'.bin' ), &
                  & MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, &
                  & mpiFileHandle, mpiError )
               DO l = nZa, nZb
                  DO k = nYa, nYb
                     mpiOffset = 2*CMPLX_DEFAULT_KIND*nX*((k-1)+nY*(l-1))
                     CALL MPI_FILE_WRITE_AT ( mpiFileHandle, mpiOffset, &
                        & Psi3a(nXa,k,l), nX, mpiCmplx, mpiStatus, &
                        & mpiError )
                  ENDDO
               ENDDO
               CALL MPI_FILE_CLOSE ( mpiFileHandle, mpiError )
               psiFileNo = psiFileNo + 1

            ELSE IF ( psiOutput == 5 ) THEN ! Compute the integrated, two-dimensional density profiles of the wave function in the 
                                            ! (x,y)-, (x,z)-, and (y,z)-planes and then write each density profile to its own 
                                            ! gnuplot (.splot) formatted file 

               ALLOCATE ( Psi2a ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc ) )
               ALLOCATE ( Psi2b ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc ) )

               Psi3b = Psi3a
               Psi2a = 0.0
               Psi2b = 0.0
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        Psi2b ( j , k ) = Psi2b ( j , k ) + ABS ( Psi3b ( j , k , l ) )**2

                     END DO

                  END DO

               END DO
               Psi2b = Psi2b * dZ
               CALL MPI_REDUCE ( Psi2b , Psi2a , ( nXb - nXa + 1 ) * ( nYb - nYa + 1 ), mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
 
               IF ( mpiRank == MPI_MASTER ) THEN

                  CALL io_write_splot_psi2 ( 'psi2xy-' , psiFileNo , nXa , nXb , nXbc , nYa , nYb , nYbc , Xa , Ya , Psi2a )

               END IF

               DEALLOCATE ( Psi2b )
               DEALLOCATE ( Psi2a )

               ALLOCATE ( Psi2a ( nXa - nXbc : nXb + nXbc , 1 - nZbc : nZ + nZbc ) )
               ALLOCATE ( Psi2b ( nXa - nXbc : nXb + nXbc , 1 - nZbc : nZ + nZbc ) )

               Psi3b = Psi3a
               Psi2a = 0.0
               Psi2b = 0.0
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        Psi2b ( j , l ) = Psi2b ( j , l ) + ABS ( Psi3b ( j , k , l ) )**2 

                     END DO

                  END DO

               END DO
               Psi2b = Psi2b * dY
               CALL MPI_REDUCE ( Psi2b , Psi2a , ( nXb - nXa + 1 ) * nZ , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
 
               IF ( mpiRank == MPI_MASTER ) THEN 

                  CALL io_write_splot_psi2 ( 'psi2xz-' , psiFileNo , nXa , nXb , nXbc , 1 , nZ , nZbc , Xa , Zc , Psi2a )

               END IF

               DEALLOCATE ( Psi2b )
               DEALLOCATE ( Psi2a )

               ALLOCATE ( Psi2a ( nYa - nYbc : nYb + nYbc , 1 - nZbc : nZ + nZbc ) )
               ALLOCATE ( Psi2b ( nYa - nYbc : nYb + nYbc , 1 - nZbc : nZ + nZbc ) )

               Psi3b = Psi3a
               Psi2a = 0.0
               Psi2b = 0.0
               DO l = nZa , nZb

                  DO k = nYa , nYb

                     DO j = nXa , nXb

                        Psi2b ( k , l ) = Psi2b ( k , l ) + ABS ( Psi3b ( j , k , l ) )**2 

                     END DO

                  END DO

               END DO
               Psi2b = Psi2b * dX 
               CALL MPI_REDUCE ( Psi2b , Psi2a , ( nYb - nYa + 1 ) * nZ , mpiReal , MPI_SUM , MPI_MASTER , MPI_COMM_WORLD , mpiError )
 
               IF ( mpiRank == MPI_MASTER ) THEN 

                  CALL io_write_splot_psi2 ( 'psi2yz-' , psiFileNo , nYa , nYb , nYbc , 1 , nZ , nZbc , Ya , Zc , Psi2a )

               END IF

               DEALLOCATE ( Psi2b )
               DEALLOCATE ( Psi2a )

               psiFileNo = psiFileNo + 1

            END IF

         END IF

         CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

!        Compute 1st stage of generalized 4th-order Runge-Kutta (GRK4L): k_1 = f ( t_n , y_n )
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOrrf , yOrrf , zOrrf , &
            & dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3a , K1 )
!        Compute the intermediate wave function for the 2nd stage of the GRK4L method: y_n + 0.5 * dT * k_1
!        Note that the intermediate wave function is only computed on the interior grid points assigned to each MPI process.
         CALL grk4_y_3d_stgx ( 1 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , &
            & Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Compute V ( x , t_n + dT / 2 ), W ( t_n + dT / 2)
         tN = t0 + ( REAL ( n ) + 0.5 ) * dT

! ----------------------------------------------------------------------------------------------------------------------------------
! What time dependence do we want to explore with Omega ( t )? 
!
!         Omega ( 1 ) = wXo ! for now, always do absolute rotation from initial state
!         Omega ( 2 ) = wYo
!         Omega ( 3 ) = wZo
!
!    ( 1 ) Spin up: Omega ( t = 0 ) = ( 0, 0 , 0 ) ----> Omega ( t ) = ( wXf , wYf , wZf ) tanh?
!
!    ( 2 ) Nutation: rocking back and forth, rotating about x-axis with a periodic twist
!
!         thetaXo = PI / 6.0 ! nutation amplitude
!         nu = 0.2 ! nutation frequency
!         thetaX = thetaXo * SIN ( nu * tN ) 
!         R = rot_rx ( thetaX )
!
!    ( 3 ) Compton flip: 
!
!         tSigma = 10.0 ! flip time
!         sigma = 1.0 ! flip rate
!         thetaX = 0.5 * PI * ( ( TANH ( sigma * ( tN - tSigma ) ) / TANH ( 0.5 * sigma * ( tF - t0 ) ) ) + 1.0 )
!         R = rot_rx ( thetaX )
!
!         Omega = MATMUL ( R , Omega )
!         wX = Omega ( 1 )
!         wY = Omega ( 2 )
!         wZ = Omega ( 3 ) 
! ----------------------------------------------------------------------------------------------------------------------------------

!        Compute 2nd stage of GRK4: 
!        k_2 = f ( t_n + 0.5 * dT , y_n + 0.5 * dT * k_1 )
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOrrf , yOrrf , zOrrf , &
            & dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3b , K2 )
!        Compute intermediate wave function for 3rd stage of GRK4: 
!        y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2
         CALL grk4_y_3d_stgx ( 2 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , &
            & Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Compute stage three of GRK4:
!        k_3 = f ( t_n + 0.5 * dT , y_n + ( 1 / 2 - 1 / lambda ) * dT * k_1 + ( 1 / lambda ) * dT * k_2
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOrrf , yOrrf , zOrrf , &
            & dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3b , K3 )
!        Compute intermediate wave function for 4th stage of GRK4: 
!        y_n + ( 1 - lambda / 2 ) * dT * k_2 + ( lambda / 2 ) * dT * k_3
         CALL grk4_y_3d_stgx ( 3 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , &
            & Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Calculate V ( x , t_n + dT ), W ( t_n + dT )
         tN = t0 + REAL ( n + 1 ) * dT

! ----------------------------------------------------------------------------------------------------------------------------------
! What time dependence do we want to explore with Omega ( t )? 
!
!         Omega ( 1 ) = wXo ! for now, always do absolute rotation from initial state
!         Omega ( 2 ) = wYo
!         Omega ( 3 ) = wZo
!
!    ( 1 ) Spin up: Omega ( t = 0 ) = ( 0, 0 , 0 ) ----> Omega ( t ) = ( wXf , wYf , wZf ) tanh?
!
!    ( 2 ) Nutation: rocking back and forth, rotating about x-axis with a periodic twist
!
!         thetaX = thetaXo * SIN ( nu * tN ) 
!         R = rot_rx ( thetaX )
!
!    ( 3 ) Compton flip:
!
!         thetaX = 0.5 * PI * ( ( TANH ( sigma * ( tN - tSigma ) ) / TANH ( 0.5 * sigma * ( tF - t0 ) ) ) + 1.0 )
!         R = rot_rx ( thetaX )
!
!         Omega = MATMUL ( R , Omega )
!         wX = Omega ( 1 )
!         wY = Omega ( 2 )
!         wZ = Omega ( 3 )  
! ----------------------------------------------------------------------------------------------------------------------------------

!        Compute the fourth stage of GRK4:
!        k_4 = f ( t_n + dT , y_n + ( 1 - lamda / 2 ) * dT * k_2 + ( lamda / 2) * dT * k_3 
         CALL grk4_f_gp_3d_rrf_cdx ( fdOrder , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , xOrrf , yOrrf , zOrrf , &
            & dX , dY , dZ , wX , wY , wZ , gS , Xa , Ya , Za , Vex3a , Psi3b , K4 )
!        Compute wave function at nth+1 time step:
!        y_{ n + 1 } = y_n + ( dT / 6 ) * [ k_1 + ( 4 - lambda ) * k_2 + lambda * k_3 + k_4 ]
         CALL grk4_y_3d_stgx ( 4 , rk4Lambda , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dTz , K1 , K2 , K3 , K4 , &
            & Psi3a , Psi3b )
!        Exchange boundary condition information among nearest neighbor processes
         CALL mpi_exchange_ghosts ( nX , nXa , nXb , nXbc , nY , nYa , nYb , nYbc , nZ , nZa , nZb , nZbc , Psi3b )
!        Update wave function each time step: y_{ n + 1 } ---> y_n        
         Psi3a = Psi3b
!        If running in ITP mode, then also renormalize condensate wave function each time step
         IF ( itpOn .EQV. .TRUE. ) THEN

             CALL evua_normalize ( MPI_MASTER , mpiReal , mpiError , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , dX , &
                & dY , dZ , Psi3a )

         END IF

      END DO

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      ! Write last wave function to a binary file using MPI/IO
      CALL MPI_FILE_OPEN ( MPI_COMM_WORLD, 'psi-501.bin', MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, mpiError )
      DO l = nZa, nZb
         DO k = nYa, nYb
            mpiOffset = 2*CMPLX_DEFAULT_KIND*nX*((k-1)+nY*(l-1))
            CALL MPI_FILE_WRITE_AT ( mpiFileHandle, mpiOffset, Psi3a(nXa,k,l), nX, mpiCmplx, mpiStatus, mpiError )
         ENDDO
      ENDDO
      CALL MPI_FILE_CLOSE ( mpiFileHandle, mpiError )

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

! --- END MAIN TIME PROPAGATION LOOP / BEGIN CLEAN UP TO STOP ----------------------------------------------------------------------

      CALL DATE_AND_TIME ( stopDate , stopTime , stopZone , StopValues )
      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     RUN STOPPED @ ', stopTime, ' ON ', stopDate, ' ... '

      END IF

      DEALLOCATE ( Psi3b )
      DEALLOCATE ( Psi3a )
      DEALLOCATE ( K4 )
      DEALLOCATE ( K3 )
      DEALLOCATE ( K2 )
      DEALLOCATE ( K1 )
      DEALLOCATE ( Vex3a )
      DEALLOCATE ( Omega )
      DEALLOCATE ( Zc )
      DEALLOCATE ( Zb )
      DEALLOCATE ( Za )
      DEALLOCATE ( Ya )
      DEALLOCATE ( Xa )
      DEALLOCATE ( MpiStatus )
      DEALLOCATE ( StopValues )
      DEALLOCATE ( StartValues )
      
      CALL MPI_FINALIZE ( mpiError )

      STOP

! --- END MAIN PROGRAM -------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE gpse_read_stdin ( fileName , fileUnit )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName
            
      INTEGER, INTENT ( IN ) :: fileUnit

      OPEN  ( UNIT = fileUnit , FILE = fileName , ACTION = 'READ' , FORM = 'FORMATTED' , STATUS = 'OLD' )

         READ  ( UNIT = fileUnit , NML = gpseIn )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         Marty Kandes, Ph.D.'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#     COPYRIGHT'
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#'     
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#         Copyright (c) 2014, 2015, 2016, 2017 Martin Charles Kandes'
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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        tF        = ', tF
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xO        = ', xO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yO        = ', yO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zO        = ', zO
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dT        = ', dT
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dX        = ', dX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dY        = ', dY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        dZ        = ', dZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xOrrf     = ', xOrrf
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yOrrf     = ', yOrrf
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zOrrf     = ', zOrrf
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wX        = ', wX
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wY        = ', wY
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZ        = ', wZ
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        gS        = ', gS
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiInput  = ', psiInput
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiOutput = ', psiOutput
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiFileNo = ', psiFileNo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        psiInit   = ', psiInit
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nXpsi     = ', nXpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nYpsi     = ', nYpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nZpsi     = ', nZpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        nRpsi     = ', nRpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        mLpsi     = ', mLpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        xOpsi     = ', xOpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        yOpsi     = ', yOpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        zOpsi     = ', zOpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        rOpsi     = ', rOpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wXpsi     = ', wXpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wYpsi     = ', wYpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wZpsi     = ', wZpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        wRpsi     = ', wRpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        pXpsi     = ', pXpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        pYpsi     = ', pYpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        pZpsi     = ', pZpsi
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexInput  = ', vexInput
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexOutput = ', vexOutput
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexFileNo = ', vexFileNo
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '#        vexInit   = ', vexInit
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
         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) '# 1 : tN , 2 : wX , 3 : wY , 4 : wZ , 5 : evuaL2Norm , 6 : evuaEa , 7 : evuaEb , &
            & 8 : evuaMuA , 9 : evuaMuB , 10 : evuaL2a , 11 : evuaL2b , 12 : evuaTx , 13 : evuaTy , 14 : evuaTz , 15 : evuaVex , &
            & 16 : evuaVmf , 17 : evuaX , 18 : evuaX2a , 19 : evuaSigXa , 20 : evuaX2b , 21 : evuaSigXb , 22 : evuaPx , &
            & 23 : evuaPx2 , 24 : evuaSigPx , 25 : evuaLxA , 26 : evuaLx2a , 27 : evuaSigLxA , 28 : evuaLxB , 29 : evuaLx2b , &
            & 30 : evuaSigLxB , 31 : evuaFx , 32 : evuaTauXa , 33 : evuaTauXb , 34 : evuaY , 35 : evuaY2a , 36 : evuaSigYa , &
            & 37 : evuaY2b , 38 : evuaSigYb , 39 : evuaPy , 40 : evuaPy2 , 41 : evuaSigPy , 42 : evuaLyA , 43 : evuaLy2a , &
            & 44 : evuaSigLyA , 45 : evuaLyB , 46 : evuaLy2b , 47 : evuaSigLyB , 48 : evuaFy , 49 : evuaTauYa , 50 : evuaTauYb , &
            & 51 : evuaZ , 52 : evuaZ2a , 53 : evuaSigZa , 54 : evuaZ2b , 55 : evuaSigZb , 56 : evuaPz , 57 : evuaPz2 , &
            & 58 : evuaSigPz , 59 : evuaLzA , 60 : evuaLz2a , 61 : evuaSigLzA , 62 : evuaLzB , 63 : evuaLz2b , 64 : evuaSigLzB , &
            & 65 : evuaFz , 66 : evuaTauZa , 67 : evuaTauZb , 68 : evuaR , 69 : evuaR2a , 70 : evuaR2b , 71 : evuaIxxA , &
            & 72 : evuaIxyA , 73 : evuaIxzA , 74 : evuaIyyA , 75 : evuaIyzA , 76 : evuaIzzA , 77 : evuaIxxB , 78 : evuaIxyB , &
            & 79 : evuaIxzB , 80 : evuaIyyB , 81 : evuaIyzB , 82 : evuaIzzB'

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

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

! ----------------------------------------------------------------------------------------------------------------------------------

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
      CALL MPI_BCAST ( nTsteps   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nTwrite   , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nX        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nXbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dNx       , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nY        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nYbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dNy       , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nZ        , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( nZbc      , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( dNz       , 1 , mpiInt      , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( t0        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
      CALL MPI_BCAST ( tF        , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
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
      CALL MPI_BCAST ( rOpsi     , 1 , mpiReal     , mpiMaster , MPI_COMM_WORLD , mpiError )
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

! ----------------------------------------------------------------------------------------------------------------------------------

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

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE mpi_copy_vex ( mpiRank , mpiSource , mpiDestination , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 ) 

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

               END DO ! possible alternative: no loop over j; call MPI_RECV ( Real3b ( nXa , nYa , l ) , nXb - nXa + 1 , ... )

            END DO

         END DO

      ELSE IF ( ( mpiRank == mpiSource ) .AND. ( mpiRank /= mpiDestination ) ) THEN

         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  CALL MPI_SSEND ( Vex3 ( j , k , l ) , 1 , mpiReal , mpiDestination , mpiRank , MPI_COMM_WORLD , mpiError )

               END DO ! possible alternative: no loop over j; call MPI_SSEND ( Real3a ( nXa , nYa , l ) , nXb - nXa + 1 , ... )

            END DO

         END DO

      END IF

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      RETURN
 
      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE mpi_copy_psi ( mpiRank , mpiSource , mpiDestination , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

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

               END DO ! possible alternative: no loop over j; call MPI_RECV ( Real3b ( nXa , nYa , l ) , nXb - nXa + 1 , ... )

            END DO

         END DO

      ELSE IF ( ( mpiRank == mpiSource ) .AND. ( mpiRank /= mpiDestination ) ) THEN

         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  CALL MPI_SSEND ( Psi3 ( j , k , l ) , 1 , mpiCmplx , mpiDestination , mpiRank , MPI_COMM_WORLD , mpiError )

               END DO ! possible alternative: no loop over j; call MPI_SSEND ( Real3a ( nXa , nYa , l ) , nXb - nXa + 1 , ... )

            END DO

         END DO

      END IF

      CALL MPI_BARRIER ( MPI_COMM_WORLD , mpiError )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

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

            CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 0 , MPI_COMM_WORLD, &
               & mpiError )

         END IF

         IF ( mpiRank /= MPI_MASTER ) THEN 

            CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 1 , MPI_COMM_WORLD , &
               & MpiStatus , mpiError )

         END IF

      ELSE 

         CALL MPI_RECV ( Psi3 ( nXa , nYa , nZa - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 0 , MPI_COMM_WORLD , &
            & MpiStatus , mpiError )

         IF ( mpiRank /= mpiProcesses - 1 ) THEN 

            CALL MPI_SEND ( Psi3 ( nXa , nYa , nZb + 1 - nZbc ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 1 , MPI_COMM_WORLD , &
               & mpiError )

         END IF

      END IF

      IF ( MODULO ( mpiRank , 2 ) == 0 ) THEN

         IF ( mpiRank /= MPI_MASTER ) THEN

            CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 2 , MPI_COMM_WORLD , mpiError )

         END IF

         IF ( mpiRank /= mpiProcesses - 1 ) THEN

            CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 3 , MPI_COMM_WORLD , &
               & MpiStatus , mpiError )

         END IF

      ELSE

         IF ( mpiRank /= mpiProcesses - 1 ) THEN

            CALL MPI_RECV ( Psi3 ( nXa , nYa , nZb + 1 ) , nX * nY * nZbc , mpiCmplx , mpiRank + 1 , 2 , MPI_COMM_WORLD , &
               & MpiStatus , mpiError )

         END IF

         CALL MPI_SEND ( Psi3 ( nXa , nYa , nZa ) , nX * nY * nZbc , mpiCmplx , mpiRank - 1 , 3 , MPI_COMM_WORLD , mpiError )

      END IF

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      END PROGRAM

! ==================================================================================================================================
