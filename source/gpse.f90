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
!     Sunday, March 16th, 2014
!
! -------------------------------------------------------------------------

      PROGRAM GPSE

      USE, INTRINSIC :: ISO_FORTRAN_ENV
      USE            :: GRID
      USE            :: IO
      USE            :: PSI
      USE            :: VEX
      
      IMPLICIT NONE

      INCLUDE 'mpif.h'

      CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.0.9'

      INTEGER, PARAMETER :: MPI_MASTER = 0

      INTEGER :: mpiDest
      INTEGER :: mpiError
      INTEGER :: mpiErrorCode
      INTEGER :: mpiErrorClass
      INTEGER :: mpiProcesses
      INTEGER :: mpiProvided
      INTEGER :: mpiRank
      INTEGER :: mpiSource
      INTEGER :: ompThreads
      INTEGER :: ompThreadID
      INTEGER :: j , k , l

      REAL :: t = 0.0

      INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIStatus

      REAL, ALLOCATABLE, DIMENSION ( :             ) :: Vex1
      REAL, ALLOCATABLE, DIMENSION ( : , :         ) :: Vex2
      REAL, ALLOCATABLE, DIMENSION ( : , : , :     ) :: Vex3
      REAL, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: Vex4
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: RePsi1
      REAL, ALLOCATABLE, DIMENSION ( : , :         ) :: RePsi2
      REAL, ALLOCATABLE, DIMENSION ( : , : , :     ) :: RePsi3
      REAL, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: RePsi4
      REAL, ALLOCATABLE, DIMENSION ( :             ) :: ImPsi1
      REAL, ALLOCATABLE, DIMENSION ( : , :         ) :: ImPsi2
      REAL, ALLOCATABLE, DIMENSION ( : , : , :     ) :: ImPsi3
      REAL, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: ImPsi4

      COMPLEX, ALLOCATABLE, DIMENSION ( :             ) :: Psi1
      COMPLEX, ALLOCATABLE, DIMENSION ( : , :         ) :: Psi2
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , :     ) :: Psi3
      COMPLEX, ALLOCATABLE, DIMENSION ( : , : , : , : ) :: Psi4

      INTEGER :: OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM

      nX = 128
      nY = 128
      nZ = 128
      xO = 0.0
      yO = 0.0
      zO = 0.0
      dX = 0.125
      dY = 0.125
      dZ = 0.125

      ALLOCATE ( MPIStatus ( MPI_STATUS_SIZE ) ) 

      CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiProvided , mpiError )
      IF ( mpiError /= MPI_SUCCESS ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'gpse: ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
         CALL MPI_ABORT ( MPI_COMM_WORLD , mpiErrorCode , mpiError )

      END IF
      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiProcesses , mpiError )
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRank , mpiError )

      WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'Hello from' , mpiRank

      IF ( mpiRank == MPI_MASTER ) THEN

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'I AM MASTER!'

         ALLOCATE ( X      ( nX                     ) )
         ALLOCATE ( Y      ( nY                     ) )
         ALLOCATE ( Z      ( nZ                     ) )
         ALLOCATE ( Vex1   ( nX * nY * nZ           ) )
         ALLOCATE ( Vex3   ( nX           , nY , nZ ) )
         ALLOCATE ( RePsi1 ( nX * nY * nZ           ) )
         ALLOCATE ( ImPsi1 ( nX * nY * nZ           ) )
         ALLOCATE ( Psi1   ( nX * nY * nZ           ) )
         ALLOCATE ( Psi3   ( nX           , nY , nZ ) )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'MASTER HAS ALLOCATED!'

         CALL regular_grid_axis ( nX , xO , dX , X )
         CALL regular_grid_axis ( nY , yO , dY , Y )
         CALL regular_grid_axis ( nZ , zO , dZ , Z )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'MASTER HAS GENERATED THE GRID!'

         DO l = 1 , nZ

            DO k = 1 , nY

               DO j = 1 , nX

                  Vex3 ( j , k , l ) = &
                     & vex_3d_shor ( X ( j ) , Y ( k ) , Z ( l ) , 0.0 , 0.0 , 0.0 , 5.0 , 1.0 , 10.0 )

                  Psi3 ( j , k , l ) = &
                     & psi_3d_se_sho_ani ( 0 , 0 , 0 , X ( j ) , Y ( k ) , Z ( l ) , 1.0 , 1.0 , 0.0 , 1.0 , 1.0 , 5.0 )

               END DO

            END DO

         END DO

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'MASTER HAS BUILT INITIAL WAVEFUNCTION!'

         Vex1 = PACK ( Vex3 , .TRUE. )
         Psi1 = PACK ( Psi3 , .TRUE. )
         RePsi1 = REAL ( Psi1 )
         ImPsi1 = AIMAG ( Psi1 )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'Hello from MASTER!'

         CALL write_vtk ( t , Vex1 , RePsi1 , ImPsi1 )

         WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'Goodbye from MASTER!'

         DEALLOCATE ( Psi3   )
         DEALLOCATE ( Psi1   )
         DEALLOCATE ( ImPsi1 )
         DEALLOCATE ( RePsi1 )
         DEALLOCATE ( Vex1   )
         DEALLOCATE ( Z      )
         DEALLOCATE ( Y      )
         DEALLOCATE ( X      )

      END IF

      WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'Goodbye from' , mpiRank

      CALL MPI_FINALIZE ( mpiError )

      DEALLOCATE ( MPIStatus )

      STOP

      END PROGRAM

! =========================================================================
