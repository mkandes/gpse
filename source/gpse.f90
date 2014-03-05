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
!     Wednesday, March 5th, 2014
!
! -------------------------------------------------------------------------

      PROGRAM GPSE

         USE, INTRINSIC :: ISO_FORTRAN_ENV
         USE GRID
         USE PSI

         IMPLICIT NONE

         INCLUDE 'mpif.h'

         CHARACTER ( LEN = * ), PARAMETER :: GPSE_VERSION_NUMBER = '0.0.8'

         INTEGER, PARAMETER :: I_AM_MASTER = 0

         INTEGER :: mpiDESTINATION
         INTEGER :: mpiERROR
         INTEGER :: mpiERRORCODE
         INTEGER :: mpiERRORCLASS
         INTEGER :: mpiNUMPROCS
         INTEGER :: mpiPROVIDED
         INTEGER :: mpiRANK
         INTEGER :: mpiSOURCE
         INTEGER :: ompNUMTHRDS
         INTEGER :: ompTHREADID

         INTEGER, ALLOCATABLE, DIMENSION ( : ) :: MPIStatus

         INTEGER :: OMP_GET_NUM_THREADS
         INTEGER :: OMP_GET_THREAD_NUM

         ALLOCATE ( MPIStatus ( MPI_STATUS_SIZE ) ) 

         CALL MPI_INIT_THREAD ( MPI_THREAD_SINGLE , mpiPROVIDED , mpiERROR )
         IF ( mpiERROR /= MPI_SUCCESS ) THEN

            WRITE ( UNIT = OUTPUT_UNIT , FMT = * ) 'gpse: ERROR - MPI_INIT_THREAD failed. Calling MPI_ABORT.'
            CALL MPI_ABORT ( MPI_COMM_WORLD , mpiERRORCODE , mpiERROR )

         END IF
         CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , mpiNUMPROCS , mpiERROR )
         CALL MPI_COMM_RANK ( MPI_COMM_WORLD , mpiRANK , mpiERROR )

         CALL MPI_FINALIZE ( mpiERROR ) 

         STOP

      END PROGRAM

! =========================================================================
