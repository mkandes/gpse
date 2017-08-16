! ==================================================================================================================================
! NAME
!
!     io [ io ] - Input-Output Module
!
! SYNOPSIS
!
!     USE :: IO
!
! DESCRIPTION  
!
!     IO is a custom Fortran module written to perform input-output operations for gpse.
!
! OPTIONS
!
! SEE ALSO
!
! BUGS
!
! HISTORY
!
! AUTHOR
!
!     Marty Kandes, Ph.D.
!     Computational & Data Science Research Specialist
!     User Services Group
!     San Diego Supercomputer Center
!     University of California, San Diego
!
! COPYRIGHT
!     
!     Copyright (c) 2014, 2015, 2016, 2017 Martin Charles Kandes
!
! LAST UPDATED
!
!     Wednesday, August 16th, 2017
!
! ----------------------------------------------------------------------------------------------------------------------------------

      MODULE IO

! --- MODULE DECLARATIONS ----------------------------------------------------------------------------------------------------------

      USE, INTRINSIC :: ISO_FORTRAN_ENV

! --- MODULE DEFINITIONS -----------------------------------------------------------------------------------------------------------
!
!     ISO_FORTRAN_ENV is the intrinsic Fortran module that provides information about the run-time environment.
!
! ----------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      PRIVATE

! --- SUBROUTINE DECLARATIONS ------------------------------------------------------------------------------------------------------

      PUBLIC :: io_read_bin_psi
      PUBLIC :: io_read_bin_vex
      PUBLIC :: io_write_bin_psi
      PUBLIC :: io_write_bin_vex
      PUBLIC :: io_write_splot_psi2
      PUBLIC :: io_write_vtk_header
      PUBLIC :: io_write_vtk_xcoordinates
      PUBLIC :: io_write_vtk_ycoordinates
      PUBLIC :: io_write_vtk_zcoordinates
      PUBLIC :: io_write_vtk_repsi
      PUBLIC :: io_write_vtk_impsi
      PUBLIC :: io_write_vtk

! --- SUBROUTINE DEFINITIONS -------------------------------------------------------------------------------------------------------
!
!     io_read_bin_psi is a PUBLIC SUBROUTINE that performs unformatted stream I/O to read a wave function from a binary file.
!
!     io_read_bin_vex is a PUBLIC SUBROUTINE that performs unformatted stream I/O to read an external potential from an binary file.
!
!     io_write_bin_psi is a PUBLIC SUBROUTINE that performs unformatted stream I/O to write a wave function to a binary file.
!
!     io_write_bin_vex is a PUBLIC SUBROUTINE that performs unformatted stream I/O to write an external potential to a binary file.
!
!     io_write_vtk_header is a PUBLIC SUBROUTINE that performs formatted stream I/O to write the standard legacy header of a vtk 
!        file for an output wave function.
!
!     io_write_vtk_xcoordinates is a PUBLIC SUBROUTINE that performs formatted stream I/O to write the x-axis coordinates of the 
!        grid points to a vtk file.
!
!     io_write_vtk_ycoordinates is a PUBLIC SUBROUTINE that performs formatted stream I/O to write the y-axis coordinates of the 
!        grid points to a vtk file.
!
!     io_write_vtk_zcoordinates is a PUBLIC SUBROUTINE that performs formatted stream I/O to write the z-axis coordinates of the 
!        grid points to a vtk file.
!
!     io_write_vtk_repsi is a PUBLIC SUBROUTINE that performs formatted stream I/O to write the real component of an output wave 
!        function to a vtk file.
!
!     io_write_vtk_impsi is a PUBLIC SUBROUTINE that performs formatted stream I/O to write the imaginary component of an output 
!        wave function to a vtk file.
!
!     io_write_vtk is a PUBLIC SUBROUTINE that performs formatted record-based I/O to write a wave function to a vtk file. This 
!         subroutine is primarily used to implement simple parallel I/O. 
!
! ----------------------------------------------------------------------------------------------------------------------------------

      CONTAINS

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_read_bin_psi ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN ) :: nXa 
      INTEGER, INTENT ( IN ) :: nXb 
      INTEGER, INTENT ( IN ) :: nXbc 
      INTEGER, INTENT ( IN ) :: nYa 
      INTEGER, INTENT ( IN ) :: nYb 
      INTEGER, INTENT ( IN ) :: nYbc 
      INTEGER, INTENT ( IN ) :: nZa 
      INTEGER, INTENT ( IN ) :: nZb 
      INTEGER, INTENT ( IN ) :: nZbc 

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( INOUT ) :: Psi3

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k , l

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'READ' , &
         & FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )

         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  READ ( UNIT = fileUnit , POS = filePosition ) Psi3 ( j , k , l )
                  INQUIRE ( UNIT = fileUnit , POS = filePosition )

               END DO

            END DO

         END DO

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_read_bin_vex ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
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

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k , l 

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'READ' , &
         & FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )

         DO l = nZa , nZb 

            DO k = nYa , nYb 

               DO j = nXa , nXb 

                  READ ( UNIT = fileUnit , POS = filePosition ) Vex3 ( j , k , l ) 
                  INQUIRE ( UNIT = fileUnit , POS = filePosition )

               END DO

            END DO

         END DO

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_bin_psi ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Psi3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
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

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k , l

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN ( UNIT = fileUnit , FILE = TRIM ( 'psi-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
            
         DO l = nZa , nZb

            DO k = nYa , nYb

               DO j = nXa , nXb

                  WRITE ( UNIT = fileUnit , POS = filePosition ) Psi3 ( j , k , l )
                  INQUIRE ( UNIT = fileUnit , POS = filePosition )

               END DO

            END DO

         END DO

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_bin_vex ( fileUnit , filePosition , nXa , nXb , nXbc , nYa , nYb , nYbc , nZa , nZb , nZbc , Vex3 )

      IMPLICIT NONE

      INTEGER, INTENT ( IN    ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN    ) :: nXa 
      INTEGER, INTENT ( IN    ) :: nXb 
      INTEGER, INTENT ( IN    ) :: nXbc 
      INTEGER, INTENT ( IN    ) :: nYa 
      INTEGER, INTENT ( IN    ) :: nYb 
      INTEGER, INTENT ( IN    ) :: nYbc 
      INTEGER, INTENT ( IN    ) :: nZa 
      INTEGER, INTENT ( IN    ) :: nZb 
      INTEGER, INTENT ( IN    ) :: nZbc 

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Vex3

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k , l 

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN ( UNIT = fileUnit , FILE = TRIM ( 'vex-'//fileUnitChar//'.bin' ) , ACCESS = 'STREAM' , ACTION = 'WRITE' , & 
         & FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' )
   
         DO l = nZa , nZb 

            DO k = nYa , nYb 

               DO j = nXa , nXb

                  WRITE ( UNIT = fileUnit , POS = filePosition ) Vex3 ( j , k , l )
                  INQUIRE ( UNIT = fileUnit , POS = filePosition )

               END DO

            END DO

         END DO

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_splot_psi2 ( fileName, fileUnit , nXa , nXb , nXbc , nYa , nYb , nYbc , X , Y , Psi2 )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Psi2

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.splot' )  , ACCESS = 'SEQUENTIAL' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

         DO k = nYa , nYb

            DO j = nXa , nXb

               WRITE ( UNIT = fileUnit , FMT = * ) Y ( k ) , X ( j ) , Psi2 ( j , k )

            END DO
            WRITE ( UNIT = fileUnit , FMT = * ) ' '

         END DO

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk_header ( fileName , fileUnit , filePosition , nX , nY , nZ , dNx , dNy , dNz )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nZ
      INTEGER, INTENT ( IN ) :: dNx
      INTEGER, INTENT ( IN ) :: dNy
      INTEGER, INTENT ( IN ) :: dNz

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

         WRITE ( UNIT = fileUnit , FMT = '(A26)' ) '# vtk DataFile Version 3.0' ! VTK file identifier and version number
         WRITE ( UNIT = fileUnit , FMT = '(A26)' ) 'STANDARD LEGACY VTK FORMAT' ! VTK file header; 256 characters maximum
         WRITE ( UNIT = fileUnit , FMT = '(A5)'  ) 'ASCII'                      ! VTK file format: ASCII or BINARY
         WRITE ( UNIT = fileUnit , FMT = '(A24)' ) 'DATASET RECTILINEAR_GRID'   ! VTK dataset structure
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX / dNx + MODULO ( nX , dNx ), &
            & nY / dNy + MODULO ( nY , dNy ), nZ / dNz + MODULO ( nZ , dNz )
         INQUIRE ( UNIT = fileUnit , POS = filePosition )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE
! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk_xcoordinates ( fileName , fileUnit , filePosition , nX , nXa , nXb , nXbc , dNx , X )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: dNx
     
      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

         WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A13,1X,I4.1,1X,A6)' ) 'X_COORDINATES' , nX / dNx + MODULO ( nX , dNx ) , 'double'
         DO j = nXa , nXb, dNx

            WRITE ( UNIT = fileUnit , FMT = * ) X ( j )

         END DO
         INQUIRE ( UNIT = fileUnit , POS = filePosition )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk_ycoordinates ( fileName , fileUnit , filePosition , nY , nYa , nYb , nYbc , dNy , Y )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files

      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: dNy

      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: k

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

         WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY / dNy + MODULO ( nY , dNy ) , 'double'
         DO k = nYa , nYb , dNy

            WRITE ( UNIT = fileUnit , FMT = * ) Y ( k )

         END DO
         INQUIRE ( UNIT = fileUnit , POS = filePosition )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk_zcoordinates ( fileName , fileUnit , filePosition , mpiSource , nZ , nZa , nZb , nZbc , dNz , Z )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ) , INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN ) :: mpiSource
      INTEGER, INTENT ( IN ) :: nZ
      INTEGER, INTENT ( IN ) :: nZa 
      INTEGER, INTENT ( IN ) :: nZb 
      INTEGER, INTENT ( IN ) :: nZbc
      INTEGER, INTENT ( IN ) :: dNz

      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: l

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

         IF ( mpiSource == 0 ) THEN ! i.e., if mpiSource == MPI_MASTER

            WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , nZ / dNz + MODULO ( nZ , dNz ) , 'double'
            DO l = nZa , nZb , dNz

               WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

            END DO

         ELSE

            WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(F23.15)' ) Z ( nZa )
            DO l = nZa + dNz , nZb , dNz

               WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

            END DO

         END IF
         INQUIRE ( UNIT = fileUnit , POS = filePosition )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk_repsi ( fileName , fileUnit , filePosition , mpiSource , nX , nXa , nXb , nXbc , dNx , nY , nYa , &
         & nYb , nYbc , dNy , nZ , nZa , nZb , nZbc , dNz , Psi3 )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN ) :: mpiSource
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: dNx
      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: dNy
      INTEGER, INTENT ( IN ) :: nZ
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc
      INTEGER, INTENT ( IN ) :: dNz

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k , l

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )
 
         IF ( mpiSource == 0 ) THEN

            WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , ( nX / dNx + MODULO ( nX , dNx ) ) * ( nY / dNy + MODULO ( nY , dNy ) ) * ( nZ / dNz + MODULO ( nZ , dNz ) ) 
            WRITE ( UNIT = fileUnit , FMT = '(A24)'          ) 'SCALARS RePsi double 1'
            WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'LOOKUP_TABLE default'
            DO l = nZa , nZb , dNz

               DO k = nYa , nYb , dNy 

                  DO j = nXa , nXb , dNx

                     WRITE ( UNIT = fileUnit , FMT = * ) REAL ( Psi3 ( j , k , l ) )

                  END DO

               END DO

            END DO

         ELSE

            DO l = nZa , nZb , dNz

               DO k = nYa , nYb , dNy

                  DO j = nXa , nXb , dNx
                        
                     WRITE ( UNIT = fileUnit , POS = filePosition , FMT = * ) REAL ( Psi3 ( j , k , l ) )
                     INQUIRE ( UNIT = fileUnit , POS = filePosition )

                  END DO

               END DO

            END DO

         END IF
         INQUIRE ( UNIT = fileUnit , POS = filePosition )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk_impsi ( fileName , fileUnit , filePosition , mpiSource , nX , nXa , nXb , nXbc , dNx , nY , nYa , &
         & nYb , nYbc , dNy , nZ , nZa , nZb , nZbc , dNz , Psi3 )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER ( KIND = 8 ), INTENT ( INOUT ) :: filePosition ! hard coded KIND = 8 to avoid integer overflow when writing large files
      INTEGER, INTENT ( IN ) :: mpiSource
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: dNx
      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: dNy
      INTEGER, INTENT ( IN ) :: nZ
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc
      INTEGER, INTENT ( IN ) :: dNz

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

      CHARACTER ( LEN = 4 ) :: fileUnitChar

      INTEGER :: j , k , l

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'.vtk' )  , ACCESS = 'STREAM' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'UNKNOWN' )

         IF ( mpiSource == 0 ) THEN

            WRITE ( UNIT = fileUnit , POS = filePosition , FMT = '(A24)' ) 'SCALARS ImPsi double 1'
            WRITE ( UNIT = fileUnit , FMT = '(A20)' ) 'LOOKUP_TABLE default'

            DO l = nZa , nZb , dNz

               DO k = nYa , nYb , dNy

                  DO j = nXa , nXb , dNx

                     WRITE ( UNIT = fileUnit , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )

                  END DO

               END DO

            END DO

         ELSE

            DO l = nZa , nZb , dNz

               DO k = nYa , nYb , dNy

                  DO j = nXa , nXb , dNx

                     WRITE ( UNIT = fileUnit , POS = filePosition , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )
                     INQUIRE ( UNIT = fileUnit , POS = filePosition )

                  END DO

               END DO

            END DO

         END IF
         INQUIRE ( UNIT = fileUnit , POS = filePosition )

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE io_write_vtk ( fileName , fileUnit , mpiRank , nX , nXa , nXb , nXbc , dNx , nY , nYa , &
         & nYb , nYbc , dNy , nZ , nZa , nZb , nZbc , dNz , X , Y , Z , Psi3 )

      IMPLICIT NONE

      CHARACTER ( LEN = * ), INTENT ( IN ) :: fileName

      INTEGER, INTENT ( IN ) :: fileUnit
      INTEGER, INTENT ( IN ) :: mpiRank
      INTEGER, INTENT ( IN ) :: nX
      INTEGER, INTENT ( IN ) :: nXa
      INTEGER, INTENT ( IN ) :: nXb
      INTEGER, INTENT ( IN ) :: nXbc
      INTEGER, INTENT ( IN ) :: dNx
      INTEGER, INTENT ( IN ) :: nY
      INTEGER, INTENT ( IN ) :: nYa
      INTEGER, INTENT ( IN ) :: nYb
      INTEGER, INTENT ( IN ) :: nYbc
      INTEGER, INTENT ( IN ) :: dNy
      INTEGER, INTENT ( IN ) :: nZ
      INTEGER, INTENT ( IN ) :: nZa
      INTEGER, INTENT ( IN ) :: nZb
      INTEGER, INTENT ( IN ) :: nZbc
      INTEGER, INTENT ( IN ) :: dNz

      REAL, DIMENSION ( nXa - nXbc : nXb + nXbc ), INTENT ( IN ) :: X
      REAL, DIMENSION ( nYa - nYbc : nYb + nYbc ), INTENT ( IN ) :: Y
      REAL, DIMENSION ( nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Z

      COMPLEX, DIMENSION ( nXa - nXbc : nXb + nXbc , nYa - nYbc : nYb + nYbc , nZa - nZbc : nZb + nZbc ), INTENT ( IN ) :: Psi3

      CHARACTER ( LEN = 4 ) :: fileUnitChar
      CHARACTER ( LEN = 4 ) :: mpiRankChar

      INTEGER :: j , k , l

      WRITE ( UNIT = fileUnitChar , FMT = '(I4.4)' ) fileUnit
      WRITE ( UNIT = mpiRankChar , FMT = '(I4.4)' ) mpiRank
      OPEN  ( UNIT = fileUnit , FILE = TRIM ( fileName//fileUnitChar//'-'//mpiRankChar//'.vtk' )  , ACCESS = 'SEQUENTIAL' , ACTION = 'WRITE' , &
         & FORM = 'FORMATTED' , STATUS = 'NEW' )

         WRITE ( UNIT = fileUnit , FMT = '(A26)' ) '# vtk DataFile Version 3.0' ! VTK file identifier and version number
         WRITE ( UNIT = fileUnit , FMT = '(A26)' ) 'STANDARD LEGACY VTK FORMAT' ! VTK file header; 256 characters maximum
         WRITE ( UNIT = fileUnit , FMT = '(A5)'  ) 'ASCII'                      ! VTK file format: ASCII or BINARY
         WRITE ( UNIT = fileUnit , FMT = '(A24)' ) 'DATASET RECTILINEAR_GRID'   ! VTK dataset structure
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I4.1,1X,I4.1,1X,I4.1)' ) 'DIMENSIONS' , nX / dNx + MODULO ( nX , dNx ), &
            & nY / dNy + MODULO ( nY , dNy ), ( nZb - nZa + 1 ) / dNz + MODULO ( nZ , dNz )
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'X_COORDINATES' , nX / dNx + MODULO ( nX , dNx ) , 'double'
         DO j = nXa , nXb, dNx

            WRITE ( UNIT = fileUnit , FMT = * ) X ( j )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Y_COORDINATES' , nY / dNy + MODULO ( nY , dNy ) , 'double'
         DO k = nYa , nYb , dNy

            WRITE ( UNIT = fileUnit , FMT = * ) Y ( k )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A13,1X,I4.1,1X,A6)' ) 'Z_COORDINATES' , ( nZb - nZa + 1 ) / dNz + MODULO ( nZb - nZb + 1 , dNz ) , 'double'
         DO l = nZa , nZb , dNz

            WRITE ( UNIT = fileUnit , FMT = * ) Z ( l )

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A10,1X,I19.1)' ) 'POINT_DATA' , ( nX / dNx + MODULO ( nX , dNx ) ) * ( nY / dNy + MODULO ( nY , dNy ) ) * ( ( nZb - nZa + 1 ) / dNz + MODULO ( nZb - nZa + 1 , dNz ) ) 
         WRITE ( UNIT = fileUnit , FMT = '(A24)'          ) 'SCALARS RePsi double 1'
         WRITE ( UNIT = fileUnit , FMT = '(A20)'          ) 'LOOKUP_TABLE default'
         DO l = nZa , nZb , dNz

            DO k = nYa , nYb , dNy 

               DO j = nXa , nXb , dNx

                  WRITE ( UNIT = fileUnit , FMT = * ) REAL ( Psi3 ( j , k , l ) )

               END DO

            END DO

         END DO
         WRITE ( UNIT = fileUnit , FMT = '(A24)' ) 'SCALARS ImPsi double 1'
         WRITE ( UNIT = fileUnit , FMT = '(A20)' ) 'LOOKUP_TABLE default'
         DO l = nZa , nZb , dNz

            DO k = nYa , nYb , dNy

               DO j = nXa , nXb , dNx

                  WRITE ( UNIT = fileUnit , FMT = * ) AIMAG ( Psi3 ( j , k , l ) )

               END DO

            END DO

         END DO

      CLOSE ( UNIT = fileUnit , STATUS = 'KEEP' )

      RETURN

      END SUBROUTINE

! ----------------------------------------------------------------------------------------------------------------------------------

      END MODULE

! ==================================================================================================================================
