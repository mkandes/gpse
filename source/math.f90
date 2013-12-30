! ==========================================================================
! NAME
!
!     math [ maTH ] - Math Module
!
! SYNOPSIS
!
! DESCRIPTION  
!
!     math is a Fortran module ... 
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
!     Copyright (c) 2013 Martin Charles Kandes
!
! LAST UPDATED
!
!     Sunday, December 29th, 2013
!
! -------------------------------------------------------------------------

      MODULE MATH

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      CHARACTER ( LEN = * ), PARAMETER :: VERSION_NUMBER = '0.0.1'

      REAL, PARAMETER, PUBLIC :: PI = 3.14159265358979323846264338327950288

      PUBLIC :: factorial
      PUBLIC :: alaguerre
      PUBLIC :: hermite
      PUBLIC :: laguerre

      CONTAINS

         INTEGER RECURSIVE FUNCTION factorial ( n ) RESULT ( nFactorial ) 
   
            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: n

         END FUNCTION

      END MODULE MATH

! =========================================================================
