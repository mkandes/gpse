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
!     Copyright (c) 2014 Martin Charles Kandes
!
! LAST UPDATED
!
!     Thursday, January 9th, 2014
!
! -------------------------------------------------------------------------

      MODULE MATH

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      CHARACTER ( LEN = * ), PARAMETER :: VERSION_NUMBER = '0.0.2'

      REAL, PARAMETER, PUBLIC :: PI = 3.14159265358979323846264338327950288

      PUBLIC :: factorial
!      PUBLIC :: alaguerre
!      PUBLIC :: hermite
!      PUBLIC :: laguerre
!     PUBLIC :: gamma

      CONTAINS

         INTEGER RECURSIVE FUNCTION factorial ( n ) RESULT ( nFactorial )
   
            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: n

            nFactorial = 0

            IF ( KIND ( n ) == INT32 ) THEN

               IF ( n == 0 ) THEN

                  nFactorial = 1

               ELSE IF ( n >= 1 .AND. n <= 12 ) THEN

                  nFactorial = n * factorial ( n - 1 )

               ELSE

                  nFactorial = -1
                  !WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'math : factorial :: ERROR - n must be an integer greater than or equal to 0, but less than or equal to 12 because KIND ( n ) = INT32'
                  STOP

               END IF

            ELSE IF ( KIND ( n ) == INT64 ) THEN

               IF ( n == 0 ) THEN

                  nFactorial = 1 

               ELSE IF ( n >= 1 .AND. n <= 20 ) THEN

                  nFactorial = n * factorial ( n - 1 ) 

               ELSE

                  nFactorial = -1
                  !WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'math : factorial :: ERROR - n must be an integer greater than or equal to 0, but less than or equal to 20 because KIND ( n ) = INT64'
                  STOP

               END IF

            ELSE

               nFactorial = -1
               !WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'math : factorial :: ERROR - KIND ( n ) != INT32 .OR. INT64'
               STOP

            END IF

            RETURN

         END FUNCTION

      END MODULE

! =========================================================================
