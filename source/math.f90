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
!     Saturday, November 22nd, 2014
!
! -------------------------------------------------------------------------

      MODULE MATH

      USE, INTRINSIC :: ISO_FORTRAN_ENV

      IMPLICIT NONE
      PRIVATE

      INTEGER                            , PARAMETER, PUBLIC :: INT_DEFAULT_KIND   = KIND ( 0 ) 
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: INT_8              = SELECTED_INT_KIND ( 2 ) 
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: INT_16             = SELECTED_INT_KIND ( 4 ) 
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: INT_32             = SELECTED_INT_KIND ( 9 ) 
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: INT_64             = SELECTED_INT_KIND ( 18 ) 
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: INT_128            = SELECTED_INT_KIND ( 38 ) 
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: REAL_DEFAULT_KIND  = KIND ( 0.0 )
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: REAL_32            = SELECTED_REAL_KIND ( 6  , 37   )    
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: REAL_64            = SELECTED_REAL_KIND ( 15 , 307  )
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: REAL_128           = SELECTED_REAL_KIND ( 33 , 4931 )
      INTEGER ( KIND = INT_DEFAULT_KIND ), PARAMETER, PUBLIC :: CMPLX_DEFAULT_KIND = KIND ( CMPLX ( 0.0 , 0.0 ) ) 

      REAL, PARAMETER, PUBLIC :: PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534

      COMPLEX, PARAMETER, PUBLIC :: I = CMPLX ( 0.0 , 1.0 )

      PUBLIC :: factorial
      PUBLIC :: alaguerre
      PUBLIC :: hermite
      PUBLIC :: laguerre
      PUBLIC :: legendre
      PUBLIC :: lambert
      PUBLIC :: bessel
      PUBLIC :: hankel
      PUBLIC :: chebyshev

      CONTAINS

         INTEGER RECURSIVE FUNCTION factorial ( n ) RESULT ( nFactorial )
   
            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: n

            nFactorial = 0

            IF ( KIND ( n ) == INT_32 ) THEN

               IF ( n == 0 ) THEN

                  nFactorial = 1

               ELSE IF ( n >= 1 .AND. n <= 12 ) THEN

                  nFactorial = n * factorial ( n - 1 )

               ELSE

                  nFactorial = -1
                  WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'math : factorial :: ERROR - n must be an integer greater than or equal to 0, but less than or equal to 12 because KIND ( n ) = INT32.'
                  STOP

               END IF

            ELSE IF ( KIND ( n ) == INT_64 ) THEN

               IF ( n == 0 ) THEN

                  nFactorial = 1 

               ELSE IF ( n >= 1 .AND. n <= 20 ) THEN

                  nFactorial = n * factorial ( n - 1 ) 

               ELSE

                  nFactorial = -1
                  WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'math : factorial :: ERROR - n must be an integer greater than or equal to 0, but less than or equal to 20 because KIND ( n ) = INT64.'
                  STOP

               END IF

            ELSE

               nFactorial = -1
               WRITE ( UNIT = ERROR_UNIT , FMT = * ) 'math : factorial :: ERROR - KIND ( n ) != INT32 .OR. INT64'
               STOP

            END IF

            RETURN

         END FUNCTION

         REAL RECURSIVE FUNCTION alaguerre ( n , k , x ) RESULT ( alaguerreNK ) 

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: n
            INTEGER, INTENT ( IN ) :: k

            REAL, INTENT ( IN ) :: x

            alaguerreNK = 0.0

            IF ( k > -1 ) THEN

               IF ( n == 0 ) THEN

                  alaguerreNK = 1.0

               ELSE IF ( n == 1 ) THEN

                  alaguerreNK = 1.0 + REAL ( k ) - x

               ELSE IF ( n >= 2 ) THEN

                  alaguerreNK = ( ( 2.0 * REAL ( n ) + REAL ( k ) - 1.0 - x ) * alaguerre ( n - 1 , k , x ) - & 
                     & ( REAL ( n ) + REAL ( k ) - 1.0 ) * alaguerre ( n - 2 , k , x ) ) / REAL ( n )

               ELSE

                  ! ERROR

               END IF
            
            ELSE

               ! ERROR 

            END IF

            RETURN 

         END FUNCTION

         REAL RECURSIVE FUNCTION hermite ( n , x ) RESULT ( hermiteN )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: n

            REAL, INTENT ( IN ) :: x

            hermiteN = 0.0

            IF ( n == 0 ) THEN

               hermiteN = 1.0

            ELSE IF ( n == 1 ) THEN

               hermiteN = 2.0 * x 

            ELSE IF ( n >= 2 ) THEN

               hermiteN = 2.0 * x * hermite ( n - 1 , x ) - 2.0 * REAL ( n - 1 ) * hermite ( n - 2 , x )    

            ELSE

               ! ERROR

            END IF

            RETURN

         END FUNCTION

         REAL RECURSIVE FUNCTION laguerre ( n , x ) RESULT ( laguerreN )

            IMPLICIT NONE

            INTEGER, INTENT ( IN ) :: n

            REAL, INTENT ( IN ) :: x

            laguerreN = 0.0

            IF ( n == 0 ) THEN

               laguerreN = 1.0

            ELSE IF ( n == 1 ) THEN

               laguerreN = 1.0 - x

            ELSE IF ( n >= 2 ) THEN

               laguerreN = ( ( 2.0 * REAL ( n ) - 1.0 - x ) * laguerre ( n - 1 , x ) - & 
                  & ( REAL ( n ) - 1.0 ) * laguerre ( n - 2 , x ) ) / REAL ( n )

            ELSE

               ! ERROR

            END IF

            RETURN

         END FUNCTION

         REAL FUNCTION legendre ( ) RESULT ( legendreN )

            IMPLICIT NONE

            legendreN = 0.0

            RETURN

         END FUNCTION

         REAL FUNCTION lambert ( ) RESULT ( lambertW )

            IMPLICIT NONE

            lambertW = 0.0

            RETURN

         END FUNCTION

         REAL FUNCTION bessel ( ) RESULT ( besselJ ) 

            IMPLICIT NONE

            besselJ = 0.0

            RETURN

         END FUNCTION

         REAL FUNCTION hankel ( ) RESULT ( hankelN )

            IMPLICIT NONE

            hankelN = 0.0

            RETURN

         END FUNCTION

         REAL FUNCTION chebyshev ( ) RESULT ( cheb )

            IMPLICIT NONE

            cheb = 0.0

            RETURN 

         END FUNCTION 

      END MODULE

! =========================================================================
