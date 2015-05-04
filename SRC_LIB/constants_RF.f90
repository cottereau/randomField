module constants_RF

    implicit none

	double precision, parameter :: PI = 3.1415926535898d0;
	double precision, parameter :: SQRT_2PI = 2.50662827463d0;
	double precision, parameter :: ZERO = 0d0;
	double precision, parameter :: TOLERANCE = 0.0000000000001
	double precision, parameter :: MIN_DOUBLE = -2.0D+307
	double precision, parameter :: MAX_DOUBLE = 2.0D+307
	integer, parameter :: ISOTROPIC = 1, &
	                      SHINOZUKA = 2, &
	                      RANDOMIZATION = 3

	integer :: TESTRANK = 0 !ONLY FOR TESTS


end module constants_RF
