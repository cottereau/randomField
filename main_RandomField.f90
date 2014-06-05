program main_RandomField

	use randomField1D
	use statistics_RF

    implicit none

    integer :: Nmc, modelDimension;
    logical :: randInit;
    double precision, dimension(:), allocatable :: xMax, kMax, corrL;
    double precision, dimension(:), allocatable :: average, stdDeviation, averageCorrL;
    integer, dimension(:), allocatable :: xNStep, kNStep;
    character (len=20) :: corrMod;
    double precision, dimension(:, :), allocatable :: randField;
    double precision :: pi = 3.1415926535898;
    integer :: i;

    write(*,*) "------------START main-----------------------------------------";

	!todo Read File --------------------------------------------------------------------------------------

	modelDimension = 1; !Number of dimensions
	randInit = .TRUE. !if it's set to false each event will have the same random numbers
	Nmc = 2; !number of Monte-Carlo experiments

	!Allocating
	allocate(xMax(Nmc));
	allocate(xNStep(Nmc));
	allocate(kMax(Nmc));
	allocate(kNStep(Nmc));
	allocate(corrL(Nmc));
	allocate(average(Nmc));
	allocate(stdDeviation(Nmc));
	allocate(averageCorrL(Nmc));

    !Variables Init (will disapear once turned all tests an "Read File" made)
    xMax   = (/(2*pi    , i=1, Nmc)/);
    xNStep = (/(4       , i=1, Nmc)/);
    kMax   = (/(2*pi    , i=1, Nmc)/);
    kNStep = (/(3       , i=1, Nmc)/);
    corrL  = (/(1       , i=1, Nmc)/); !correlation length
    corrMod = "gaussian";

	!Selecting problem dimension-----------------------------------------------------------------------
    if(modelDimension.eq.1) then
    	call createRandomField1D(Nmc, randInit, xMax, kMax, xNStep, kNStep, corrL, corrMod, randField);
    	call set_Statistics1D(randField, average, stdDeviation, averageCorrL)
    else if(modelDimension.eq.2) then
    	!call createRandomField2D
    else if(modelDimension.eq.3) then
    	!call createRandomField3D
    else
    	stop "model dimension not accepted"
    endif

	call Disp2D(randField, "randField");
	call Disp1D(average, "average");
	call Disp1D(stdDeviation, "stdDeviation");
	call Disp1D(averageCorrL, "averageCorrL");

	!todo Write Result File--------------------------------------------------------------------------------------


	!Deallocating
	deallocate(xMax);
	deallocate(xNStep);
	deallocate(kMax);
	deallocate(kNStep);
	deallocate(corrL);

	deallocate(randField);

	deallocate(average);
	deallocate(stdDeviation);
	deallocate(averageCorrL);

    write(*,*) "------------END main-----------------------";
contains

end program main_RandomField
