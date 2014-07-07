program main_RandomField

	use randomFieldND
	use statistics_RF
	use writeResultFile_RF

    implicit none

    !In each quantity above we have: quantity (number of Monte Carlo Experiments, dimension)

    !INPUTS
    integer                                          :: Nmc, modelDimension;
    logical                                          :: randInit;
    character (len=30)                               :: fileName;
    character (len=15), dimension(:),    allocatable :: corrMod; !Only one corrMod for each Monte Carlo Experiment (is the same in every dimension)
    double precision,   dimension(:, :), allocatable :: corrL, xMax, xPeriod;

	!OUTPUTS (in this case variables to be filled in in the proces)
	double precision, dimension(:, :), allocatable :: kMax;
	integer, dimension(:, :),          allocatable :: xNStep, kNStep;
    double precision, dimension(:, :), allocatable :: randField;
    double precision, dimension(:),    allocatable :: average, stdDeviation, averageCorrL;

	!Local variables
    double precision :: pi = 3.1415926535898;
    double precision :: startTime = 0, endTime = 0;
    integer :: i;

	call CPU_TIME(startTime)

	write(*,*) "";
    write(*,*) "------------START main-----------------------------------------";
    write(*,*) "";

	!todo Read File --------------------------------------------------------------------------------------

	Nmc = 1000; !number of Monte-Carlo experiments
	modelDimension = 2; !Number of dimensions
	randInit = .TRUE. !if it's set to false each event will have the same random numbers
	fileName = "TestResult2DModel"

    !Allocating
	allocate(corrMod(Nmc));
	allocate(corrL  (modelDimension, Nmc));
	allocate(xMax   (modelDimension, Nmc));
	allocate(xPeriod(modelDimension, Nmc));


	write(*,*) ">>>>>>>>> Variables initialization: corrL, xMax and xPeriod";
	!ONLY FOR TESTS (should disapear once turned all tests an "Read File" made)----------------------
	    corrMod(:)   = "gaussian";

	    corrL  (1, :) = (/(3,       i=1, Nmc)/);
	    xMax   (1, :) = (/(2*pi,    i=1, Nmc)/);
	    xPeriod(1, :) = 10*corrL(1, :);

	if(modelDimension > 1) then
	    corrL  (2, :) = (/(1,       i=1, Nmc)/);
	    xMax   (2, :) = (/(2*pi,    i=1, Nmc)/);
		xPeriod(2, :) = 10*corrL(2, :);
	end if

	if(modelDimension > 2) then
	    corrL  (3, :) = (/(1,       i=1, Nmc)/);
	    xMax   (3, :) = (/(2*pi,    i=1, Nmc)/);
		xPeriod(3, :) = 10*corrL(3, :);
	end if

	if(modelDimension > 3) then
	    stop "model dimension not accepted"
	end if

!	call Disp2D(corrL,   "corrL"  );
!	call Disp2D(xMax,    "xMax"   );
!	call Disp2D(xPeriod, "xPeriod");

	call createRandomFieldND(Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
    						     kMax, xNStep, kNStep, randField);

	call set_StatisticsND   (randField, xNStep, xMax, average, stdDeviation, averageCorrL)

	call CPU_TIME(endTime)

	call writeResultFileND  (corrMod, corrL, xMax, xPeriod,   &
    						 kMax, xNStep, kNStep, randField, &
    						 average, stdDeviation, averageCorrL, fileName, endTime - startTime)

	write(*,*) ">>>>>>>>> Showing Results";
	write(*,*) "";
	write(*,*) "     Global Average = ", calculateAverage(average);
	write(*,*) "Global StdDeviation = ", calculateAverage(stdDeviation);
	write(*,*) " Correlation Length = ", averageCorrL;
	write(*,*) "     Total Time (s) = ", endTime - startTime;

	!Deallocating
		deallocate(corrMod);
		deallocate(corrL);
		deallocate(xMax);
		deallocate(xPeriod);

		deallocate(kMax);
		deallocate(xNStep);
		deallocate(kNStep);
		deallocate(randField);

		deallocate(average);
		deallocate(stdDeviation);
		deallocate(averageCorrL);

	write(*,*) "";
    write(*,*) "------------END main-----------------------";
	write(*,*) "";

end program main_RandomField
