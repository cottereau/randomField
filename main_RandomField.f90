program main_RandomField

	use randomField1D
	use randomFieldND
	use statistics_RF
	use writeResultFile_RF

    implicit none

    !In each quantity above we have: quantity (number of Monte Carlo Experiments, dimension)

    !INPUTS
    integer :: Nmc, modelDimension;
    logical :: randInit;
    character (len=30) :: fileName;
    character (len=15), dimension(:), allocatable :: corrMod; !Only one corrMod for each Monte Carlo Experiment (is the same in every dimension)
    double precision, dimension(:, :), allocatable :: corrL, xMax, xPeriod;

	!OUTPUTS (in this case variables to be filled in in the proces)
	double precision, dimension(:, :), allocatable :: kMax;
	integer, dimension(:, :), allocatable          :: xNStep, kNStep;
    double precision, dimension(:, :), allocatable :: randField;
    double precision, dimension(:, :), allocatable :: average, stdDeviation, averageCorrL;

	!Local variables
    double precision :: pi = 3.1415926535898;
    integer :: i, truc;

    write(*,*) "------------START main-----------------------------------------";

	!todo Read File --------------------------------------------------------------------------------------

	Nmc = 3; !number of Monte-Carlo experiments
	modelDimension = 1; !Number of dimensions
	randInit = .TRUE. !if it's set to false each event will have the same random numbers
	fileName = "TestResultNDModel"
	truc = 2

	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	!>>>>>>>>>>>>>> 1D >>>>>>>>>>>>>>>>>>>>>>>>>
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(modelDimension.eq.1) then

    	!Allocating
		allocate(corrMod(Nmc));
		allocate(corrL  (modelDimension, Nmc));
		allocate(xMax   (modelDimension, Nmc));
		allocate(xPeriod(modelDimension, Nmc));

		allocate(average     (1, Nmc));
		allocate(stdDeviation(1, Nmc));
		allocate(averageCorrL(1, Nmc));

	    !Variables Init (should disapear once turned all tests an "Read File" made)
	    corrMod(:)   = "gaussian";
	    corrL  (1, :) = (/(0.5,     i=1, Nmc)/);
	    xMax   (1, :) = (/(2*pi,    i=1, Nmc)/);
	    xPeriod(1, :) = (/(5,       i=1, Nmc)/);

		call Disp2D(corrL, "corrL");
		call Disp2D(xMax, "xMax");
		call Disp2D(xPeriod, "xPeriod");

    	call createRandomFieldND(Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
    						     kMax, xNStep, kNStep, randField);
	   	call set_StatisticsND(randField, xNStep, kNStep, average, stdDeviation, averageCorrL)
	   	call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
    						   kMax, xNStep, kNStep, randField, &
    						   average, stdDeviation, averageCorrL, fileName)

!	  !>>>>>>>>>>>OLD MODEL START
!	   	!Allocating
!		allocate(corrMod(Nmc));
!		allocate(corrL(Nmc,modelDimension));
!		allocate(xMax(Nmc,modelDimension));
!		allocate(xPeriod(Nmc,modelDimension));
!
!		allocate(average(modelDimension, Nmc));
!		allocate(stdDeviation(modelDimension, Nmc));
!		allocate(averageCorrL(modelDimension, Nmc));
!
!	    !Variables Init (should disapear once turned all tests an "Read File" made)
!	    corrMod(:)   = "gaussian";
!	    corrL(:,1)   = (/(0.5    , i=1, Nmc)/); !correlation length
!	    xMax(:,1)    = (/(2*pi   , i=1, Nmc)/);
!	    xPeriod(:,1) = (/(5      , i=1, Nmc)/);
!
!	   	call createRandomField1D(Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
!	   						     kMax, xNStep, kNStep, randField);
!	   	call set_Statistics1D(randField, xNStep, kNStep, average, stdDeviation, averageCorrL)
!
!	   	call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
!	   						   kMax, xNStep, kNStep, randField, &
!	   						   average, stdDeviation, averageCorrL, fileName)
!	  !>>>>>>>>>>>OLD MODEL END


	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	!>>>>>>>>>>>>>> 2D >>>>>>>>>>>>>>>>>>>>>>>>>
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else if(modelDimension.eq.2) then

    	!Allocating
		allocate(corrMod(Nmc));
		allocate(corrL  (modelDimension, Nmc));
		allocate(xMax   (modelDimension, Nmc));
		allocate(xPeriod(modelDimension, Nmc));

		allocate(average     (modelDimension, 1));
		allocate(stdDeviation(modelDimension, 1));
		allocate(averageCorrL(modelDimension, 1));

	    !Variables Init (should disapear once turned all tests an "Read File" made)
	    corrMod(:)   = "gaussian";
	    corrL(:,1)   = (/(1,       i=1, Nmc)/);
	    corrL(:,2)   = (/(2,       i=1, Nmc)/);
	    xMax(:,1)    = (/(2*pi,    i=1, Nmc)/);
	    xMax(:,2)    = (/(2*pi,    i=1, Nmc)/);
	    xPeriod(:,1) = (/(5,       i=1, Nmc)/);
		xPeriod(:,2) = (/(5,       i=1, Nmc)/);

		call Disp2D(corrL, "corrL");
		call Disp2D(xMax, "xMax");
		call Disp2D(xPeriod, "xPeriod");

    	call createRandomFieldND(Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
    						     kMax, xNStep, kNStep, randField);
	   	call set_StatisticsND(randField, xNStep, kNStep, average, stdDeviation, averageCorrL)
	   	call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
    						   kMax, xNStep, kNStep, randField, &
    						   average, stdDeviation, averageCorrL, fileName)

	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	!>>>>>>>>>>>>>> 3D >>>>>>>>>>>>>>>>>>>>>>>>>
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else if(modelDimension.eq.3) then

    	!Allocating
		allocate(corrMod(Nmc));
		allocate(corrL  (modelDimension, Nmc));
		allocate(xMax   (modelDimension, Nmc));
		allocate(xPeriod(modelDimension, Nmc));

		allocate(average     (1, Nmc));
		allocate(stdDeviation(1, Nmc));
		allocate(averageCorrL(1, Nmc));

	    !Variables Init (should disapear once turned all tests an "Read File" made)
	    corrMod(:)   = "gaussian";
	    corrL  (1, :) = (/(1,       i=1, Nmc)/);
	    corrL  (2, :) = (/(2,       i=1, Nmc)/);
	    corrL  (3, :) = (/(1,       i=1, Nmc)/);
	    xMax   (1, :) = (/(2*pi,    i=1, Nmc)/);
	    xMax   (2, :) = (/(2*pi,    i=1, Nmc)/);
	    xMax   (3, :) = (/(2*pi,    i=1, Nmc)/);
	    xPeriod(1, :) = (/(5,       i=1, Nmc)/);
		xPeriod(2, :) = (/(5,       i=1, Nmc)/);
		xPeriod(3, :) = (/(5,       i=1, Nmc)/);

		call Disp2D(corrL, "corrL");
		call Disp2D(xMax, "xMax");
		call Disp2D(xPeriod, "xPeriod");

    	call createRandomFieldND(Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
    						     kMax, xNStep, kNStep, randField);
	   	call set_StatisticsND(randField, xNStep, kNStep, average, stdDeviation, averageCorrL)
	   	call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
    						   kMax, xNStep, kNStep, randField, &
    						   average, stdDeviation, averageCorrL, fileName)
    else
    	stop "model dimension not accepted"
    endif

	call Disp2D(randField, "randField");
	call Disp2D(average, "average");
	call Disp2D(stdDeviation, "stdDeviation");
	call Disp2D(averageCorrL, "averageCorrL");

	!todo Write Result File--------------------------------------------------------------------------------------


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

    write(*,*) "------------END main-----------------------";
contains

end program main_RandomField
