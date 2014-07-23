program main_RandomField

	use randomFieldND
	use statistics_RF
	use writeResultFile_RF
	use readFile_RF
	use mpi

    implicit none

    !In each quantity above we have: quantity (number of Monte Carlo Experiments, dimension)

    !INPUTS
    integer                                        :: Nmc, modelDimension;
    character (len=30)                             :: inputName, outputName;
    character (len=15)                             :: corrMod;
    double precision,   dimension(:),  allocatable :: corrL, xMax, xPeriod;

	!OUTPUTS (in this case variables to be filled in in the proces)
	double precision, dimension(:),    allocatable :: kMax;
	integer, dimension(:),             allocatable :: xNStep, kNStep;
    double precision, dimension(:, :), allocatable :: randField;
    double precision, dimension(:),    allocatable :: average, stdDeviation, averageCorrL;

	!Local variables
    double precision :: pi = 3.1415926535898;
    double precision :: startTime = 0, endTime = 0;
    integer :: i;
    character(len=30), dimension(:,:), allocatable :: dataTable

	call CPU_TIME(startTime)

	write(*,*) "";
    write(*,*) "------------START main-----------------------------------------";
    write(*,*) "";

	write(*,*) ">>>>>>>>> Reading file input";

	inputName  = "input01"
	outputName = "Test1"

	call set_DataTable(inputName, dataTable)
	call Disp2Dchar (dataTable, "dataTable");

	call read_DataTable(dataTable, "Nmc",  Nmc)
	call read_DataTable(dataTable, "corrMod", corrMod)
	call read_DataTable(dataTable, "corrL", corrL)
	call read_DataTable(dataTable, "xMax", xMax)
	call read_DataTable(dataTable, "xPeriod", xPeriod)

	deallocate(dataTable)

	modelDimension = size(xMax)

	write(*,*) ">>>>>>>>INPUT Data"
	write(*,*) "Nmc            = ", Nmc !number of Monte-Carlo experiments
	write(*,*) "corrMod        = ", corrMod
	write(*,*) "corrL          = ", corrL
	write(*,*) "xMax           = ", xMax
	write(*,*) "xPeriod        = ", xPeriod
	write(*,*) "modelDimension = ", modelDimension


	if(Nmc < 1) stop "Number of events should be a positive integer"
	if(modelDimension < 1) stop "modelDimension should be a positive integer"

!	call Disp2D(corrL,   "corrL"  );
!	call Disp2D(xMax,    "xMax"   );
!	call Disp2D(xPeriod, "xPeriod");

	call createRandomFieldND(Nmc, corrMod, corrL, xMax, xPeriod, &
    						     kMax, xNStep, kNStep, randField);

	call set_StatisticsND(randField, xNStep, xMax, average, stdDeviation, averageCorrL)

	call CPU_TIME(endTime)

	call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
    						 kMax, xNStep, kNStep, randField, &
    						 average, stdDeviation, averageCorrL, outputName, endTime - startTime)

    call writeResultHDF5(xMax, kMax, xNStep, randField, outputName)

	call CPU_TIME(endTime)

	write(*,*) ">>>>>>>>> Showing Results";
	write(*,*) "";
	write(*,*) "     Global Average = ", calculateAverage(average);
	write(*,*) "Global StdDeviation = ", calculateAverage(stdDeviation);
	if (Nmc == 1) write(*,*) " Correlation Length not computed (Nmc = 1) = ", averageCorrL;
	if (Nmc >  1) write(*,*) " Correlation Length = ", averageCorrL;
	write(*,*) "     Total Time (s) = ", endTime - startTime;

	!Deallocating
		if(allocated(corrL))   deallocate(corrL);
		if(allocated(xMax))    deallocate(xMax);
		if(allocated(xPeriod)) deallocate(xPeriod);

		if(allocated(kMax))      deallocate(kMax);
		if(allocated(xNStep))    deallocate(xNStep);
		if(allocated(kNStep))    deallocate(kNStep);
		if(allocated(randField)) deallocate(randField);

		if(allocated(average))      deallocate(average);
		if(allocated(stdDeviation)) deallocate(stdDeviation);
		if(allocated(averageCorrL)) deallocate(averageCorrL);

	write(*,*) "";
    write(*,*) "------------END main-----------------------";
	write(*,*) "";

end program main_RandomField
