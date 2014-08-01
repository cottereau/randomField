program main_RandomField

	use randomFieldND
	use statistics_RF
	use writeResultFile_RF
	use readFile_RF
	use displayCarvalhol
	use mpi

    implicit none

    !In each quantity above we have: quantity (number of Monte Carlo Experiments, dimension)

    !INPUTS
    integer                                        :: Nmc, nDim;
    character (len=30), parameter                  :: inputName = "input01", outputName = "Test1";
    character (len=15)                             :: corrMod;
    double precision,   dimension(:),  allocatable :: corrL, xMax, xPeriod;

	!OUTPUTS (in this case variables to be filled in in the proces)
	double precision, dimension(:),    allocatable :: kMax;
	integer         , dimension(:),    allocatable :: xNStep, kNStep;
    double precision, dimension(:, :), allocatable :: randField, totalRandField;
    double precision, dimension(:),    allocatable :: evntAvg, evntStdDev, procCorrL;
    double precision, dimension(:),    allocatable :: ptAvg, ptStdDev, globalCorrL;
    double precision                               :: globalAvg, globalStdDev;
    double precision, dimension(:),    allocatable :: compAvg, compStdDev, compCorrL;

	!Local variables
    integer          :: i, xNStepTotal;
    integer          :: code, rang, error, nb_procs;
    double precision :: pi = 3.1415926535898;
    character(len=30), dimension(:,:), allocatable :: dataTable
    double precision , dimension(:)  , allocatable :: sumRF, sumRFsquare, &
                                                      totalSumRF, totalSumRFsquare;

	call MPI_INIT(code)

	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

	if(rang == 0) then
	    write(*,*) "------------START MPI processing-----------------------------------------";
		write(*,*) "";
	    write(*,*) "------------START main-----------------------------------------";
	    write(*,*) "";
		write(*,*) "Number of procs = ", nb_procs
		write(*,*) ">>>>>>>>> Reading file input";
	end if

	!Reading Input
	call set_DataTable(inputName, dataTable)
	if(rang == 0) call Disp2Dchar (dataTable, "dataTable");

	call read_DataTable(dataTable, "Nmc",  Nmc)
	call read_DataTable(dataTable, "corrMod", corrMod)
	call read_DataTable(dataTable, "corrL", corrL)
	call read_DataTable(dataTable, "xMax", xMax)
	call read_DataTable(dataTable, "xPeriod", xPeriod)

	deallocate(dataTable)

	nDim        = size(xMax)

	if(rang == 0) then
		write(*,*) ">>>>>>>>INPUT Read Data"
		write(*,*) "Nmc            = ", Nmc !number of Monte-Carlo experiments
		write(*,*) "corrMod        = ", corrMod
		write(*,*) "corrL          = ", corrL
		write(*,*) "xMax           = ", xMax
		write(*,*) "xPeriod        = ", xPeriod
		write(*,*) "nDim           = ", nDim
	end if


	if(Nmc < 1) then
		write(*,*) "Number of events should be a positive integer"
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if
	if(nDim < 1) then
		write(*,*) "nDim should be a positive integer"
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if

	!Calculating a part of the random field in each processor
	if(rang == 0) write(*,*) ">>>>>>>>> Calculating Random field";
	call createRandomFieldND(Nmc, corrMod, corrL, xMax, xPeriod,   &
    						     kMax, xNStep, kNStep, randField);

	xNStepTotal = product(xNStep) !should be after "createRandomFieldND"

	write (*,*) "Proc = ", rang
	!call DispCarvalhol(randField, "randField")

	!Calculating Statistics
	if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics";
	call set_StatisticsND_MPI(randField, xNStep, xMax, rang,                  &
							  sumRF, sumRFsquare, ptAvg, ptStdDev, procCorrL, &
							  totalSumRF, totalSumRFsquare,                   &
							  evntAvg, evntStdDev,                            &
							  globalAvg, globalStdDev, globalCorrL)

	!Writing results and comparison statistics (the whole random field in one processor)
	call gather_RF_MPI(randField, totalRandField, xNStepTotal, rang) !Gathering all random fields in proc 0

	if(rang == 0) then
		!call DispCarvalhol(totalRandField, "totalRandField")
!		call set_StatisticsND(totalRandField, xNStep, xMax,        &
!					          compAvg, compStdDev, compCorrL)

		!call DispCarvalhol(totalRandField, "totalRandField")

!		call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
!						 kMax, xNStep, kNStep, totalRandField, &
!						 compAvg, compStdDev, compCorrL, outputName)

!		call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
!	    						 kMax, xNStep, kNStep, totalRandField, &
!	    						 compAvg, compStdDev, compCorrL, outputName)
!
!	    call writeResultHDF5(xMax, kMax, xNStep, totalRandField, outputName)



!		write(*,*) ">>>>>>>>> Showing Results";
!		write(*,*) "";
!		write(*,*) "     Global Average = ", calculateAverage(average);
!		write(*,*) "Global StdDeviation = ", calculateAverage(stdDeviation);
!		if (Nmc == 1) write(*,*) " Correlation Length not computed (Nmc = 1) = ", averageCorrL;
!		if (Nmc >  1) write(*,*) " Correlation Length = ", averageCorrL;
!		write(*,*) "     Total Time (s) = ", endTime - startTime;
	end if

	!Deallocating
		if(allocated(corrL))   deallocate(corrL);
		if(allocated(xMax))    deallocate(xMax);
		if(allocated(xPeriod)) deallocate(xPeriod);

		if(allocated(kMax))      deallocate(kMax);
		if(allocated(xNStep))    deallocate(xNStep);
		if(allocated(kNStep))    deallocate(kNStep);
		if(allocated(randField)) deallocate(randField);

		if(allocated(evntAvg))      deallocate(evntAvg);
		if(allocated(evntStdDev))   deallocate(evntStdDev);
		if(allocated(procCorrL))    deallocate(procCorrL);
		if(allocated(ptAvg))        deallocate(ptAvg);
		if(allocated(ptStdDev))     deallocate(ptStdDev);
		if(allocated(globalCorrL))  deallocate(globalCorrL);
		if(allocated(sumRF))        deallocate(sumRF)
		if(allocated(sumRFsquare))  deallocate(sumRFsquare)
		if(allocated(totalSumRF))        deallocate(totalSumRF)
		if(allocated(totalSumRFsquare))  deallocate(totalSumRFsquare)
		if(allocated(totalRandField))    deallocate(totalRandField);

	if(rang == 0) then
		write(*,*) "";
	    write(*,*) "------------END main-----------------------";
		write(*,*) "";
	end if

	call MPI_FINALIZE(code)

end program main_RandomField
