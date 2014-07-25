program main_RandomField

	use randomFieldND
	use statistics_RF
	use writeResultFile_RF
	use readFile_RF

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
    double precision, dimension(:, :), allocatable :: randField, totalRandField, tempRandField;
    double precision, dimension(:),    allocatable :: average, stdDeviation, averageCorrL;

	!Local variables
    double precision :: pi = 3.1415926535898;
    double precision :: startTime = 0, endTime = 0;
    integer :: i, xTotalnStep;
    integer :: code, rang, error, nb_procs, dblBitSize, newExtent;
    integer :: nombre_bloc, longueur_bloc, pas, ancien_type, type_RF, type_Temp
    character(len=30), dimension(:,:), allocatable :: dataTable

    !Test variables (to be deleted)
!    integer, parameter :: nb_valeurs=64, sqrtNb = 8;
!    integer            :: longueur_tranche, j, k
!    double precision, dimension(sqrtNb,sqrtNb) :: donnees
!    double precision, dimension(:,:), allocatable :: valeurs

	call CPU_TIME(startTime)

	call MPI_INIT(code)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

	if(rang == 0) then
		write(*,*) "";
	    write(*,*) "------------START main-----------------------------------------";
	    write(*,*) "";

		write(*,*) "Number of procs = ", nb_procs

		write(*,*) ">>>>>>>>> Reading file input";
	end if

	inputName  = "input01"
	outputName = "Test1"

	call set_DataTable(inputName, dataTable)
	if(rang == 0) call Disp2Dchar (dataTable, "dataTable");

	call read_DataTable(dataTable, "Nmc",  Nmc)
	call read_DataTable(dataTable, "corrMod", corrMod)
	call read_DataTable(dataTable, "corrL", corrL)
	call read_DataTable(dataTable, "xMax", xMax)
	call read_DataTable(dataTable, "xPeriod", xPeriod)

	deallocate(dataTable)

	modelDimension = size(xMax)

	if(rang == 0) then
		write(*,*) ">>>>>>>>INPUT Data"
		write(*,*) "Nmc            = ", Nmc !number of Monte-Carlo experiments
		write(*,*) "corrMod        = ", corrMod
		write(*,*) "corrL          = ", corrL
		write(*,*) "xMax           = ", xMax
		write(*,*) "xPeriod        = ", xPeriod
		write(*,*) "modelDimension = ", modelDimension
	end if


	if(Nmc < 1) then
		write(*,*) "Number of events should be a positive integer"
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if
	if(modelDimension < 1) then
		write(*,*) "modelDimension should be a positive integer"
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if

!	call Disp2D(corrL,   "corrL"  );
!	call Disp2D(xMax,    "xMax"   );
!	call Disp2D(xPeriod, "xPeriod");

	call createRandomFieldND(Nmc, corrMod, corrL, xMax, xPeriod, &
    						     kMax, xNStep, kNStep, randField);

!	    write(*,*) "Proc ", rang, " xStart = ", xStart
!		write(*,*) "Proc ", rang, "   xEnd = ", xEnd
!	    write(*,*) "Proc ", rang, "sizeUnif", sizeUnif
!	    write(*,*) "Proc ", rang, "sizeLoc", sizeLoc
!		call DispCarvalhol(randField, "randField")

	if(rang == 0) then
		allocate (tempRandField(nb_procs*size(randField,1), size(randField,2)))
		totalRandField = 0d0
	end if

	nombre_bloc = size(randField,2)
	longueur_bloc = size(randField,1)
	pas = size(randField,1)*nb_procs
	ancien_type = MPI_DOUBLE_PRECISION

!	if(rang == 0) then
!		write(*,*) "nombre_bloc = ",nombre_bloc
!		write(*,*) "longueur_bloc = ",longueur_bloc
!		write(*,*) "pas = ",pas
!	end if

	call MPI_TYPE_VECTOR(nombre_bloc, longueur_bloc, pas, ancien_type, type_Temp, code) !type_Temp and code are outputs
	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblBitSize, code)
	newExtent = size(randField,1)*dblBitSize
	call MPI_TYPE_CREATE_RESIZED(type_Temp, 0, newExtent, type_RF, code) !Changing the starting point to the next "slice"
	call MPI_TYPE_COMMIT(type_RF, code)

!    write(*,*) "Proc ", rang, "of", nb_procs
!	call DispCarvalhol(randField, "randField")


	call MPI_GATHER(randField,      size(randField), MPI_DOUBLE_PRECISION, &
	                tempRandField, 1,              type_RF,               &
	                0,              MPI_COMM_WORLD,  code)

	!if(rang == 0) call DispCarvalhol(tempRandField, "tempRandField")

	call MPI_TYPE_FREE(type_RF, code)

	!START Test Gather--------------------------------------------------------

!	call MPI_TYPE_VECTOR(4, 4, 8, MPI_DOUBLE_PRECISION, type_Temp, code)
!	call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblBitSize, code)
!	newExtent = 4*dblBitSize
!	call MPI_TYPE_CREATE_RESIZED(type_Temp, 0, newExtent, type_RF, code)
!	call MPI_TYPE_COMMIT(type_RF, code)
!
!	longueur_tranche = nb_valeurs/(nb_procs)
!	k = longueur_tranche/4
!	allocate(valeurs(k, k))
!	do j = 1, k
!		valeurs(:, j) = (/((i+j), i=1, k)/)
!	end do
!	write(*,*) "Proc ", rang, "of", nb_procs
!	call DispCarvalhol(valeurs, "valeurs")
!	!call MPI_BARRIER(MPI_COMM_WORLD, code)
!	call MPI_GATHER(valeurs, longueur_tranche, MPI_DOUBLE_PRECISION,    &
!	                donnees, 1, type_RF, &
!	                0, MPI_COMM_WORLD, code)
!	if(rang == 0) call DispCarvalhol(donnees, "donnees")
!	call MPI_TYPE_FREE(type_RF, code)

	!END Test Gather--------------------------------------------------------

	if(rang == 0) then
		!STILL TO PARALELIZE
		xTotalnStep = product(xNStep)
		allocate (totalRandField(xTotalnStep, Nmc))
		totalRandField = tempRandField(1:xTotalnStep,:)
		deallocate(tempRandField)

		!call DispCarvalhol(totalRandField)

		call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
						 kMax, xNStep, kNStep, totalRandField, &
						 average, stdDeviation, averageCorrL, outputName, endTime - startTime)

		call set_StatisticsND(totalRandField, xNStep, xMax, average, stdDeviation, averageCorrL)

		call CPU_TIME(endTime)

		call writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
	    						 kMax, xNStep, kNStep, randField, &
	    						 average, stdDeviation, averageCorrL, outputName, endTime - startTime)

	    call writeResultHDF5(xMax, kMax, xNStep, totalRandField, outputName)

		call CPU_TIME(endTime)

		write(*,*) ">>>>>>>>> Showing Results";
		write(*,*) "";
		write(*,*) "     Global Average = ", calculateAverage(average);
		write(*,*) "Global StdDeviation = ", calculateAverage(stdDeviation);
		if (Nmc == 1) write(*,*) " Correlation Length not computed (Nmc = 1) = ", averageCorrL;
		if (Nmc >  1) write(*,*) " Correlation Length = ", averageCorrL;
		write(*,*) "     Total Time (s) = ", endTime - startTime;
	end if

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

	if(rang == 0) then
		write(*,*) "";
	    write(*,*) "------------END main-----------------------";
		write(*,*) "";
	end if

	call MPI_FINALIZE(code)

end program main_RandomField
