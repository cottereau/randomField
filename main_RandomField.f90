program main_RandomField

	use randomFieldND
	use statistics_RF
	use obsolete_RF
	use writeResultFile_RF
	use readFile_RF
	use displayCarvalhol
	use mpi
	use test_func_RF

    implicit none

    !INPUTS
    integer                                        :: Nmc, nDim;
    character (len=30), parameter                  :: mainInput        = "input_main"
    character (len=50), parameter                  :: sampleSpecFolder = "../sampleSpec"
    character (len=50), parameter                  :: meshFolder       = "../mesh"
    character (len=50), parameter                  :: outputFolder     = "../output"
    character (len=50)                             :: sampleSpecName, meshName, outputName;
    character (len=15)                             :: corrMod, margiFirst, meshType, meshMod;
    double precision,   dimension(:),  allocatable :: corrL;
    double precision                               :: fieldAvg, fieldVar;
    logical                                        :: calculateStat = .true.
    integer         , parameter                    :: pointsPerCorrL = 10;
    integer         , parameter                    :: method = 1; !1 for Victor, 2 for Shinozuka

	!OUTPUTS (in this case variables to be filled in in the proces)
	integer         , dimension(:)   , allocatable :: xNStep, kNStep;
    integer         , dimension(:, :), allocatable :: all_xNStep;
	double precision, dimension(:)   , allocatable :: kMax, xMax, xMin, xStep;
    double precision, dimension(:, :), allocatable :: randField, xPoints;
    double precision, dimension(:, :), allocatable :: all_RandField, all_xPoints;
    double precision, dimension(:, :), allocatable :: all_xMin, all_xMax;
    double precision, dimension(:)   , allocatable :: evntAvg, evntStdDev, procCorrL, globalCorrL;
    double precision, dimension(:)   , allocatable :: compEvntAvg, compEvntStdDev, compGlobCorrL;
    double precision                               :: globalAvg, globalStdDev, compGlobAvg, compGlobStdDev;
!    double precision  , dimension(:) , allocatable :: tests_average_Gauss, tests_stdDev_Gauss, &
!    												  tests_average_Trans, tests_stdDev_Trans, &
!    												  changeProp;

	!LOCAL VARIABLES
    logical            :: structured = .false., calcComp = .false.
    integer            :: i, baseStep, testSize;
    integer            :: code, rang, error, nb_procs;
    double precision   :: pi = 3.1415926535898;
    character(len=30)  :: rangChar;
    character(len=200) :: path
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    character(len=30) , dimension(:,:), allocatable :: dataTable;
    integer            :: nIterations
    double precision   :: mult_Step, add_Step;
    logical            :: variate_Nmc = .false., variate_xStep = .false.
    double precision  , dimension(:)  , allocatable :: sumRF, sumRFsquare, &
                                                       totalSumRF, totalSumRFsquare;

	call MPI_INIT(code)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

	if(rang == 0) then
	    write(*,*) "------------START main-----------------------------------------";
		write(*,*) "";
	    write(*,*) "------------START MPI processing-----------------------------------------";
		write(*,*) "Number of procs = ", nb_procs
	    write(*,*) "";
		write(*,*) ">>>>>>>>> Reading file input";
	end if

	!Allocating
	allocate(HDF5Name(1))

	!Reading Filenames Input
	call set_DataTable("../"//mainInput, dataTable)
	!if(rang == 0) call DispCarvalhol (dataTable, "input_main");
	call read_DataTable(dataTable, "sampleSpecName" , sampleSpecName)
	call read_DataTable(dataTable, "meshName"       , meshName)
	call read_DataTable(dataTable, "outputName"     , outputName)
	deallocate(dataTable)

	!Reading Input
	path = trim(adjustL(sampleSpecFolder))//"/"//sampleSpecName
	call set_DataTable(path, dataTable)
	!if(rang == 0) call DispCarvalhol (dataTable, sampleSpecName);
	call read_DataTable(dataTable, "Nmc"        ,  Nmc)
	call read_DataTable(dataTable, "corrMod"    , corrMod)
	call read_DataTable(dataTable, "margiFirst" , margiFirst)
	call read_DataTable(dataTable, "corrL"      , corrL)
	call read_DataTable(dataTable, "fieldAvg"   , fieldAvg)
	call read_DataTable(dataTable, "fieldVar"   , fieldVar)
	deallocate(dataTable)

	nDim = size(corrL)

	!Reading Mesh
	path = trim(adjustL(meshFolder))//"/"//meshName
	call set_DataTable(path, dataTable)
	!if(rang == 0) call DispCarvalhol (dataTable, meshName);
	call read_DataTable(dataTable, "meshType", meshType)
	call read_DataTable(dataTable, "meshMod", meshMod)
	if (meshType == "structured" .or. meshType == "unstructured") then
		select case (meshMod)
			case("manual")
				write(rangChar,'(I7)') rang
				call read_DataTable(dataTable, "Max"//trim(adjustl(rangChar)), xMax)
				call read_DataTable(dataTable, "Min"//trim(adjustl(rangChar)), xMin)
			case("automatic")
				baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
				if (nb_procs == baseStep**nDim) then
					call read_DataTable(dataTable, "Max", xMax)
					call read_DataTable(dataTable, "Min", xMin)
					call read_DataTable(dataTable, "Step", xStep)
					call set_ProcMaxMin(xMin, xMax, rang, nb_procs)
				else
					write(*,*) "In meshing automatic mode the number of processors should be an integer like (N)**nDim, where N is an integer and nDim is the number of dimensions"
					call MPI_ABORT(MPI_COMM_WORLD, error, code)
				end if
		end select
		nDim = size(xMax)
	end if
	deallocate(dataTable)

	!Showing read inputs
	if(rang == 0) then
		write(*,*) ""
		write(*,*) ">>>>>>>>INPUT Read Data"
		write(*,*) "##IN/OUT"
		write(*,*) "sampleSpecName = ", sampleSpecName
		write(*,*) "outputName     = ", outputName
		write(*,*) ""
		write(*,*) "##MESH"
		write(*,*) "meshName       = ", meshName
		write(*,*) "meshType       = ", meshType
		write(*,*) "meshMod        = ", meshMod
		write(*,*) "xMin           = ", xMin
		write(*,*) "xMax           = ", xMax
		write(*,*) "xStep          = ", xStep
		write(*,*) ""
		write(*,*) "##EVENTS"
		write(*,*) "nDim           = ", nDim
		write(*,*) "Nmc            = ", Nmc !number of Monte-Carlo experiments
		write(*,*) "corrMod        = ", corrMod
		write(*,*) "margiFirst     = ", margiFirst
		write(*,*) "corrL          = ", corrL
		write(*,*) "fieldAvg       = ", fieldAvg
		write(*,*) "fieldVar       = ", fieldVar
	end if

	!Input Validation
	if(Nmc < 1) then
		write(*,*) ""
		write(*,*) "ERROR - Number of events should be a positive integer"
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if
	if(nDim < 1) then
		write(*,*) ""
		write(*,*) "ERROR - nDim should be a positive integer"
		write(*,*) "nDim           = ", nDim
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if
	do i = 1, nDim
		if(corrL(i) <= 0.0d0) then
			write(*,*) ""
			write(*,*) "ERROR - corrL should be a positive number greater than 0.0"
			write(*,*) "corrL ", i, "of proc ", rang, " is ", corrL(i)
			call MPI_ABORT(MPI_COMM_WORLD, error, code)
		end if
	end do

	path = trim(adjustL(outputFolder))//"/"//outputName
	!write(*,*) "Path before = ", path

	!/////////////TESTING BLOCK
    variate_Nmc   = .true.
    variate_xStep = .false.
    nIterations   = 10
    mult_Step     = 2
    add_Step      = 0

	call func_test_001_convergence(xMin, xMax, xStep,               &
								   corrL, corrMod, pointsPerCorrL,  &
							 	   margiFirst, fieldAvg, fieldVar,  &
							 	   Nmc, method,                     &
							 	   variate_Nmc, variate_xStep,      &
							 	   nIterations,                     &
							 	   mult_Step, add_Step,             &
							 	   path)

!	!///////////// START COVERGENCE TEST BLOCK
!	!Creating Standard Gaussian Field
!	testSize = 20
!	allocate(tests_average_Gauss(testSize), tests_stdDev_Gauss(testSize))
!	allocate(tests_average_Trans(testSize), tests_stdDev_Trans(testSize))
!	allocate(changeProp(testSize))
!
!	!Put here what is contant over interactions
!	if(rang == 0) write(*,*) ">>>>>>>>> Creating xPoints";
!	call set_XPoints(corrL, xMin, xMax, xStep, xPoints, pointsPerCorrL)
!	call dispCarvalhol(transpose(xPoints(:, 1:20)), "transpose(xPoints(:, 1:20))")
!
!	!Put here the starting point of what will change over interactions
!	allocate(randField(size(xPoints,2), Nmc))
!	randField = 0.0
!	call dispCarvalhol(randField(1:20, :), "randField(1:20, :) before loop")
!
!	do i = 1, testSize
!
!		changeProp(i) = dble(Nmc)
!
!		if(rang == 0) write(*,*) ">>>>>>>>> Creating Standard Gaussian Random field (unstructured)";
!		call createStandardGaussianFieldUnstruct(xPoints,                                        &
!		                                         corrMod, margiFirst, corrL, fieldAvg, fieldVar, &
!		                                         Nmc, method, randField);
!		if(calculateStat) then
!			if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics Before Transformation (unstructured)";
!			call set_Statistics_MPI(randField, xPoints, rang,             &
!							 		evntAvg, evntStdDev, procCorrL,       &
!							 		tests_average_Gauss(i), tests_stdDev_Gauss(i), globalCorrL)
!		end if
!
!		if(rang == 0) write(*,*) ">>>>>>>>> Transforming Random field (unstructured)";
!		call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)
!		if(calculateStat) then
!			if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics After Transformation (unstructured)";
!			call set_Statistics_MPI(randField, xPoints, rang,             &
!							 		evntAvg, evntStdDev, procCorrL,       &
!							 		tests_average_Trans(i), tests_stdDev_Trans(i), globalCorrL)
!		end if
!
!		!Put here what change in each interaction
!		Nmc = 2*Nmc
!		deallocate(randField)
!		allocate(randField(size(xPoints,2), Nmc))
!		!!!!!!!!!!!!!!
!
!	end do
!
!	!Printing
!	write(*,'(A4, A10, A1, 4A20)') "Iter", "Nmc", ")", "Avg-Gauss", "StdDev-Gauss", &
!		 											  "Avg-Trans", "StdDev-Trans"
!	do i = 1, testSize
!
!		write(*,'(I4, F10.0, A1, 4F20.8)') i,     changeProp(i) , ")", tests_average_Gauss(i), tests_stdDev_Gauss(i), &
!														   tests_average_Trans(i), tests_stdDev_Trans(i)
!	end do
!
!	call dispCarvalhol(tests_average_Gauss, "tests_average_Gauss")
!	call dispCarvalhol(tests_stdDev_Gauss, "tests_stdDev_Gauss")
!	call dispCarvalhol(tests_average_Trans, "tests_average_Trans")
!	call dispCarvalhol(tests_stdDev_Trans, "tests_stdDev_Trans")
!
!	deallocate(tests_average_Gauss, tests_stdDev_Gauss)
!	deallocate(tests_average_Trans, tests_stdDev_Trans)
!	deallocate(changeProp)
!	!///////////// END COVERGENCE TEST BLOCK

	!--------------------------------------------------------------------------------------
	!------------------------------UNSTRUCTURED--------------------------------------------
	!--------------------------------------------------------------------------------------

	!if(meshType == "unstructured") then
	if(.false.) then

		!Creating x points for each processor (in the future this would be one of the inputs)
		write(*,*) ""
		if(rang == 0) write(*,*) ">>>>>>>>> Creating xPoints for unstructured mesh";
		!allocate(xNStep(nDim))
		!xNStep  = pointsPerCorrL*ceiling((xMax-xMin)/corrL);
		!call set_XPoints(corrL, xMin, xMax, xPoints, pointsPerCorrL)
		!call dispCarvalhol(transpose(xPoints(:, 1:20)), "transpose(xPoints(:, 1:20))")

		!CHANGE WITH TRANSPOSITION
		!allocate(randField(size(xPoints,2), Nmc))

		!Calculating a part of the random field in each processor
		!if(rang == 0) write(*,*) ">>>>>>>>> Calculating Random field (unstructured)";
		!call createStandardGaussianFieldUnstructShinozuka(xPoints, corrL, corrMod, Nmc, randField)
		call createStandardGaussianFieldUnstructVictor(xPoints, corrL, corrMod, Nmc, randField);
		if(calculateStat) then
			if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics Before Transformation (unstructured)";
			call set_Statistics_MPI(randField, xPoints, rang,             &
							 		evntAvg, evntStdDev, procCorrL,       &
							 		globalAvg, globalStdDev, globalCorrL)
			if(rang == 0) then
				write(*,*) "    STD Avg = ", 0.0d0;
				write(*,*) "    MPI Avg = ", globalAvg;
				write(*,*) " "
				write(*,*) "    STD Var = ", 1.0d0;
				write(*,*) "    MPI Var = ", globalStdDev**2;
				write(*,*) " "
			end if
		end if
		call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)
		!call dispCarvalhol(xPoints(1:20, :), "xPoints(1:20, :)")
		!call dispCarvalhol(randField(1:20, :), "randField(1:20, :)")
				!Calculating Statistics
		if(calculateStat) then
			if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics After Transformation (unstructured)";
			call set_Statistics_MPI(randField, xPoints, rang,             &
							 		evntAvg, evntStdDev, procCorrL,       &
							 		globalAvg, globalStdDev, globalCorrL)
			if(rang == 0) then
				write(*,*) "  INPUT Avg = ", fieldAvg
				write(*,*) "    MPI Avg = ", globalAvg;
				write(*,*) " "
				write(*,*) "  INPUT Var = ", fieldVar
				write(*,*) "    MPI Var = ", globalStdDev**2;
				write(*,*) " "
			end if
		end if


		if(rang == 0) write(*,*) ">>>>>>>>> Writing Result File";
		path = trim(adjustL(outputFolder))//"/"//outputName
		!call write_ResultHDF5(xPoints, randField, path, rang) !HDF5 of this processor
		if(rang == 0) write(*,*) ">>>>>>>>> Write hdf5";
		call write_ResultHDF5(xPoints, randField, outputName, rang, outputFolder, &
    						  MPI_COMM_WORLD, ["proc"], [rang], HDF5Name(1))
		if(rang == 0) write(*,*) ">>>>>>>>> Write XMF";
		call writeXMF_RF_MPI(Nmc, HDF5Name, [size(randField,1)], nDim, outputName, rang, outputFolder, &
    						 MPI_COMM_WORLD, ".")



!		!Writing results and comparison statistics (To be deleted)
!		if(rang == 0) write(*,*) ">>>>>>>>> Calculating Comparison Statistics (unstructured)";
!		call set_allRandField(randField, xPoints, rang,    &
!						      all_RandField, all_xPoints)
!		!if(rang == 0) call DispCarvalhol(all_xPoints, "all_xPoints")
!		!if(rang == 0) call DispCarvalhol(all_RandField, "all_RandField")
!		if(rang == 0) then
!			call set_CompStatistics(all_RandField, all_xPoints,                   &
!			                        compEvntAvg, compEvntStdDev,                  &
!			                        compGlobAvg, compGlobStdDev, compGlobCorrL)
!		end if
!
!		if (rang == 0 .and. nDim == 1) then
!			path = trim(adjustL(outputFolder))//"/"//outputName
!			call write_MatlabTable (all_RandField, path)
!		end if

!	!--------------------------------------------------------------------------------------
!	!--------------------------------STRUCTURED (need to be updated)--------------------------------------------
!	!--------------------------------------------------------------------------------------
!
!	else if(meshType == "structured") then
!
!	else
!		write(*,*) "ERROR - meshtype not accepted"
!		write(*,*) "meshtype = ", meshtype
!		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if

!	if(rang == 0 .and. calculateStat) then
!		calcComp = .false.
!		write(*,*) "Showing Statistics------------------------------- "
!		write(*,*) ""
!		if(Nmc < 11) then
!			write(*,*) ">>BY EVENT"
!			call DispCarvalhol(evntAvg, "    MPI Avg = ")
!			if (calcComp) call DispCarvalhol(compEvntAvg, "   Comp Avg = ")
!			call DispCarvalhol(evntStdDev**2, "    MPI Var = ")
!			if (calcComp) call DispCarvalhol(compEvntStdDev**2, "   Comp Var = ")
!		end if
!
!		write(*,*) ">>GLOBAL"
!		write(*,*) "  INPUT Avg = ", fieldAvg
!		if (calcComp) write(*,*) "   Comp Avg = ", compGlobAvg;
!		write(*,*) "    MPI Avg = ", globalAvg;
!		write(*,*) " "
!		write(*,*) "  INPUT Var = ", fieldVar
!		if (calcComp) write(*,*) "   Comp Var = ", compGlobStdDev**2;
!		write(*,*) "    MPI Var = ", globalStdDev**2;
!		write(*,*) " "
!		!write(*,*) "INPUT CorrL = ", corrL
!		!if (calcComp) write(*,*) " Comp CorrL = ", compGlobCorrL;
!		!write(*,*) "  MPI CorrL = ", globalCorrL;
!		!call DispCarvalhol(globalCorrL, "  MPI CorrL ")
!		!call DispCarvalhol(compGlobCorrL, " Comp CorrL = ")
!	end if


	!Deallocating
		if(allocated(corrL))   deallocate(corrL);
		if(allocated(xMax))    deallocate(xMax);
		if(allocated(xMin))    deallocate(xMin);
		if(allocated(kMax))      deallocate(kMax);
		if(allocated(xNStep))    deallocate(xNStep);
		if(allocated(kNStep))    deallocate(kNStep);
		if(allocated(xPoints))   deallocate(xPoints);
		if(allocated(randField)) deallocate(randField);

		if(allocated(evntAvg))      deallocate(evntAvg);
		if(allocated(evntStdDev))   deallocate(evntStdDev);
		if(allocated(procCorrL))    deallocate(procCorrL);
		if(allocated(globalCorrL))  deallocate(globalCorrL);
		if(allocated(sumRF))        deallocate(sumRF)
		if(allocated(sumRFsquare))  deallocate(sumRFsquare)
		if(allocated(totalSumRF))        deallocate(totalSumRF)
		if(allocated(totalSumRFsquare))  deallocate(totalSumRFsquare)
		if(allocated(all_RandField))     deallocate(all_RandField);
		if(allocated(all_xPoints))       deallocate(all_xPoints);
		if(allocated(all_xNStep)) deallocate(all_xNStep)
		if(allocated(all_xMin))   deallocate(all_xMin)
		if(allocated(all_xMax))   deallocate(all_xMax)

		if(allocated(compEvntAvg))    deallocate(compEvntAvg);
		if(allocated(compEvntStdDev)) deallocate(compEvntStdDev);
		if(allocated(compGlobCorrL))  deallocate(compGlobCorrL);

	if(rang == 0) then
		write(*,*) "";
	    write(*,*) "------------END main-----------------------";
		write(*,*) "";
	end if

	call MPI_FINALIZE(code)

end program main_RandomField
