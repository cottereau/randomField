program main_RandomField

	use mpi
	use constants_RF
	use readFile_RF
	use randomFieldND
	use test_func_RF
	use statistics_RF
	use mesh_RF
	use writeResultFile_RF
	use displayCarvalhol
	use charFunctions
	use write_Log_File
	use systemUt_RF
	use common_variables_RF
	use type_RF
	use type_MESH

    implicit none

    !INPUTS
    integer                                        :: Nmc, nDim;
    character (len=30), parameter                  :: mesh_input = "mesh_input"
    character (len=30), parameter                  :: gen_input  = "gen_input"
    character (len=30), parameter                  :: test_input = "test_input"
    character (len=15)                             :: corrMod, margiFirst, meshType, meshMod;
    double precision,   dimension(:),  allocatable :: corrL;
    double precision                               :: fieldAvg, fieldVar;
    logical                                        :: calculateStat = .true.
    integer                                        :: method !1 for Isotropic, 2 for Shinozuka
    logical :: step_variate, nmc_variate, corrL_variate
    logical :: step_speedUp, nmc_speedUp, corrL_fix_pointsPerCorrL
    integer :: step_nIter, nmc_nIter, corrL_nIter
    integer :: nmc_initial, corrL_pointsPerCorrL
    integer :: compiler = 2 !1 for gfortran and 2 for ifort
    integer :: seedStart = 45

    double precision, dimension(:), allocatable :: step_mult, step_add, step_initial
    double precision, dimension(:), allocatable :: corrL_mult, corrL_add, corrL_initial
    integer         , dimension(:), allocatable :: step_div
    double precision :: nmc_mult, nmc_add

	!OUTPUTS (in this case variables to be filled in in the proces)
	integer         , dimension(:)   , allocatable :: xNStep;
    double precision, dimension(:, :), allocatable :: randField, xPoints;
    double precision, dimension(:)   , allocatable :: xMax, xMin, xStep;

	!LOCAL VARIABLES
    logical            :: file_exist
    integer            :: i, baseStep, xNTotal, nIter, all_xNTotal;
    integer            :: comm, code, rang, error, nb_procs;
    character(len=30)  :: rangChar;
    character(len=200) :: path
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    character(len=30) , dimension(:,:), allocatable :: dataTable;
    double precision  , dimension(:),   allocatable :: xMinGlob, xMaxGlob
    integer, dimension(:), allocatable :: seed

    character(len=200) :: single_path
    double precision :: t1, t2, t3, t4;
    double precision :: all_t1, all_t2, all_t3, all_t4;
    double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                      avg_Trans, stdDev_Trans
    double precision  , dimension(:,:) , allocatable :: avg_Gauss_evnt, stdDev_Gauss_evnt, &
                                                        avg_Trans_evnt, stdDev_Trans_evnt
    double precision  , dimension(:,:) , allocatable :: avg_Gauss_point, stdDev_Gauss_point, &
                                                        avg_Trans_point, stdDev_Trans_point

    type(RF) :: testRF
    type(MESH) :: testMESH

    !Initializing MPI
    !call init_communication(MPI_COMM_WORLD)
    call MPI_INIT(code)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

    !write(*,*) "rang = ",rang
    !write(*,*) "nb_procs = ",nb_procs

    comm = MPI_COMM_WORLD

    testRF%rang = 54321
    testMESH%rang = 54321

    write(*,*) "testRF%rang = ", testRF%rang
    write(*,*) "testMESH%rang = ", testMESH%rang

    if(rang == 0) write(*,*)  "-> Initialize MPI_communications"
    !Initializing folders
    if(rang == 0) write(*,*)  "-> Initialize Folders"
    call init_basic_folders()
    !Initializing logFiles
    if(rang == 0) write(*,*)  "-> Initialize logFiles"
    call init_log_file(stringNumb_join(string_vec_join([results_path,"/",log_folder_name,"/",log_filename]), rang), rang)

    if(rang == 0)then
        write(*,*)
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "************************                             ***********************"
        write(*,*) "************************  RANDOM FIELD LIBRARY TEST  ***********************"
        write(*,*) "************************                             ***********************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*)
    end if

	if(rang == 0) write(*,*) "  Number of procs = ", nb_procs
	if(rang == 0) write(*,*) ""

    !Reading Inputs
    if(rang == 0) write(*,*)  "-> Reading inputs"
    write(get_fileId(),*) "-> Reading inputs"
    call read_input()

    !Making test folders
    call init_test_folders()

    !Initial allocation
    call allocate_init()

    !Generating random fields
    if(rang == 0) write(*,*)  "-> Generating random fields"
    write(get_fileId(),*) "-> Generating random fields"
    write(get_fileId(),*) "   It works"

    if(step_variate) then
        if(step_speedUp) then
            if(rang == 0) write(*,*)  "-> Running xStep variation tests (speed)"
            write(get_fileId(),*) "-> Running xStep variation tests (speed)"
            call func_test_001_speed_xStep(xMin, xMax,                    &
                                          corrL, corrMod,                 &
                                          margiFirst, fieldAvg, fieldVar, &
                                          Nmc, method,                    &
                                          step_initial, step_nIter,       &
                                          step_mult, step_add,            &
                                          xMinGlob, xMaxGlob)
        else
            if(rang == 0) write(*,*)  "-> Running xStep variation tests (convergence)"
            write(get_fileId(),*) "-> Running xStep variation tests (convergence)"
            call func_test_003_conv_xStep(xMin, xMax,                     &
                                          corrL, corrMod,                 &
                                          margiFirst, fieldAvg, fieldVar, &
                                          Nmc, method,                    &
                                          step_initial, step_nIter,       &
                                          step_div,                       &
                                          xMinGlob, xMaxGlob)
        end if
    end if
    if(nmc_variate) then
        call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        if(nmc_speedUp) then
            if(rang == 0) write(*,*)  "-> Running Nmc variation tests"
            write(get_fileId(),*) "-> Running Nmc variation tests"
            call func_test_002_speed_Nmc(xMin, xMax,                     &
                                        corrL, corrMod,                  &
                                        margiFirst, fieldAvg, fieldVar,  &
                                        xStep, method,                   &
                                        nmc_initial, nmc_nIter,          &
                                        nmc_mult, nmc_add,               &
                                        xMinGlob, xMaxGlob)
        else
            if(rang == 0) write(*,*)  "-> Running Nmc variation tests (convergence)"
            write(get_fileId(),*) "-> Running Nmc variation tests (convergence)"
            call func_test_004_conv_Nmc(xMin, xMax,                     &
                                       corrL, corrMod,                  &
                                       margiFirst, fieldAvg, fieldVar,  &
                                       xStep, method,                   &
                                       nmc_initial, nmc_nIter,          &
                                       nmc_mult, nmc_add,               &
                                       xMinGlob, xMaxGlob)
        end if
        call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
    end if

    if(corrL_variate) then
        if(rang == 0) write(*,*)  "-> Running corrL variation tests"
        write(get_fileId(),*) "-> Running corrL variation tests"
        if(corrL_fix_pointsPerCorrL) then
            call func_test_005_conv_corrL(xMin, xMax,                     &
                                          xStep, corrMod,                 &
                                          margiFirst, fieldAvg, fieldVar, &
                                          Nmc, method,                    &
                                          corrL_initial, corrL_nIter,     &
                                          corrL_mult, corrL_add, corrL_pointsPerCorrL, &
                                          xMinGlob, xMaxGlob)
        else
            call func_test_005_conv_corrL(xMin, xMax,                     &
                                      xStep, corrMod,                 &
                                      margiFirst, fieldAvg, fieldVar, &
                                      Nmc, method,                    &
                                      corrL_initial, corrL_nIter,     &
                                      corrL_mult, corrL_add)
        end if
    end if

    !SINGLE REALIZATION
    if(.not. (step_variate.or. nmc_variate .or. corrL_variate)) then
            if(rang == 0) write(*,*)  "-> Single realization"
            write(get_fileId(),*) "-> Single realization"
        call single_realization()
    end if

	!Deallocating
	call deallocate_all()


	if(rang == 0) then
	    write(*,*) ""
		write(*,*) "---------------------------------------------------------------------";
	    write(*,*) "-----------------END RANDOM FIELD LIBRARY TEST-----------------------";
		write(*,*) "---------------------------------------------------------------------";
        write(*,*) ""
	end if

    !Finalizing MPI
	call end_communication()

        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
    contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_communication(comm_local)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local

            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rang, code)
            call MPI_COMM_SIZE(comm_local, nb_procs, code)

            write(*,*) "rang = ",rang
            write(*,*) "nb_procs = ",nb_procs

            comm = comm_local

        end subroutine init_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_basic_folders()
            implicit none
            !LOCAL
            logical :: dirExists
            integer, dimension(8) :: date_time
            character(len=10), dimension(3) :: strings

            !date_time_label
            call date_and_time(strings(1), strings(2), strings(3), date_time)
            results_folder_name = strings(1)(3:8)//"_"//strings(2)(1:6)//"_res"

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            write(*,*) "results_folder_name = ", trim(results_folder_name)

            !log_folder_name     = strings(1)(3:8)//"_"//strings(2)(1:6)//"_log"
            log_folder_name     = trim(adjustL(results_folder_name))//"/log"


            !write(*,*) "results_folder_name = ", results_folder_name

            !results folder check
            call create_folder(results_path, ".", rang, comm, compiler)

            !results_date_time folder creation
            call delete_folder(results_folder_name, results_path, rang, comm, compiler)
            call create_folder(results_folder_name, results_path, rang, comm, compiler)

            !log folder creation
            call delete_folder(log_folder_name, results_path, rang, comm, compiler)
            call create_folder(log_folder_name, results_path, rang, comm, compiler)

        end subroutine init_basic_folders

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_test_folders()
            implicit none

            character(len=200) :: path

            if(step_variate) then
                call create_folder(xStep_folder_name, &
                                   string_vec_join([results_path,"/",results_folder_name]), &
                                   rang, comm, compiler)
                path = string_vec_join([results_path,"/",results_folder_name, "/", xStep_folder_name])
                call create_folder("xmf", path, rang, comm, compiler)
                call create_folder("h5", path, rang, comm, compiler)
            end if
            if(corrL_variate) then
                call create_folder(corrL_folder_name, &
                                   string_vec_join([results_path,"/",results_folder_name]), &
                                   rang, comm, compiler)
                path = string_vec_join([results_path,"/",results_folder_name, "/", corrL_folder_name])
                call create_folder("xmf", path, rang, comm, compiler)
                call create_folder("h5", path, rang, comm, compiler)
            end if
            if(nmc_variate) then
                call create_folder(Nmc_folder_name, &
                                   string_vec_join([results_path,"/",results_folder_name]), &
                                   rang, comm, compiler)
                path = string_vec_join([results_path,"/",results_folder_name, "/", Nmc_folder_name])
                call create_folder("xmf", path, rang, comm, compiler)
                call create_folder("h5", path, rang, comm, compiler)
            end if
            if(.not. (step_variate.or. nmc_variate .or. corrL_variate)) then
                path = string_vec_join([results_path,"/",results_folder_name])
                call create_folder("xmf", path, rang, comm, compiler)
                call create_folder("h5", path, rang, comm, compiler)
            end if
        end subroutine init_test_folders

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine allocate_init()
        end subroutine allocate_init

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_input()

            implicit none

            !double precision, allocatable, dimension(:) :: xMaxIni, xMinIni

            !Reading Generation Statistis Input
            write(get_fileId(),*) "    Reading Generation Statistics Input"
            !path = trim(adjustL(sampleSpecFolder))//"/"//sampleSpecName
            path = gen_input
            path = adjustL(path)
            write(get_fileId(),*) "        file: ", trim(path)
            inquire(file=path, exist = file_exist)
            if(.not. file_exist) then
                write(get_fileId(),*) "ERROR - The file ", path, " was not found"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if

            call MPI_BARRIER (comm ,code)

            call set_DataTable(path, dataTable)
            !if(rang == 0) call DispCarvalhol (dataTable, path);
            call read_DataTable(dataTable, "Nmc"        ,  Nmc)
            call read_DataTable(dataTable, "corrMod"    , corrMod)
            call read_DataTable(dataTable, "margiFirst" , margiFirst)
            call read_DataTable(dataTable, "corrL"      , corrL)
            call read_DataTable(dataTable, "fieldAvg"   , fieldAvg)
            call read_DataTable(dataTable, "fieldVar"   , fieldVar)
            call read_DataTable(dataTable, "method"     , method)
            deallocate(dataTable)

            !Allocation
            nDim = size(corrL)
            allocate (xMax(nDim))
            allocate (xMin(nDim))
            allocate (xMaxGlob(nDim))
            allocate (xMinGlob(nDim))
            allocate (xStep(nDim))

            !Reading Mesh
            write(get_fileId(),*) "    Reading Mesh Input"
            path = mesh_input
            path = adjustL(path)
            write(get_fileId(),*) "        file: ", trim(path)
            inquire(file=path, exist = file_exist)
            if(.not. file_exist) then
                write(*,*) "ERROR - The file ", path, " was not found"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if

            call MPI_BARRIER (comm ,code)

            call set_DataTable(path, dataTable)
            !if(rang == 0) call DispCarvalhol (dataTable, path);
            call read_DataTable(dataTable, "meshType", meshType)
            call read_DataTable(dataTable, "meshMod", meshMod)

            if (meshType == "structured" .or. meshType == "unstructured") then
               select case (meshMod)
                   case("manual")
                       write(rangChar,'(I7)') rang
                       call read_DataTable(dataTable, "Max"//trim(adjustl(rangChar)), xMax)
                       call read_DataTable(dataTable, "Min"//trim(adjustl(rangChar)), xMin)
                   case("automatic")
                        !allocate(xMaxIni(nDim))
                        !allocate(xMinIni(nDim))
                        baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
                        !if (nb_procs == baseStep**nDim) then
                            call read_DataTable(dataTable, "Max", xMaxGlob)
                            call read_DataTable(dataTable, "Min", xMinGlob)
                            call read_DataTable(dataTable, "Step", xStep)
                            call set_Local_Extremes_Mesh (xMin, xMax, xMinGlob, xMaxGlob, rang, nb_procs)

                        !else
                        !    write(*,*) "In meshing automatic mode the number of processors should be an integer like (N)**nDim, where N is an integer and nDim is the number of dimensions"
                        !    call MPI_ABORT(MPI_COMM_WORLD, error, code)
                        !end if
                end select
            else
                write(*,*) "meshType not accepted: ", meshType
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            deallocate(dataTable)

            !Turning xStep and corrL to the nearest (rounding down) step possible
            write(get_fileId(),*) "    Correcting xStep"
            call recalculate_xStep(xMin, xMax, xStep)
            !write(get_fileId(),*) "    Correcting corrL"
            !call recalculate_corrL(xMin, xMax, corrL)


            !Looking for test_input
            write(get_fileId(),*) "    Checking for Test Input"
            step_variate  = .false.
            nmc_variate   = .false.
            corrL_variate = .false.
            path = test_input
            path = adjustL(path)
            write(get_fileId(),*) "        file: ", trim(path)
            inquire(file=path, exist = file_exist)
            if(file_exist) then
                call set_DataTable(path, dataTable)
                !if(rang == 0) call DispCarvalhol (dataTable, path);
                call read_DataTable(dataTable, "step_variate", step_variate)
                call read_DataTable(dataTable, "nmc_variate", nmc_variate)
                call read_DataTable(dataTable, "corrL_variate", corrL_variate)

                if(step_variate) then
                    allocate(step_mult(nDim))
                    allocate(step_add(nDim))
                    allocate(step_initial(nDim))
                    call read_DataTable(dataTable, "step_speedUp", step_speedUp)
                    call read_DataTable(dataTable, "step_initial", step_initial)
                    call read_DataTable(dataTable, "step_nIter", step_nIter)
                    if(step_speedUp) then
                        call read_DataTable(dataTable, "step_mult", step_mult)
                        call read_DataTable(dataTable, "step_add", step_add)
                    else
                        call read_DataTable(dataTable, "step_div", step_div)
                    end if

                    call recalculate_xStep(xMin, xMax, step_initial)
                end if

                if(corrL_variate) then
                    allocate(corrL_mult(nDim))
                    allocate(corrL_add(nDim))
                    allocate(corrL_initial(nDim))
                    call read_DataTable(dataTable, "corrL_fix_pointsPerCorrL", corrL_fix_pointsPerCorrL)
                    call read_DataTable(dataTable, "corrL_initial", corrL_initial)
                    call read_DataTable(dataTable, "corrL_nIter", corrL_nIter)
                    call read_DataTable(dataTable, "corrL_mult", corrL_mult)
                    call read_DataTable(dataTable, "corrL_add", corrL_add)
                    if(corrL_fix_pointsPerCorrL) then
                        call read_DataTable(dataTable,"corrL_pointsPerCorrL", corrL_pointsPerCorrL)
                        call recalculate_corrL(xMin, xMax, corrL_initial)
                    end if
                end if

                if(nmc_variate) then
                    call read_DataTable(dataTable, "nmc_speedUp", nmc_speedUp)
                    call read_DataTable(dataTable, "nmc_initial", nmc_initial)
                    call read_DataTable(dataTable, "nmc_mult", nmc_mult)
                    call read_DataTable(dataTable, "nmc_add", nmc_add)
                    call read_DataTable(dataTable, "nmc_nIter", nmc_nIter)
                end if

                deallocate(dataTable)
            else
                if(rang == 0) write(*,*) "The file ", path, " doesn't exist (no tests will be performed)"
            end if

            !Showing read inputs
            if(rang == 0) then
                write(*,*) ""
                write(*,*) "     MESH"
                write(*,*) "     |meshType = ", meshType
                write(*,*) "     |meshMod  = ", meshMod
                write(*,*) "     |xMin     = ", xMin
                write(*,*) "     |xMax     = ", xMax
                write(*,*) "     |xStep    = ", xStep
                write(*,*) ""
                write(*,*) "     EVENTS"
                write(*,*) "     |nDim       = ", nDim
                write(*,*) "     |Nmc        = ", Nmc !number of Monte-Carlo experiments
                write(*,*) "     |corrMod    = ", corrMod
                write(*,*) "     |margiFirst = ", margiFirst
                write(*,*) "     |corrL      = ", corrL
                write(*,*) "     |fieldAvg   = ", fieldAvg
                write(*,*) "     |fieldVar   = ", fieldVar
                write(*,*) "     |method     = ", method
                write(*,*) ""
                if(step_variate) then
                    write(*,*) "     TESTS STEP"
                    write(*,*) "     |step_initial = ", step_initial
                    write(*,*) "     |step_nIter   = ", step_nIter
                    write(*,*) "     |step_speedUp = ", step_speedUp
                    if(step_speedUp) then
                        write(*,*) "     |step_mult    = ", step_mult
                        write(*,*) "     |step_add     = ", step_add
                    else
                        write(*,*) "     |step_div     = ", step_div
                    end if
                end if
                if(corrL_variate) then
                    write(*,*) "     TESTS corrL"
                    write(*,*) "     |corrL_initial            = ", corrL_initial
                    write(*,*) "     |corrL_nIter              = ", corrL_nIter
                    write(*,*) "     |corrL_mult               = ", corrL_mult
                    write(*,*) "     |corrL_add                = ", corrL_add
                    write(*,*) "     |corrL_fix_pointsPerCorrL = ", corrL_fix_pointsPerCorrL
                    if(corrL_fix_pointsPerCorrL) then
                        write(*,*)"     |corrL_pointsPerCorrL     = ", corrL_pointsPerCorrL
                    end if
                end if
                if(nmc_variate) then
                    write(*,*) "     TESTS NMC"
                    write(*,*) "     |nmc_initial = ", nmc_initial
                    write(*,*) "     |nmc_mult    = ", nmc_mult
                    write(*,*) "     |nmc_add     = ", nmc_add
                    write(*,*) "     |nmc_nIter   = ", nmc_nIter
                end if
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

        end subroutine read_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine single_realization()

            call calculate_random_seed(seed, seedStart)
            call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
            call set_XPoints(xMin, xMax, xStep, xPoints)
            call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)

            xNTotal = size(xPoints, 2)
            nIter = 1

            !Discovering the total number of points in all procs
            call MPI_ALLREDUCE (xNTotal, all_xNTotal,1,MPI_INTEGER, &
                                MPI_SUM,comm,code)

            write(get_fileId(),*) "xMin = ", xMin
            write(get_fileId(),*) "xMax = ", xMax
            write(get_fileId(),*) "xMinGlob = ", xMinGlob
            write(get_fileId(),*) "xMaxGlob = ", xMaxGlob
            write(get_fileId(),*) "xStep = ", xStep
            write(get_fileId(),*) "xNTotal = ", xNTotal
            write(get_fileId(),*) "nIter = ", nIter
            write(get_fileId(),*) "Nmc = ", Nmc

            allocate(randField(size(xPoints,2), Nmc))
            allocate(avg_Gauss(nIter), stdDev_Gauss(nIter))
            allocate(avg_Trans(nIter), stdDev_Trans(nIter))
            allocate(avg_Gauss_evnt(nIter, Nmc), stdDev_Gauss_evnt(nIter, Nmc))
            allocate(avg_Trans_evnt(nIter, Nmc), stdDev_Trans_evnt(nIter, Nmc))
            allocate(avg_Gauss_point(all_xNTotal, nIter), stdDev_Gauss_point(all_xNTotal, nIter))
            allocate(avg_Trans_point(all_xNTotal, nIter), stdDev_Trans_point(all_xNTotal, nIter))

    !        write(get_fileId(),*) "shape(avg_Gauss_point) = ", shape(avg_Gauss_point)
    !        write(get_fileId(),*) "shape(stdDev_Gauss_point) = ", shape(stdDev_Gauss_point)
!
!            i = 1 !Iteration number
!
!            single_path = string_vec_join([results_path,"/",results_folder_name])
!            write(get_fileId(),*) "single_path = ", single_path
!
!            write(get_fileId(),*) "-> Creating Standard Gaussian Random field";
!            t1 = MPI_Wtime();
!            write(*,*) "t1 = ", t1
!            call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
!                                                      randField, method, seed)
!            t2 = MPI_Wtime();
!            write(*,*) "t2 - t1 = ", t2 - t1
!            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
!            call write_Mono_XMF_h5(xPoints, randField, "gauss_", rang, single_path, &
!                                                MPI_COMM_WORLD, ["proc_"], [rang], 0)
!
!            write(get_fileId(),*) "-> Calculating Statistics Before Transformation";

    !        write(*,*) "BEFORE TRANSFORMATION "
    !        write(*,*) "Average = ", sum(randField)/(xNTotal*Nmc)
    !        write(*,*) "Variance = ", sum(randField**2)/(xNTotal*Nmc) - (sum(randField)/(xNTotal*Nmc))**2


!            call calculate_average_and_stdVar_MPI(randField,                                  &
!                                                  avg_Gauss(i), stdDev_Gauss(i),              &
!                                                  comm,                                       &
!                                                  avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:))
    !                                              avg_Gauss_point(:,i), stdDev_Gauss_point(:,i))
!            t3 = MPI_Wtime();
!            call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)
!            t4 = MPI_Wtime();
!
!            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
!                call write_Mono_XMF_h5(xPoints, randField, "trans_", rang, single_path, &
!                                                MPI_COMM_WORLD, ["_proc_"], [rang], 0)
!
!            write(get_fileId(),*) "-> Calculating Statistics After Transformation";

    !        write(*,*) "AFTER TRANSFORMATION "
    !        write(*,*) "Average = ", sum(randField)/(xNTotal*Nmc)
    !        write(*,*) "Variance = ", sum(randField**2)/(xNTotal*Nmc) - (sum(randField)/(xNTotal*Nmc))**2

!            call calculate_average_and_stdVar_MPI(randField,                                   &
!                                                  avg_Trans(i), stdDev_Trans(i),               &
!                                                  comm,                                        &
!                                                  avg_Trans_evnt(i,:), stdDev_Trans_evnt(i,:))
!                                                  !avg_Trans_point(:,i), stdDev_Trans_point(:,i))
!
!            !Discovering the global time
!            call MPI_ALLREDUCE (t1, all_t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            call MPI_ALLREDUCE (t4, all_t4, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!
!            call write_generation_spec(xMinGlob, xMaxGlob, xStep,      &
!                                       corrL, corrMod,                 &
!                                       margiFirst, fieldAvg, fieldVar, &
!                                       Nmc, method, seed,              &
!                                       rang, single_path, "singleGen", &
!                                       [all_t1,all_t2,all_t3,all_t4])
!                                      !avg_Gauss(i), stdDev_Gauss(i),                 &
!                                      !avg_Trans(i), stdDev_Trans(i),                 &
!                                      !avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:),   &
!                                      !avg_Trans_evnt(i,:), stdDev_Trans_evnt(i,:))
!                                      !avg_Gauss_point(:,i), stdDev_Gauss_point(:,i), &
!                                      !avg_Trans_point(:,i), stdDev_Trans_point(:,i))

!            if(rang == 0) then
!                write(*,*) "-----STATISTICS-----"
!                write(*,*) "     Gauss By Event (", all_xNTotal, " points per event)"
!                write(*,*) "          avg_Gauss_evnt FIRST = ", avg_Gauss_evnt(i, 1)
!                write(*,*) "          avg_Gauss_evnt MIDDLE = ", avg_Gauss_evnt(i, (1+size(avg_Gauss_evnt,2))/2)
!                write(*,*) "          avg_Gauss_evnt LAST  = ", avg_Gauss_evnt(i, size(avg_Gauss_evnt,2))
!                write(*,*) "          stdDev_Gauss_evnt FIRST = ", stdDev_Gauss_evnt(i, 1)
!                write(*,*) "          stdDev_Gauss_evnt MIDDLE = ", stdDev_Gauss_evnt(i, (1+size(stdDev_Gauss_evnt,2))/2)
!                write(*,*) "          stdDev_Gauss_evnt LAST  = ", stdDev_Gauss_evnt(i, size(stdDev_Gauss_evnt,2))
!
!                write(*,*) "     Gauss By Point (", Nmc, " events per point)"
!                write(*,*) "          avg_Gauss_point FIRST = ", avg_Gauss_point(1, i)
!                write(*,*) "          avg_Gauss_point MIDDLE = ", avg_Gauss_point((1+size(avg_Gauss_point,1))/2, i)
!                write(*,*) "          avg_Gauss_point LAST  = ", avg_Gauss_point(size(avg_Gauss_point,1), i)
!                write(*,*) "          stdDev_Gauss_point FIRST = ", stdDev_Gauss_point(1, i)
!                write(*,*) "          stdDev_Gauss_point MIDDLE = ", stdDev_Gauss_point((1+size(stdDev_Gauss_point,1))/2, i)
!                write(*,*) "          stdDev_Gauss_point LAST  = ", stdDev_Gauss_point(size(stdDev_Gauss_point,1)/2, i)
!
!                write(*,*) ""
!                write(*,*) "     GLOBAL "
!                write(*,*) ""
!                write(*,'(A16, 4A20)') " )", "Avg-Gauss", "StdDev-Gauss", &
!                    "Avg-Trans", "StdDev-Trans"
!
!                write(*,'(A16, 4F20.8)') " )",              &
!                        avg_Gauss(i), stdDev_Gauss(i), &
!                        avg_Trans(i), stdDev_Trans(i)
!                write(*,'(A16, 4F20.8)') "Reference      )", &
!                    0.0d0, 1.0d0,  &
!                    fieldAvg, fieldVar
!            end if

            deallocate(avg_Gauss, stdDev_Gauss)
            deallocate(avg_Trans, stdDev_Trans)
            deallocate(avg_Gauss_evnt, stdDev_Gauss_evnt)
            deallocate(avg_Trans_evnt, stdDev_Trans_evnt)
            deallocate(avg_Gauss_point, stdDev_Gauss_point)
            deallocate(avg_Trans_point, stdDev_Trans_point)
        end subroutine single_realization

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine end_communication()
            call finalize_log_file(rang)
            call MPI_FINALIZE(code)
        end subroutine end_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine deallocate_all()
            if(allocated(corrL))   deallocate(corrL);
            if(allocated(xMax))    deallocate(xMax);
            if(allocated(xMin))    deallocate(xMin);
            if(allocated(xMaxGlob)) deallocate(xMaxGlob)
            if(allocated(xMinGlob)) deallocate(xMinGlob)
            if(allocated(xStep))   deallocate(xStep);
            if(allocated(xPoints))   deallocate(xPoints);
            if(allocated(randField)) deallocate(randField);
            if(allocated(seed)) deallocate(seed);

            if(allocated(step_mult)) deallocate(step_mult);
            if(allocated(step_add))  deallocate(step_add);
            if(allocated(step_div))  deallocate(step_div);

        end subroutine deallocate_all

end program main_RandomField
