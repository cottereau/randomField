program main_RandomField

	use mpi
	use constants_RF
	use readFile_RF
	use randomFieldND
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
    integer                       :: nDim, Nmc;
    character (len=30), parameter :: mesh_input = "mesh_input"
    character (len=30), parameter :: gen_input  = "gen_input"
    character (len=30), parameter :: test_input = "test_input"
    logical :: step_variate, nmc_variate, corrL_variate
    logical :: step_speedUp, nmc_speedUp, corrL_fix_pointsPerCorrL
    integer :: step_nIter, nmc_nIter, corrL_nIter
    integer :: nmc_initial, corrL_pointsPerCorrL
    integer :: compiler = 2 !1 for gfortran and 2 for ifort
    logical :: writeFiles = .false.
    logical :: sameFolder = .false.

    double precision, dimension(:), allocatable :: step_mult, step_add, step_initial
    double precision, dimension(:), allocatable :: corrL_mult, corrL_add, corrL_initial
    integer         , dimension(:), allocatable :: step_div
    double precision :: nmc_mult, nmc_add

	!LOCAL VARIABLES
    logical            :: file_exist
    integer            :: i, baseStep, nIter, all_xNTotal
    integer            :: comm, code, rang, error, nb_procs;
    character(len=30)  :: rangChar;
    character(len=200) :: path
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    character(len=30) , dimension(:,:), allocatable :: dataTable;
    integer, dimension(:), allocatable :: seed
    logical, dimension(:), allocatable :: periods


    double precision :: t1, t2, t3, t4, t5, t6;
    double precision :: all_t1, all_t2, all_t3, all_t4, all_t5, all_t6;
    double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                      avg_Trans, stdDev_Trans
    double precision  , dimension(:,:) , allocatable :: avg_Gauss_evnt, stdDev_Gauss_evnt, &
                                                        avg_Trans_evnt, stdDev_Trans_evnt
    double precision  , dimension(:,:) , allocatable :: avg_Gauss_point, stdDev_Gauss_point, &
                                                        avg_Trans_point, stdDev_Trans_point

    type(RF)   :: RDF
    type(MESH) :: MSH
    !type(TEST) :: TST



    !Initializing MPI
    call init_communication(MPI_COMM_WORLD)

    if(rang == 0) write(*,*) "-> MPI_communications started"
    if(rang == 0) write(*,*) "         nb_procs = ", nb_procs

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

    !Initial allocation
    call allocate_init()

    !Generating random fields
    if(rang == 0) write(*,*)  "-> Generating random fields"
    write(get_fileId(),*) "-> Generating random fields"

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

            comm = get_fileId(rang)
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

            if(sameFolder) results_folder_name = "res" !ONLY FOR TESTS

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            if(rang == 0) write(*,*) "results_folder_name = ", trim(results_folder_name)

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(sameFolder) log_folder_name     = ".." !ONLY FOR TESTS

            !results folder creation
            call create_folder(results_path, ".", rang, comm, compiler)

            !results_date_time folder creation
            call delete_folder(results_folder_name, results_path, rang, comm, compiler)
            call create_folder(results_folder_name, results_path, rang, comm, compiler)

            !log folder creation
            call delete_folder(log_folder_name, results_path, rang, comm, compiler)
            call create_folder(log_folder_name, results_path, rang, comm, compiler)

            !create xmf and h5 folders
            if(writeFiles) then
                path = string_vec_join([results_path,"/",results_folder_name])
                call create_folder("xmf", path, rang, comm, compiler)
                call create_folder("h5", path, rang, comm, compiler)
            end if


        end subroutine init_basic_folders


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

            integer ::independent

            !Reading Mesh------------------------------------------------------------------
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
            call read_DataTable(dataTable, "nDim", nDim)
            call read_DataTable(dataTable, "meshType", MSH%meshType)
            call read_DataTable(dataTable, "meshMod", MSH%meshMod)
            call init_MESH(MSH, nDim, comm, rang, nb_procs)

            if (MSH%meshType == "structured" .or. MSH%meshType == "unstructured") then
               select case (MSH%meshMod)
                   case("manual")
                       write(rangChar,'(I7)') rang
                       call read_DataTable(dataTable, "Max"//trim(adjustl(rangChar)), MSH%xMax)
                       call read_DataTable(dataTable, "Min"//trim(adjustl(rangChar)), MSH%xMin)
                   case("automatic")
                        baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
                        call read_DataTable(dataTable, "Max", MSH%xMaxGlob)
                        call read_DataTable(dataTable, "Min", MSH%xMinGlob)
                        call read_DataTable(dataTable, "Step", MSH%xStep)
                end select
            else
                write(*,*) "meshType not accepted: ", MSH%meshType
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if

            deallocate(dataTable)

            !Reading Generation Input---------------------------------------------------------------
            write(get_fileId(),*) "    Reading Generation Input"
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
            call read_DataTable(dataTable, "nDim", nDim)
            call read_DataTable(dataTable, "Nmc", Nmc)
            call init_RF(RDF, nDim, Nmc, comm, rang, nb_procs)
            call read_DataTable(dataTable, "corrMod"    , RDF%corrMod)
            call read_DataTable(dataTable, "margiFirst" , RDF%margiFirst)
            call read_DataTable(dataTable, "fieldAvg"   , RDF%fieldAvg)
            call read_DataTable(dataTable, "fieldVar"   , RDF%fieldVar)
            call read_DataTable(dataTable, "method"     , RDF%method)
            call read_DataTable(dataTable, "seedStart"  , RDF%seedStart)
            call read_DataTable(dataTable, "corrL"      , RDF%corrL)
            call read_DataTable(dataTable, "independent", independent)

            if(independent == 1) then
                RDF%independent = .true.
                MSH%independent = .true.

                call read_DataTable(dataTable, "overlap", MSH%overlap)
            else
                RDF%independent = .false.
                MSH%independent = .false.
            end if

            if(RDF%nb_procs == 1) then
                write(get_fileId(),*) "WARNING!! Independent generation in a single processor."
                write(get_fileId(),*) " "
                write(get_fileId(),*) "--OLD values--"
                write(get_fileId(),*) "RDF%independent = ", RDF%independent
                write(get_fileId(),*) "MSH%independent = ", MSH%independent
                write(get_fileId(),*) "MSH%overlap     = ", MSH%overlap
                if(RDF%independent) MSH%overlap = -2.0D0
                RDF%independent = .false.
                MSH%independent = .false.
                write(get_fileId(),*) " "
                write(get_fileId(),*) "--NEW values (changed)--"
                write(get_fileId(),*) "RDF%independent = ", RDF%independent
                write(get_fileId(),*) "MSH%independent = ", MSH%independent
                write(get_fileId(),*) "MSH%overlap     = ", MSH%overlap
            end if

            deallocate(dataTable)

            !Input Validation
            if(RDF%Nmc < 1) then
               write(*,*) ""
               write(*,*) "ERROR - Number of events should be a positive integer"
               write(*,*) "RDF%Nmc = ", RDF%Nmc
               call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            if((RDF%nDim < 1) .or. (MSH%nDim < 1) .or. (RDF%nDim /= MSH%nDim)) then
               write(*,*) ""
               write(*,*) "ERROR - nDim should be a positive integer and should be the same in both mesh and generation files"
               write(*,*) "RDF%nDim = ", RDF%nDim
               write(*,*) "MSH%nDim = ", MSH%nDim
               call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            do i = 1, RDF%nDim
               if(RDF%corrL(i) <= 0.0d0) then
                   write(*,*) ""
                   write(*,*) "ERROR - corrL should be a positive number greater than 0.0"
                   write(*,*) "corrL ", i, "of proc ", RDF%rang, " is ", RDF%corrL(i)
                   call MPI_ABORT(MPI_COMM_WORLD, error, code)
               end if
            end do

        end subroutine read_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine define_topography()

            allocate(periods(MSH%nDim))
            periods(:) = .false.

            write(get_fileId(),*) "-> set_procPerDim"

            call set_procPerDim (MSH)
            write(get_fileId(),*) "-> MPI_CART_CREATE"
            call MPI_CART_CREATE (MSH%comm, MSH%nDim, MSH%procPerDim, periods, .false., MSH%topComm, code)
            write(get_fileId(),*) "-> MPI_CART_COORDS"
            call MPI_CART_COORDS (MSH%topComm, MSH%rang, MSH%nDim, MSH%coords, code)

            if(RDF%independent) then
                write(get_fileId(),*) "-> redefine_Global_Extremes"
                call redefine_Global_Extremes (MSH, RDF, pointsPerCorrL = 10)
                write(get_fileId(),*) "-> set_Local_Extremes_From_Coords"
                call set_Local_Extremes_From_Coords (MSH)
                write(get_fileId(),*) "-> set_neighbours"
                call set_neighbours (MSH)
                write(get_fileId(),*) "-> redefine_Overlap"
                call redefine_Overlap (MSH, RDF)
                write(get_fileId(),*) "-> redefine_extremes"
                call redefine_extremes (MSH, RDF%corrL)
            else
                call set_Local_Extremes_Mesh (MSH%xMin, MSH%xMax, MSH%xMinGlob, MSH%xMaxGlob, MSH%rang, MSH%nb_procs)
            end if

            call show_MESH(MSH, "MSH", unit_in = get_fileId())

            deallocate(periods)

        end subroutine define_topography

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine single_realization()

            implicit none
            double precision, dimension(:), allocatable :: seedStartVec

            write(get_fileId(),*) "Defining Topography"
            call define_topography()
            write(get_fileId(),*) "Initializing Random Seed"
            call calculate_random_seed(RDF%seed, RDF%seedStart)
            call init_random_seed(RDF%seed)

            if(MSH%independent) then
                !Define independent seed in each proc
                call calculate_random_seed(RDF%seed, RDF%seedStart+RDF%rang)
                call init_random_seed(RDF%seed)
                !Building xPoints
                write(get_fileId(),*) "Setting xPoints (independent)"
                call set_XPoints_independent(MSH, RDF, RDF%xPoints_Local)
            else
                !Building xPoints
                call set_XPoints(MSH, RDF, RDF%xPoints_Local)
            end if


            call allocate_randField(RDF, RDF%randField_Local)

            single_path = string_vec_join([results_path,"/",results_folder_name])
            write(get_fileId(),*) "single_path = ", single_path

            !Discovering the total number of points in all procs
            write(get_fileId(),*) "Discovering total number of points (MPI_ALLREDUCE)"
            call MPI_ALLREDUCE (RDF%xNTotal, all_xNTotal,1,MPI_INTEGER, &
                                MPI_SUM,comm,code)

            write(get_fileId(),*) "Discovering t1 (MPI_ALLREDUCE)"
            t1 = MPI_Wtime();
            call MPI_ALLREDUCE (t1, all_t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
            if(RDF%rang == 0) write(*,*) "Time Zero = ", all_t1

            write(get_fileId(),*) "Generating Random Field"
            call create_RF_Unstruct_Init (RDF, MSH)
            call show_RF(RDF, "RDF", unit_in = get_fileId())

            t2 = MPI_Wtime();
            call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
            if(RDF%rang == 0) write(*,*) "Generation Time = ", all_t2 - all_t1
            all_t3 = -1.0D0

            if(writeFiles) then
                write(get_fileId(),*) "-> Writing XMF and hdf5 files";
                call write_Mono_XMF_h5(RDF%xPoints, RDF%randField, "trans_", RDF%rang, single_path, &
                                                    MPI_COMM_WORLD, ["_proc_"], [RDF%rang], 0)
                t3 = MPI_Wtime();
                call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
                if(RDF%rang == 0) write(*,*) "Writing Files Time = ", all_t3 - all_t2
            end if

            call write_generation_spec(MSH, RDF, single_path, "singleGen", &
                                       [all_t1,all_t2,all_t3])


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

            call finalize_MESH(MSH)
            call finalize_RF(RDF)

            if(allocated(step_mult)) deallocate(step_mult);
            if(allocated(step_add))  deallocate(step_add);
            if(allocated(step_div))  deallocate(step_div);

        end subroutine deallocate_all

end program main_RandomField

!!TRASH

!    if(step_variate) then
!        if(step_speedUp) then
!            if(rang == 0) write(*,*)  "-> Running xStep variation tests (speed)"
!            write(get_fileId(),*) "-> Running xStep variation tests (speed)"
!            call func_test_001_speed_xStep(xMin, xMax,                    &
!                                          corrL, corrMod,                 &
!                                          margiFirst, fieldAvg, fieldVar, &
!                                          Nmc, method,                    &
!                                          step_initial, step_nIter,       &
!                                          step_mult, step_add,            &
!                                          xMinGlob, xMaxGlob)
!        else
!            if(rang == 0) write(*,*)  "-> Running xStep variation tests (convergence)"
!            write(get_fileId(),*) "-> Running xStep variation tests (convergence)"
!            call func_test_003_conv_xStep(xMin, xMax,                     &
!                                          corrL, corrMod,                 &
!                                          margiFirst, fieldAvg, fieldVar, &
!                                          Nmc, method,                    &
!                                          step_initial, step_nIter,       &
!                                          step_div,                       &
!                                          xMinGlob, xMaxGlob)
!        end if
!    end if
!    if(nmc_variate) then
!        !call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
!        if(nmc_speedUp) then
!            if(rang == 0) write(*,*)  "-> Running Nmc variation tests"
!            write(get_fileId(),*) "-> Running Nmc variation tests"
!            call func_test_002_speed_Nmc(xMin, xMax,                     &
!                                        corrL, corrMod,                  &
!                                        margiFirst, fieldAvg, fieldVar,  &
!                                        xStep, method,                   &
!                                        nmc_initial, nmc_nIter,          &
!                                        nmc_mult, nmc_add,               &
!                                        xMinGlob, xMaxGlob)
!        else
!            if(rang == 0) write(*,*)  "-> Running Nmc variation tests (convergence)"
!            write(get_fileId(),*) "-> Running Nmc variation tests (convergence)"
!            call func_test_004_conv_Nmc(xMin, xMax,                     &
!                                       corrL, corrMod,                  &
!                                       margiFirst, fieldAvg, fieldVar,  &
!                                       xStep, method,                   &
!                                       nmc_initial, nmc_nIter,          &
!                                       nmc_mult, nmc_add,               &
!                                       xMinGlob, xMaxGlob)
!        end if
!        !call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
!    end if
!
!    if(corrL_variate) then
!        if(rang == 0) write(*,*)  "-> Running corrL variation tests"
!        write(get_fileId(),*) "-> Running corrL variation tests"
!        if(corrL_fix_pointsPerCorrL) then
!            call func_test_005_conv_corrL(xMin, xMax,                     &
!                                          xStep, corrMod,                 &
!                                          margiFirst, fieldAvg, fieldVar, &
!                                          Nmc, method,                    &
!                                          corrL_initial, corrL_nIter,     &
!                                          corrL_mult, corrL_add, corrL_pointsPerCorrL, &
!                                          xMinGlob, xMaxGlob)
!        else
!            call func_test_005_conv_corrL(xMin, xMax,                     &
!                                      xStep, corrMod,                 &
!                                      margiFirst, fieldAvg, fieldVar, &
!                                      Nmc, method,                    &
!                                      corrL_initial, corrL_nIter,     &
!                                      corrL_mult, corrL_add)
!        end if
!    end if

!        !---------------------------------------------------------------------------------
!        !---------------------------------------------------------------------------------
!        !---------------------------------------------------------------------------------
!        !---------------------------------------------------------------------------------
!        subroutine init_test_folders()
!            implicit none
!
!            character(len=200) :: path
!
!            if(step_variate) then
!                call create_folder(xStep_folder_name, &
!                                   string_vec_join([results_path,"/",results_folder_name]), &
!                                   rang, comm, compiler)
!                path = string_vec_join([results_path,"/",results_folder_name, "/", xStep_folder_name])
!                call create_folder("xmf", path, rang, comm, compiler)
!                call create_folder("h5", path, rang, comm, compiler)
!            end if
!            if(corrL_variate) then
!                call create_folder(corrL_folder_name, &
!                                   string_vec_join([results_path,"/",results_folder_name]), &
!                                   rang, comm, compiler)
!                path = string_vec_join([results_path,"/",results_folder_name, "/", corrL_folder_name])
!                call create_folder("xmf", path, rang, comm, compiler)
!                call create_folder("h5", path, rang, comm, compiler)
!            end if
!            if(nmc_variate) then
!                call create_folder(Nmc_folder_name, &
!                                   string_vec_join([results_path,"/",results_folder_name]), &
!                                   rang, comm, compiler)
!                path = string_vec_join([results_path,"/",results_folder_name, "/", Nmc_folder_name])
!                call create_folder("xmf", path, rang, comm, compiler)
!                call create_folder("h5", path, rang, comm, compiler)
!            end if
!            if(.not. (step_variate.or. nmc_variate .or. corrL_variate)) then
!                path = string_vec_join([results_path,"/",results_folder_name])
!                call create_folder("xmf", path, rang, comm, compiler)
!                call create_folder("h5", path, rang, comm, compiler)
!            end if
!        end subroutine init_test_folders





!            !nIter = 1
!
!            !Discovering the total number of points in all procs
!            call MPI_ALLREDUCE (RDF%xNTotal, all_xNTotal,1,MPI_INTEGER, &
!                                MPI_SUM,comm,code)
!
!            !t1 = MPI_Wtime();
!            !call MPI_ALLREDUCE (t1, all_t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            !if(RDF%rang == 0) write(*,*) "Time Zero = ", all_t1
!
!            call create_RF_Unstruct_Init (RDF, MSH)
!
!            !t2 = MPI_Wtime();
!            !call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            !if(RDF%rang == 0) write(*,*) "Generation Time = ", all_t2 - all_t1
!
!
!            i = 1 !Iteration number
!
!            !---------------------------------------------------------------------------------
!            write(get_fileId(),*) "-> Creating Standard Gaussian Random field";
!
!            t1 = MPI_Wtime();
!            call MPI_ALLREDUCE (t1, all_t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            if(RDF%rang == 0) write(*,*) "Time Zero = ", all_t1
!
!            call create_Std_Gaussian_Field_Unstruct (RDF%xPoints, RDF%corrL, RDF%corrMod, RDF%Nmc,  &
!                                                     RDF%randField, RDF%method, RDF%seed)
!
!            t2 = MPI_Wtime();
!            call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            if(RDF%rang == 0) write(*,*) "Generation Time = ", all_t2 - all_t1
!
!
            !write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            !call write_Mono_XMF_h5(RDF%xPoints, RDF%randField, "gauss_", RDF%rang, single_path, &
            !                                    MPI_COMM_WORLD, ["proc_"], [RDF%rang], 0)

            !t3 = MPI_Wtime();
            !call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
            !if(RDF%rang == 0) write(*,*) "Writing Files Time = ", all_t3 - all_t2

            !call write_generation_spec(MSH, RDF, single_path, "singleGen", &
            !                           [all_t1,all_t2,all_t3])
!
!            t3 = MPI_Wtime();
!            call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            if(RDF%rang == 0) write(*,*) "Writing Files Time = ", all_t3 - all_t2
!
!            if(RDF%rang == 0) call show_MESH(MSH)
            !if(RDF%rang == 0) call show_RF(RDF)
            !---------------------------------------------------------------------------------

!            write(get_fileId(),*) "-> Calculating Statistics Before Transformation";

    !        write(*,*) "BEFORE TRANSFORMATION "
    !        write(*,*) "Average = ", sum(randField)/(xNTotal*Nmc)
    !        write(*,*) "Variance = ", sum(randField**2)/(xNTotal*Nmc) - (sum(randField)/(xNTotal*Nmc))**2


!            call calculate_average_and_stdVar_MPI(randField,                                  &
!                                                  avg_Gauss(i), stdDev_Gauss(i),              &
!                                                  comm,                                       &
!                                                  avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:))
    !                                              avg_Gauss_point(:,i), stdDev_Gauss_point(:,i))

            !---------------------------------------------------------------------------------
!            write(get_fileId(),*) "-> Transforming Random field";
!
!            t4 = MPI_Wtime();
!            call MPI_ALLREDUCE (t4, all_t4, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            if(RDF%rang == 0) write(*,*) "Time from zero till here = ", all_t4 - all_t1
!
!            !call multiVariateTransformation (RDF%margiFirst, RDF%fieldAvg, RDF%fieldVar, RDF%randField)
!
!            t5 = MPI_Wtime();
!            call MPI_ALLREDUCE (t5, all_t5, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            if(RDF%rang == 0) write(*,*) "Transformation Time = ", all_t5 - all_t4

!            !START Only For Plot 4 proc 2D
!            if(RDF%rang == 0) then
!                RDF%xPoints(1,:) = RDF%xPoints(1,:) -2.0D0
!                RDF%xPoints(2,:) = RDF%xPoints(2,:) -2.0D0
!            else if(RDF%rang == 1) then
!                RDF%xPoints(1,:) = RDF%xPoints(1,:) -2.0D0
!                RDF%xPoints(2,:) = RDF%xPoints(2,:) +2.0D0
!            else if(RDF%rang == 2) then
!                RDF%xPoints(1,:) = RDF%xPoints(1,:) +2.0D0
!                RDF%xPoints(2,:) = RDF%xPoints(2,:) -2.0D0
!            else if(RDF%rang == 3) then
!                RDF%xPoints(1,:) = RDF%xPoints(1,:) +2.0D0
!                RDF%xPoints(2,:) = RDF%xPoints(2,:) +2.0D0
!            end if
!            !END Only For Plot
!
!            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
!            call write_Mono_XMF_h5(RDF%xPoints, RDF%randField, "trans_", RDF%rang, single_path, &
!                                                MPI_COMM_WORLD, ["_proc_"], [RDF%rang], 0)
!
!            t6 = MPI_Wtime();
!            call MPI_ALLREDUCE (t6, all_t6, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
!            if(RDF%rang == 0) write(*,*) "Writing Files Time = ", all_t6 - all_t5
!
!            if(RDF%rang == 0) call show_RF(RDF)
            !call dispCarvalhol(RDF%randField, "RDF%randField")
            !---------------------------------------------------------------------------------
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
!

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

!            deallocate(avg_Gauss, stdDev_Gauss)
!            deallocate(avg_Trans, stdDev_Trans)
!            deallocate(avg_Gauss_evnt, stdDev_Gauss_evnt)
!            deallocate(avg_Trans_evnt, stdDev_Trans_evnt)
!            deallocate(avg_Gauss_point, stdDev_Gauss_point)
!            deallocate(avg_Trans_point, stdDev_Trans_point)











!            !Looking for test_input
!            write(get_fileId(),*) "    Checking for Test Input"
!            step_variate  = .false.
!            nmc_variate   = .false.
!            corrL_variate = .false.
!            path = test_input
!            path = adjustL(path)
!            write(get_fileId(),*) "        file: ", trim(path)
!            inquire(file=path, exist = file_exist)
!            if(file_exist) then
!                call set_DataTable(path, dataTable)
!                !if(rang == 0) call DispCarvalhol (dataTable, path);
!                call read_DataTable(dataTable, "step_variate", step_variate)
!                call read_DataTable(dataTable, "nmc_variate", nmc_variate)
!                call read_DataTable(dataTable, "corrL_variate", corrL_variate)
!
!                if(step_variate) then
!                    allocate(step_mult(nDim))
!                    allocate(step_add(nDim))
!                    allocate(step_initial(nDim))
!                    call read_DataTable(dataTable, "step_speedUp", step_speedUp)
!                    call read_DataTable(dataTable, "step_initial", step_initial)
!                    call read_DataTable(dataTable, "step_nIter", step_nIter)
!                    if(step_speedUp) then
!                        call read_DataTable(dataTable, "step_mult", step_mult)
!                        call read_DataTable(dataTable, "step_add", step_add)
!                    else
!                        call read_DataTable(dataTable, "step_div", step_div)
!                    end if
!
!                    !call recalculate_xStep(xMin, xMax, step_initial)
!                end if
!
!                if(corrL_variate) then
!                    allocate(corrL_mult(nDim))
!                    allocate(corrL_add(nDim))
!                    allocate(corrL_initial(nDim))
!                    call read_DataTable(dataTable, "corrL_fix_pointsPerCorrL", corrL_fix_pointsPerCorrL)
!                    call read_DataTable(dataTable, "corrL_initial", corrL_initial)
!                    call read_DataTable(dataTable, "corrL_nIter", corrL_nIter)
!                    call read_DataTable(dataTable, "corrL_mult", corrL_mult)
!                    call read_DataTable(dataTable, "corrL_add", corrL_add)
!                    if(corrL_fix_pointsPerCorrL) then
!                        call read_DataTable(dataTable,"corrL_pointsPerCorrL", corrL_pointsPerCorrL)
!                        !call recalculate_corrL(xMin, xMax, corrL_initial)
!                    end if
!                end if
!
!                if(nmc_variate) then
!                    call read_DataTable(dataTable, "nmc_speedUp", nmc_speedUp)
!                    call read_DataTable(dataTable, "nmc_initial", nmc_initial)
!                    call read_DataTable(dataTable, "nmc_mult", nmc_mult)
!                    call read_DataTable(dataTable, "nmc_add", nmc_add)
!                    call read_DataTable(dataTable, "nmc_nIter", nmc_nIter)
!                end if
!
!                deallocate(dataTable)
!            else
!                if(rang == 0) write(*,*) "The file ", path, " doesn't exist (no tests will be performed)"
!            end if



