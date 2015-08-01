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
	use readUNV_RF

    implicit none

    !INPUTS
    integer                       :: nDim, Nmc;
    character (len=30), parameter :: mesh_input = "mesh_input"
    character (len=30), parameter :: gen_input  = "gen_input"
    character (len=30), parameter :: test_input = "test_input"
    character (len=30), parameter :: unv_input  = "unv_files/Maroua1.unv "!"unv_files/Luciano_Cube.unv"
    logical :: step_variate, nmc_variate, corrL_variate
    logical :: step_speedUp, nmc_speedUp, corrL_fix_pointsPerCorrL
    integer :: step_nIter, nmc_nIter, corrL_nIter
    integer :: nmc_initial, corrL_pointsPerCorrL
    integer :: compiler = 2 !1 for gfortran and 2 for ifort
    logical :: writeFiles = .true.
    logical :: sameFolder = .true.
    logical :: explodedView = .false.
    integer :: outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc

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
    double precision, dimension(:,:), allocatable, target :: coordList
    integer         , dimension(:,:), allocatable :: connectList
    logical :: monotype


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
            !call read_DataTable(dataTable, "meshType", MSH%meshType)
            call read_DataTable(dataTable, "meshMod", MSH%meshMod)
            call init_MESH(MSH, nDim, comm, rang, nb_procs)

            select case (MSH%meshMod)
                case("unv")
                    write(*,*) "Read input UNV"
                    call readUNV(unv_input, nDim, coordList, connectList, monotype, rang, nb_procs, comm)
                    !call dispCarvalhol(coordList, "coordList", "(F20.5)",unit_in = get_fileId())
                    !call dispCarvalhol(connectList, "connectList",unit_in = get_fileId())
                case("automatic")
                    baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
                    call read_DataTable(dataTable, "Min", MSH%xMinGlob)
                    call read_DataTable(dataTable, "Max", MSH%xMaxGlob)
                    call read_DataTable(dataTable, "Step", MSH%xStep)
                    !write(rangChar,'(I7)') rang
                    !call read_DataTable(dataTable, "Max"//trim(adjustl(rangChar)), MSH%xMax)
                    !call read_DataTable(dataTable, "Min"//trim(adjustl(rangChar)), MSH%xMin)
                case default
                    write(*,*) "meshMod not accepted: ", MSH%meshMod
                    call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end select

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


            if (MSH%meshMod == "unv") then
                RDF%xPoints => coordList
                write(get_fileId(),*) "-> defining_UNV_extremes"
                call get_Global_Extremes_Mesh(RDF%xPoints, MSH%xMinGlob, MSH%xMaxGlob, RDF%comm)
                RDF%xNTotal = MSH%xNTotal
                RDF%xMinBound = MSH%xMinGlob
                RDF%xMaxBound = MSH%xMaxGlob
                MSH%xMinBound = MSH%xMinGlob
                MSH%xMaxBound = MSH%xMaxGlob
                MSH%xMin    = minval(RDF%xPoints, 2)
                MSH%xMax    = maxval(RDF%xPoints, 2)
                MSH%xNTotal = size(RDF%xPoints, 2)
                RDF%xNTotal = MSH%xNTotal

            else
                write(get_fileId(),*) "-> set_procPerDim"
                call set_procPerDim (MSH)
                write(get_fileId(),*) "-> MPI_CART_CREATE"
                call MPI_CART_CREATE (MSH%comm, MSH%nDim, MSH%procPerDim, periods, .false., MSH%topComm, code)
                write(get_fileId(),*) "-> MPI_CART_COORDS"
                call MPI_CART_COORDS (MSH%topComm, MSH%rang, MSH%nDim, MSH%coords, code)
                write(get_fileId(),*) "-> redefine_Global_Extremes"
                write(get_fileId(),*) "     BEFORE:"
                write(get_fileId(),*) "     MSH%xStep    = ", MSH%xStep
                write(get_fileId(),*) "     MSH%xMinGlob = ", MSH%xMinGlob
                write(get_fileId(),*) "     MSH%xMaxGlob = ", MSH%xMaxGlob
                call redefine_Global_Extremes (MSH, RDF, pointsPerCorrL = 10)
                write(get_fileId(),*) "     AFTER:"
                write(get_fileId(),*) "     MSH%xStep    = ", MSH%xStep
                write(get_fileId(),*) "     MSH%xMinGlob = ", MSH%xMinGlob
                write(get_fileId(),*) "     MSH%xMaxGlob = ", MSH%xMaxGlob
                write(get_fileId(),*) " "
                write(get_fileId(),*) "-> set_Local_Extremes_From_Coords"
                call set_Local_Extremes_From_Coords (MSH)
                write(get_fileId(),*) "-> Getting Global Matrix Reference"
                call get_XPoints_globCoords(RDF, MSH)
                write(get_fileId(),*) "     RDF%origin = ", RDF%origin
                write(get_fileId(),*) "     RDF%stride = ", RDF%stride
                if(RDF%independent) then
                    write(get_fileId(),*) "-> set_neighbours"
                    call set_neighbours (MSH)
                    write(get_fileId(),*) "-> redefine_Overlap"
                    call redefine_Overlap (MSH, RDF)
                    write(get_fileId(),*) "-> redefine_extremes_for_overlap"
                    call redefine_extremes_for_overlap (MSH, RDF%corrL)
                end if

            end if

            deallocate(periods)

        end subroutine define_topography

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine single_realization()

            implicit none
            double precision, dimension(:), allocatable :: seedStartVec

            write(get_fileId(),*) "-> Defining Topography"
            call define_topography()

            write(get_fileId(),*) "-> Initializing Random Seed"
            if(MSH%independent) then
                !Define independent seed in each proc
                call calculate_random_seed(RDF%seed, RDF%seedStart+RDF%rang)
                call init_random_seed(RDF%seed)

            else
                call calculate_random_seed(RDF%seed, RDF%seedStart)
                call init_random_seed(RDF%seed)
            end if

            write(get_fileId(),*) "-> Setting xPoints"
            call set_XPoints(MSH, RDF, RDF%xPoints_Local)
            write(get_fileId(),*) "     shape(RDF%xPoints)    = ", shape(RDF%xPoints)
            write(get_fileId(),*) "     maxval(RDF%xPoints,2) = ", maxval(RDF%xPoints,2)
            write(get_fileId(),*) "     minval(RDF%xPoints,2) = ", minval(RDF%xPoints,2)

            if(size(RDF%xPoints,2) <50) call dispCarvalhol(transpose(RDF%xPoints), "transpose(RDF%xPoints)", "(F20.5)",unit_in = get_fileId())

            call allocate_randField(RDF, RDF%randField_Local)

            write(get_fileId(),*) "-> Setting folder path"
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

            t2 = MPI_Wtime();
            call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,comm,code)
            if(RDF%rang == 0) write(*,*) "Generation Time = ", all_t2 - all_t1
            all_t3 = -1.0D0

            if(explodedView .and. RDF%independent) then
                do i = 1, RDF%nDim
                    RDF%xPoints(i,:) = RDF%xPoints(i,:) + 1.5*MSH%coords(i)*RDF%corrL(i)*MSH%overlap
                end do
            end if

            if(writeFiles) then
                write(get_fileId(),*) "-> Writing XMF and hdf5 files";

                if(outputStyle==1) then
                    write(get_fileId(),*) "   (Parallel)";
                    call write_Mono_XMF_h5(RDF, MSH, connectList, monotype, "trans_", RDF%rang, single_path, &
                                                    MPI_COMM_WORLD, ["_All"], [0], 0, style=outputStyle)
                else
                    write(get_fileId(),*) "   (Per Proc)";
                    call write_Mono_XMF_h5(RDF, MSH, connectList, monotype, "trans_", RDF%rang, single_path, &
                                                    MPI_COMM_WORLD, ["_All"], [RDF%rang], 0, style=outputStyle)
                end if
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

            if(allocated(coordList))   deallocate(coordList)
            if(allocated(connectList)) deallocate(connectList)

        end subroutine deallocate_all

end program main_RandomField
