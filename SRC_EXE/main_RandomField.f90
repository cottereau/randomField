program main_RandomField

	use mpi
	use constants_RF
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
	use type_inputRF
	use calls_RF
	use sample_RF



    implicit none

    !INPUTS
    integer :: nDim, Nmc;
    integer :: compiler = 2 !1 for gfortran and 2 for ifort
    logical :: writeFiles = .true.
    logical :: writeDataSet = .true.
    logical :: sameFolder = .true.
    integer :: outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc

	!LOCAL VARIABLES
    logical            :: file_exist
    integer            :: i, baseStep, nIter
    integer            :: rang
    character(len=30)  :: rangChar;
    character(len=200) :: path, logFilePath
    integer, dimension(:), allocatable :: seed
    double precision, dimension(5) :: times, all_times

    !DEVEL
    !integer, dimension(2) :: nFields = [2,1] !Number of independent fields in each dimension
    integer               :: fieldNumber
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    !double precision, dimension(:), allocatable :: procStart
    type(MESH)            :: globMSH
    integer               :: group, groupComm, groupMax
    integer               :: code, nTotalProcs
    logical               :: validProc
    double precision, dimension(:), allocatable :: stepProc, procExtent, overlap, meshStep
    double precision, dimension(:,:), allocatable :: subdivisionCoords
    type(IPT_RF)  :: IPT

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, IPT%comm, IPT%rang, IPT%nb_procs)
    rang = IPT%rang
    IPT%writeDataSet = writeDataSet

    times(1) = MPI_Wtime() !Initial Time

    if(IPT%rang == 0)then
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

    if(rang == 0) write(*,*) "-> MPI_communications started"
    if(rang == 0) write(*,*) "         running on: "
    if(rang == 0) call system("pwd")
    if(rang == 0) write(*,*) "         nb_procs    = ", IPT%nb_procs
    if(rang == 0) write(*,*) "         outputStyle = ", outputStyle

    !Initializing folders
    if(rang == 0) write(*,*)  "-> Initialize Folders"
    call init_basic_folders(IPT%comm)

#ifdef MAKELOG
    if(rang == 0) write(*,*) "IFDEF MAKELOG DEFINED"

    !Initializing logFiles
    if(rang == 0) write(*,*)  "-> Initialize logFiles"
    logFilePath = trim(adjustL(&
                      string_join_many(results_path,"/",log_folder_name,"/",log_filename)))
    if(rang == 0) write(*,*)  " logFilePath = ", trim(adjustL(logFilePath)), "<RANK>"
    call init_log_file(trim(adjustL(logFilePath)), rang, IPT%log_ID)
#else
    if(rang == 0) write(*,*) "IFDEF MAKELOG NOT DEFINED"
#endif

    !READING INPUTS--------------------------------------------------
    !----------------------------------------------------------------
    if(rang == 0) write(*,*)  "-> Reading inputs"
    call wLog("-> Reading inputs")
    !Reading Mesh---------------------------------------------------
    if(rang == 0) write(*,*)  "     -> Reading Mesh Input"
    call wLog("     -> Reading Mesh Input")
    path = mesh_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_mesh_input(path, IPT)
    !Reading Generation Input---------------------------------------
    if(rang == 0) write(*,*)  "     -> Reading Generation Input"
    call wLog("     -> Reading Generation Input")
    path = gen_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_generation_input(path, IPT)
    !Validating Inputs----------------------------------------------
    call wLog("    Validating Inputs")
    call validate_input(IPT)
    if(rang == 0) call show_IPT_RF(IPT)
#ifdef MAKELOG
    call show_IPT_RF(IPT, forLog_in=.true.)
#endif

    times(2) = MPI_Wtime() !Reading Inputs


    !Initial allocation---------------------------------------------
    call allocate_init()

    if(IPT%nProcPerField > IPT%nb_procs) then
        write(*,*) "ERROR, IPT%nProcPerField cannot be superior to IPT%nb_procs"
        stop(" ")
    end if

    !Changing xMaxGlob and xMinGlob according to localization level
    if(rang == 0) write(*,*) "-> REDEFINE xMaxGlob and xMinGlob----------------------------------------"
    call wLog("-> REDEFINE xMaxGlob and xMinGlob----------------------------------------")
    call redefineIPTlimits(IPT, IPT%xMinGlob, IPT%xMaxGlob, IPT%localizationLevel)
    if(rang == 0) write(*,*) "IPT%xMinGlob"
    if(rang == 0) write(*,*) IPT%xMinGlob
    if(rang == 0) write(*,*) "IPT%xMaxGlob"
    if(rang == 0) write(*,*) IPT%xMaxGlob
    if(rang == 0) write(*,*) "IPT%localizationLevel"
    if(rang == 0) write(*,*) IPT%localizationLevel
    call wLog("IPT%xMinGlob")
    call wLog(IPT%xMinGlob)
    call wLog("IPT%xMaxGlob")
    call wLog(IPT%xMaxGlob)
    call wLog("IPT%localizationLevel")
    call wLog(IPT%localizationLevel)

    !Building Subdivisions
    if(rang == 0) write(*,*) " "
    if(rang == 0) write(*,*) "-> DIVIDING----------------------------------------"
    call wLog("-> DIVIDING----------------------------------------")
    allocate(stepProc(IPT%nDim_mesh))
    allocate(procExtent(IPT%nDim_mesh))
    allocate(overlap(IPT%nDim_mesh))
    allocate(subdivisionCoords(IPT%nDim_mesh, product(IPT%nFields**IPT%localizationLevel)))
    call build_subdivisions(IPT, globMSH, groupMax, &
                            group, groupComm, stepProc, procExtent, overlap)
    call setGrid(subdivisionCoords, globMSH%xMinGlob, stepProc, IPT%nFields**IPT%localizationLevel, inverse=.true.)
    !if(rang == 0) call DispCarvalhol(subdivisionCoords, "subdivisionCoords")
    if(rang == 0) write(*,*) "Max Coord = ", subdivisionCoords(:, size(subdivisionCoords,2)) + stepProc

    times(3) = MPI_Wtime() !Organizing Collective Writing

    !Making all realizations
    if(.true.)then
    if(rang == 0) write(*,*) " "
    if(rang == 0) write(*,*) "-> SAMPLING----------------------------------------"
    call wLog("-> SAMPLING----------------------------------------")
    allocate(HDF5Name(product(IPT%nFields**IPT%localizationLevel)))
    do i = 1, product(IPT%nFields**IPT%localizationLevel)
        if(rang == 0) write(*,*) " "
        if(rang == 0) write(*,*) " "
        if(rang == 0) write(*,*)  "-> Making Field ", i
        call wLog("-> Making Field")
        fieldNumber = i;
        if(mod(i, groupMax) == group) then
            !call wLog("Proc")
            !call wLog(rang)
            !call wLog("dealing with field")
            !call wLog(fieldNumber)
            !call wLog("     Trying communication")
            call MPI_BARRIER(groupComm, code)
            call single_realization(IPT, globMSH, writeFiles, outputStyle, sameFolder, &
                                    groupComm, fieldNumber, subdivisionCoords(:,i), stepProc, HDF5Name(i))

        end if
    end do
    end if

    call finalize_MESH(globMSH)

    times(4) = MPI_Wtime() !Generation Time

    !Combining realizations (localization)
    if(product(IPT%nFields) /= 1) then
        if(rang == 0) write(*,*) " "
        if(rang == 0) write(*,*) "-> COMBINING----------------------------------------"
        call wLog("-> COMBINING----------------------------------------")
        call combine_subdivisions(IPT, writeFiles, outputStyle, sameFolder, stepProc, procExtent, overlap, times(1))
    else
        if(rang == 0) write(*,*) "-> NOTHING TO BE COMBINED----------------------------------------"
    end if

    times(5) = MPI_Wtime() !Localization Time

    call MPI_ALLREDUCE (times, all_times, size(times), MPI_DOUBLE_PRECISION, MPI_SUM, IPT%comm,code)

    if(rang == 0) write(*,*) ""
    if(rang == 0) write(*,*) ""
    if(rang == 0) write(*,*) "AVERAGE TIMES (WALL)------------------------ "
    if(rang == 0) write(*,*) "Reading Inputs   = ", (all_times(2) - all_times(1))/dble(IPT%nb_procs)
    if(rang == 0) write(*,*) "Pre Organization = ", (all_times(3) - all_times(2))/dble(IPT%nb_procs)
    if(rang == 0) write(*,*) "Generation       = ", (all_times(4) - all_times(3))/dble(IPT%nb_procs)
    if(rang == 0) write(*,*) "Localization     = ", (all_times(5) - all_times(4))/dble(IPT%nb_procs)
    if(rang == 0) write(*,*) ""

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
        subroutine init_communication(comm_local, comm, rang, nb_procs)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local
            !OUTPUT
            integer, intent(out) :: comm, rang, nb_procs
            !LOCAL
            integer :: code

            !call MPI_Init_thread(MPI_THREAD_MULTIPLE, &provided
            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rang, code)
            call MPI_COMM_SIZE(comm_local, nb_procs, code)

            comm = comm_local

        end subroutine init_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_basic_folders(comm)
            implicit none
            !INPUT
            integer, intent(in) :: comm
            !LOCAL
            logical :: dirExists
            integer, dimension(8) :: date_time
            integer :: code
            character(len=10), dimension(3) :: strings

            !date_time_label
            call date_and_time(strings(1), strings(2), strings(3), date_time)
            results_folder_name = strings(1)(3:8)//"_"//strings(2)(1:6)//"_res"

            if(sameFolder) results_folder_name = "res"

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(sameFolder) log_folder_name     = ".."

            call create_folder(log_folder_name, results_path, rang, comm, compiler)

            call wLog("-> Setting folder path")
            single_path = string_join_many(results_path,"/",results_folder_name)
            call wLog("     single_path = "//trim(single_path))

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
        subroutine end_communication()
            !LOCAL
            integer :: code

            call finalize_log_file()
            call MPI_FINALIZE(code)
        end subroutine end_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine deallocate_all()

            !if(allocated(offset)) deallocate(offset)
            !if(allocated(locDims)) deallocate(locDims)
            !if(allocated(procStart)) deallocate(procStart)
            !if(allocated(periods)) deallocate(periods)
            if(allocated(overlap)) deallocate(overlap)
            if(allocated(HDF5Name)) deallocate(HDF5Name)
            if(allocated(stepProc)) deallocate(stepProc)
            if(allocated(subdivisionCoords)) deallocate(subdivisionCoords)

            call finalize_IPT_RF(IPT)
            call finalize_MESH(globMSH)
            !call finalize_RF(globRDF)

        end subroutine deallocate_all

end program main_RandomField
