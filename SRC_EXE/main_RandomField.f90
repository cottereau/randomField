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


	!LOCAL VARIABLES
    character(len=200) :: path, logFilePath
    double precision, dimension(7) :: times, all_times
    integer               :: code

    !INPUT VARIABLES
    type(IPT_RF)  :: IPT_Temp !Only for practicity purposes
    type(IPT_RF)  :: IPT !The one that shoud be initialized when calling from an external program

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, IPT_Temp%comm, IPT_Temp%rang, IPT_Temp%nb_procs)

    !Options
    IPT_Temp%writeDataSet = .true.
    IPT_Temp%writeUNVinterpolation = .true.
    IPT_Temp%sameFolder = .true.
    IPT_Temp%outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc
    IPT_Temp%delete_intermediate_files = .false.
    IPT_Temp%sampleFields = .true.
    !Ps: IPT_Temp%ignoreTillLocLevel Defined after Reading Inputs

    times(1) = MPI_Wtime() !Initial Time

    if(IPT_Temp%rang == 0)then
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

    if(IPT_Temp%rang == 0) write(*,*) "-> MPI_communications started"
    if(IPT_Temp%rang == 0) write(*,*) "         running on: "
    if(IPT_Temp%rang == 0) call system("pwd")
    if(IPT_Temp%rang == 0) write(*,*) "         nb_procs    = ", IPT_Temp%nb_procs
    if(IPT_Temp%rang == 0) write(*,*) "         outputStyle = ", IPT_Temp%outputStyle

    !Initializing folders
    if(IPT_Temp%rang == 0) write(*,*)  "-> Initialize Folders"
    call init_basic_folders(IPT_Temp%comm)

#ifdef MAKELOG
    if(IPT_Temp%rang == 0) write(*,*) "IFDEF MAKELOG DEFINED"

    !Initializing logFiles
    if(IPT_Temp%rang == 0) write(*,*)  "-> Initialize logFiles"
    logFilePath = trim(adjustL(&
                      string_join_many(results_path,"/",log_folder_name,"/",log_filename)))
    if(IPT_Temp%rang == 0) write(*,*)  " logFilePath = ", trim(adjustL(logFilePath)), "<RANK>"
    call init_log_file(trim(adjustL(logFilePath)), IPT_Temp%rang, IPT_Temp%log_ID, IPT_Temp%nb_procs)
#else
    if(IPT_Temp%rang == 0) write(*,*) "IFDEF MAKELOG NOT DEFINED"
#endif

    !READING INPUTS--------------------------------------------------
    !----------------------------------------------------------------
    if(IPT_Temp%rang == 0) write(*,*)  "-> Reading inputs"
    call wLog("-> Reading inputs")
    !Reading Mesh---------------------------------------------------
    if(IPT_Temp%rang == 0) write(*,*)  "     -> Reading Mesh Input"
    call wLog("     -> Reading Mesh Input")
    path = mesh_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_mesh_input(path, IPT_Temp)
    !Reading Generation Input---------------------------------------
    if(IPT_Temp%rang == 0) write(*,*)  "     -> Reading Generation Input"
    call wLog("     -> Reading Generation Input")
    path = gen_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_generation_input(path, IPT_Temp)

    IPT_Temp%ignoreTillLocLevel = IPT_Temp%localizationLevel - 1 !<1 doesn't affetct the behaviour of the program (for restarts)

    !Validating Inputs----------------------------------------------
    if(IPT_Temp%rang == 0) write(*,*)  "     -> Validating Input (IPT_Temp)"
    call wLog("    Validating Inputs (IPT_Temp)")
    call validate_input(IPT_Temp)
    if(IPT_Temp%rang == 0) call show_IPT_RF(IPT_Temp, "IPT_Temp")
#ifdef MAKELOG
    call show_IPT_RF(IPT_Temp, forLog_in=.true.)
#endif

    times(2) = MPI_Wtime() !Reading Inputs


    !Initial allocation---------------------------------------------
    call allocate_init()

    if(IPT_Temp%nProcPerField > IPT_Temp%nb_procs) then
        write(*,*) "ERROR, IPT_Temp%nProcPerField cannot be superior to IPT_Temp%nb_procs"
        stop(" ")
    end if

    !Initialize Inputs
    if(IPT_Temp%rang == 0) write(*,*)  " "
    if(IPT_Temp%rang == 0) write(*,*)  " "
    if(IPT_Temp%rang == 0) write(*,*)  "-> Initializing Input (IPT)"
    call init_IPT_RF(&
        IPT, &
        log_ID = IPT_Temp%log_ID, &
        comm = IPT_Temp%comm, &
        rang = IPT_Temp%rang, &
        nb_procs = IPT_Temp%nb_procs, &
        nDim = IPT_Temp%nDim_gen, &
        meshMod = IPT_Temp%meshMod, &
        xMaxGlob_in = IPT_Temp%xMaxGlob_in, &
        xMinGlob_in = IPT_Temp%xMinGlob_in, &
        pointsPerCorrL = IPT_Temp%pointsPerCorrL, &
        procPerDim = IPT_Temp%procPerDim, &
        coordList_local = IPT_Temp%coordList_local, &
        connectList_local = IPT_Temp%connectList_local, &
        monotype = IPT_Temp%monotype, &
        unv = IPT_Temp%unv, &
        unv_path = IPT_Temp%unv_path, &
        fieldAvg = IPT_Temp%fieldAvg, &
        fieldVar = IPT_Temp%fieldVar, &
        corrL_in = IPT_Temp%corrL_in, &
        overlap_in = IPT_Temp%overlap_in, &
        corrMod = IPT_Temp%corrMod, &
        margiFirst = IPT_Temp%margiFirst, &
        method = IPT_Temp%method, &
        Nmc = IPT_Temp%Nmc, &
        seedStart = IPT_Temp%seedStart, &
        nProcPerField = IPT_Temp%nProcPerField, &
        nFields = IPT_Temp%nFields, &
        localizationLevel = IPT_Temp%localizationLevel, &
        writeDataSet = IPT_Temp%writeDataSet, &
        sameFolder = IPT_Temp%sameFolder, &
        outputStyle = IPT_Temp%outputStyle, &
        delete_intermediate_files = IPT_Temp%delete_intermediate_files, &
        ignoreTillLocLevel = IPT_Temp%ignoreTillLocLevel, &
        sampleFields = IPT_Temp%sampleFields, &
        writeUNVinterpolation = IPT_Temp%writeUNVinterpolation)

    !Generating random fields
    call make_random_field(IPT, times)

    call MPI_ALLREDUCE (times, all_times, size(times), MPI_DOUBLE_PRECISION, MPI_SUM, IPT%comm,code)

    if(IPT_Temp%rang == 0) write(*,*) ""
    if(IPT_Temp%rang == 0) write(*,*) ""
    if(IPT_Temp%rang == 0) write(*,*) "AVERAGE TIMES (WALL)------------------------ "
    if(IPT_Temp%rang == 0) write(*,*) "Reading Inputs   = ", (all_times(2) - all_times(1))/dble(IPT_Temp%nb_procs)
    if(IPT_Temp%rang == 0) write(*,*) "Pre Organization = ", (all_times(3) - all_times(2))/dble(IPT_Temp%nb_procs)
    if(IPT_Temp%rang == 0) write(*,*) "Generation       = ", (all_times(4) - all_times(3))/dble(IPT_Temp%nb_procs)
    if(IPT_Temp%rang == 0) write(*,*) "Localization Int = ", (all_times(5) - all_times(4))/dble(IPT_Temp%nb_procs)
    if(IPT_Temp%rang == 0) write(*,*) "Localization Ext = ", (all_times(6) - all_times(5))/dble(IPT_Temp%nb_procs)
    if(IPT_Temp%rang == 0) write(*,*) "Writing Files    = ", (all_times(7) - all_times(6))/dble(IPT_Temp%nb_procs)
    if(IPT_Temp%rang == 0) write(*,*) ""

	!Deallocating
    call deallocate_all()

    if(IPT_Temp%rang == 0) then
        write(*,*) ""
        write(*,*) "---------------------------------------------------------------------";
        write(*,*) "-----------------END RANDOM FIELD LIBRARY TEST-----------------------";
        write(*,*) "---------------------------------------------------------------------";
        write(*,*) ""
    end if

    call finalize_IPT_RF(IPT_Temp)

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
            integer, dimension(8) :: date_time
            integer :: code
            character(len=10), dimension(3) :: strings

            !date_time_label
            call date_and_time(strings(1), strings(2), strings(3), date_time)
            results_folder_name = strings(1)(3:8)//"_"//strings(2)(1:6)//"_res"

            if(IPT_Temp%sameFolder) results_folder_name = "res"

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(IPT_Temp%sameFolder) log_folder_name     = ".."

            call create_folder(log_folder_name, results_path, IPT_Temp%rang, comm)

            if(IPT_Temp%rang == 0) write(*,*) "-> Setting folder path"
            single_path = string_join_many(results_path,"/",results_folder_name)
            if(IPT_Temp%rang == 0) write(*,*) "     single_path = "//trim(single_path)

            !create xmf and h5 folders
            !if(writeFiles) then
            path = string_vec_join([results_path,"/",results_folder_name])
            call create_folder("xmf", path, IPT_Temp%rang, comm)
            call create_folder("h5", path, IPT_Temp%rang, comm)
            !end if

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

            call finalize_IPT_RF(IPT)

        end subroutine deallocate_all

end program main_RandomField
