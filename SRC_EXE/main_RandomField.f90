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
    !integer            :: i
    !integer            :: rang
    !double precision, dimension(:), allocatable :: gen_times, temp_gen_times
    character(len=200) :: path, logFilePath

    double precision, dimension(5) :: times, all_times

    !DEVEL
    !integer               :: fieldNumber
    !character(len=110), dimension(:)  , allocatable :: HDF5Name
    !type(MESH)            :: globMSH
    !integer               :: group, groupComm, groupMax
    integer               :: code

    !double precision, dimension(:), allocatable :: stepProc, procExtent, overlap
    !double precision, dimension(:,:), allocatable :: subdivisionCoords
    !double precision :: t_bef, t_aft
    type(IPT_RF)  :: IPT

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, IPT%comm, IPT%rang, IPT%nb_procs)
    !rang = IPT%rang
    IPT%writeDataSet = .true.
    IPT%writeUNVinterpolation = .true.
    IPT%sameFolder = .true.
    IPT%outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc
    IPT%delete_intermediate_files = .false.
    IPT%ignoreTillLocLevel = 0 !<1 doesn't affetct the behaviour of the program (for restarts)
    IPT%sampleFields = .true.

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

    if(IPT%rang == 0) write(*,*) "-> MPI_communications started"
    if(IPT%rang == 0) write(*,*) "         running on: "
    if(IPT%rang == 0) call system("pwd")
    if(IPT%rang == 0) write(*,*) "         nb_procs    = ", IPT%nb_procs
    if(IPT%rang == 0) write(*,*) "         outputStyle = ", IPT%outputStyle

    !Initializing folders
    if(IPT%rang == 0) write(*,*)  "-> Initialize Folders"
    call init_basic_folders(IPT%comm)

#ifdef MAKELOG
    if(IPT%rang == 0) write(*,*) "IFDEF MAKELOG DEFINED"

    !Initializing logFiles
    if(IPT%rang == 0) write(*,*)  "-> Initialize logFiles"
    logFilePath = trim(adjustL(&
                      string_join_many(results_path,"/",log_folder_name,"/",log_filename)))
    if(IPT%rang == 0) write(*,*)  " logFilePath = ", trim(adjustL(logFilePath)), "<RANK>"
    call init_log_file(trim(adjustL(logFilePath)), IPT%rang, IPT%log_ID, IPT%nb_procs)
#else
    if(IPT%rang == 0) write(*,*) "IFDEF MAKELOG NOT DEFINED"
#endif

    !READING INPUTS--------------------------------------------------
    !----------------------------------------------------------------
    if(IPT%rang == 0) write(*,*)  "-> Reading inputs"
    call wLog("-> Reading inputs")
    !Reading Mesh---------------------------------------------------
    if(IPT%rang == 0) write(*,*)  "     -> Reading Mesh Input"
    call wLog("     -> Reading Mesh Input")
    path = mesh_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_mesh_input(path, IPT)
    !Reading Generation Input---------------------------------------
    if(IPT%rang == 0) write(*,*)  "     -> Reading Generation Input"
    call wLog("     -> Reading Generation Input")
    path = gen_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_generation_input(path, IPT)
    !Validating Inputs----------------------------------------------
    call wLog("    Validating Inputs")
    call validate_input(IPT)
    if(IPT%rang == 0) call show_IPT_RF(IPT)
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
    if(IPT%rang == 0) write(*,*) "-> REDEFINE xMaxGlob and xMinGlob----------------------------------------"
    call wLog("-> REDEFINE xMaxGlob and xMinGlob----------------------------------------")
    call redefineIPTlimits(IPT, IPT%xMinGlob, IPT%xMaxGlob, IPT%localizationLevel)
    if(IPT%rang == 0) write(*,*) "IPT%xMinGlob"
    if(IPT%rang == 0) write(*,*) IPT%xMinGlob
    if(IPT%rang == 0) write(*,*) "IPT%xMaxGlob"
    if(IPT%rang == 0) write(*,*) IPT%xMaxGlob
    if(IPT%rang == 0) write(*,*) "IPT%localizationLevel"
    if(IPT%rang == 0) write(*,*) IPT%localizationLevel
    call wLog("IPT%xMinGlob")
    call wLog(IPT%xMinGlob)
    call wLog("IPT%xMaxGlob")
    call wLog(IPT%xMaxGlob)
    call wLog("IPT%localizationLevel")
    call wLog(IPT%localizationLevel)

    !Generating random fields
    call make_random_field(IPT, times, product(IPT%nFields**IPT%localizationLevel))

    call MPI_ALLREDUCE (times, all_times, size(times), MPI_DOUBLE_PRECISION, MPI_SUM, IPT%comm,code)

    if(IPT%rang == 0) write(*,*) ""
    if(IPT%rang == 0) write(*,*) ""
    if(IPT%rang == 0) write(*,*) "AVERAGE TIMES (WALL)------------------------ "
    if(IPT%rang == 0) write(*,*) "Reading Inputs   = ", (all_times(2) - all_times(1))/dble(IPT%nb_procs)
    if(IPT%rang == 0) write(*,*) "Pre Organization = ", (all_times(3) - all_times(2))/dble(IPT%nb_procs)
    if(IPT%rang == 0) write(*,*) "Generation       = ", (all_times(4) - all_times(3))/dble(IPT%nb_procs)
    if(IPT%rang == 0) write(*,*) "Localization     = ", (all_times(5) - all_times(4))/dble(IPT%nb_procs)
    if(IPT%rang == 0) write(*,*) ""

	!Deallocating
    call deallocate_all()

    if(IPT%rang == 0) then
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
            integer, dimension(8) :: date_time
            integer :: code
            character(len=10), dimension(3) :: strings

            !date_time_label
            call date_and_time(strings(1), strings(2), strings(3), date_time)
            results_folder_name = strings(1)(3:8)//"_"//strings(2)(1:6)//"_res"

            if(IPT%sameFolder) results_folder_name = "res"

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(IPT%sameFolder) log_folder_name     = ".."

            call create_folder(log_folder_name, results_path, IPT%rang, comm)

            if(IPT%rang == 0) write(*,*) "-> Setting folder path"
            single_path = string_join_many(results_path,"/",results_folder_name)
            if(IPT%rang == 0) write(*,*) "     single_path = "//trim(single_path)

            !create xmf and h5 folders
            !if(writeFiles) then
            path = string_vec_join([results_path,"/",results_folder_name])
            call create_folder("xmf", path, IPT%rang, comm)
            call create_folder("h5", path, IPT%rang, comm)
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
            !call finalize_RF(globRDF)

        end subroutine deallocate_all

end program main_RandomField
