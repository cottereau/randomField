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
    logical :: sameFolder = .true.
    integer :: outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc

	!LOCAL VARIABLES
    logical            :: file_exist
    integer            :: i, baseStep, nIter
    integer            :: rang
    character(len=30)  :: rangChar;
    character(len=200) :: path, logFilePath
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    integer, dimension(:), allocatable :: seed

    type(IPT_RF)  :: IPT

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, IPT%comm, IPT%rang, IPT%nb_procs)
    rang = IPT%rang

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
    call wLog("    Reading Mesh Input")
    path = mesh_input
    path = adjustL(path)
    call wLog("        file: "//trim(path))
    call read_mesh_input(path, IPT)
    !Reading Generation Input---------------------------------------
    call wLog("    Reading Generation Input")
    path = gen_input
    path = adjustL(path)
    call wLog("        file: "//trim(path))
    call read_generation_input(path, IPT)
    !Validating Inputs----------------------------------------------
    call wLog("    Validating Inputs")
    call validate_input(IPT)
#ifdef MAKELOG
    call show_IPT_RF(IPT, forLog_in=.true.)
#endif
    !Initial allocation---------------------------------------------
    call allocate_init()

    !SINGLE REALIZATION
    if(rang == 0) write(*,*)  "-> Single realization"
    call wLog("-> Single realization")
    call single_realization(IPT, writeFiles, outputStyle, sameFolder)

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

            if(sameFolder) results_folder_name = "res" !ONLY FOR TESTS

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(sameFolder) log_folder_name     = ".."

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
