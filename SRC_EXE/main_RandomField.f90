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
    character(len=200) :: path, logFilePath, randFieldFilePath
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    character(len=100) :: BBoxPartFileName
    integer, dimension(:), allocatable :: seed

    !DEVEL
    integer, dimension(2) :: nFields = [2,1] !Number of independent fields in each dimension
    integer               :: nProcPerField = 3, fieldNumber
    logical, dimension(2) :: periods
    integer(HSIZE_T), dimension(2) :: offset, locDims
    double precision, dimension(2) :: procStart
    type(MESH)            :: globMSH
    type(RF)              :: globRDF
    integer               :: group, groupComm, groupNbProcs, groupRang, groupMax
    integer               :: code, nTotalProcs
    logical               :: validProc

    character(len=50) :: attr_Name, dset="samples"
    integer :: hdferr
    integer(HID_T) :: file_id, attr_id, space_id, dset_id, mem_id
    integer(HSIZE_T), dimension(2) :: locShape
    integer(HSIZE_T), dimension(2) :: zero2D

    type(IPT_RF)  :: IPT

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, IPT%comm, IPT%rang, IPT%nb_procs)
    rang = IPT%rang
    IPT%writeDataSet = writeDataSet

    !

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
    !call wLog("        file: "//trim(path))
    call read_mesh_input(path, IPT)
    !Reading Generation Input---------------------------------------
    call wLog("    Reading Generation Input")
    path = gen_input
    path = adjustL(path)
    !call wLog("        file: "//trim(path))
    call read_generation_input(path, IPT)
    !Validating Inputs----------------------------------------------
    call wLog("    Validating Inputs")
    call validate_input(IPT)
#ifdef MAKELOG
    call show_IPT_RF(IPT, forLog_in=.true.)
#endif
    if(rang == 0) call show_IPT_RF(IPT)
    !Initial allocation---------------------------------------------
    call allocate_init()

    !DIVISING FIELD
    periods(:) = .false.
    globMSH%nDim = 2;
    fieldNumber = 1;

    call init_MESH(globMSH, IPT, IPT%comm, product(nFields), IPT%rang)

    call wLog("-> set_procPerDim")
    !call set_procPerDim (globMSH, globMSH%procPerDim)
    globMSH%procPerDim  = nFields
    globMSH%independent = .true.
    call wLog("     globMSH%procPerDim")
    call wLog(globMSH%procPerDim)
    call wLog("-> round_basic_inputs")
    call round_basic_inputs(globMSH, globMSH%xStep, globMSH%overlap)
    call wLog("-> set_global_extremes")
    call set_global_extremes(globMSH, globMSH%xMaxGlob, globMSH%xMinGlob, globMSH%procExtent, procStart)
    call wLog("     globMSH%procExtent = ")
    call wLog(globMSH%procExtent)

    !DEFINING GROUPS AND COMMUNICATORS
    group    = IPT%rang/nProcPerField
    groupMax = IPT%nb_procs/nProcPerField

    call MPI_COMM_SPLIT(IPT%comm, group, IPT%rang, groupComm, code)
    call MPI_COMM_SIZE(groupComm, groupNbProcs, code)
    call MPI_COMM_RANK(groupComm, groupRang, code)

    call wLog("     groupMax = ")
    call wLog(groupMax)
    call wLog("     group = ")
    call wLog(group)
    call wLog("     groupComm = ")
    call wLog(groupComm)
    call wLog("     groupNbProcs = ")
    call wLog(groupNbProcs)
    call wLog("     groupRang = ")
    call wLog(groupRang)

    !Making all realizations
    IPT%nb_procs = groupNbProcs
    allocate(HDF5Name(product(nFields)))
    do i = 1, product(nFields)
        if(rang == 0) write(*,*)  "-> Single realization"
        call wLog("-> Single realization")
        fieldNumber = i;
        if(mod(i, groupMax) == group) then
            call wLog("Proc")
            call wLog(rang)
            call wLog("dealing with field")
            call wLog(fieldNumber)
            call wLog("     Trying communication")
            call MPI_BARRIER(groupComm, code)
            call single_realization(IPT, globMSH, writeFiles, outputStyle, sameFolder, &
                                    nProcPerField, groupComm, fieldNumber, HDF5Name(i))

        end if
    end do

    call finalize_MESH(globMSH)
    call wLog("     Waiting for the other procs")
    call MPI_BARRIER(groupComm, code)

    !Gluing fields together
    call wLog("     PUTTING FIELDS TOGETHER")
    call MPI_COMM_SIZE(IPT%comm, IPT%nb_procs, code)
    call MPI_COMM_RANK(IPT%comm, IPT%rang, code)

    if(product(nFields) > IPT%nb_procs) stop("Too little processors for this number of fields")

    group = 0
    if(IPT%rang < product(nFields)) then
        group = 1
    end if
    call wLog("     group =")
    call wLog(group)

    call MPI_COMM_SPLIT(IPT%comm, group, IPT%rang, groupComm, code)

    IPT%comm       = groupComm
    IPT%nb_procs   = product(nFields)
    IPT%procPerDim = nFields

    if(group == 1) then

        call init_MESH(globMSH, IPT, IPT%comm, IPT%nb_procs, IPT%rang)
        call init_RF(globRDF, IPT, IPT%comm, IPT%nb_procs, IPT%rang)
        globMSH%procPerDim  = nFields
        globMSH%independent = .true.
        call wLog("     globMSH%procPerDim")
        call wLog(globMSH%procPerDim)
        call wLog("-> round_basic_inputs")
        call round_basic_inputs(globMSH, globMSH%xStep, globMSH%overlap)
        call wLog("-> set_global_extremes")
        call set_global_extremes(globMSH, globMSH%xMaxGlob, globMSH%xMinGlob, globMSH%procExtent, procStart)
        call wLog("     globMSH%procExtent = ")
        call wLog(globMSH%procExtent)

        call wLog("-> set_communications_topology")
        call set_communications_topology(globMSH)
        call wLog("-> round_basic_inputs")
        call round_basic_inputs(globMSH, globMSH%xStep, globMSH%overlap)
        call wLog("-> set_global_extremes")
        call set_global_extremes(globMSH, globMSH%xMaxGlob, globMSH%xMinGlob, globMSH%procExtent, globMSH%procStart)
        call wLog("     globMSH%procStart = ")
        call wLog(globMSH%procStart)
        call wLog("     globMSH%procExtent = ")
        call wLog(globMSH%procExtent)
        call wLog("-> set_local_bounding_box")
        call set_local_bounding_box(globMSH,&
                                    globMSH%xMinBound, globMSH%xMaxBound, &
                                    globMSH%xNStep, globMSH%xNTotal, globMSH%origin, validProc)

        !call set_validProcs_comm(validProc, groupComm, globMSH%rang, &
        !                         globMSH%validProc, globRDF%validProc, globMSH%comm, globRDF%comm, &
        !                         globMSH%nb_procs, globRDF%nb_procs, globMSH%rang, globRDF%rang)

        call wLog("-> Setting xPoints")
        call set_xPoints(globMSH, globRDF, globRDF%xPoints_Local)
        call wLog("      maxval(RDF%xPoints,2) = ")
        call wLog(maxval(globRDF%xPoints,2))
        call wLog( "      minval(RDF%xPoints,2) = ")
        call wLog(minval(globRDF%xPoints,2))

        call wLog("-> Reading Random Field")
        BBoxPartFileName = string_join_many(BBoxFileName,"_part",numb2String(IPT%rang+1,5))
        randFieldFilePath = string_join_many(single_path,"/h5/",BBoxPartFileName,".h5")
        call wLog("     Reading on:")
        call wLog(randFieldFilePath)
        call wLog("     Allocating random field")
        call allocate_randField(globRDF, globMSH%xNStep, globRDF%randField_Local)
        call wLog("     shape(globRDF%randField_Local)")
        call wLog(shape(globRDF%randField_Local))
        call wLog("     Read dataset")
        !MATRIX
        call h5open_f(hdferr) ! Initialize FORTRAN interface.
        call h5fopen_f(trim(randFieldFilePath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
        if(hdferr /= 0) stop("ERROR OPENING FILE inside read_RF_h5_File_Table")
        call h5dopen_f(file_id, trim(dset), dset_id, hdferr)! Open Dataset
        locShape = shape(globRDF%randField_Local)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, globRDF%randField_Local, locShape,  hdferr)
        call h5dclose_f(dset_id, hdferr) !Close Dataset
        call h5fclose_f(file_id, hdferr) ! Close the file.
        call h5close_f(hdferr) ! Close FORTRAN interface.
!
!        if(RDF%nb_procs > 1) then
!            call wLog("")
!            call wLog("GENERATING OVERLAP")
!            call wLog("-------------------------------")
!            if(RDF%rang == 0) write(*,*)"GENERATING OVERLAP"
!            if(RDF%rang == 0) write(*,*) "-------------------------------"
!            call wLog("")
!
!            !RDF%randField = 1.0 ! For Tests
!            if(RDF%rang == 0) write(*,*) "    ->Applying Weighting Functions"
!            call wLog("    ->Applying Weighting Functions on Field")
!            call applyWeightingFunctions_OnMatrix(globRDF, globMSH, partitionType)
!            if(RDF%rang == 0) write(*,*) "    ->addNeighboursFields"
!            call wLog("    ->addNeighboursFields")
!            call addNeighboursFields(globRDF, globMSH)
!        end if


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

            if(sameFolder) results_folder_name = "res"

            call MPI_BARRIER (comm ,code) !Necessary because each proc can have a different time
            call MPI_BCAST (results_folder_name, 100, MPI_CHARACTER, 0, comm, code)

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(sameFolder) log_folder_name     = ".."

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
            call finalize_MESH(globMSH)
            call finalize_RF(globRDF)

        end subroutine deallocate_all

end program main_RandomField
