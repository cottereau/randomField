module sample_RF

    use displayCarvalhol
    use write_Log_File
    use math_RF
    use constants_RF
    use mpi
    use writeResultFile_RF
    use type_RF
    use type_MESH
    use type_inputRF
    use localization_RF
    use common_variables_RF
    use randomFieldND
    use mesh_RF
    use topography_RF
    !use calls_RF

    implicit none

contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------

        subroutine redefineIPTlimits(IPT, xMinGlob, xMaxGlob, localizationLevel)

            implicit none
            !INPUT
            type(IPT_RF), intent(in)  :: IPT
            !OUTPUT
            double precision, dimension(:), intent(out) :: xMinGlob
            double precision, dimension(:), intent(out) :: xMaxGlob
            integer, intent(inout) :: localizationLevel

            !LOCAL
            type(MESH) :: globMSH

            double precision, dimension(IPT%nDim_mesh) :: xRangeTotal, xMaxTotal, xMinTotal, xRangeGlob
            double precision, dimension(IPT%nDim_mesh) :: overlap
            logical :: locLevelOK

            call wLog("-> Redefining input xMinGlob, xMinGlob based on localization level")
            locLevelOK = .false.


            call init_MESH(globMSH, IPT, IPT%comm, IPT%rang)

            call wLog("->  set_procPerDim")
            call wLog("     globMSH%procPerDim")
            call wLog(globMSH%procPerDim)
            call wLog("-> round_basic_inputs")
            call round_basic_inputs(globMSH, globMSH%xStep, globMSH%overlap)
            call wLog("-> set_global_extremes")
            call set_global_extremes(globMSH, globMSH%xMaxGlob, globMSH%xMinGlob, globMSH%procExtent, &
                                     globMSH%procStart)


            xMaxTotal   = globMSH%xMaxGlob
            xMinTotal   = globMSH%xMinGlob
            xRangeTotal = xMaxTotal - xMinTotal
            overlap     = globMSH%overlap

            call wLog("xMaxTotal = ")
            call wLog(xMaxTotal)
            call wLog("xMinTotal = ")
            call wLog(xMinTotal)
            call wLog("xRangeTotal = ")
            call wLog(xRangeTotal)
            call finalize_MESH(globMSH)

            do while(.not. locLevelOK)

                xRangeGlob = (xRangeTotal + overlap*dble((IPT%nFields**(IPT%localizationLevel-1)) - 1))&
                             /dble(IPT%nFields**(IPT%localizationLevel-1))
                call wLog("xRangeGlob = ")
                call wLog(xRangeGlob)

                if(all(xRangeGlob < 0)) then
                    if(IPT%rang == 0) then
                        write(*,*) "Localization level not adapted to mesh requirements, will be reduced"
                        write(*,*) "OLD localizationLevel = ", localizationLevel
                    end if
                    call wLog("Localization level not adapted to mesh requirements, will be reduced")
                    call wLog("OLD localizationLevel = ")
                    call wLog(localizationLevel)

                    localizationLevel = localizationLevel -1

                    if(IPT%rang == 0) write(*,*) "NEW localizationLevel = ", localizationLevel
                    call wLog("NEW localizationLevel = ")
                    call wLog(localizationLevel)

                    if(localizationLevel < 0) then
                        locLevelOK = .true.
                        write(*,*) "ERROR NEGATIVE LOCALIZATION ON RANK = ", IPT%rang
                        call wLog("ERROR NEGATIVE LOCALIZATION ON RANK = ")
                        call wLog(IPT%rang)
                        call wLog(localizationLevel)
                        stop ("ERROR inside redefineIPTlimits, localization level smaller than 0")
                    end if
                else
                    locLevelOK = .true.
                    xMinGlob   = xMinTotal
                    xMaxGlob   = xMinTotal + xRangeGlob


                    call wLog("xMinGlob")
                    call wLog(xMinGlob)
                    call wLog("xMaxGlob")
                    call wLog(xMaxGlob)

                end if

            end do

        end subroutine redefineIPTlimits

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------

        subroutine build_subdivisions(IPT, globMSH, groupMax, group, groupComm, stepProc, procExtent, overlap)

            implicit none
            !INPUT
            type(IPT_RF), intent(in)  :: IPT
            !OUTPUT
            type(MESH), intent(out)  :: globMSH
            integer, intent(out) :: group, groupMax, groupComm
            double precision, dimension(:), intent(out) :: stepProc
            double precision, dimension(:), intent(out) :: procExtent
            double precision, dimension(:), intent(out) :: overlap

            !LOCAL
            integer :: code

            call init_MESH(globMSH, IPT, IPT%comm, IPT%rang)

            call wLog("->  set_procPerDim")
            !globMSH%independent = .true.
            call wLog("     globMSH%procPerDim")
            call wLog(globMSH%procPerDim)
            call wLog("-> round_basic_inputs")
            call round_basic_inputs(globMSH, globMSH%xStep, globMSH%overlap)
            overlap = globMSH%overlap

            call wLog("-> set_global_extremes")
            call set_global_extremes(globMSH, globMSH%xMaxGlob, globMSH%xMinGlob, globMSH%procExtent, &
                                     globMSH%procStart, stepProc)
            call wLog("     globMSH%procExtent = ")
            call wLog(globMSH%procExtent)

            procExtent = globMSH%procExtent

            !DEFINING GROUPS AND COMMUNICATORS
            group    = IPT%rang/IPT%nProcPerField
            groupMax = IPT%nb_procs/IPT%nProcPerField

            call MPI_COMM_SPLIT(IPT%comm, group, IPT%rang, groupComm, code)

            call wLog("     group = ")
            call wLog(group)
            call wLog("     groupMax = ")
            call wLog(groupMax)
            call wLog("     groupComm = ")
            call wLog(groupComm)

        end subroutine build_subdivisions

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------

        subroutine single_realization(IPT, globMSH, &
                                      fieldComm, fieldNumber, subdivisionStart, stepProc, h5fullPath)

            implicit none
            !INPUT
            type(IPT_RF), intent(in) :: IPT
            type(MESH)  , intent(in) :: globMSH
            integer, intent(in) :: fieldComm, fieldNumber
            double precision, dimension(:), intent(in) :: subdivisionStart, stepProc

            !LOCAL
            type(RF)      :: RDF
            type(MESH)    :: MSH

            double precision :: t1, t2, t3
            double precision :: all_t1, all_t2, all_t3
            integer :: code
            integer :: all_xNTotal
            integer :: i
            logical :: validProc
            integer :: rang
            integer :: newNbProcs, newRang
            character(len=200) :: BBoxPartFileName, BBoxPath
            character(len=*)   :: h5fullPath

            integer, dimension(IPT%nDim_gen) :: globCoord


            validProc = .true.
            call MPI_COMM_RANK(fieldComm, newRang, code)
            call MPI_COMM_SIZE(fieldComm, newNbProcs, code)

            call init_MESH(MSH, IPT, fieldComm, newRang, newNbProcs)
            call init_RF(RDF, IPT, fieldComm, newNbProcs, newRang)

            !Copy from globMSH
            MSH%xStep    = globMSH%xStep
            MSH%overlap  = globMSH%overlap
            MSH%xMinGlob = subdivisionStart
            MSH%xMaxGlob = MSH%xMinGlob + globMSH%procExtent
            MSH%procExtent = globMSH%procExtent

            MSH%procPerDim(:) = 1
            MSH%procPerDim(MSH%nDim) = newNbProcs
            MSH%coords = 0
            MSH%coords(MSH%nDim) = newRang
            globCoord = nint((subdivisionStart - globMSH%xMinGlob)/stepProc)
            call wLog("  globCoord = ")
            call wLog(globCoord)

            call wLog("-> set_local_bounding_box")

            call set_local_bounding_box(MSH,&
                                        MSH%xMinBound, MSH%xMaxBound, &
                                        MSH%xNStep, MSH%xNTotal, MSH%origin, validProc, &
                                        localization = .false.)
            rang = MSH%rang
            call set_validProcs_comm(validProc, fieldComm, rang, &
                                     MSH%validProc, RDF%validProc, MSH%comm, RDF%comm, &
                                     MSH%nb_procs, RDF%nb_procs, MSH%rang, RDF%rang)

            call wLog("MSH%xNTotal = ")
            call wLog(MSH%xNTotal)

            if(validProc) then

                call wLog("-> Initializing Random Seed")

                if(RDF%seedStart >= 0) then
                    !Deterministic Seed
                    RDF%seedStart = RDF%seedStart + fieldNumber
                    call calculate_random_seed(RDF%seed, RDF%seedStart)
                else
                    !Stochastic Seed
                    if(RDF%rang == 0) call calculate_random_seed(RDF%seed, RDF%seedStart)
                    call MPI_BCAST (RDF%seed, size(RDF%seed), MPI_INTEGER, 0, RDF%comm, code)
                    RDF%seed = RDF%seed + fieldNumber
                end if

                call init_random_seed(RDF%seed)

                call wLog("      RDF%seed = ")
                call wLog(RDF%seed)
                call wLog(" ")

                call wLog("-> Setting xPoints")
                call set_xPoints(MSH, RDF, RDF%xPoints_Local)
                call wLog("      maxval(RDF%xPoints,2) = ")
                call wLog(maxval(RDF%xPoints,2))
                call wLog( "      minval(RDF%xPoints,2) = ")
                call wLog(minval(RDF%xPoints,2))

                !i = size(RDF%xPoints,2)
                !if(i>50) i = 50
                !call dispCarvalhol(transpose(RDF%xPoints(:,1:i)), "transpose(RDF%xPoints)", "(F20.5)",unit_in = RDF%log_ID)

                !Discovering the total number of points in all procs
                call MPI_ALLREDUCE (MSH%xNTotal, all_xNTotal,1,MPI_INTEGER, &
                                    MPI_SUM,MSH%comm,code)
                !Getting Initial Time
                t1 = MPI_Wtime();
                call MPI_ALLREDUCE (t1, all_t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MSH%comm,code)
                !if(RDF%rang == 0) write(*,*) "Time Zero = ", all_t1

                call wLog("-> Generating Random Field")
                call wLog("     Allocating random field")
                call allocate_randField(RDF, MSH%xNStep, RDF%randField_Local)
                call wLog("     shape(RDF%randField)")
                call wLog(shape(RDF%randField))
                call wLog("     Calculating sample")
                call create_RF_Unstruct_Init (RDF, MSH)

                call wLog("      maxval(RDF%randField,1) = ")
                call wLog(maxval(RDF%randField,1))
                call wLog( "      minval(RDF%randField,1) = ")
                call wLog(minval(RDF%randField,1))


                t2 = MPI_Wtime();
                call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MSH%comm,code)
                RDF%gen_CPU_Time = all_t2 - all_t1
                !if(RDF%rang == 0) write(*,*) "Generation CPU Time = ", all_t2 - all_t1
                call wLog ("    Generation CPU Time (s)")
                call wLog (RDF%gen_CPU_Time)
                all_t3 = -1.0D0

                !RDF%randField = 1.0D0 ! For Tests

                call wLog("-> Writing XMF and hdf5 files");

                if(.true.) then
                    call wLog("IPT%outputStyle");
                    call wLog(IPT%outputStyle);
                    if(RDF%rang == 0) write(*,*) "-> Writing 'Bounding Box' XMF and hdf5 files"
                    BBoxPartFileName = string_join_many(BBoxFileName,"_L00")
                    do i = 1, size(globCoord)
                        BBoxPartFileName = string_join_many(BBoxPartFileName,"_",numb2String(globCoord(i)+1,3))
                    end do

                    if(IPT%outputStyle==1) then
                        call wLog("   (Parallel)");
                        call wLog("minval(RDF%randField,1) =")
                        call wLog(minval(RDF%randField,1))
                        call wLog("maxval(RDF%randField,1) =")
                        call wLog(maxval(RDF%randField,1))
                        call write_Mono_XMF_h5(RDF, MSH, BBoxPartFileName, RDF%rang, single_path, &
                                               MSH%comm, ["_Part"], [fieldNumber], fieldNumber, style=IPT%outputStyle, &
                                               HDF5FullPath = BBoxPath, writeDataSet = IPT%writeDataSet, localization=.false.)
                    else
                        call wLog("   (Per Proc)");
                        call write_Mono_XMF_h5(RDF, MSH, BBoxPartFileName, RDF%rang, single_path, &
                                               MSH%comm, ["_Part"], [RDF%rang], 0, style=IPT%outputStyle, &
                                               HDF5FullPath = BBoxPath, writeDataSet = IPT%writeDataSet, localization=.false.)

                    end if

                    t3 = MPI_Wtime();
                    call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MSH%comm,code)
                    call wLog ("    Writing Files CPU Time (s)")
                    call wLog (all_t3 - all_t2)

                    h5fullPath = BBoxPath
                end if
            end if

            if(IPT%rang == 0) call write_stat_input("./stat_input", BBoxPath)

            call finalize_MESH(MSH)
            call finalize_RF(RDF)

        end subroutine single_realization

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine combine_subdivisions(IPT, outputStyle, &
                                        stepProc, procExtent, overlap, &
                                        t_ref, t_prep, t_gen, gen_times, nGenGroups, &
                                        delete_intermediate_files, ignoreTillLocLevel, &
                                        finalFieldPath)
            implicit none
            !INPUT
            type(IPT_RF), intent(in)  :: IPT
            integer, intent(in) :: outputStyle!1: parallel hdf5, 2: hdf5 per proc
            double precision, dimension(:), intent(in) :: stepProc, procExtent, overlap
            double precision, intent(in) :: t_ref, t_prep, t_gen
            double precision, dimension(:), intent(in) :: gen_times
            integer, intent(in) :: nGenGroups
            logical, intent(in) :: delete_intermediate_files
            integer, intent(in) :: ignoreTillLocLevel

            !OUTPUT
            character(len=200) :: finalFieldPath

            !LOCAL
            logical :: writeFiles = .true.
            type(MESH) :: globMSH
            type(RF)   :: globRDF
            integer :: group, fieldNumber, groupComm, groupMax
            integer :: code
            character(len=200) :: randFieldFilePath, XMFFilePath
            character(len=200) :: BBoxPartFileName
            logical :: validProc
            integer :: partitionType=1
            integer :: locLevel, locIter
            integer(HSIZE_T), dimension(2) :: locShape
            type(IPT_RF) :: newIPT
            double precision, dimension(IPT%nDim_mesh, product(IPT%nFields)) :: subdivisionCoords
            double precision, dimension(:, :), allocatable :: offsetCoords
            logical :: last_round
            integer, dimension(IPT%nDim_mesh) :: procCoord, locStep
            integer, dimension(IPT%nDim_mesh) :: totalFieldPerDim, nFields_level, nFields_level_before
            integer :: i, validProcGroup, validProcComm, prodNFields
            double precision :: t_loc, t_trans
            double precision :: all_t_loc, all_t_prep, all_t_gen, all_t_trans
            integer(kind=8) :: xNTotal
            integer, dimension(IPT%nDim_mesh) :: minPosProc, maxPosProc
            logical :: writeDataSet

            !H5 LOCAL
            character(len=50) :: dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, dset_id
            double precision, dimension(IPT%nDim_mesh) :: ones

            !Gluing fields together

            call MPI_BARRIER(IPT%comm, code)

            ones = 1.0D0
            call setGrid(subdivisionCoords, 0*ones, ones, IPT%nFields, inverse=.true.)

            if(product(IPT%nFields) > IPT%nb_procs) stop("Too little processors for this number of fields")

            prodNFields = product(IPT%nFields)
            group    = IPT%rang/prodNFields
            groupMax = IPT%nb_procs/prodNFields

            call wLog("     IPT%rang =")
            call wLog(IPT%rang)
            call wLog("     group =")
            call wLog(group)
            call wLog(" IPT%nb_procs =")
            call wLog(IPT%nb_procs)


            validProcGroup = 0
            if(group < groupMax) validProcGroup = 1

            call MPI_COMM_SPLIT(IPT%comm, validProcGroup, IPT%rang, validProcComm, code)
            call MPI_COMM_SPLIT(IPT%comm, group, IPT%rang, groupComm, code)

            call wLog("TESTING---------------")
            call wLog("IPT%comm")
            call wLog(IPT%comm)
            call wLog("group")
            call wLog(group)
            call wLog("IPT%rang")
            call wLog(IPT%rang)
            call wLog("groupComm")
            call wLog(groupComm)

            if(validProcGroup == 1) then

                call init_IPT_RF(newIPT, IPT%nDim_mesh, IPT%log_ID, IPT%rang)
                call copy_IPT_RF(newIPT, IPT)

                newIPT%overlap = overlap
                !newIPT%independent = .true.
                newIPT%comm = groupComm
                call MPI_COMM_RANK(groupComm, newIPT%rang, code)
                call MPI_COMM_SIZE(groupComm, newIPT%nb_procs, code)
                totalFieldPerDim = IPT%nFields**(IPT%localizationLevel)
                !if (IPT%rang == 0) write(*,*) "totalFieldPerDim = ", totalFieldPerDim

                do locLevel = 1, IPT%localizationLevel

                    if(IPT%rang == 0) write(*,*)  "  locLevel = ", locLevel, "-------------"
                    call wLog("     locLevel -----------------------------------")
                    call wLog(locLevel)

                    if(locLevel <= ignoreTillLocLevel) then
                        if(IPT%rang == 0) write(*,*)  "Ignoring This Localization Level"
                        call wLog("Ignoring This Localization Level")
                        cycle
                    end if

                    if(allocated(offsetCoords)) deallocate(offsetCoords)
                    nFields_level = IPT%nFields**(locLevel)
                    nFields_level_before = IPT%nFields**(locLevel-1)
                    locStep = totalFieldPerDim/nFields_level
                    if (IPT%rang == 0) write(*,*) "    locStep = ", locStep
                    !if (IPT%rang == 0) write(*,*) " "
                    if(allocated(offsetCoords)) deallocate(offsetCoords)
                    allocate(offsetCoords(IPT%nDim_mesh, product(locStep)))
                    call setGrid(offsetCoords, 0*ones, dble(nFields_level)*stepProc, locStep, inverse=.true.)


                    !if (IPT%rang == 0) call DispCarvalhol(offsetCoords, "offsetCoords")

                    do locIter = 1, size(offsetCoords,2)

                        if(mod(locIter-1, groupMax) /= group) cycle

                        if(IPT%rang == 0) write(*,*)  "         locIter = ", locIter
                        call wLog("     Waiting for the other procs")
                        call wLog("  locIter = ")
                        call wLog(locIter)
                        call MPI_BARRIER(groupComm, code)


                        call wLog(" subdivisionCoords(:,IPT%rang+1)")
                        call wLog(subdivisionCoords(:,newIPT%rang+1))
                        call wLog(" (offsetCoords(:,locIter)/(stepProc*nFields_level_before))")
                        call wLog((offsetCoords(:,locIter)/(stepProc*nFields_level_before)))

                        procCoord = nint(subdivisionCoords(:,newIPT%rang+1)) + &
                                    nint(offsetCoords(:,locIter)/(stepProc*nFields_level_before))
                        call wLog("     procCoord -----------------------------------")
                        call wLog(procCoord)

                        newIPT%xMinGlob = IPT%xMinGlob + (offsetCoords(:,locIter))
                        newIPT%xMaxGlob = newIPT%xMinGlob &
                                          + dble(nFields_level)*procExtent &
                                          - dble(nFields_level-1)*(procExtent-stepProc)

                        call wLog("     newIPT%xMinGlob")
                        call wLog(newIPT%xMinGlob)
                        call wLog("     newIPT%xMaxGlob")
                        call wLog(newIPT%xMaxGlob)

                        call init_MESH(globMSH, newIPT, groupComm, newIPT%rang)
                        call init_RF(globRDF, newIPT, groupComm, product(IPT%nFields), newIPT%rang)

                        call wLog("newIPT%nb_procs")
                        call wLog(newIPT%nb_procs)
                        call wLog("RDF%nb_procs")
                        call wLog(globRDF%nb_procs)

                        call wLog("-> set_communications_topology")
                        call set_communications_topology(globMSH, globMSH%coords, globMSH%neigh, &
                                                         globMSH%neighShift, globMSH%considerNeighbour, &
                                                         globMSH%mappingFromShift, globMSH%op_neigh)

                        call wLog("-> round_basic_inputs")
                        call round_basic_inputs(globMSH, globMSH%xStep, globMSH%overlap)

                        call wLog("-> set_global_extremes")
                        call set_global_extremes(globMSH, globMSH%xMaxGlob, globMSH%xMinGlob, globMSH%procExtent, globMSH%procStart)


                        call wLog("     globMSH%procExtent = ")
                        call wLog(globMSH%procExtent)
                        call wLog("     globMSH%xMaxGlob = ")
                        call wLog(globMSH%xMaxGlob)
                        call wLog("     globMSH%procStart = ")
                        call wLog(globMSH%procStart)
                        call wLog("     globMSH%procExtent = ")
                        call wLog(globMSH%procExtent)
                        call wLog("-> set_local_bounding_box")
                        call set_local_bounding_box(globMSH,&
                                                    globMSH%xMinBound, globMSH%xMaxBound, &
                                                    globMSH%xNStep, globMSH%xNTotal, globMSH%origin, validProc, localization = .true.)
                        !call wLog("     globMSH%xMinBound = ")
                        !call wLog(globMSH%xMinBound)
                        !call wLog("     globMSH%xMaxBound = ")
                        !call wLog(globMSH%xMaxBound)
                        !call wLog("     globMSH%xNStep = ")
                        !call wLog(globMSH%xNStep)
                        !call wLog("     globMSH%xNTotal = ")
                        !call wLog(globMSH%xNTotal)
                        !call wLog("     globMSH%origin = ")
                        !call wLog(globMSH%origin)


                        call wLog("-> set_overlap_geometry")
                        call set_overlap_geometry (globMSH, globMSH%xMinInt, globMSH%xMaxInt, globMSH%xMinExt, globMSH%xMaxExt, &
                                                   globMSH%xMaxNeigh, globMSH%xMinNeigh, globMSH%xOrNeigh,                      &
                                                   globMSH%nOvlpMax, globMSH%nOvlpPoints)

                        call wLog("-> Setting xPoints")
                        call set_xPoints(globMSH, globRDF, globRDF%xPoints_Local)
                        call wLog("      maxval(RDF%xPoints,2) = ")
                        call wLog(maxval(globRDF%xPoints,2))
                        call wLog( "      minval(RDF%xPoints,2) = ")
                        call wLog(minval(globRDF%xPoints,2))

                        if(globRDF%rang == 0) write(*,*)"READING RANDOM FIELD IN H5 FILE"
                        if(globRDF%rang == 0) write(*,*) "-------------------------------"

                        call wLog("-> Reading Random Field")

                        BBoxPartFileName = string_join_many(BBoxFileName,"_L",numb2String(locLevel-1,2))
                        do i = 1, size(procCoord)
                            BBoxPartFileName = string_join_many(BBoxPartFileName,"_",numb2String(procCoord(i)+1,3))
                        end do
                        randFieldFilePath = string_join_many(single_path,"/h5/",BBoxPartFileName,".h5")
                        XMFFilePath = string_join_many(single_path,"/xmf/",BBoxPartFileName,".xmf")

                        call wLog("     single_path =")
                        call wLog(single_path)
                        call wLog("     BBoxPartFileName =")
                        call wLog(BBoxPartFileName)
                        call wLog("     File:")
                        call wLog(randFieldFilePath)

                        call wLog("     Allocating random field")
                        call allocate_randField(globRDF, globMSH%xNStep, globRDF%randField_Local)
                        call wLog("     shape(globRDF%randField_Local)")
                        call wLog(shape(globRDF%randField_Local))
                        call wLog("     Read dataset")

                        !HDF5
                        call h5open_f(hdferr) ! Initialize FORTRAN interface.
                        call h5fopen_f(trim(randFieldFilePath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
                        if(hdferr /= 0) then
                            write(*,*) "Inside rang ", IPT%rang
                            stop("ERROR OPENING FILE inside read_RF_h5_File_Table")
                        end if
                        call h5dopen_f(file_id, trim(dset), dset_id, hdferr)! Open Dataset
                        locShape = shape(globRDF%randField_Local)
                        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, globRDF%randField_Local, locShape,  hdferr)
                        call h5dclose_f(dset_id, hdferr) !Close Dataset
                        call h5fclose_f(file_id, hdferr) ! Close the file.
                        call h5close_f(hdferr) ! Close FORTRAN interface.

                        if(delete_intermediate_files) call system("rm "//trim(randFieldFilePath))
                        if(delete_intermediate_files) call system("rm "//trim(XMFFilePath))
                        !call DispCarvalhol(globRDF%randField_Local, "globRDF%randField_Local")

                        if(globRDF%nb_procs > 1) then
                            call wLog("")
                            call wLog("GENERATING OVERLAP")
                            call wLog("-------------------------------")
                            if(globRDF%rang == 0) write(*,*)"GENERATING OVERLAP"
                            if(globRDF%rang == 0) write(*,*) "-------------------------------"
                            call wLog("")

                            if(globRDF%rang == 0) write(*,*) "    ->Applying Weighting Functions"
                            call wLog("    ->Applying Weighting Functions on Field")
                            call applyWeightingFunctions_OnMatrix(globRDF, globMSH, partitionType)
                            if(globRDF%rang == 0) write(*,*) "    ->addNeighboursFields"
                            call wLog("    ->addNeighboursFields")
                            !call addNeighboursFields(globRDF, globMSH)
                            call addNeighboursFieldsV2(globRDF, globMSH)
                        end if



                        if(size(offsetCoords,2) == 1) then
                            last_round = .true.
                            t_loc = MPI_Wtime(); !Localization Time
                            !call MPI_ALLREDUCE (t_loc, all_t_loc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, globMSH%comm, code)
                            all_t_prep = (t_prep - t_ref) * dble(IPT%nb_procs)
                            all_t_gen  = (t_gen - t_prep) * dble(IPT%nb_procs)
                            all_t_loc  = (t_loc - t_gen) * dble(IPT%nb_procs)
                            !Things to do only on the final field
                            if(globRDF%rang == 0) write(*,*)"ADJUSTEMENTS ON THE FINAL FIELD"
                            call wLog("-> Adjustement on the Final Field")
                            minPosProc = find_xNStep(globMSH%xMinGlob, globMSH%xMinInt , globMSH%xStep) - globMSH%origin + 1
                            maxPosProc = find_xNStep(globMSH%xMinGlob, globMSH%xMaxBound , globMSH%xStep) - globMSH%origin + 1
                            call wLog("     Normalizing sample")
                            xNTotal = product(nint((globMSH%xMaxGlob-globMSH%xMinGlob)/globMSH%xStep,8) +1)
                            call normalize_randField(globRDF, xNTotal, globRDF%randField, minPosProc, maxPosProc)
                            call multiVariateTransformation (globRDF%margiFirst, globRDF%fieldAvg, globRDF%fieldVar, globRDF%randField)
                            t_trans = MPI_Wtime(); !Transformation Time
                            all_t_trans = (t_trans - t_loc) * dble(IPT%nb_procs)
                            globRDF%prep_CPU_Time = all_t_prep
                            globRDF%gen_CPU_Time = all_t_gen
                            globRDF%loc_CPU_Time = all_t_loc
                            globRDF%trans_CPU_Time = all_t_trans
                            if(globRDF%rang == 0) call write_generation_spec(globMSH, globRDF, single_path, "singleGen", &
                                                       IPT%nFields, IPT%localizationLevel, &
                                                       gen_times, nGenGroups, IPT%nb_procs)
                        end if

                        call wLog("-> Writing XMF and hdf5 files FOR GLOBAL FIELD")
                        if(writeFiles) then
                            call wLog("outputStyle");
                            call wLog(outputStyle);
                            if(globRDF%rang == 0) write(*,*) "-> Writing 'Bounding Box' XMF and hdf5 files"
                            BBoxPartFileName = string_join_many(BBoxFileName,"_L",numb2String(locLevel,2))
                            do i = 1, size(procCoord)
                                 BBoxPartFileName = string_join_many(BBoxPartFileName,"_", &
                                                   numb2String( &
                                                   nint(offsetCoords(i,locIter)/(stepProc(i)*nFields_level(i)))+1,3))
                            end do

                            call wLog("  fileName");
                            call wLog(BBoxPartFileName);
                            randFieldFilePath = string_join_many(single_path,"/h5/",BBoxPartFileName,".h5")

                            writeDataSet = .true.
                            if(last_round) then
                                finalFieldPath = randFieldFilePath
                                if(IPT%rang == 0) then
                                    call write_stat_input("./stat_input", randFieldFilePath)
                                    writeDataSet = IPT%writeDataSet
                                end if
                            end if

                            if(outputStyle==1) then
                                call wLog("   (Parallel)");
                                call wLog("minval(RDF%randField,1) =")
                                call wLog(minval(globRDF%randField,1))
                                call wLog("maxval(RDF%randField,1) =")
                                call wLog(maxval(globRDF%randField,1))
                                call write_Mono_XMF_h5(globRDF, globMSH, BBoxPartFileName, globRDF%rang, single_path, &
                                                       globMSH%comm, ["_ALL"], [0], fieldNumber, style=outputStyle, &
                                                       writeDataSet = IPT%writeDataSet, localization=.true.)
                            else
                                call wLog("   (Per Proc)");
                                call write_Mono_XMF_h5(globRDF, globMSH, BBoxPartFileName, globRDF%rang, single_path, &
                                                       globMSH%comm, ["_procOnlyShape"], [IPT%rang], 0, style=outputStyle, &
                                                       writeDataSet = IPT%writeDataSet, localization=.true.)
                            end if

                            if(.not. writeDataSet) then
                                call system("rm -r "//string_join_many(single_path,"/h5"))
                                call system("rm -r "//string_join_many(single_path,"/xmf"))
                            end if

                        end if


                        call finalize_MESH(globMSH)
                        call finalize_RF(globRDF)

                    end do

                    if(IPT%rang == 0) write(*,*)  "  Ending Localization Level ", locLevel
                    call wLog("     Syncing Procs")
                    call MPI_BARRIER(validProcComm, code) !We need all the procs to finish before going to next localizationLevel

                end do

                call finalize_IPT_RF(newIPT)

            end if

            if(allocated(offsetCoords)) deallocate(offsetCoords)

        end subroutine combine_subdivisions

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine interpolateToMesh(BBoxFileName, coordList, UNV_randField, rang)
            implicit none

            !INPUT
            character(len=*), intent(in)      :: BBoxFileName
            double precision, dimension(:,:), intent(in) :: coordList
            integer, intent(in) :: rang
            !OUTPUT
            double precision, dimension(:,:), intent(out) :: UNV_randField
            !LOCAL
            integer :: nDim, Nmc
            !logical :: independent
            character(len=50) :: attr_Name, dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, space_id, dset_id, mem_id
            integer(HSIZE_T), dimension(size(coordList,1)) :: offset, locDims
            integer, dimension(size(coordList,1)) :: xNStep, coordPosInt
            integer, dimension(size(coordList,1), 2**size(coordList,1)) :: neighCoord
            double precision, dimension(size(coordList,1)) :: coordPos
            double precision, dimension(size(coordList,1)) :: distance
            double precision, dimension(size(coordList,1)) :: xMinGlob, xMaxGlob, xStep
            double precision, dimension(:,:), allocatable, target :: BB_randField
            double precision, dimension(size(coordList,1)) :: xMin_Loc_UNV, xMax_Loc_UNV
            integer, dimension(size(coordList,1)) :: minPos, maxPos, extent
            integer(HSIZE_T), dimension(2) :: locShape, zero2D
            integer :: i, j
            double precision, dimension(:,:)    , pointer :: BB_2D
            double precision, dimension(:,:,:)  , pointer :: BB_3D
            double precision :: weight


            call h5open_f(hdferr) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(BBoxFileName), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
            if(hdferr /= 0) stop("ERROR OPENING FILE")
            !write(*,*) "hdferr = ", hdferr

            !READING SCALARS----------------------------
            !BOOL
            !attr_name = "independent"
            !call read_h5attr_bool(file_id, trim(adjustL(attr_name)), independent)

            !INTEGERS
            attr_name = "nDim"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), nDim)
            attr_name = "Nmc"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), Nmc)

            !DOUBLE VEC
            attr_name = "xMinGlob"
            call read_h5attr_real_vec(file_id, attr_name, xMinGlob)
            attr_name = "xMaxGlob"
            call read_h5attr_real_vec(file_id, attr_name, xMaxGlob)
            attr_name = "xStep"
            call read_h5attr_real_vec(file_id, attr_name, xStep)

            xNStep = nint((xMaxGlob-xMinGlob)/xStep) +1

            !DEFINE LOCAL BOUNDING BOX
            do i = 1, nDim
                xMin_Loc_UNV(i) = minval(coordList(i,:))
                xMax_Loc_UNV(i) = maxval(coordList(i,:))
            end do

            minPos = floor((xMin_Loc_UNV-xMinGlob)/xStep) + 1
            maxPos = ceiling((xMax_Loc_UNV-xMinGlob)/xStep) + 1
            where(minPos < 1) minPos = 1
            where(maxPos > xNStep) maxPos = xNStep

            extent = maxPos - minPos + 1

            call wLog("xMin_Loc_UNV")
            call wLog(xMin_Loc_UNV)
            call wLog("xMax_Loc_UNV")
            call wLog(xMax_Loc_UNV)
            call wLog("minPos BB")
            call wLog(minPos)
            call wLog("maxPos BB")
            call wLog(maxPos)
            call wLog("extent BB")
            call wLog(extent)

            allocate(BB_randField(product(extent),1))

            !READING MATRIX BLOCK----------------------------------------------------------------
            locDims = extent
            call h5dopen_f(file_id, trim(dset), dset_id, hdferr)! Open Dataset
            call h5dget_space_f(dset_id, space_id, hdferr) !Open Dataset Space
            !call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr) !Measure Dataset Space

            !allocate(STA%randField(product(dims),1))
            offset = minPos-1
            locShape = shape(BB_randField)
            zero2D = 0
            call wLog(" locShape = ")
            call wLog(int(locShape))
            call wLog(" offset   = ")
            call wLog(int(offset))
            call wLog(" locDims  = ")
            call wLog(int(locDims))
            !For hyperslab lecture

            !IN
            call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, locDims, hdferr) !Select Hyperslab IN

            !OUT
            call h5screate_simple_f(2, locShape, mem_id, hdferr) !Create memory dataspace
            call h5sselect_hyperslab_f(mem_id, H5S_SELECT_SET_F, zero2D, locShape, hdferr) !Select Hyperslab OUT
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, BB_randField, locShape, hdferr, mem_id, space_id) !Read Dataset Hyperslab


            call h5dclose_f(dset_id, hdferr) !Close Dataset
            call h5sclose_f(space_id, hdferr) !Close Dataset Space

            call h5fclose_f(file_id, hdferr) ! Close the file.
            call h5close_f(hdferr) ! Close FORTRAN interface.


            !INTERPOLATION
            if(nDim == 2) then
                BB_2D(minPos(1):maxPos(1), minPos(2):maxPos(2)) => BB_randField
                neighCoord(:,1) = [0, 0]
                neighCoord(:,2) = [1, 0]
                neighCoord(:,3) = [0, 1]
                neighCoord(:,4) = [1, 1]

            else if(nDim == 3) then
                BB_3D(minPos(1):maxPos(1), minPos(2):maxPos(2), minPos(3):maxPos(3)) => BB_randField
                neighCoord(:,1) = [0, 0, 0]
                neighCoord(:,2) = [1, 0, 0]
                neighCoord(:,3) = [0, 1, 0]
                neighCoord(:,4) = [0, 0, 1]
                neighCoord(:,5) = [1, 1, 0]
                neighCoord(:,6) = [0, 1, 1]
                neighCoord(:,7) = [1, 0, 1]
                neighCoord(:,8) = [1, 1, 1]

            end if

            UNV_randField(:,:) = 0

            do i = 1, size(coordList, 2)

                !MAPPING COORD In BB
                coordPos = ((coordList(:,i)-xMinGlob)/xStep) + 1.0D0
                coordPosInt = floor(coordPos)
                where(coordPosInt == maxPos) coordPosInt = coordPosInt - 1 !Dealing with points on the positive border
                if(any(coordPosInt<0)) stop("coordPosInt smaller than 1")

                !Applying Values
                do j = 1, size(neighCoord, 2)

                    distance(:) = abs(coordPos - dble(coordPosInt+neighCoord(:,j)))
                    weight      = product(1.0D0 - distance)
                    !weight      = 1.0D0 - sqrt(sum(distance**2))/sqrt(dble(nDim))

                    !write(*,*) "weight = ", weight

                    if(any(coordPosInt(:)+neighCoord(:,j) > maxPos)) then
                        call wLog("Error in rang ")
                        call wLog(rang)
                        call wLog("   coordPos = ")
                        call wLog(coordPos)
                        call wLog("          j = ")
                        call wLog(j)
                        call wLog("coordPosInt(:)+neighCoord(:,j) = ")
                        call wLog(coordPosInt(:)+neighCoord(:,j))
                        call wLog("maxPos = ")
                        call wLog(maxPos)
                        !stop(" ERROR! UNV TRIED POSITION OUT OF RANGE")

                    end if

                    if(any(coordPosInt(:)+neighCoord(:,j) < minPos)) then
                        call wLog("Error in rang ")
                        call wLog(rang)
                        call wLog("   coordPos = ")
                        call wLog(coordPos)
                        call wLog("          j = ")
                        call wLog(j)
                        call wLog("coordPosInt(:)+neighCoord(:,j) = ")
                        call wLog(coordPosInt(:)+neighCoord(:,j))
                        call wLog("minPos = ")
                        call wLog(minPos)
                        !stop(" ERROR! UNV TRIED POSITION OUT OF RANGE")
                    end if

                    if(nDim == 2) then
                        UNV_randField(i,1) = UNV_randField(i,1) +                  &
                                             (                                     &
                                             BB_2D(coordPosInt(1)+neighCoord(1,j), &
                                                   coordPosInt(2)+neighCoord(2,j)) &
                                             * weight                              &
                                             )
                    else if (nDim == 3) then
                        UNV_randField(i,1) = UNV_randField(i,1) +                  &
                                             (                                     &
                                             BB_3D(coordPosInt(1)+neighCoord(1,j), &
                                                   coordPosInt(2)+neighCoord(2,j), &
                                                   coordPosInt(3)+neighCoord(3,j)) &
                                             * weight                              &
                                             )
                    end if

                end do


!                if(nDim == 2) then
!                    UNV_randField(i,1) = BB_2D(coordPosInt(1), coordPosInt(2))
!                else if (nDim == 3) then
!                    UNV_randField(i,1) = BB_3D(coordPosInt(1), coordPosInt(2), coordPosInt(3))
!                end if

!                coordPosInt = nint(coordPos)
!
!                if(nDim == 2) then
!                    UNV_randField(i,1) = BB_2D(coordPosInt(1), coordPosInt(2))
!                else if (nDim == 3) then
!                    UNV_randField(i,1) = BB_3D(coordPosInt(1), coordPosInt(2), coordPosInt(3))
!                end if

            end do




            if (allocated(BB_randField))   deallocate(BB_randField)

        end subroutine interpolateToMesh


end module sample_RF
