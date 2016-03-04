module calls_RF

    use displayCarvalhol
    use write_Log_File
    use math_RF
    use constants_RF
    use mpi
    use writeResultFile_RF
    use type_RF
    use localization_RF
    use type_MESH
    use common_variables_RF
    use randomFieldND
    use mesh_RF
    use type_inputRF
    use sample_RF

    implicit none

    interface createRandomField
       module procedure create_RF_Unstruct_Init
    end interface createRandomField

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine make_random_field (IPT, times, nTotalFields)
        !INPUT
        type(IPT_RF)  :: IPT
        double precision, dimension(:), intent(inout) :: times
        integer, intent(in)               :: nTotalFields
        !LOCAL
        type(MESH)            :: globMSH
        double precision, dimension(:,:), allocatable :: UNV_randField
        integer               :: gen_group, gen_Comm, gen_groupMax
        integer               :: loc_group, loc_Comm, loc_groupMax
        integer               :: code
        double precision, dimension(IPT%nDim_mesh) :: stepProc, procExtent, overlap, gen_GroupRange
        !double precision, dimension(IPT%nDim_mesh, product(IPT%nFields)) :: offsetId
        double precision, dimension(IPT%nDim_mesh, nTotalFields) :: subdivisionCoords
        integer         , dimension(IPT%nDim_mesh, nTotalFields) :: subdivisionId
        double precision      :: t_bef, t_aft
        integer               :: fieldNumber
        character(len=200) :: BBoxPath, MONO_FileName
        character(len=200), dimension(:), allocatable :: MONO_FileNames
        character(len=110), dimension(nTotalFields) :: HDF5Name
        double precision, dimension(nTotalFields) :: gen_times, temp_gen_times
        integer :: i, d, countFields, j
        integer :: nSamplesInProc, nSamplesInAllProc, rest, sum_SamplesInProc
        integer :: loc_nbProcs, gen_nbProcs, gen_rang, loc_rang
        double precision, dimension(:,:,:), allocatable, target :: randField_inProc
        double precision, dimension(:,:), allocatable, target :: randField_Group
        double precision, dimension(:,:), allocatable :: randField_Local
        double precision, dimension(:), allocatable ::unityPartition
        integer(kind=8) :: xNTotal_Proc, xNTotal_Group
        double precision, dimension(:,:), allocatable :: xMinFiles, xMaxFiles
        double precision, dimension(IPT%nDim_gen) :: ones, xMin_Group
        integer, dimension(IPT%nDim_gen) :: xNStep_Proc, xNStep_Group
        integer, dimension(IPT%nDim_gen) :: locStep, minP, maxP
        double precision, dimension(:, :), pointer :: RF_2D_Proc, RF_2D_Group
        double precision, dimension(:, :, :), pointer :: RF_3D_Proc, RF_3D_Group


        !Building Subdivisions
        if(IPT%rang == 0) write(*,*) " "
        if(IPT%rang == 0) write(*,*) "-> DIVIDING----------------------------------------"
        call wLog("-> DIVIDING----------------------------------------")
        !allocate(subdivisionCoords(IPT%nDim_mesh, product(IPT%nFields**IPT%localizationLevel)))
        !call build_subdivisions(IPT, globMSH, gen_groupMax, &
        !                        gen_group, gen_Comm, stepProc, procExtent, overlap)

        call build_subdivisions(IPT, globMSH, &
                                stepProc, procExtent, overlap, &
                                gen_groupMax, gen_group, gen_Comm, gen_nbProcs, &
                                gen_rang, &
                                loc_groupMax, loc_group, loc_Comm, loc_nbProcs, &
                                loc_rang)
        call set_communications_topology(globMSH, globMSH%coords, globMSH%neigh, &
                                         globMSH%neighShift, globMSH%considerNeighbour, &
                                         globMSH%mappingFromShift, globMSH%op_neigh, gen_group)

        !Init
        ones = 1.0D0
        locStep = IPT%nFields**(IPT%localizationLevel)
        !call setGrid(subdivisionId, 0.0D0*ones, ones, locStep, inverse=.true.)
        call setGrid(subdivisionCoords, globMSH%xMinGlob, stepProc, locStep, inverse=.true.)
        do i = 1, size(subdivisionCoords, 2)
            subdivisionId(:,i) = nint((subdivisionCoords(:,i)-globMSH%xMinGlob)/stepProc)
        end do
        !if(rang == 0) call DispCarvalhol(subdivisionCoords, "subdivisionCoords")
        !if(IPT%rang == 0) call DispCarvalhol(subdivisionId, "subdivisionId")
        if(IPT%rang == 0) write(*,*) "Max Coord = ", subdivisionCoords(:, size(subdivisionCoords,2)) + stepProc
        call wLog("globMSH%coords = ")
        call wLog(globMSH%coords)
        call show_MESH(globMSH, "Dopped Neighbours", forLog_in=.true.)

        !Discovering number of fields in each proc

        nSamplesInProc = 0
        xNTotal_Proc   = 0

        if(gen_rang == 0) then
            nSamplesInAllProc = int(nTotalFields/(loc_nbProcs))
            rest              = nTotalFields - (nSamplesInAllProc*loc_nbProcs)
            nSamplesInProc    = nSamplesInAllProc
            if(IPT%rang < rest) nSamplesInProc = nSamplesInProc + 1


            !if(IPT%rang == 0) write(*,*) " groupMax     = ", groupMax

            if(IPT%rang == 0) write(*,*) " procExtent = ", procExtent
            !if(IPT%rang == 0) write(*,*) " find_xNStep(xMaxExt=procExtent, xStep=globMSH%xStep) = ", &
            !find_xNStep(xMaxExt=procExtent, xStep=globMSH%xStep)

            xNStep_Proc = find_xNStep(xMaxExt=procExtent, xStep=globMSH%xStep)
            xNTotal_Proc = product(int(xNStep_Proc, 8))

            allocate(randField_inProc(xNTotal_Proc, IPT%Nmc, nSamplesInProc))

            if(.not. IPT%delete_intermediate_files) allocate(MONO_FileNames(nSamplesInProc))
            allocate(xMinFiles(IPT%nDim_gen , nSamplesInProc))
            allocate(xMaxFiles(IPT%nDim_gen, nSamplesInProc))

        end if

        if(IPT%rang == 0) write(*,*) " nTotalFields      = ", nTotalFields
        call MPI_ALLREDUCE (nSamplesInProc, sum_SamplesInProc, 1, MPI_INTEGER, MPI_SUM, IPT%comm, code)
        if(IPT%rang == 0) write(*,*) " sum_SamplesInProc = ", sum_SamplesInProc

        call wLog("nSamplesInProc = ")
        call wLog(nSamplesInProc)
        call wLog("xNTotal_Proc = ")
        call wLog(xNTotal_Proc)


        call MPI_BARRIER(IPT%comm, code)
        times(3) = MPI_Wtime() !Organizing Collective Writing

        !Making all realizations
        gen_times(:) = 0.0D0
        countFields  = 0
        !if(.false)then
        if(IPT%sampleFields)then
            if(IPT%rang == 0) write(*,*) " "
            if(IPT%rang == 0) write(*,*) "-> SAMPLING----------------------------------------"
            call wLog("-> SAMPLING----------------------------------------")
            do i = 1, product(IPT%nFields**IPT%localizationLevel)
            !do i = 1, 1 !FOR TESTS
                !if(mod(i, gen_groupMax) == gen_group) then
                if(all(subdivisionId(:,i)/(IPT%nFields) == globMSH%coords)) then
                    if(gen_rang == 0) write(*,*)  "-> Group ", gen_group, " making Field ", i
                    call wLog("-> Making Field")
                    call wLog(i)
                    fieldNumber = i;
                    countFields = countFields + 1
                    !call wLog("Proc")
                    !call wLog(rang)
                    !call wLog("dealing with field")
                    !call wLog(fieldNumber)
                    !call wLog("     Trying communication")
                    t_bef = MPI_Wtime()
                    !write(*,*) "Before communicator"
                    call MPI_BARRIER(gen_Comm, code)
                    allocate(randField_Local(IPT%rang+3,1))
                    randField_Local = IPT%rang
                    !write(*,*) "After communicator"
                    call single_realization(IPT, globMSH, &
                                            gen_Comm, fieldNumber, subdivisionCoords(:,i), &
                                            stepProc, randField_Local, HDF5Name(i))
                    call wLog("Gathering Sample")
                    call gather_sample(randField_Local, randField_inProc(:,:,countFields), &
                                       gen_rang, gen_nbProcs, gen_comm)

                    if(.not. IPT%delete_intermediate_files .and. gen_rang == 0) then
                        call wLog("Writing intermediate generation file")
                        xMinFiles(:, countFields) = subdivisionCoords(:,i)
                        xMaxFiles(:, countFields) = xMinFiles(:, countFields) + globMSH%procExtent
                        MONO_FileNames(countFields) = "MONO"
                        do d = 1, IPT%nDim_gen
                            MONO_FileNames(countFields) = &
                            string_join_many(MONO_FileNames(countFields),"_",numb2String(subdivisionId(d,i)+1,3))
                        end do
                        !MONO_FileName = string_join_many(MONO_FileName,"_group_",numb2String(gen_group,3))

                        call write_MONO_proc_result(xMinFiles(:, countFields), xMaxFiles(:, countFields), &
                                                    globMSH%xStep, globMSH%nDim, &
                                                    randField_inProc(:,1,countFields), &
                                                    MONO_FileNames(countFields), single_path)
                    end if
                    if(allocated(randField_Local)) deallocate(randField_Local)
                    t_aft = MPI_Wtime()
                    temp_gen_times(i) = t_aft-t_bef
                    !write(*,*) "After single"
                end if
            end do
        end if

        !Localization Inside Group
        gen_GroupRange = (dble(IPT%nFields)**(dble(IPT%localizationLevel-1)))*stepProc &
                         + globMSH%corrL*globMSH%overlap - 2.0D0*globMSH%xStep
        xNStep_Group  = find_xNStep(xMaxExt=gen_GroupRange, xStep=globMSH%xStep)
        xNTotal_Group = product(int(xNStep_Group,8))
        xMin_Group    = globMSH%xMinGlob + &
                        ((dble(IPT%nFields)**(dble(IPT%localizationLevel-1)))*stepProc)*globMSH%coords
        call wLog("gen_GroupRange = ")
        call wLog(gen_GroupRange)
        call wLog("xNStep_Group = ")
        call wLog(xNStep_Group)

        allocate(randField_Group(xNTotal_Group, IPT%Nmc))

        if(IPT%nDim_gen == 2) then
            RF_2D_Group(1:xNStep_Group(1),1:xNStep_Group(2)) => randField_Group
        else if(IPT%nDim_gen == 3) then
            RF_3D_Group(1:xNStep_Group(1),1:xNStep_Group(2),1:xNStep_Group(3)) => randField_Group
        end if

        if(gen_rang == 0) then
            allocate(unityPartition(xNTotal_Proc))

            call wLog("Generating Partition of Unity")
            call generateUnityPartition_Matrix(xNStep_Proc, globMSH%overlap, globMSH%corrL, globMSH%xStep,&
                                         1, unityPartition, globMSH%nDim)
            call write_MONO_proc_result(globMSH%procExtent*0, globMSH%procExtent, &
                                        globMSH%xStep, globMSH%nDim, &
                                        unityPartition, "UnityPartition", single_path)

            do i = 1, nSamplesInProc
                !Multiplication
                do j = 1, IPT%Nmc
                    randField_inProc(:,j,i) = randField_inProc(:,j,i)*unityPartition
                    if(.not. IPT%delete_intermediate_files) then
                        MONO_FileName = MONO_FileNames(i)
                        MONO_FileName = string_join_many("Loc_",MONO_FileName)

                        if(gen_rang == 0) call write_MONO_proc_result(xMinFiles(:, i), xMaxFiles(:, i), &
                                                                      globMSH%xStep, globMSH%nDim, &
                                                                      randField_inProc(:,j,i), MONO_FileName, &
                                                                      single_path)

                    end if
                end do

                !Sum
                minP = find_xNStep(xMin_Group, xMinFiles(:, i), globMSH%xStep)
                maxP = minP + find_xNStep(xMaxExt=procExtent, xStep=globMSH%xStep)
                if(IPT%nDim_gen == 2) then
                    RF_2D_Proc(1:xNStep_Proc(1),1:xNStep_Proc(2)) => randField_inProc(:,1,i)

                    RF_2D_Group(minP(1):maxP(1),minP(2):maxP(2)) = RF_2D_Proc &
                                                       + RF_2D_Group(minP(1):maxP(1),minP(2):maxP(2))
                else if(IPT%nDim_gen == 3) then
                    RF_3D_Proc(1:xNStep_Proc(1),1:xNStep_Proc(2),1:xNStep_Proc(3)) => randField_inProc(:,1,i)

                    RF_3D_Group(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)) = RF_3D_Proc &
                                              + RF_3D_Group(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
                end if

            end do


        end if




!
!        call MPI_BARRIER(IPT%comm, code)
!        times(4) = MPI_Wtime() !Generation Time
!
!        call MPI_ALLREDUCE (temp_gen_times, gen_times, size(gen_times), MPI_DOUBLE_PRECISION, MPI_SUM, IPT%comm,code)

!        !Combining realizations (localization)
!        if(.true.) then
!            if(IPT%rang == 0) write(*,*) " "
!            if(IPT%rang == 0) write(*,*) "-> COMBINING----------------------------------------"
!            call wLog("-> COMBINING----------------------------------------")
!            call combine_subdivisions(IPT, IPT%outputStyle, stepProc, procExtent, &
!                                      overlap, times(1), times(3), times(4), gen_times(:), &
!                                      IPT%delete_intermediate_files, IPT%ignoreTillLocLevel, &
!                                      loc_group, loc_Comm, loc_groupMax, BBoxPath)
!        end if
!
!        times(5) = MPI_Wtime() !Localization Time
!
!
!        !Writing Interpolation File
!        call MPI_BARRIER(IPT%comm, code)
!        if(IPT%unv .and. IPT%writeUNVinterpolation .and. IPT%outputStyle == 1) then
!            if(IPT%rang == 0) write(*,*) "-> Writing 'UNV' XMF and hdf5 files for"
!            if(IPT%rang == 0) write(*,*) trim(adjustL(IPT%unv_path))
!            allocate(UNV_randField(size(IPT%coordList,2),1))
!            if(IPT%rang == 0) write(*,*) "  Source:"
!            if(IPT%rang == 0) write(*,*) BBoxPath
!            if(IPT%rang == 0) write(*,*) "-> INTERPOLATING TO GIVEN MESH----------------------------------------"
!            call wLog("-> INTERPOLATING TO GIVEN MESH----------------------------------------")
!            call interpolateToMesh(BBoxPath, IPT%coordList, UNV_randField, IPT%rang)
!            call write_UNV_XMF_h5(UNV_randField, IPT%coordList, IPT%connectList, &
!                                  "UNV_", IPT%rang, single_path, &
!                                  IPT%comm, 0)
!            if(allocated(UNV_randField)) deallocate(UNV_randField)
!        end if

        call finalize_MESH(globMSH)
        if(allocated(xMinFiles))        deallocate(xMinFiles)
        if(allocated(xMaxFiles))        deallocate(xMaxFiles)
        if(allocated(unityPartition))   deallocate(unityPartition)
        if(allocated(randField_inProc)) deallocate(randField_inProc)
        if(allocated(randField_Local))  deallocate(randField_Local)
        if(allocated(randField_Group))  deallocate(randField_Group)
        if(associated(RF_2D_Proc))  nullify(RF_2D_Proc)
        if(associated(RF_2D_Group)) nullify(RF_2D_Group)
        if(associated(RF_3D_Proc))  nullify(RF_3D_Proc)
        if(associated(RF_3D_Group)) nullify(RF_3D_Group)

    end subroutine make_random_field
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine create_RF_Unstruct_noInit (xPoints, corrL, corrMod, Nmc,   &
!                                          randField, method, seedStart,   &
!                                          margiFirst, fieldAvg, fieldVar, &
!                                          comm, rang, nb_procs, calculate, MSH)
!        !INPUT
!        double precision, dimension(1:, 1:), intent(in), target :: xPoints;
!        double precision, dimension(1:)    , intent(in) :: corrL;
!        integer                            , intent(in) :: corrMod;
!        integer                            , intent(in) :: Nmc;
!        integer                            , intent(in) :: method
!        integer                            , intent(in) :: seedStart
!        integer                            , intent(in) :: margiFirst;
!        double precision                   , intent(in) :: fieldAvg
!        double precision                   , intent(in) :: fieldVar;
!        integer                            , intent(in) :: comm, rang, nb_procs
!        logical, dimension(1:), optional   , intent(in) :: calculate
!        type(MESH), intent(inout) :: MSH
!
!        !OUTPUT
!        double precision, dimension(:, :), intent(out), target :: randField;
!
!        !LOCAL
!        type(RF) :: RDF
!
!        write(*,*) "Inside create_RF_Unstruct_noInit"
!
!        !Initializing RF
!        call init_RF(RDF, size(corrL), Nmc, comm, rang, nb_procs)
!        RDF%xPoints   => xPoints
!        RDF%randField => randField
!        RDF%xNTotal    = size(RDF%xPoints, 2)
!        RDF%corrL      = corrL
!        RDF%corrMod    = corrMod
!        RDF%Nmc        = Nmc
!        RDF%method     = method
!        RDF%seedStart  = seedStart
!        RDF%margiFirst = margiFirst
!        RDF%fieldAvg   = fieldAvg
!        RDF%fieldVar   = fieldVar
!        if(present(calculate)) RDF%calculate  = calculate
!
!        call create_RF_Unstruct_Init(RDF, MSH)
!
!    end subroutine create_RF_Unstruct_noInit

end module calls_RF
