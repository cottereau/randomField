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
        integer               :: group, groupComm, groupMax
        integer               :: code
        double precision, dimension(IPT%nDim_mesh) :: stepProc, procExtent, overlap
        double precision, dimension(IPT%nDim_mesh, nTotalFields) :: subdivisionCoords
        double precision      :: t_bef, t_aft
        integer               :: fieldNumber
        character(len=110), dimension(nTotalFields) :: HDF5Name
        double precision, dimension(nTotalFields) :: gen_times, temp_gen_times
        integer :: i


        !Building Subdivisions
        if(IPT%rang == 0) write(*,*) " "
        if(IPT%rang == 0) write(*,*) "-> DIVIDING----------------------------------------"
        call wLog("-> DIVIDING----------------------------------------")
        !allocate(subdivisionCoords(IPT%nDim_mesh, product(IPT%nFields**IPT%localizationLevel)))
        call build_subdivisions(IPT, globMSH, groupMax, &
                                group, groupComm, stepProc, procExtent, overlap)
        call setGrid(subdivisionCoords, globMSH%xMinGlob, stepProc, IPT%nFields**IPT%localizationLevel, inverse=.true.)
        !if(rang == 0) call DispCarvalhol(subdivisionCoords, "subdivisionCoords")
        if(IPT%rang == 0) write(*,*) "Max Coord = ", subdivisionCoords(:, size(subdivisionCoords,2)) + stepProc

        call MPI_BARRIER(IPT%comm, code)
        times(3) = MPI_Wtime() !Organizing Collective Writing

        !Making all realizations
        gen_times(:) = 0.0D0
        if(IPT%sampleFields)then
            if(IPT%rang == 0) write(*,*) " "
            if(IPT%rang == 0) write(*,*) "-> SAMPLING----------------------------------------"
            call wLog("-> SAMPLING----------------------------------------")
            do i = 1, product(IPT%nFields**IPT%localizationLevel)
                if(mod(i, groupMax) == group) then
                    if(mod(IPT%rang,IPT%nProcPerField) == 0) write(*,*)  "-> Group ", group, " making Field ", i
                    call wLog("-> Making Field")
                    call wLog(i)
                    fieldNumber = i;
                    !call wLog("Proc")
                    !call wLog(rang)
                    !call wLog("dealing with field")
                    !call wLog(fieldNumber)
                    !call wLog("     Trying communication")
                    t_bef = MPI_Wtime()
                    call MPI_BARRIER(groupComm, code)
                    call single_realization(IPT, globMSH, IPT%outputStyle, &
                                            groupComm, fieldNumber, subdivisionCoords(:,i), stepProc, HDF5Name(i))
                    t_aft = MPI_Wtime()
                    temp_gen_times(i) = t_aft-t_bef

                end if
            end do
        end if

        call finalize_MESH(globMSH)

        call MPI_BARRIER(IPT%comm, code)
        times(4) = MPI_Wtime() !Generation Time

        call MPI_ALLREDUCE (temp_gen_times, gen_times, size(gen_times), MPI_DOUBLE_PRECISION, MPI_SUM, IPT%comm,code)
        !if(allocated(temp_gen_times)) deallocate(temp_gen_times)

        !Combining realizations (localization)
        if(.true.) then
            if(IPT%rang == 0) write(*,*) " "
            if(IPT%rang == 0) write(*,*) "-> COMBINING----------------------------------------"
            call wLog("-> COMBINING----------------------------------------")
            call combine_subdivisions(IPT, IPT%outputStyle, stepProc, procExtent, &
                                      overlap, times(1), times(3), times(4), gen_times(:), &
                                      groupMax, IPT%delete_intermediate_files, IPT%ignoreTillLocLevel)
        end if

        times(5) = MPI_Wtime() !Localization Time

        call finalize_MESH(globMSH)

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
