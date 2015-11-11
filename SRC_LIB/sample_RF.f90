module sample_RF

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
    use topography_RF
    use calls_RF

    implicit none

contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------

        subroutine single_realization(IPT)

            implicit none
            !INPUT
            type(IPT_RF), intent(in) :: IPT
            !LOCAL
            type(RF)      :: RDF
            type(MESH)    :: MSH
            double precision, dimension(IPT%nDim_gen) :: procExtent, procStart
            double precision, dimension(:), allocatable :: seedStartVec
            double precision :: tLoc1, tLoc2
            double precision :: t1, t2, t3, t4, t5, t6;
            double precision :: all_t1, all_t2, all_t3, all_t4, all_t5, all_t6;
            logical :: writeFiles = .true.
            logical :: sameFolder = .true.
            integer :: outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc
            integer :: code
            integer :: all_xNTotal
            integer :: i
            !integer, dimension(:,:), allocatable :: meshGrid
            !integer, dimension(4) :: stepTest = [4, 5, 6, 7]


            !allocate(meshGrid(size(stepTest), product(stepTest)))

            !call dispCarvalhol(transpose(meshGridInt(stepTest)), "transpose(meshGridInt(stepTest))")
            !deallocate(meshGrid)

            call init_MESH(MSH, IPT)
            call init_RF(RDF, IPT)

            call wLog("-> set_communications_topology")
            call set_communications_topology(MSH)
            call wLog("-> round_basic_inputs")
            call round_basic_inputs(MSH, MSH%xStep, MSH%overlap)
            call wLog("-> set_global_extremes")
            call set_global_extremes(MSH, MSH%xMaxGlob, MSH%xMinGlob, procExtent, procStart)
            call wLog("     procStart = ")
            call wLog(procStart)
            call wLog("     procExtent = ")
            call wLog(procExtent)
            call set_local_bounding_box(MSH, procStart, procExtent,&
                                        MSH%xMinBound, MSH%xMaxBound, MSH%xNStep)
            call set_overlap_geometry (MSH, MSH%xMinInt, MSH%xMaxInt, MSH%xMinExt, MSH%xMaxExt, &
                                       MSH%xMaxNeigh, MSH%xMinNeigh, MSH%xOrNeigh)


            call wLog("-> Initializing Random Seed")
            if(MSH%independent) then
                !Define independent seed in each proc
                call calculate_random_seed(RDF%seed, RDF%seedStart+RDF%rang)
                call init_random_seed(RDF%seed)
            else
                call calculate_random_seed(RDF%seed, RDF%seedStart)
                call init_random_seed(RDF%seed)
            end if
            call wLog("      RDF%seed = ")
            call wLog(RDF%seed)
            call wLog(" ")

            if(RDF%method /= FFT .or. outputStyle == 2) then
                call wLog("-> Setting xPoints")
                call set_XPoints(MSH, RDF, RDF%xPoints_Local)
                call wLog("      maxval(RDF%xPoints,2) = ")
                call wLog(maxval(RDF%xPoints,2))
                call wLog( "      minval(RDF%xPoints,2) = ")
                call wLog(minval(RDF%xPoints,2))

                !i = size(RDF%xPoints,2)
                !if(i>50) i = 50
                !call dispCarvalhol(transpose(RDF%xPoints(:,1:i)), "transpose(RDF%xPoints)", "(F20.5)",unit_in = RDF%log_ID)
            end if

            call wLog("-> Allocating random field")
            call allocate_randField(RDF, MSH%xNStep, RDF%randField_Local)
            call wLog("     shape(RDF%randField)")
            call wLog(shape(RDF%randField))

            call wLog("-> Setting folder path")
            single_path = string_vec_join([results_path,"/",results_folder_name])
            call wLog("     single_path = "//trim(single_path))

            !Discovering the total number of points in all procs
            call MPI_ALLREDUCE (RDF%xNTotal, all_xNTotal,1,MPI_INTEGER, &
                                MPI_SUM,IPT%comm,code)
            !Getting Initial Time
            t1 = MPI_Wtime();
            call MPI_ALLREDUCE (t1, all_t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,IPT%comm,code)
            !if(RDF%rang == 0) write(*,*) "Time Zero = ", all_t1
!
!            call wLog("-> Generating Random Field")
!            call create_RF_Unstruct_Init (RDF, MSH)
!
!            if(outputStyle == 1 .and. MSH%meshMod == msh_AUTO .and. RDF%independent) then
!                if(RDF%method /= FFT) then
!                    call wLog(" ")
!                    call wLog("-> Reordering Random Field")
!                    if(RDF%rang == 0) write(*,*) "-> Reordering Random Field"
!                    tLoc1 = MPI_Wtime()
!                    call reorderRandomFieldStruct(RDF, MSH)
!                    tLoc2 = MPI_Wtime()
!                    call wLog("       time (s)")
!                    call wLog(tLoc2 - tLoc1)
!                end if
!           end if
!
!            !i = size(RDF%xPoints,2)
!            !if(i>50) i = 50
!            !call dispCarvalhol(transpose(RDF%xPoints(:,1:i)), "transpose(RDF%xPoints)", "(F20.5)",unit_in = RDF%log_ID)
!
!
!            !i = size(RDF%xPoints,2)
!            !if(i>50) i = 50
!            !call dispCarvalhol(RDF%xPoints(1:i,:), "xPoints", "(F20.5)",unit_in = RDF%log_ID)
!
!            !call show_RF(RDF, forLog_in = .true.)
!
!            !call multiVariateTransformation (RDF%margiFirst, RDF%fieldAvg, RDF%fieldVar, RDF%randField)
!
!            t2 = MPI_Wtime();
!            call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,IPT%comm,code)
!            RDF%gen_CPU_Time = all_t2 - all_t1
!            if(RDF%rang == 0) write(*,*) "Generation CPU Time = ", all_t2 - all_t1
!            call wLog ("    Generation CPU Time (s)")
!            call wLog (RDF%gen_CPU_Time)
!            all_t3 = -1.0D0
!
!            call wLog("-> Writing XMF and hdf5 files");
!
!            if(writeFiles) then
!                call wLog("outputStyle");
!                call wLog(outputStyle);
!                if(RDF%rang == 0) write(*,*) "-> Writing XMF and hdf5 files"
!
!                if(outputStyle==1) then
!                    call wLog("   (Parallel)");
!                    call write_Mono_XMF_h5(RDF, MSH, IPT%connectList, IPT%monotype, "trans_", RDF%rang, single_path, &
!                                                    MPI_COMM_WORLD, ["_All"], [0], 0, style=outputStyle)
!                else
!                    call wLog("   (Per Proc)");
!                    call write_Mono_XMF_h5(RDF, MSH, IPT%connectList, IPT%monotype, "trans_", RDF%rang, single_path, &
!                                                    MPI_COMM_WORLD, ["_All"], [RDF%rang], 0, style=outputStyle)
!                end if
!                t3 = MPI_Wtime();
!                call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,IPT%comm,code)
!                if(RDF%rang == 0) write(*,*) "Writing Files Time = ", all_t3 - all_t2
!                call wLog ("    Writing Files CPU Time (s)")
!                call wLog (all_t3 - all_t2)
!            end if
!
!            call write_generation_spec(MSH, RDF, single_path, "singleGen", &
!                                       [all_t1,all_t2,all_t3])

            call finalize_MESH(MSH)
            call finalize_RF(RDF)

        end subroutine single_realization

end module sample_RF
