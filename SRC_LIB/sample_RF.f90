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
    use calls_RF

    implicit none

contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------

        subroutine single_realization(IPT, globMSH, writeFiles, outputStyle, sameFolder, nProcPerField, fieldComm, fieldNumber, h5fullPath)

            implicit none
            !INPUT
            type(IPT_RF), intent(in) :: IPT
            type(MESH)  , intent(in) :: globMSH
            logical, intent(in) :: writeFiles
            logical, intent(in) :: sameFolder
            integer, intent(in) :: outputStyle!1: parallel hdf5, 2: hdf5 per proc
            integer              , intent(in) :: nProcPerField
            integer, intent(in) :: fieldComm, fieldNumber

            !LOCAL
            type(RF)      :: RDF
            type(MESH)    :: MSH
            double precision, dimension(IPT%nDim_gen) :: procExtent, procStart
            double precision, dimension(:), allocatable :: seedStartVec
            double precision, dimension(:,:), allocatable :: UNV_randField
            double precision :: tLoc1, tLoc2
            double precision :: t1, t2, t3, t4, t5, t6;
            double precision :: all_t1, all_t2, all_t3, all_t4, all_t5, all_t6;
            integer :: code
            integer :: all_xNTotal
            integer :: i
            logical :: validProc
            integer :: newComm, newNbProcs, newRang
            character(len=100) :: BBoxPartFileName, BBoxPath
            character(len=*)   :: h5fullPath
            integer(kind=8) :: xNTotal


            validProc = .true.
            call MPI_COMM_RANK(fieldComm, newRang, code)
            call MPI_COMM_SIZE(fieldComm, newNbProcs, code)

            !call wLog("-> validateProc")
            !call validateProc(validProc, IPT, newComm, newNbProcs, newRang)

            !if(validProc) then

                call init_MESH(MSH, IPT, fieldComm, newNbProcs, newRang)
                call init_RF(RDF, IPT, fieldComm, newNbProcs, newRang)

                !Copy from globMSH
                MSH%xStep   = globMSH%xStep
                MSH%overlap = globMSH%overlap
                MSH%xMinGlob = 0.0D0
                MSH%xMaxGlob = globMSH%procExtent
                MSH%procExtent = globMSH%procExtent
                MSH%independent = .false.
                RDF%independent = .false.
                MSH%procPerDim = 1
                MSH%procPerDim(MSH%nDim) = newNbProcs
                MSH%coords = 0
                MSH%coords(MSH%nDim) = newRang

                !call wLog("-> set_communications_topology")
                !call set_communications_topology(MSH)
                !call wLog("-> round_basic_inputs")
                !call round_basic_inputs(MSH, MSH%xStep, MSH%overlap)
                !call wLog("-> set_global_extremes")
                !call set_global_extremes(MSH, MSH%xMaxGlob, MSH%xMinGlob, procExtent, procStart)
                !call wLog("     procStart = ")
                !call wLog(procStart)
                !call wLog("     procExtent = ")
                !call wLog(procExtent)
                !call wLog("-> set_local_bounding_box")

                call set_local_bounding_box(MSH,&
                                            MSH%xMinBound, MSH%xMaxBound, &
                                            MSH%xNStep, MSH%xNTotal, MSH%origin, validProc)

                call set_validProcs_comm(validProc, fieldComm, MSH%rang, &
                                         MSH%validProc, RDF%validProc, MSH%comm, RDF%comm, &
                                         MSH%nb_procs, RDF%nb_procs, MSH%rang, RDF%rang)

                if(validProc) then
                    !call MPI_BARRIER(MSH%comm, code)

    !                    call wLog("-> set_overlap_geometry")
    !                    call set_overlap_geometry (MSH, MSH%xMinInt, MSH%xMaxInt, MSH%xMinExt, MSH%xMaxExt, &
    !                                               MSH%xMaxNeigh, MSH%xMinNeigh, MSH%xOrNeigh, MSH%nOvlpPoints)


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

                    call wLog("-> Setting folder path")
                    single_path = string_join_many(results_path,"/",results_folder_name)
                    call wLog("     single_path = "//trim(single_path))

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

                    !call wLog("     Normalizing sample")
                    !xNTotal = product(nint((MSH%xMaxGlob-MSH%xMinGlob)/MSH%xStep) +1)
                    !call normalize_randField(RDF, xNTotal, RDF%randField)

                    !i = size(RDF%xPoints,2)
                    !if(i>50) i = 50
                    !call dispCarvalhol(transpose(RDF%xPoints(:,1:i)), "transpose(RDF%xPoints)", "(F20.5)",unit_in = RDF%log_ID)


                    !i = size(RDF%xPoints,2)
                    !if(i>50) i = 50
                    !call dispCarvalhol(RDF%xPoints(1:i,:), "xPoints", "(F20.5)",unit_in = RDF%log_ID)

                    !call show_RF(RDF, forLog_in = .true.)

                    !call multiVariateTransformation (RDF%margiFirst, RDF%fieldAvg, RDF%fieldVar, RDF%randField)

                    t2 = MPI_Wtime();
                    call MPI_ALLREDUCE (t2, all_t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MSH%comm,code)
                    RDF%gen_CPU_Time = all_t2 - all_t1
                    if(RDF%rang == 0) write(*,*) "Generation CPU Time = ", all_t2 - all_t1
                    call wLog ("    Generation CPU Time (s)")
                    call wLog (RDF%gen_CPU_Time)
                    all_t3 = -1.0D0

                    call wLog("-> Writing XMF and hdf5 files");

                    if(writeFiles) then
                        call wLog("outputStyle");
                        call wLog(outputStyle);
                        if(RDF%rang == 0) write(*,*) "-> Writing 'Bounding Box' XMF and hdf5 files"
                        BBoxPartFileName = string_join_many(BBoxFileName,"_part",numb2String(fieldNumber,5))

                        if(outputStyle==1) then
                            call wLog("   (Parallel)");
                            call wLog("minval(RDF%randField,1) =")
                            call wLog(minval(RDF%randField,1))
                            call wLog("maxval(RDF%randField,1) =")
                            call wLog(maxval(RDF%randField,1))
                            call write_Mono_XMF_h5(RDF, MSH, IPT%connectList, IPT%monotype, BBoxPartFileName, RDF%rang, single_path, &
                                                   MSH%comm, ["_Part"], [fieldNumber], fieldNumber, style=outputStyle, meshMod = msh_AUTO, &
                                                   HDF5FullPath = BBoxPath, writeDataSet = IPT%writeDataSet)
                        else
                            call wLog("   (Per Proc)");
                            call write_Mono_XMF_h5(RDF, MSH, IPT%connectList, IPT%monotype, BBoxPartFileName, RDF%rang, single_path, &
                                                   MSH%comm, ["_Part"], [RDF%rang], 0, style=outputStyle, meshMod = msh_AUTO, &
                                                   HDF5FullPath = BBoxPath, writeDataSet = IPT%writeDataSet)

                        end if

                        t3 = MPI_Wtime();
                        call MPI_ALLREDUCE (t3, all_t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM,MSH%comm,code)
                        if(RDF%rang == 0) write(*,*) "Writing Files Time = ", all_t3 - all_t2
                        call wLog ("    Writing Files CPU Time (s)")
                        call wLog (all_t3 - all_t2)

                        h5fullPath = BBoxPath
                    end if
!
!                    call write_generation_spec(MSH, RDF, single_path, "singleGen", &
!                                               [all_t1,all_t2,all_t3])
                end if
!            !end if
!
!            call MPI_BARRIER(IPT%comm, code) !Waiting for the generation to finish
!
!
!
!            if(IPT%unv .and. writeFiles .and. outputStyle == 1) then
!                if(IPT%rang == 0) write(*,*) "-> Writing 'UNV' XMF and hdf5 files for"
!                if(IPT%rang == 0) write(*,*) IPT%unv_path
!                allocate(UNV_randField(size(IPT%coordList,2),1))
!                if(RDF%rang == 0) write(*,*) "  Source:"
!                if(RDF%rang == 0) write(*,*) BBoxPath
!                call interpolateToUNV(BBoxPath, IPT%coordList, UNV_randField, IPT%rang)
!                call write_UNV_XMF_h5(UNV_randField, IPT%coordList, IPT%connectList, IPT%monotype, &
!                                      "UNV_", RDF%rang, single_path, &
!                                      MSH%comm, ["_All_UNV"], [RDF%rang], 0)
!            end if
!
!            if(IPT%rang == 0) call write_stat_input("./stat_input", BBoxPath)

            call finalize_MESH(MSH)
            call finalize_RF(RDF)
            if(allocated(UNV_randField)) deallocate(UNV_randField)

        end subroutine single_realization

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine interpolateToUNV(BBoxFileName, coordList, UNV_randField, rang)
            implicit none

            !INPUT
            character(len=*), intent(in)      :: BBoxFileName
            double precision, dimension(:,:), intent(in) :: coordList
            integer, intent(in) :: rang
            !OUTPUT
            double precision, dimension(:,:), intent(out) :: UNV_randField
            !LOCAL
            integer :: nDim, Nmc, method, corrMod, margiFirst
            logical :: independent
            character(len=50) :: attr_Name, dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, attr_id, space_id, dset_id, mem_id
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
            integer :: i, j, d
            double precision, dimension(:,:)    , pointer :: BB_2D
            double precision, dimension(:,:,:)  , pointer :: BB_3D
            double precision :: weight


            call h5open_f(hdferr) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(BBoxFileName), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
            if(hdferr /= 0) stop("ERROR OPENING FILE")
            !write(*,*) "hdferr = ", hdferr

            !READING SCALARS----------------------------
            !BOOL
            attr_name = "independent"
            call read_h5attr_bool(file_id, trim(adjustL(attr_name)), independent)

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

        end subroutine interpolateToUNV

end module sample_RF
