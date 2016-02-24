program main_Stat

    use mpi
    use hdf5
    use constants_RF
    use hdf5_RF
    use statistics_RF
    use displayCarvalhol
    use write_Log_File
    use readFile_RF
    use systemUt_RF
    use common_variables_RF
    use charFunctions
    use type_STAT
    use math_RF


    implicit none

    !INPUTS
    integer :: code, rang, nb_procs, comm;

    type(STAT) :: STA

    character(len=200) :: resPath
    character(len=200), parameter :: meshPath = "./mesh_input"
    character(len=200), parameter :: genPath = "./gen_input"
    character(len=200), parameter :: statPath = "./stat_input"
    character(len=200), parameter :: outPath = "./results/res/singleGen"
    integer :: nFiles

    !LOCAL
    integer :: i

    comm = MPI_COMM_WORLD
    call init(comm)
    STA%comm     = comm
    STA%rang     = rang
    STA%nb_procs = nb_procs

    if(STA%rang == 0) then
        write(*,*) "How many files? "
        read(*,*) nFiles
        write(*,*) "nFiles = ", nFiles
    end if

    call MPI_BCAST (nFiles, 1, MPI_INTEGER, 0, STA%comm, code)

    if(STA%rang == 0) then
        write(*,*) "  -----------------------------------------------"
        write(*,*) "  -------------CALCULATING STATISTICS------------"
        write(*,*) "  -----------------------------------------------"
        write(*,*) " "
    end if

    !write(*,*) "STA%rang = ", STA%rang

    do i = 1, nFiles
        if(STA%rang == 0) then
            write(*,*) "File ", i
            read(*,*) resPath
            write(*,*) "path = ", resPath
        end if

        call MPI_BCAST (resPath, len(resPath), MPI_CHARACTER, 0, STA%comm, code)

        !write(*,*) "path = ", resPath

        if(STA%rang == 0) write(*,*) "-> Reading HDF5 Attributes"
        call read_RF_h5_File_Attributes()
        if(STA%rang == 0) write(*,*) "-> Setting Local Range"
        call set_Local_Range()
        call set_Sk_Dir()
        if(STA%rang == 0) write(*,*) "-> Reading HDF5 Tables"
        call read_RF_h5_File_Table()
        if(STA%rang == 0) write(*,*) "-> Calculating Average and StdVar"
        call calculate_average_and_stdVar_MPI(STA)
        if(STA%rang == 0) write(*,*) "-> Recontructing Spectrum"
        call rebuild_Sk(STA)
        if(STA%rang == 0) write(*,*) "-> Calculating Correlation Length"
        call rebuild_corrL(STA, STA%corrL_out)
        if(STA%rang == 0) write(*,*) "-> Writing Statistics on File"
        if(STA%rang == 0) call write_StatisticsOnH5(STA, resPath)

        call show_STAT(STA, "Calculated Statistics", 6)

        call finalize_STAT(STA)

    end do

    call finalize()

    contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init(comm_local)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local

            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rang, code)
            call MPI_COMM_SIZE(comm_local, nb_procs, code)

        end subroutine init

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_RF_h5_File_Attributes()

            !LOCAL
            integer :: nDim, Nmc, method, corrMod, margiFirst
            logical :: independent
            character(len=50) :: attr_Name
            integer :: hdferr
            integer(HID_T) :: file_id


            if(STA%rang == 0) write(*,*) " Searching for file: ",resPath

            call h5open_f(hdferr) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(resPath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
            if(hdferr /= 0) stop("ERROR OPENING FILE")
            !write(*,*) "hdferr = ", hdferr

            !READING SCALARS----------------------------
            !BOOL
            !attr_name = "independent"
            !call read_h5attr_bool(file_id, trim(adjustL(attr_name)), STA%independent)

            !INTEGERS
            attr_name = "nDim"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), nDim)
            attr_name = "Nmc"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), Nmc)
            attr_name = "method"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), method)
            attr_name = "corrMod"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), corrMod)
            attr_name = "margiFirst"
            call read_h5attr_int(file_id, trim(adjustL(attr_name)), margiFirst)

            !Allocating Vectors
            call init_STAT(STA, nDim, Nmc, method, corrMod, margiFirst, independent)

            !DOUBLE VEC
            attr_name = "xMinGlob"
            call read_h5attr_real_vec(file_id, attr_name, STA%xMinGlob)
            attr_name = "xMaxGlob"
            call read_h5attr_real_vec(file_id, attr_name, STA%xMaxGlob)
            attr_name = "xStep"
            call read_h5attr_real_vec(file_id, attr_name, STA%xStep)
            attr_name = "corrL"
            call read_h5attr_real_vec(file_id, attr_name, STA%corrL)
            attr_name = "overlap"
            call read_h5attr_real_vec(file_id, attr_name, STA%overlap)

            STA%xNStep = nint((STA%xMaxGlob-STA%xMinGlob)/STA%xStep,8) +1

            call h5fclose_f(file_id, hdferr) ! Close the file.
            call h5close_f(hdferr) ! Close FORTRAN interface.

        end subroutine read_RF_h5_File_Attributes

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine write_StatisticsOnH5(STA, resPath)
            implicit none
            !INPUT
            type(STAT) :: STA
            character(len=*), intent(in) :: resPath
            !LOCAL
            character(len=200) :: attr_name
            integer(HID_T)  :: file_id       !File identifier
            integer :: error
            integer :: i
            logical :: attr_exists

            write(*,*) "Writing statistics on: ", trim(resPath)
            call h5open_f(error) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(resPath), H5F_ACC_RDWR_F, file_id, error) !Open File

            !DOUBLE
            !write(*,*) "DOUBLE"
            attr_name = "globalAvg"
            call h5aexists_by_name_f(file_id, ".", trim(adjustL(attr_name)), attr_exists, error)
            !write(*,*) "attr_exists 1 = ", attr_exists
            if(attr_exists) then
                call h5adelete_f(file_id, trim(adjustL(attr_name)), error)
            end if
            call write_h5attr_real(file_id, trim(adjustL(attr_name)), STA%globalAvg)

            attr_name = "globalStdDev"
            call h5aexists_by_name_f(file_id, ".", trim(adjustL(attr_name)), attr_exists, error)
            !write(*,*) "attr_exists 2 = ", attr_exists
            if(attr_exists) then
                call h5adelete_f(file_id, trim(adjustL(attr_name)), error)
            end if
            call write_h5attr_real(file_id, trim(adjustL(attr_name)), STA%globalStdDev)


            !DOUBLE VEC
            !write(*,*) "DOUBLE VEC"
            attr_name = "evntAvg"
            call h5aexists_by_name_f(file_id, ".", trim(adjustL(attr_name)), attr_exists, error)
            !write(*,*) "attr_exists 3 = ", attr_exists
            if(attr_exists) then
                call h5adelete_f(file_id, trim(adjustL(attr_name)), error)
            end if
            call write_h5attr_real_vec(file_id, attr_name, STA%evntAvg)


            attr_name = "evntStdDev"
            call h5aexists_by_name_f(file_id, ".", trim(adjustL(attr_name)), attr_exists, error)
            !write(*,*) "attr_exists 4 = ", attr_exists
            if(attr_exists) then
                call h5adelete_f(file_id, trim(adjustL(attr_name)), error)
            end if
            call write_h5attr_real_vec(file_id, attr_name, STA%evntStdDev)


            attr_name = "corrL_out"
            call h5aexists_by_name_f(file_id, ".", trim(adjustL(attr_name)), attr_exists, error)
            !write(*,*) "attr_exists 5 = ", attr_exists
            if(attr_exists) then
                call h5adelete_f(file_id, trim(adjustL(attr_name)), error)
            end if
            call write_h5attr_real_vec(file_id, attr_name, STA%corrL_out)

            !write(*,*) "Sk"
            do i = 1, STA%nDim
                attr_name = stringNumb_join("Sk_",i)
                call h5aexists_by_name_f(file_id, ".", trim(adjustL(attr_name)), attr_exists, error)
                !write(*,*) "attr_exists 6 7 8= ", attr_exists
                if(attr_exists) then
                    call h5adelete_f(file_id, trim(adjustL(attr_name)), error)
                end if
                call write_h5attr_real_vec(file_id, attr_name, STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2)))

            end do

            call h5fclose_f(file_id, error)! Close the file.
            call h5close_f(error) ! Close FORTRAN interface

        end subroutine write_StatisticsOnH5

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine set_Local_Range()
            implicit none

            !LOCAL
            integer :: topComm
            logical, dimension(STA%nDim) :: periods
            integer :: code

            periods(:) = .false.

            call set_procPerDim_2 (STA, STA%procPerDim)
            call MPI_CART_CREATE (STA%comm, STA%nDim, STA%procPerDim, periods, .false., topComm, code)
            call MPI_CART_COORDS (topComm, STA%rang, STA%nDim, STA%coords, code)

            STA%localRange(:,1) = STA%coords * STA%xNStep/STA%procPerDim + 1
            STA%localRange(:,2) = (STA%coords + 1) * STA%xNStep/STA%procPerDim + 1

            where(STA%coords == STA%procPerDim-1)
                STA%localRange(:,2) = STA%xNStep
            elsewhere
                STA%localRange(:,2) = STA%localRange(:,2) - 1
            end where

            STA%xNStep_Loc = STA%localRange(:,2) - STA%localRange(:,1) + 1

        end subroutine set_Local_Range

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_procPerDim_2 (STA, procPerDim)
        implicit none

        !INPUT
        type(STAT), intent(in) :: STA

        !OUTPUT
        integer, dimension(:) :: procPerDim

        !LOCAL VARIABLES
        integer :: i, j;
        double  precision :: procRootDim, logProc2;

        if(size(procPerDim)/=STA%nDim) then
            write(*,*) "Error inside 'set_procPerDim_2', dimensions are not compatible"
            write(*,*) "size(procPerDim) = ", size(procPerDim)
            write(*,*) "MSH%nDim         = ", STA%nDim
            stop(" ")
        end if


        procRootDim = dble(STA%nb_procs)**(1/dble(STA%nDim))
        logProc2   = log(dble(STA%nb_procs))/log(2.0D0)

        if (areEqual(procRootDim, dble(nint(procRootDim)))) then
            if(STA%rang == 0) write(*,*) "    Exact Division"
            procPerDim(:) = nint(dble(STA%nb_procs)**(1.0d0/STA%nDim))

        else if(areEqual(logProc2, dble(nint(logProc2)))) then
            if(STA%rang == 0) write(*,*) "    Power of two"

            procPerDim(:) = 1
            if(STA%nb_procs /= 1) then
                do j = 1, nint(logProc2)
                    i = cyclicMod(j, STA%nDim)
                    procPerDim(i) = procPerDim(i)*2
                end do
            end if

        else
            stop ("ERROR, inside 'set_procPerDim_2', no mesh division algorithm for this number of procs")
        end if

        if(STA%rang == 0) write(*,*) "    procPerDim = ", procPerDim

    end subroutine set_procPerDim_2

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine set_Sk_Dir()
            implicit none

            !LOCAL
            integer :: i

            allocate(STA%Sk_Dir(sum(STA%xNStep_Loc)))

            do i = 1, STA%nDim

                if(i == 1) then
                    STA%Sk_Ind(i,1) = 1
                else
                    STA%Sk_Ind(i,1) = sum(STA%xNStep_Loc(1:i-1)) + 1
                end if

                STA%Sk_Ind(i,2) = sum(STA%xNStep_Loc(1:i))

            end do

        end subroutine set_Sk_Dir

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_RF_h5_File_Table()

            !LOCAL
            character(len=50) :: dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, space_id, dset_id, mem_id
            !integer(HSIZE_T), dimension(STA%nDim) :: dims, maxdims
            integer(HSIZE_T), dimension(STA%nDim) :: offset, locDims
            integer(HSIZE_T), dimension(2) :: locShape
            integer(HSIZE_T), dimension(2) :: zero2D
            !double precision, dimension(:,:), allocatable :: locRF

            locDims(:) = STA%localRange(:,2) - STA%localRange(:,1) + 1
            !write(*,*) " locDims = ", locDims

            if(STA%rang == 0) write(*,*) " Searching for file: ",resPath
            call h5open_f(hdferr) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(resPath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
            if(hdferr /= 0) stop("ERROR OPENING FILE inside read_RF_h5_File_Table")
            !write(*,*) "hdferr = ", hdferr

            !Allocating Vectors

            !MATRIX
            call h5dopen_f(file_id, trim(dset), dset_id, hdferr)! Open Dataset
            call h5dget_space_f(dset_id, space_id, hdferr) !Open Dataset Space
            !call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr) !Measure Dataset Space

            !allocate(STA%randField(product(dims),1))
            allocate(STA%randField(product(locDims),1))
            offset = STA%localRange(:,1)-1
            locShape = shape(STA%randField)
            zero2D = 0
            !write(*,*) " locShape = ", locShape
            !write(*,*) " offset   = ", offset
            !write(*,*) " locDims  = ", locDims
            !For hyperslab lecture
            !IN
            !if(STA%independent) then
                call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, locDims, hdferr) !Select Hyperslab IN
            !else
            !    write(*,*) "Developement to non-independent INPUT"
            !end if
            !call h5sselect_hyperslab_f(space_id, H5S_SELECT_OR, offset, locDims, hdferr) !Add Selection to he hyperslab IN

            !OUT
            call h5screate_simple_f(2, locShape, mem_id, hdferr) !Create memory dataspace
            call h5sselect_hyperslab_f(mem_id, H5S_SELECT_SET_F, zero2D, locShape, hdferr) !Select Hyperslab OUT
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, STA%randField, locShape, hdferr, mem_id, space_id) !Read Dataset Hyperslab

            !call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, STA%randField, dims, hdferr) !Read Dataset
            call h5dclose_f(dset_id, hdferr) !Close Dataset
            call h5sclose_f(space_id, hdferr) !Close Dataset Space

            !write(*,*) "dims    = ", dims
            !write(*,*) "maxdims = ", maxdims
            !if(STA%rang == 1) write(*,*) "STA%randField(:,1)            = ", STA%randField(:,1)
            !write(*,*) "STA%randField(1:10,1)    = ", STA%randField(1:10,1)
            !write(*,*) "STA%randField(160:170,1) = ", STA%randField(160:170,1)

            call h5fclose_f(file_id, hdferr) ! Close the file.
            call h5close_f(hdferr) ! Close FORTRAN interface.

            !if(allocated(locRF)) deallocate(locRF)

        end subroutine read_RF_h5_File_Table

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine finalize()
            implicit none
            call finalize_STAT(STA)
            call MPI_FINALIZE(code)
        end subroutine finalize

end program main_Stat
