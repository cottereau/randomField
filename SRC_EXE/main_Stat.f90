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
    integer :: code, commLocal, rang, nb_procs, comm;
    integer :: nDim, Nmc, method, corrMod, margiFirst
    type(STAT) :: STA
    character(len=70) :: testeChar, testInt
    character(len=200), parameter :: resPath = "./results/res/h5/samples-ALLprocStruct.h5"
    character(len=200), parameter :: meshPath = "./mesh_input"
    character(len=200), parameter :: genPath = "./gen_input"
    character(len=200), parameter :: outPath = "./results/res/singleGen"

    comm = MPI_COMM_WORLD

    call init(comm)
    if(STA%rang == 0) then
        write(*,*) "  -----------------------------------------------"
        write(*,*) "  -------------CALCULATING STATISTICS------------"
        write(*,*) "  -----------------------------------------------"
        write(*,*) " "
    end if
    !call read_RF_h5_File_Attributes()
    !call set_Local_Range()
    !call set_Sk_Dir()
    !call read_RF_h5_File_Table()
    !call calculate_average_and_stdVar_MPI(STA)
    !call rebuild_Sk(STA)
    !call show_STAT(STA, "HEY STAT", 6)

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

            STA%comm     = comm
            STA%rang     = rang
            STA%nb_procs = nb_procs

        end subroutine init

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_RF_h5_File_Attributes()

            !LOCAL
            integer :: nDim, Nmc, method, corrMod, margiFirst
            logical :: independent
            character(len=50) :: attr_Name, dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, attr_id, space_id, dset_id
            double precision, dimension(:), allocatable :: locData

            if(STA%rang == 0) write(*,*) " Searching for file: ",resPath

            call h5open_f(hdferr) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(resPath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
            if(hdferr /= 0) stop("ERROR OPENING FILE")
            !write(*,*) "hdferr = ", hdferr

            !READING SCALARS----------------------------
            !BOOL
            attr_name = "independent"
            call read_h5attr_bool(file_id, trim(adjustL(attr_name)), STA%independent)

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
        subroutine set_Local_Range()
            implicit none

            !LOCAL
            integer :: topComm
            logical, dimension(STA%nDim) :: periods
            integer :: code

            periods(:) = .false.

            write(*,*) "-> Setting Local Range"
            call set_procPerDim (STA%nb_procs, STA%nDim, STA%procPerDim)
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
            character(len=50) :: attr_Name, dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, attr_id, space_id, dset_id, mem_id
            !integer(HSIZE_T), dimension(STA%nDim) :: dims, maxdims
            integer(HSIZE_T), dimension(STA%nDim) :: offset, locDims
            integer(HSIZE_T), dimension(2) :: locShape
            integer(HSIZE_T), dimension(2) :: zero2D
            !double precision, dimension(:,:), allocatable :: locRF

            locDims(:) = STA%localRange(:,2) - STA%localRange(:,1) + 1
            write(*,*) " locDims = ", locDims

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
            write(*,*) " locShape = ", locShape
            write(*,*) " offset   = ", offset
            write(*,*) " locDims  = ", locDims
            !For hyperslab lecture
            !IN
            if(STA%independent) then
                call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, locDims, hdferr) !Select Hyperslab IN
            else
                write(*,*) "Developement to non-independent INPUT"
            end if
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
