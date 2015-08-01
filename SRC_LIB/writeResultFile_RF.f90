module writeResultFile_RF

    use displayCarvalhol
    use math_RF
    use hdf5
    use mpi
    use constants_RF
    use write_Log_File
    use type_RF
    use type_MESH

contains


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_Mono_XMF_h5(RDF, MSH, connectList, monotype, fileName, rang, folderPath, &
                                            communicator, labelsH5, indexesH5, indexXMF, style)

        implicit none

        !INPUTS
        type(RF), intent(in)   :: RDF;
        type(MESH), intent(in) :: MSH;
        integer         , dimension(1:,1:), intent(in) :: connectList;
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator, indexXMF
        character(len=*), dimension(1:)   , intent(in) :: labelsH5
        integer         , dimension(1:)   , intent(in) :: indexesH5
        logical, intent(in) :: monotype
        integer, intent(in) :: style !1 for parallel h5 writing
                                     !2 for sequential per processor h5 writing
                                     !3 for gathered monoprocesor h5 writing
        !LOCAL
        character(len=110), dimension(3) :: HDF5Names
        integer, dimension(3) :: xSz, ySz
        integer, dimension(MSH%nDim) :: sizeSamples
        character(len=110) :: XMFName, HDF5Name, fileName2


        !!!!!!!!!!!!HDF5
        write(get_fileId(),*) "-> Writing h5 file in", trim(adjustL(folderPath))//"/h5";

        select case (style)
            case(1)
                select case (MSH%meshMod)
                    case("unv")
                        fileName2 = "samples"
                        call write_pHDF5_Unstr(double_Data=RDF%randField, &
                                               fileName=fileName2, &
                                               rang=rang, &
                                               folderPath=trim(adjustL(folderPath))//"/h5", &
                                               communicator=communicator, &
                                               HDF5Name=HDF5Names(1), xSz=xSz(1), ySz=ySz(1), transp = .true.)
                        fileName2 = "nodes"
                        call write_pHDF5_Unstr(double_Data=RDF%xPoints, &
                                               fileName=fileName2, &
                                               rang=rang, &
                                               folderPath=trim(adjustL(folderPath))//"/h5", &
                                               communicator=communicator, &
                                               HDF5Name=HDF5Names(2), xSz=xSz(2), ySz=ySz(2), transp = .false.)
                        fileName2 = "connect"
                        call write_pHDF5_Unstr(integer_Data=connectList-1, &
                                               fileName=fileName2, &
                                               rang=rang, &
                                               folderPath=trim(adjustL(folderPath))//"/h5", &
                                               communicator=communicator, &
                                               HDF5Name=HDF5Names(3), xSz=xSz(3), ySz=ySz(3), transp = .false.)
                    case("automatic")
                        fileName2 = "samples"
                        call write_pHDF5_Str(  MSH=MSH, &
                                               RDF=RDF, &
                                               fileName=fileName2, &
                                               rang=rang, &
                                               folderPath=trim(adjustL(folderPath))//"/h5", &
                                               communicator=communicator, &
                                               HDF5Name=HDF5Name)
                    case default
                        stop("In hdf5 style 1 writing - meshMod not implemented")
                end select
            case(2)
                call write_HDF5_Unstr_per_proc(RDF%xPoints, RDF%randField, fileName, rang, trim(adjustL(folderPath))//"/h5", &
                                   communicator, labelsH5, indexesH5, HDF5Name)
            case default
                stop("hdf5 writing style not implemented")
        end select


        !!!!!!!!!!!!XMF
        write(get_fileId(),*) "-> Writing XMF file in", trim(adjustL(folderPath))//"/xmf";
        XMFName = stringNumb_join(trim(adjustL(fileName))//"it_", indexXMF)

        select case (style)
            case(1)
                select case (MSH%meshMod)
                    case("unv")
                        call write_pHDF5_Unstr_XMF(HDF5Names, xSz, ySz, XMFName, &
                                               rang, trim(adjustL(folderPath))//"/xmf", &
                                               communicator, "../h5")
                    case("automatic")
                        call write_pHDF5_Str_XMF(HDF5Name, MSH, fileName, &
                                                 rang, trim(adjustL(folderPath))//"/xmf", &
                                                 communicator, "../h5")
                    case default
                        stop("In XMF style 1 writing - meshMod not implemented")
                end select


            case(2)
                call write_HDF5_Unstr_per_proc_XMF(1, [HDF5name], [size(RDF%xPoints,2)], [.true.], size(RDF%xPoints,1), XMFName, &
                                    rang, trim(adjustL(folderPath))//"/xmf", &
                                    communicator, "../h5")
            case default
                stop("hdf5/XMF writing style not implemented")
        end select

    end subroutine write_Mono_XMF_h5

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_pHDF5_Unstr(double_Data, integer_Data, fileName, rang, folderPath, &
                                 communicator, HDF5Name, xSz, ySz, transp)
        implicit none

        !INPUTS
        double precision, dimension(1:,1:), intent(in), optional :: double_Data;
        integer         , dimension(1:,1:), intent(in), optional :: integer_Data;
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator
        logical                           , intent(in) :: transp

        !OUTPUTS
        character(len=110) , optional  , intent(out) ::HDF5Name
        integer, intent(out) :: xSz, ySz

        !HDF5 VARIABLES
        character(len=110)             :: fileHDF5Name, fullPath !File name
        character(len=30)              :: eventName, coordName;        !Dataset names
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: dspace_id     !Dataspace identifier
        integer(HID_T)                 :: memspace      ! Dataspace identifier in memory
        integer(HID_T)                 :: plist_id      ! Property list identifier
        integer(HID_T)                 :: filespace     ! Dataspace identifier in file
        integer                        :: rank = 2      !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(2) :: dims !Dataset dimensions
        integer                        :: error !Error flag
        integer                        :: info, code
        integer(HSIZE_T), dimension(2) :: count
        integer(HSSIZE_T), dimension(2) :: offset

        !LOCAL VARIABLES
        integer :: yDim, xDim, i
        integer :: nb_procs
        character (len=12) :: numberStr, rangStr;
        integer, dimension(:), allocatable :: all_n_Dim
        character(LEN=8) :: dsetname



        write(get_fileId(),*) "------------START Writing result HDF5 file (MPI)-----------------------";
        write(get_fileId(),*) "fileName         = ", fileName
        write(get_fileId(),*) "folderPath       = ", folderPath

        if(present(double_Data) .eqv. present(integer_Data)) then
            write(*,*) "present(double_Data)  = ", present(double_Data)
            write(*,*) "present(integer_Data) = ", present(integer_Data)
            stop ("Too many ot too few arguments inside write_pHDF5_Unstr")
        end if

        if(present(integer_Data)) then
            xDim = size(integer_Data, 2)
            yDim = size(integer_Data , 1)
        else
            xDim = size(double_Data, 2)
            yDim = size(double_Data , 1)
        end if

        info = MPI_INFO_NULL

        !Discovering needed global information
        call MPI_COMM_SIZE(communicator, nb_procs, code)
        allocate(all_n_Dim(nb_procs))
        if(transp) then
            call MPI_ALLGATHER(yDim    , 1, MPI_INTEGER,     &
                               all_n_Dim, 1, MPI_INTEGER,     &
                               communicator    , code)
        else
            call MPI_ALLGATHER(xDim    , 1, MPI_INTEGER,     &
                               all_n_Dim, 1, MPI_INTEGER,     &
                               communicator    , code)
        end if

        !Creating file name
        write(*,*) "xDim = ", xDim
        write(*,*) "yDim = ", yDim
        write(*,*) "all_n_Dim = ", all_n_Dim
        dsetname = trim(adjustL(fileName)) ! Dataset name
        fileHDF5Name = trim(fileName)//"-ALLproc.h5"
        fullPath     = string_join(folderPath,"/"//fileHDF5Name)
        write(get_fileId(),*) "' fileHDF5Name = ", fileHDF5Name

        !PREPARING ENVIROMENT
        dims = [yDim, sum(all_n_Dim)]
        if(transp) dims = [xDim, sum(all_n_Dim)]
        ySz = dims(1)
        xSz = dims(2)
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) !creates property list (plist_id)
        call h5pset_fapl_mpio_f(plist_id, communicator, info, error) !sets property list
        call h5fcreate_f(fullPath, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id) !create file with the property list (file_id)
        call h5pclose_f(plist_id, error) !closes property list (plist_id)
        call h5screate_simple_f(rank, dims, filespace, error) ! Create the dataspace, the size of the whole table (filespace).

        !DEFINITION OF THE TYPE OF THE TABLE
        if(present(integer_Data)) then
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, filespace, &
                          dset_id, error) ! Appropriates the dataspace to the dataset (dset_id)
        else
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
                          dset_id, error) ! Appropriates the dataspace to the dataset (dset_id)
        end if
        call h5sclose_f(filespace, error) ! Closes the dataspace (it was already used (filespace)

        !CHOOSING SPACE IN MEMORY FOR THIS PROC
        count = [yDim, xDim] !Dimensions to write
        if(transp) count = [xDim, yDim]
        offset = [0, sum(all_n_Dim(1:rang))] !Lines Offset to start writing
        call h5screate_simple_f(rank, count, memspace, error)  !Initialize memspace

        ! Select hyperslab in the file.
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)

        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        !! Write the dataset collectively.
        if(present(integer_Data)) then
            if(transp) then
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, transpose(integer_Data), dims, error, &
                        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, integer_Data, dims, error, &
                        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            end if
        else
            if(transp) then
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, transpose(double_Data), dims, error, &
                        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, double_Data, dims, error, &
                        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            end if
        end if
        ! Close dataspaces.
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)

        ! Close the dataset and property list.
        call h5dclose_f(dset_id, error)
        call h5pclose_f(plist_id, error)

        ! Close the file.
        call h5fclose_f(file_id, error)

        ! Close FORTRAN predefined datatypes.
        call h5close_f(error)

        !write(*,*) "fileHDF5Name = ", fileHDF5Name
        if(present(HDF5Name)) then
            HDF5Name = trim(adjustL(fileHDF5Name))
            HDF5Name = adjustL(HDF5Name)
            write(get_fileId(),*) "'inside write HDF5' output -- HDF5Name = ", HDF5Name
        end if

        write(get_fileId(),*) "------------END Writing result HDF5 file (MPI)-----------------------";

    end subroutine write_pHDF5_Unstr


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_pHDF5_Unstr_XMF(HDF5nameList, xSz, ySz, fileName, &
                                     rang, folderPath, &
                                     communicator, HDF5relativePath)
        implicit none

        !INPUTS
        character(len=*), dimension(:), intent(in) :: HDF5nameList
        integer, dimension(:), intent(in) :: xSz, ySz
        character(len=*)                  , intent(in) :: filename;
        integer, intent(in) :: rang
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator
        character(len=*)                  , intent(in) :: HDF5relativePath

        !LOCAL VARIABLES
        integer             :: Nmc, i, j, file, nDim;
        character (len=110) :: fileXMFName, fullPathXMF, HDF5path;
        character (len=35)  :: eventName, meshName;

        write(get_fileId(),*) "------------START Writing result XMF file-----------------------";

        Nmc  = ySz(1)
        nDim = ySz(2)


        if(rang == 0 .and. (nDim == 1 .or. nDim == 2 .or. nDim == 3)) then

            fileXMFName = string_join(fileName,".xmf")
            fullPathXMF = string_join(folderPath, "/"//fileXMFName)
            HDF5path = string_join(HDF5relativePath, "/")
            write(get_fileId(),*) "fileXMFName = ", fileXMFName
            write(get_fileId(),*) "fullPathXMF = ", fullPathXMF
            write(get_fileId(),*) "HDF5path    = ", HDF5path

            !Building file
            file=21;
            !write(*,*) "Flag 2"
            open (unit = file , file = fullPathXMF, action = 'write')

            write (file,'(A)'      )'<?xml version="1.0" ?>'
            write (file,'(A)'      )'<Xdmf Version="2.0">'
            write (file,'(A)'      )' <Domain>'
            i = 1
            write (file,'(5A)'     )'   <DataItem Name="samples" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(xSz(i))), ' ',trim(numb2String(ySz(i)))  ,'">'
            write (file,'(4A)'     )'        ',trim(HDF5path),trim(HDF5nameList(i)),':/samples'
            write (file,'(A)'      )'   </DataItem>'
            i = 2
            write (file,'(5A)'     )'   <DataItem Name="nodes" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(xSz(i))), ' ',trim(numb2String(ySz(i)))  ,'">'
            write (file,'(4A)'     )'        ',trim(HDF5path),trim(HDF5nameList(i)),':/nodes'
            write (file,'(A)'      )'   </DataItem>'
            i = 3
            write (file,'(5A)'     )'   <DataItem Name="connect" Format="HDF" DataType="Int" Dimensions="',trim(numb2String(xSz(i))), ' ',trim(numb2String(ySz(i)))  ,'">'
            write (file,'(4A)'     )'        ',trim(HDF5path),trim(HDF5nameList(i)),':/connect'
            write (file,'(A)'      )'   </DataItem>'

            write (file,'(A)'      )'  <Grid CollectionType="Spatial" GridType="Collection">' !Opens the Collection

            write (file,'(A)'     )'   <Grid Name="Group1">'
            write (file,'(3A)'    )'     <Topology Type="Hexahedron" NumberOfElements="',trim(numb2String(xSz(i))),'">'
            write (file,'(A)'     )'       <DataItem Reference="XML">'
            write (file,'(A)'     )'         /Xdmf/Domain/DataItem[@Name="connect"]'
            write (file,'(A)'     )'       </DataItem>'
            write (file,'(A)'     )'     </Topology>'
            if(nDim == 1) write (file,'(A)'      )'     <Geometry GeometryType="X">'
            if(nDim == 2) write (file,'(A)'      )'     <Geometry GeometryType="XY">'
            if(nDim == 3) write (file,'(A)'      )'     <Geometry GeometryType="XYZ">'
            write (file,'(A)'     )'       <DataItem Reference="XML">'
            write (file,'(A)'     )'         /Xdmf/Domain/DataItem[@Name="nodes"]'
            write (file,'(A)'     )'       </DataItem>'
            write (file,'(A)'     )'     </Geometry>'
            !TODO
            do i = 1, 1
                write (file,'(3A)'     )'     <Attribute Name="RF1" Center="Node" AttributeType="Scalar">'
                write (file,'(A)'     )'       <DataItem Reference="XML">'
                write (file,'(A)'     )'         /Xdmf/Domain/DataItem[@Name="samples"]'
                write (file,'(A)'     )'       </DataItem>'
                write (file,'(A)'      )'     </Attribute>'
            end do

            write (file,'(A)'      )'   </Grid>' !END Writing the data of one subdomain
            write (file,'(A)'      )'  </Grid>' !Closes the Collection
            write (file,'(A)'      )' </Domain>'
            write (file,'(A)'      )'</Xdmf>'

            close(file)
        else if(rang == 0) then
            write(*,*) "No XDMF Model prepared for this dimension"
            write(*,*) "nDim = ", nDim
            write(*,*) "The file won't be created"
        end if

        write(get_fileId(),*) "------------END Writing result XMF file-----------------------";

    end subroutine write_pHDF5_Unstr_XMF

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_HDF5_Unstr_per_proc(xPoints, randField, fileName, rang, folderPath, &
        communicator, labels, indexes, HDF5Name)
        implicit none

        !INPUTS
        double precision,  dimension(1:,1:), intent(in) :: xPoints, randField;
        character (len=*)                , intent(in) :: filename;
        integer                          , intent(in) :: rang;
        character(len=*)                 , intent(in) :: folderPath
        integer                          , intent(in) :: communicator
        character(len=*) , dimension(1:), optional  , intent(in) :: labels
        integer          , dimension(1:), optional  , intent(in) :: indexes

        !OUTPUTS
        character(len=110) , optional  , intent(out) ::HDF5Name

        !HDF5 VARIABLES
        character(len=110)             :: fileHDF5Name, fullPath !File name
        character(len=30)              :: eventName, coordName;        !Dataset names
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: dspace_id     !Dataspace identifier
        integer                        :: rank = 2      !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(2) :: dims !Dataset dimensions
        integer                        :: error !Error flag

        !LOCAL VARIABLES
        integer :: nDim, Nmc, nPoints, i
        integer :: effectComm
        character (len=12) :: numberStr, rangStr;
        !double precision, dimension(:,:), allocatable :: grid_data

        !        if(rang == 0) then
        !                write(*,*) "";
        !               write(*,*) "------------START Writing result HDF5 file (MPI)-----------------------";
        !               write(*,*) "";
        !          end if

          !if(rang == 0) then
              !write (*,*) "lbound(xPoints) = ", lbound(xPoints)
              !write (*,*) "lbound(randField) = ", lbound(randField)
              !write (*,*) "size(xPoints,1) = ", size(xPoints,1)
              !write (*,*) "size(xPoints,2) = ", size(xPoints,2)
              !write (*,*) "xPoints(1, 1) = ", xPoints(1, 1)
              !write (*,*) "xPoints(1:10, :) = ", xPoints(1:10, :)
              !call dispCarvalhol(xPoints(:,:)  , "xPoints(:,:)"  , "F30.5")
              !call dispCarvalhol(randField(:,:), "randField(:,:)", "F30.5")
          !end if

        effectComm = communicator
        nDim       = size(xPoints , 1)
        nPoints    = size(randField, 1)
        Nmc        = size(randField, 2)

        !Creating file name
        !write(*,*) "fileName = ", fileName
        !write(*,*) "labels = ", labels
        !write(*,*) "indexes = ", indexes
        !write(*,*) "folderPath = ", folderPath

        if(.not. present(labels)) then
            write(rangStr,'(I8)'  ) rang
            rangStr      = adjustL(rangStr)
            fileHDF5Name = trim(fileName)//"-proc"//trim(rangStr)//".h5"
        else
            fileHDF5Name = fileName
            do i = 1, size(labels)
                fileHDF5Name =  string_join(fileHDF5Name,stringNumb_join(labels(i), indexes(i)))
            end do
        end if

        !write(*,*) "fileHDF5Name = ", fileHDF5Name

        fileHDF5Name = string_join(fileHDF5Name,".h5")
        fullPath     = string_join(folderPath,"/"//fileHDF5Name)

        !write(*,*) "'inside write HDF5' -- fileHDF5Name = ", fileHDF5Name

        if (nDim > 3) then
            write(*,*) "Dimension exceeds 3, HDF file won't be created"
        else

            !if(rang == 0) write(*,*) ">>>>>>>>> Opening file";
            call h5open_f(error) ! Initialize FORTRAN interface.
            call h5fcreate_f(fullPath, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.

            !if(rang == 0) write(*,*) ">>>>>>>>> Creating Coordinates dataset 'XYZ table'";

            dims = shape(xPoints)
            write(coordName,'(A)') "XYZ"
            call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
            call h5dcreate_f(file_id, coordName, H5T_NATIVE_DOUBLE,  &
                dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                xPoints, dims, error) ! Write the dataset.
            call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
            call h5sclose_f(dspace_id, error) ! Terminate access to the data space.

            !if(rang == 0) write(*,*) ">>>>>>>>> Creating Quantities dataset 'random field'";
            dims(1) = size(randField,1)
            dims(2) = 1 !One random field in each dataset

            do i = 1, Nmc
                write(numberStr,'(I8)'  ) i
                numberStr = adjustL(numberStr)
                write(eventName,'(2A)') "RF_", trim(numberStr)
                !write(*,*) "i = ", i
                !if(Nmc < 11) write(*,*) "eventName = ", eventName

                call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
                call h5dcreate_f(file_id, eventName, H5T_NATIVE_DOUBLE,  &
                    dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField(:,i), dims, error) ! Write the dataset.
                call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
                call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
            end do

            !if(rang == 0) write(*,*) ">>>>>>>>> Closing file";
            call h5fclose_f(file_id, error) ! Close the file.
            call h5close_f(error) ! Close FORTRAN interface.

        end if

        !write(*,*) "fileHDF5Name = ", fileHDF5Name
        if(present(HDF5Name)) then
            HDF5Name = trim(adjustL(fileHDF5Name))
            HDF5Name = adjustL(HDF5Name)
            !write(*,*) "'inside write HDF5' -- HDF5Name = "
            !write(*,*) HDF5Name
        end if

    end subroutine write_HDF5_Unstr_per_proc

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_HDF5_Unstr_per_proc_XMF(nSamples, HDF5nameList, nPointList, mask, nDim, fileName, rang, folderPath, &
                                              communicator, HDF5relativePath, attName, byProc)
        implicit none

        !INPUTS
        integer                           , intent(in) :: nSamples, nDim;
        logical         , dimension(:)    , intent(in) :: mask
        character(len=*), dimension(1:)   , intent(in) :: HDF5nameList;
        integer         , dimension(1:)   , intent(in) :: nPointList;
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator
        character(len=*), optional        , intent(in) :: HDF5relativePath
        character(len=*), dimension(1:)   , optional  , intent(in) :: attName;
        logical                           , optional  , intent(in) :: byProc !To divise the mesh by proc and not by subdomain

        !LOCAL VARIABLES
        integer             :: Nmc, i, j, file, nb_procs, code;
        integer             :: effectComm
        character (len=110) :: fileXMFName, fullPathXMF, effecHDF5path;
        character (len=35)  :: eventName, meshName;
        character (len=50) , dimension(:), allocatable :: effectAttName;
        integer            , dimension(:), allocatable :: all_nPointList
        character (len=110), dimension(:), allocatable :: all_HDF5nameList
        logical            , dimension(:), allocatable :: all_mask

        effectComm = communicator


        write(get_fileId(),*) "------------START Writing result XMF file-----------------------";

        call MPI_COMM_SIZE(effectComm, nb_procs, code)

        !Common parameters
        Nmc         = nSamples
        write(get_fileId(),*) "HDF5nameList = ", HDF5nameList

        if(rang == 0) then
            allocate(all_nPointList(nb_procs*size(HDF5nameList)))
            allocate(all_HDF5nameList(nb_procs*size(HDF5nameList)))
            allocate(all_mask(nb_procs*size(mask)))
            allocate(effectAttName(Nmc))
        end if

        call MPI_GATHER(nPointList    , size(nPointList), MPI_INTEGER,     &
                        all_nPointList, size(nPointList), MPI_INTEGER,     &
                        0         , effectComm    , code)

        call MPI_GATHER(HDF5nameList    , len(HDF5nameList)*size(HDF5nameList), MPI_CHARACTER,     &
                        all_HDF5nameList, len(HDF5nameList)*size(HDF5nameList), MPI_CHARACTER,     &
                        0              , effectComm       , code)

        call MPI_GATHER(mask    , size(mask), MPI_LOGICAL,     &
                        all_mask, size(mask), MPI_LOGICAL,     &
                        0      , effectComm       , code)

        if(rang == 0) write(get_fileId(),*) "all_nPointList   = ", all_nPointList
        if(rang == 0) write(get_fileId(),*) "all_HDF5nameList = ", all_HDF5nameList
        if(rang == 0) write(get_fileId(),*) "all_mask = ", all_mask

        if(rang == 0 .and. (nDim == 1 .or. nDim == 2 .or. nDim == 3)) then

            fileXMFName = string_join(fileName,".xmf")
            fullPathXMF = string_join(folderPath, "/"//fileXMFName)
            !if(rang == 0) write(*,*) "all_nPointList in rang 0 = ", all_nPointList
            !if(rang == 0) write(*,*) "all_HDF5nameList in rang 0 = ", all_HDF5nameList
            write(get_fileId(),*) "fileXMFName = ", fileXMFName
            write(get_fileId(),*) "fullPathXMF = ", fullPathXMF
            !Optional inputs
            !write(*,*) "treating attName"
            if(present(attName)) then
                if (size(attName)==Nmc) then
                    effectAttName = attName
                end if
            else
                do i = 1, Nmc
                    effectAttName(i) = stringNumb_join("RF_", i)
                end do
            end if

            !write(*,*) "treating HDF5relativePath"
            if(present(HDF5relativePath)) then
                effecHDF5path = string_join(HDF5relativePath, "/")
            else
                effecHDF5path = "./"
            end if

            write(get_fileId(),*) "effecHDF5path   = ", effecHDF5path

            !Building file
            file=21;
            !write(*,*) "Flag 2"
            open (unit = file , file = fullPathXMF, action = 'write')

            write (file,'(A)'      )'<?xml version="1.0" ?>'
            write (file,'(A)'      )'<Xdmf Version="2.0">'
            write (file,'(A)'      )' <Domain>'
            write (file,'(A)'      )'  <Grid CollectionType="Spatial" GridType="Collection">' !Opens the Collection

            do j = 1, size(all_HDF5nameList)
                if(all_mask(j)) then
                    !write(*,*) trim(all_HDF5nameList(j))
                    !write(numberStr,'(I)'  ) j
                    !numberStr = adjustL(numberStr)
                    !write(meshName,'(2A)' ) "meshRF_", trim(adjustL(numberStr))
                    meshName = trim(stringNumb_join("SubD_", mod(j-1,size(HDF5nameList))))//&
                               trim(stringNumb_join("Proc_", int((j-1)/size(HDF5nameList))))
                    if(present(byProc)) then
                        if(byProc) meshName = stringNumb_join("Proc_", int((j-1)/size(HDF5nameList)))
                        !write(*,*) "                       j-1 = ",j-1
                        !write(*,*) "        size(HDF5nameList) = ",size(HDF5nameList)
                        !write(*,*) "int((j-1)/size(HDF5nameList) = ",int((j-1)/size(HDF5nameList))
                    end if
                        !write(*,*) "size(all_HDF5nameList) = ", size(all_HDF5nameList)
                        !write(*,*) "trim(effecHDF5path) = ", trim(effecHDF5path)
                        !write(*,*) "trim(all_HDF5nameList(j)) = ", trim(all_HDF5nameList(j))
                        !write(*,*) "'           ',trim(effecHDF5path),trim(all_HDF5nameList(j)),':/XYZ'= ", '           ',trim(effecHDF5path),trim(all_HDF5nameList(j)),':/XYZ'
                    write (file,'(3A)'     )'   <Grid Name="',trim(meshName),'" GridType="Uniform">' !START Writing the data of one subdomain
                    write (file,'(3A)'     )'    <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',trim(numb2String(all_nPointList(j))),'">'
                    write (file,'(A)'      )'    </Topology>'
                    if(nDim == 1) write (file,'(A)'      )'     <Geometry GeometryType="X">'
                    if(nDim == 2) write (file,'(A)'      )'     <Geometry GeometryType="XY">'
                    if(nDim == 3) write (file,'(A)'      )'     <Geometry GeometryType="XYZ">'
                    write (file,'(5A)'     )'      <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(all_nPointList(j))), ' ',trim(numb2String(nDim))  ,'">'
                    write (file,'(4A)'     )'           ',trim(effecHDF5path),trim(all_HDF5nameList(j)),':/XYZ'
                    write (file,'(A)'      )'      </DataItem>'
                    write (file,'(A)'      )'    </Geometry>'

                    do i = 1, Nmc
                        write (file,'(3A)'     )'     <Attribute Name="',trim(effectAttName(i)),'" Center="Node" AttributeType="Scalar">'
                        write (file,'(3A)'     )'      <DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(all_nPointList(j))),'">'
                        write (file,'(5A)'     )'          ',trim(effecHDF5path),trim(all_HDF5nameList(j)),":/", trim(stringNumb_join("RF_", i))
                        write (file,'(A)'      )'       </DataItem>'
                        write (file,'(A)'      )'     </Attribute>'
                    end do

                    write (file,'(A)'      )'   </Grid>' !END Writing the data of one subdomain
                end if
            end do !END Loop Over HDF5 names

            write (file,'(A)'      )'  </Grid>' !Closes the Collection
            write (file,'(A)'      )' </Domain>'
            write (file,'(A)'      )'</Xdmf>'

            close(file)
        else if(rang == 0) then
            write(*,*) "No XDMF Model prepared for this dimension"
            write(*,*) "nDim = ", nDim
            write(*,*) "The file won't be created"
        end if

        if(rang == 0) then
            deallocate(all_nPointList)
            deallocate(all_HDF5nameList)
            deallocate(effectAttName)

            write(get_fileId(),*) "------------END Writing result XMF file-----------------------";

        end if

    end subroutine write_HDF5_Unstr_per_proc_XMF

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_pHDF5_Str(MSH, RDF, fileName, rang, folderPath, &
                                 communicator, HDF5Name)
        implicit none

        !INPUTS
        type(MESH), intent(in) :: MSH
        type(RF)  , intent(in) :: RDF
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator

        !OUTPUTS
        character(len=110) , optional  , intent(out) ::HDF5Name

        !HDF5 VARIABLES
        character(len=110)             :: fileHDF5Name, fullPath !File name
        character(len=30)              :: eventName, coordName;        !Dataset names
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: dspace_id     !Dataspace identifier
        integer(HID_T)                 :: memspace      ! Dataspace identifier in memory
        integer(HID_T)                 :: plist_id      ! Property list identifier
        integer(HID_T)                 :: filespace     ! Dataspace identifier in file
        integer                        :: rank !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(MSH%nDim) :: dims !Dataset dimensions
        integer                        :: error !Error flag
        integer                        :: info, code
        integer(HSIZE_T) , dimension(MSH%nDim) :: count
        integer(HSSIZE_T), dimension(MSH%nDim) :: offset

        !LOCAL VARIABLES
        integer :: yDim, xDim, i
        integer :: nb_procs
        character (len=12) :: numberStr, rangStr;
        integer, dimension(MSH%nDim) :: total_xNStep
        character(LEN=8) :: dsetname



        write(get_fileId(),*) "------------START Writing result HDF5 file (MPI)-----------------------";
        write(get_fileId(),*) "fileName         = ", fileName
        write(get_fileId(),*) "folderPath       = ", folderPath


        info = MPI_INFO_NULL
        rank = MSH%nDim

        !Discovering needed global information
        write(get_fileId(),*) "RDF%xNStep = ", RDF%xNStep
        write(get_fileId(),*) "MSH%nDim   = ", MSH%nDim
        total_xNStep = nint((MSH%xMaxGlob - MSH%xMinGlob)/MSH%xStep) + 1
        !call MPI_ALLREDUCE (RDF%xNStep, total_xNStep, MSH%nDim, MPI_INTEGER, MPI_SUM, communicator, code)


        !Creating file name
        dsetname = trim(adjustL(fileName)) ! Dataset name
        fileHDF5Name = trim(fileName)//"-ALLprocStruct.h5"
        fullPath     = string_join(folderPath,"/"//fileHDF5Name)
        write(get_fileId(),*) "' fileHDF5Name = ", fileHDF5Name

        !PREPARING ENVIROMENT
        dims = total_xNStep
        !dims = dims(size(dims):1:-1) !transpose
        write(get_fileId(),*) "dims = ", dims

        !if(transp) dims = [xDim, sum(all_n_Dim)]
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) !creates property list (plist_id)
        call h5pset_fapl_mpio_f(plist_id, communicator, info, error) !sets property list
        call h5fcreate_f(fullPath, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id) !create file with the property list (file_id)
        call h5pclose_f(plist_id, error) !closes property list (plist_id)
        call h5screate_simple_f(rank, dims, filespace, error) ! Create the dataspace, the size of the whole table (filespace).

        !DEFINITION OF THE TYPE OF THE TABLE

        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
                          dset_id, error) ! Appropriates the dataspace to the dataset (dset_id)
        call h5sclose_f(filespace, error) ! Closes the dataspace (it was already used (filespace)

        !CHOOSING SPACE IN MEMORY FOR THIS PROC
        count = RDF%xNStep
        !count = count(size(count):1:-1) !transpose
        write(get_fileId(),*) "count = ", count
        call h5screate_simple_f(rank, count, memspace, error)  !Initialize memspace

        ! Select hyperslab in the file.
        offset = RDF%origin - 1!Lines Offset to start writing
        write(get_fileId(),*) "offset = ", offset
        !offset = [0, sum(all_n_Dim(1:rang))] !Lines Offset to start writing
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)

        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Write dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, RDF%randField, dims, error, &
                file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

        ! Close dataspaces.
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)

        ! Close the dataset and property list.
        call h5dclose_f(dset_id, error)
        call h5pclose_f(plist_id, error)

        ! Close the file.
        call h5fclose_f(file_id, error)

        ! Close FORTRAN predefined datatypes.
        call h5close_f(error)

        !write(*,*) "fileHDF5Name = ", fileHDF5Name
        if(present(HDF5Name)) then
            HDF5Name = trim(adjustL(fileHDF5Name))
            HDF5Name = adjustL(HDF5Name)
            write(get_fileId(),*) "'inside write HDF5' output -- HDF5Name = ", HDF5Name
        end if

        write(get_fileId(),*) "------------END Writing result HDF5 file (MPI)-----------------------";

    end subroutine write_pHDF5_Str

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_pHDF5_Str_XMF(HDF5nameList, MSH, fileName, &
                                     rang, folderPath, &
                                     communicator, HDF5relativePath)

        implicit none

        !INPUTS
        type(MESH), intent(in) :: MSH
        character(len=*), intent(in) :: HDF5nameList
        character(len=*)              , intent(in) :: filename;
        integer, intent(in) :: rang
        character(len=*)              , intent(in) :: folderPath
        integer                       , intent(in) :: communicator
        character(len=*)              , intent(in) :: HDF5relativePath

        !LOCAL VARIABLES
        integer             :: i, j, file, nDim;
        character (len=110) :: fileXMFName, fullPathXMF, HDF5path, dimText;
        character (len=35)  :: eventName, meshName;
        integer, dimension(MSH%nDim) :: total_xNStep


        write(get_fileId(),*) "------------START Writing result XMF file-----------------------";

        nDim = MSH%nDim

        if(rang == 0 .and. (nDim == 1 .or. nDim == 2 .or. nDim == 3)) then

            fileXMFName = string_join(fileName,".xmf")
            fullPathXMF = string_join(folderPath, "/"//fileXMFName)
            HDF5path = string_join(HDF5relativePath, "/")
            write(get_fileId(),*) "fileXMFName = ", fileXMFName
            write(get_fileId(),*) "fullPathXMF = ", fullPathXMF
            write(get_fileId(),*) "HDF5path    = ", HDF5path

            total_xNStep = nint((MSH%xMaxGlob - MSH%xMinGlob)/MSH%xStep) + 1

            dimText = ""
            do i = size(total_xNStep), 1, -1
                dimText = trim(dimText)//" "//trim(numb2String(total_xNStep(i)))
            end do
            dimText = trim(adjustL(dimText))

            write(get_fileId(),*) "dimText    = ", dimText

            !Building file
            file=21;
            !write(*,*) "Flag 2"
            open (unit = file , file = fullPathXMF, action = 'write')

            write (file,'(A)'      )'<?xml version="1.0" ?>'
            write (file,'(A)'      )'<Xdmf Version="2.0">'
            write (file,'(A)'      )' <Domain>'
            i = 1
            write (file,'(3A)'     )'   <DataItem Name="samples" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(dimText),'">'
            write (file,'(4A)'     )'        ',trim(HDF5path),trim(HDF5nameList),':/samples'
            write (file,'(A)'      )'   </DataItem>'
            write (file,'(A)'      )'  <Grid GridType="Collection" CollectionType="Spatial">' !Opens the Collection

            write (file,'(A)'     )'   <Grid Name="Group1">'
            write (file,'(3A)'    )'     <Topology TopologyType="3DCoRectMesh" Dimensions="',trim(dimText),'"/>'
            if(nDim == 1) write (file,'(A)'      )'     <Geometry GeometryType="ORIGIN_DX">'
            if(nDim == 2) write (file,'(A)'      )'     <Geometry GeometryType="ORIGIN_DXDY">'
            if(nDim == 3) write (file,'(A)'      )'     <Geometry GeometryType="ORIGIN_DXDYDZ">'
            write (file,'(3A)'     )'   <DataItem Name="origin" Format="XML" DataType="Float" Precision="8" Dimensions="',trim(numb2String(nDim)),'">'
            if(nDim == 1) write (file,'(A,F25.10)'      )'    ', MSH%xMinGlob(1)
            if(nDim == 2) write (file,'(A,F25.10)'      )'    ', MSH%xMinGlob(2), ' ', MSH%xMinGlob(1)
            if(nDim == 3) write (file,'(A,F25.10)'      )'    ', MSH%xMinGlob(3), ' ', MSH%xMinGlob(2), ' ', MSH%xMinGlob(1)
            write (file,'(A)'      )'   </DataItem>'
            write (file,'(3A)'     )'   <DataItem Name="step" Format="XML" DataType="Float" Precision="8" Dimensions="',trim(numb2String(nDim)),'">'
            if(nDim == 1) write (file,'(A,F25.10)'      )'    ', MSH%xStep(1)
            if(nDim == 2) write (file,'(A,F25.10)'      )'    ', MSH%xStep(2), ' ', MSH%xStep(1)
            if(nDim == 3) write (file,'(A,F25.10)'      )'    ', MSH%xStep(3), ' ', MSH%xStep(2), ' ', MSH%xStep(1)
            write (file,'(A)'      )'   </DataItem>'
            write (file,'(A)'     )'     </Geometry>'

            !TODO
            do i = 1, 1
                write (file,'(3A)'     )'     <Attribute Name="RF1" Center="Node" AttributeType="Scalar">'
                write (file,'(A)'     )'       <DataItem Reference="XML">'
                write (file,'(A)'     )'         /Xdmf/Domain/DataItem[@Name="samples"]'
                write (file,'(A)'     )'       </DataItem>'
                write (file,'(A)'      )'     </Attribute>'
            end do

            write (file,'(A)'      )'   </Grid>' !END Writing the data of one subdomain
            write (file,'(A)'      )'  </Grid>' !Closes the Collection
            write (file,'(A)'      )' </Domain>'
            write (file,'(A)'      )'</Xdmf>'

            close(file)
        else if(rang == 0) then
            write(*,*) "No XDMF Model prepared for this dimension"
            write(*,*) "nDim = ", nDim
            write(*,*) "The file won't be created"
        end if

        write(get_fileId(),*) "------------END Writing result XMF file-----------------------";

    end subroutine write_pHDF5_Str_XMF

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_generation_spec(MSH, RDF, folderpath, name, timeVec)

        implicit none

        !INPUTS
        type(RF)   :: RDF
        type(MESH) :: MSH
        character (len=*)             , intent(in) :: folderpath, name;
        double precision, dimension(:), intent(in), optional :: timeVec

        !LOCAL
        integer :: fileId, code
        character (len=40) :: doubleFmt
        integer(kind=8) :: sum_xNTotal, sum_kNTotal

        call MPI_REDUCE (RDF%xNTotal,sum_xNTotal,1,MPI_INTEGER8,MPI_SUM,0,RDF%comm,code)
        call MPI_REDUCE (RDF%kNTotal,sum_kNTotal,1,MPI_INTEGER8,MPI_SUM,0,RDF%comm,code)

        if(RDF%rang == 0) then
            fileId = 15

            write(doubleFmt, fmt="(A,i1,A)") "(", RDF%nDim, "F30.15)";

            open (unit = fileId , file = string_vec_join([folderpath, "/", name]), action = 'write')

            write(fileId,*) "FILE:", name
            write(fileId,*) "--nb_procs-----------------------"
            write(fileId,fmt = "(I20)") RDF%nb_procs
            write(fileId,*) "--nDim-----------------------"
            write(fileId,fmt = "(I20)") RDF%nDim
            write(fileId,*) "--xMinGlob-----------------------"
            write(fileId,fmt = doubleFmt) MSH%xMinGlob
            write(fileId,*) "--xMaxGlob-----------------------"
            write(fileId,fmt = doubleFmt) MSH%xMaxGlob
            write(fileId,*) "--xMin_Loc-----------------------"
            write(fileId,fmt = doubleFmt) MSH%xMin
            write(fileId,*) "--xMax_Loc-----------------------"
            write(fileId,fmt = doubleFmt) MSH%xMax
            write(fileId,*) "--kMax_Loc-----------------------"
            write(fileId,fmt = doubleFmt) RDF%kMax
            write(fileId,*) "--xNTotal_Loc-----------------------"
            write(fileId,fmt = "(I20)") RDF%xNTotal
            write(fileId,*) "--kNTotal_Loc-----------------------"
            write(fileId,fmt = "(I20)") RDF%kNTotal
            write(fileId,*) "--xStep-----------------------"
            write(fileId,fmt = doubleFmt) MSH%xStep
            write(fileId,*) "--sum_xNTotal-----------------------"
            write(fileId,fmt = "(I20)") sum_xNTotal
            write(fileId,*) "--sum_kNTotal-----------------------"
            write(fileId,fmt = "(I20)") sum_kNTotal
            write(fileId,*) "--corrL-----------------------"
            write(fileId,fmt = doubleFmt) RDF%corrL
            write(fileId,*) "--corrMod-----------------------"
            write(fileId,*) RDF%corrMod
            write(fileId,*) "--margiFirst-----------------------"
            write(fileId,*) RDF%margiFirst
            write(fileId,*) "--fieldAvg-----------------------"
            write(fileId,fmt = doubleFmt) RDF%fieldAvg
            write(fileId,*) "--fieldVar-----------------------"
            write(fileId,fmt = doubleFmt) RDF%fieldVar
            write(fileId,*) "--Nmc-----------------------"
            write(fileId,fmt = "(I20)") RDF%Nmc
            write(fileId,*) "--method-----------------------"
            if(RDF%method == ISOTROPIC) write(fileId,*) "ISOTROPIC"
            if(RDF%method == SHINOZUKA) write(fileId,*) "SHINOZUKA"
            if(RDF%method == RANDOMIZATION) write(fileId,*) "RANDOMIZATION"
            write(fileId,*) "--independent-----------------------"
            if(MSH%overlap == -2.0D0) then
                write(fileId,*) .true. !Exception for monoproc cases
            else
                write(fileId,*) RDF%independent
            end if
            write(fileId,*) "--overlap-----------------------"
            write(fileId,*) MSH%overlap
            write(fileId,*) "--Seed-----------------------"
            write(fileId,fmt = "(I20)") RDF%seed
            if(present(timeVec)) then
                write(fileId,*) "--timeVec-----------------------"
                write(fileId,fmt = "(F30.15)") timeVec
            end if

            close(fileId)

        end if


    end subroutine write_generation_spec

end module writeResultFile_RF
