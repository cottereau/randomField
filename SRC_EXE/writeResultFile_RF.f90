module writeResultFile_RF

    use displayCarvalhol
    use statistics_RF
    use math_RF
    use hdf5
    use mpi
    use constants_RF
    use write_Log_File
    use type_RF
    use type_MESH
    use type_TEST

contains


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_ResultHDF5Unstruct_MPI(xPoints, randField, fileName, rang, folderPath, &
                                            communicator, labels, indexes, HDF5Name)
        implicit none

        !INPUTS
        double precision, dimension(1:,1:), intent(in) :: xPoints, randField;
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator
        character(len=*), dimension(1:), optional, intent(in) :: labels
        integer         , dimension(1:), optional, intent(in) :: indexes

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


        write(get_fileId(),*) "------------START Writing result HDF5 file (MPI)-----------------------";
        !write(get_fileId(),*) "lbound(xPoints) = ", lbound(xPoints)
        !write(get_fileId(),*) "lbound(randField) = ", lbound(randField)
        write(get_fileId(),*) "shape(xPoints)   = ", shape(xPoints)
        write(get_fileId(),*) "shape(randField) = ", shape(randField)
        write(get_fileId(),*) "fileName         = ", fileName
        write(get_fileId(),*) "folderPath       = ", folderPath


        effectComm = communicator
        nDim       = size(xPoints , 1)
        nPoints    = size(randField, 1)
        Nmc        = size(randField, 2)

        !Creating file name
        !write(*,*) "labels = ", labels
        !write(*,*) "indexes = ", indexes

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

        write(get_fileId(),*) "' fileHDF5Name = ", fileHDF5Name

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
            write(get_fileId(),*) "'inside write HDF5' output -- HDF5Name = ", HDF5Name
        end if

        write(get_fileId(),*) "------------END Writing result HDF5 file (MPI)-----------------------";

    end subroutine write_ResultHDF5Unstruct_MPI


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine writeXMF_RF_MPI(nSamples, HDF5nameList, nPointList, mask, nDim, fileName, rang, folderPath, &
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

    end subroutine writeXMF_RF_MPI

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
        kMax, xNStep, kNStep, randField, &
        average, stdDeviation, averageCorrL, fileName, time)

        implicit none
        !INPUT
        character (len=*)                             :: corrMod;
        double precision,  dimension(:),   intent(in) :: corrL, xMax, xPeriod, kMax;
        integer,           dimension(:),   intent(in) :: xNStep, kNStep;
        double precision,  dimension(:,:), intent(in) :: randField;
        double precision,  dimension(:),   intent(in) :: average, stdDeviation, averageCorrL;
        character (len=*)              ,   intent(in) :: filename;
        double precision,  optional,       intent(in) :: time;

        !OUTPUT - Result File

        !LOCAL VARIABLES
        integer            :: i, j, file, nLines, nColumns, nDim;
        character (len=40) :: titleFmt, stringFmt, doubleFmt, intFmt, dimFmt;
        character (len=50) :: path
        character (len=20), dimension(:), allocatable :: eventLabel;

        write(*,*) "";
        write(*,*) "------------START Writing result file-----------------------";
        write(*,*) "";

        path = trim(adjustL(filename));

        nLines   = size(randField,1);
        nColumns = size(randField,2);
        nDim     = size(xMax,     1);
        titleFmt = "A18"

        allocate(eventLabel(nColumns))

        write(stringFmt, fmt="(A,i10,A)") "(t25,", nColumns, "A25)";
        write(intFmt,    fmt="(A,i10,A)") "(t25,", nColumns, "I25)";
        write(doubleFmt, fmt="(A,i10,A)") "(t25,", nColumns, "F25.16)";
        write(dimFmt,    fmt="(A,i10,A)") "(t25,", nDim, "F25.16)";

        do i = 1, nColumns
            write(eventLabel(i), fmt="(A10,i10)") "Event", i;
        end do

        !>>>>>>>>>>>>>>>>>>>>
        file=10;
        open (unit = file , file = path, action = 'write')

        write (file, *) "DIMENSION = ", size(xMax,1);
        write (file, *) "";

        write (file,"("//titleFmt//","//stringFmt//")") "", eventLabel;

        write (file, *) "INPUTS";
        write (file,"("//titleFmt//","//stringFmt//")") "corrMod", adjustr(corrMod);
        write (file,"("//titleFmt//","//doubleFmt//")") "corrL",   corrL;
        write (file,"("//titleFmt//","//doubleFmt//")") "xMax",    xMax;
        write (file,"("//titleFmt//","//doubleFmt//")") "xPeriod", xPeriod;

        write (file, *) "";
        write (file, *) "COMPUTED";
        write (file,"("//titleFmt//","//doubleFmt//")") "kMax",   kMax;
        write (file,"("//titleFmt//","//intFmt   //")") "xNStep", xNStep;
        write (file,"("//titleFmt//","//intFmt   //")") "kNStep", kNStep;

        write(*,*) ">>>>>>>>> Writing Random Field";
        write (file, *) "";
        write (file, *) "RANDOM FIELD";
        write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)),i=1, size(randField,1));;

        write(*,*) ">>>>>>>>> Writing Statistics";
        write (file, *) "";
        write (file, *) "STATISTICS - POINT BY POINT";
        write (file,"("//titleFmt//","//stringFmt//","//stringFmt//")") "", "average", "stdDeviation";
        do i = 1, size(average)
            write (file,"("//titleFmt//","//"(t25,2F25.16)"//")") "", average(i), stdDeviation(i);
        end do

        !        write (file, *) "STATISTICS - GLOBAL";
        !        write (file, *) "";
        !        write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "average = ",   calculateAverage(average)
        !        write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "stdDeviation = ", calculateAverage(stdDeviation);
        !        write (file,"("//titleFmt//","//dimFmt//")") "avgCorrL = ", averageCorrL;

        if(present(time)) write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "Total time (s) = ", time;
        close(file)

    end subroutine writeResultFileND

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_Mono_XMF_h5(xPoints, randField, fileName, rang, folderPath, &
                                            communicator, labelsH5, indexesH5, indexXMF)

        implicit none

        !INPUTS
        double precision, dimension(1:,1:), intent(in) :: xPoints, randField;
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator, indexXMF
        character(len=*), dimension(1:)   , intent(in) :: labelsH5
        integer         , dimension(1:)   , intent(in) :: indexesH5

        !LOCAL
        character(len=110) :: HDF5Name, XMFName

        write(get_fileId(),*) "-> Writing h5 file in", trim(adjustL(folderPath))//"/h5";

        call write_ResultHDF5Unstruct_MPI(xPoints, randField, fileName, rang, trim(adjustL(folderPath))//"/h5", &
                                          MPI_COMM_WORLD, labelsH5, indexesH5, HDF5Name)


        write(get_fileId(),*) "-> Writing XMF file in", trim(adjustL(folderPath))//"/xmf";
        XMFName = stringNumb_join(trim(adjustL(fileName))//"it_", indexXMF)

        call writeXMF_RF_MPI(1, [HDF5name], [size(xPoints,2)], [.true.], size(xPoints,1), XMFName, &
                             rang, trim(adjustL(folderPath))//"/xmf", &
                             MPI_COMM_WORLD, "../h5")

    end subroutine write_Mono_XMF_h5

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_MatlabTable(randField, filename)

        implicit none
        !INPUT
        double precision,  dimension(:,:), intent(in) :: randField;
        character (len=*)              ,   intent(in) :: filename;

        !OUTPUT - Result File

        !LOCAL VARIABLES
        integer            :: i, j, file, nColumns;
        character (len=40) :: doubleFmt;
        character (len=50) :: path

        write(*,*) "";
        write(*,*) "------------START Writing result file-----------------------";
        write(*,*) "";

        path = trim(adjustL(filename));

        nColumns = size(randField,2);

        write(doubleFmt, fmt="(A,i10,A)") "(", nColumns, "F25.16)";

        !>>>>>>>>>>>>>>>>>>>>
        file=11;
        open (unit = file , file = trim(path)//"MATLAB", action = 'write')
        write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)), i=1, size(randField,1));
        close(file)

    end subroutine write_MatlabTable

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

        !call show_RF(RDF, "RDF")

        !call MPI_REDUCE (RDF%xNTotal,sum_xNTotal,1,MPI_INT,MPI_SUM,0,RDF%comm,code)
        !call MPI_REDUCE (RDF%kNTotal,sum_kNTotal,1,MPI_INT,MPI_SUM,0,RDF%comm,code)

        call MPI_REDUCE (RDF%xNTotal,sum_xNTotal,1,MPI_INTEGER8,MPI_SUM,0,RDF%comm,code)
        call MPI_REDUCE (RDF%kNTotal,sum_kNTotal,1,MPI_INTEGER8,MPI_SUM,0,RDF%comm,code)

        !write(*,*) "RDF%rang = ", RDF%rang, "", "RDF%xNTotal = ", RDF%xNTotal, "RDF%kNTotal = ", RDF%kNTotal

        if(RDF%rang == 0) then
            fileId = 15

            write(doubleFmt, fmt="(A,i1,A)") "(", RDF%nDim, "F30.15)";
            !write(*,*) "doubleFmt = ", doubleFmt

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
            write(fileId,*) RDF%independent
            write(fileId,*) "--overlap-----------------------"
            write(fileId,*) MSH%overlap
            write(fileId,*) "--Seed-----------------------"
            write(fileId,fmt = "(I20)") RDF%seed
            if(present(timeVec)) then
                write(fileId,*) "--timeVec-----------------------"
                write(fileId,fmt = "(F30.15)") timeVec
            end if

!            if(present(avg_Gauss)) then
!                write(fileId,*) "--avg_Gauss-----------------------"
!                write(fileId,fmt = "(F30.15)") avg_Gauss
!            end if
!            if(present(stdDev_Gauss)) then
!                write(fileId,*) "--stdDev_Gauss-----------------------"
!                write(fileId,fmt = "(F30.15)") stdDev_Gauss
!            end if
!            if(present(avg_Trans)) then
!                write(fileId,*) "--avg_Trans-----------------------"
!                write(fileId,fmt = "(F30.15)") avg_Trans
!            end if
!            if(present(stdDev_Trans)) then
!                write(fileId,*) "--stdDev_Trans-----------------------"
!                write(fileId,fmt = "(F30.15)") stdDev_Trans
!            end if
!
!            if(present(avg_Gauss_evnt)) then
!                write(fileId,*) "--avg_Gauss_evnt-----------------------"
!                write(fileId,fmt = "(F30.15)") avg_Gauss_evnt
!            end if
!            if(present(stdDev_Gauss_evnt)) then
!                write(fileId,*) "--stdDev_Gauss_evnt-----------------------"
!                write(fileId,fmt = "(F30.15)") stdDev_Gauss_evnt
!            end if
!            if(present(avg_Trans_evnt)) then
!                write(fileId,*) "--avg_Trans_evnt-----------------------"
!                write(fileId,fmt = "(F30.15)") avg_Trans_evnt
!            end if
!            if(present(stdDev_Trans_evnt)) then
!                write(fileId,*) "--stdDev_Trans_evnt-----------------------"
!                write(fileId,fmt = "(F30.15)") stdDev_Trans_evnt
!            end if
!
!            if(present(avg_Gauss_point)) then
!                write(fileId,*) "--avg_Gauss_point-----------------------"
!                write(fileId,fmt = "(F30.15)") avg_Gauss_point
!            end if
!            if(present(stdDev_Gauss_point)) then
!                write(fileId,*) "--stdDev_Gauss_point-----------------------"
!                write(fileId,fmt = "(F30.15)") stdDev_Gauss_point
!            end if
!            if(present(avg_Trans_point)) then
!                write(fileId,*) "--avg_Trans_point-----------------------"
!                write(fileId,fmt = "(F30.15)") avg_Trans_point
!            end if
!            if(present(stdDev_Trans_point)) then
!                write(fileId,*) "--stdDev_Trans_point-----------------------"
!                write(fileId,fmt = "(F30.15)") stdDev_Trans_point
!            end if

            close(fileId)

        end if


    end subroutine write_generation_spec

end module writeResultFile_RF
