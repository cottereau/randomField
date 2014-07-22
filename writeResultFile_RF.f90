module writeResultFile_RF

    use displayCarvalhol
    use statistics_RF, only : calculateAverage
    use math_RF
    use hdf5

contains


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
        character (len=*), optional,       intent(in) :: filename;
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

!        if(present(filename)) path = "./Results/"//trim(filename);
!        if(.not.present(filename)) path = "./Results/Result0.txt";

        if(present(filename)) path = trim(filename);
        if(.not.present(filename)) path = "unnamedResult";

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

        write (file, *) "STATISTICS - GLOBAL";
        write (file, *) "";
        write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "average = ",   calculateAverage(average)
        write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "stdDeviation = ", calculateAverage(stdDeviation);
        write (file,"("//titleFmt//","//dimFmt//")") "avgCorrL = ", averageCorrL;

        if(present(time)) write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "Total time (s) = ", time;
        close(file)

        !>>>>>>>>>>>>>>>>>>>>
        file=11;
        open (unit = file , file = trim(path)//"MATLAB", action = 'write')

        write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)), i=1, size(randField,1));;

        close(file)

        deallocate(eventLabel)

        write(*,*) "";
        write(*,*) "------------END Writing result file-----------------------";
        write(*,*) "";

    end subroutine writeResultFileND

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine writeResultHDF5(xMax, kMax, xNStep, randField, fileName)
        implicit none

        !INPUTS
        double precision,  dimension(:),   intent(in) :: xMax, kMax;
        integer,           dimension(:),   intent(in) :: xNStep;
        double precision,  dimension(:,:), intent(in) :: randField;
        character (len=*),                 intent(in) :: filename;

        !HDF5 VARIABLES
        character(len=35)              :: fileHDF5Name !File name
        character(len=4),  parameter   :: dsetXYZ = "XYZ", dsetRF  = "RF"  !Dataset name
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: dspace_id     !Dataspace identifier
        integer                        :: rank = 2      !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(2) :: dims !Dataset dimensions
        integer                        :: error !Error flag

		!LOCAL VARIABLES
        double precision, dimension(:,:), allocatable ::grid_data
        integer :: nDim, Nmc, nPoints, i
        character (len=35) :: eventName;
        character (len=12) :: numberStr;

        write(*,*) "";
        write(*,*) "------------START Writing result HDF5 file-----------------------";
        write(*,*) "";

        fileHDF5Name = trim(fileName)//".h5"
        nDim         = size(xNStep , 1)
        nPoints      = size(randField, 1)
        Nmc          = size(randField, 2)

        write(*,*) "fileHDF5Name", fileHDF5Name

        if (nDim > 3) then
        	write(*,*) "Dimension exceeds 3, HDF file won't be created"

! START - TO BE MADE VERSION (for every Nmc > 1)
!        else
!
!			write(*,*) ">>>>>>>>> Opening file";
!	        call h5open_f(error) ! Initialize FORTRAN interface.
!	        call h5fcreate_f(fileHDF5Name, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.
!
!			write(*,*) ">>>>>>>>> Creating Coordinates dataset 'XYZ table'";
!
!			allocate (grid_data(3, nPoints)) !3 lines to put X, Y and Z
!	        grid_data = 0D0
!
!	        dims = shape(grid_data)
!	        write(*,*) "dims   = ", dims
!
!	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
!	        call h5dcreate_f(file_id, dsetXYZ, H5T_NATIVE_DOUBLE,  &
!	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
!
!
!
!			do i = 1, nPoints
!				call get_Permutation(i, xMax, xNStep, grid_data(1:nDim,i));
!			end do
!
!	 		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid_data, dims, error) ! Write the dataset.
!	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
!	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
!
!			deallocate (grid_data)
!
!			write(*,*) ">>>>>>>>> Creating Quantities dataset 'random field'";
!	        dims(1) = size(randField,1)
!	        dims(2) = 1 !One random field in each dataset
!			write(*,*) "dims   = ", dims;
!
!			do i = 1, Nmc
!				write(numberStr,'(I)'  ) i
!				numberStr = adjustL(numberStr)
!				write(eventName,'(A,A)') "RF_", numberStr
!				write(*,*) eventName
!
!		        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
!		        call h5dcreate_f(file_id, eventName, H5T_NATIVE_DOUBLE,  &
!		                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
!		        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField(:,i), dims, error) ! Write the dataset.
!		        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
!		        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
!		        write(*,*) "error = ", error
!		    end do
!
!			write(*,*) ">>>>>>>>> Closing file";
!	        call h5fclose_f(file_id, error) ! Close the file.
!	        call h5close_f(error) ! Close FORTRAN interface.
!
!	        call writeXMF_RF(randField, fileHDF5Name, fileName)
!
!        end if
! END - TO BE MADE VERSION (for every Nmc > 1)

! START - VERSION TESTED FOR Nmc = 1
        else

			write(*,*) ">>>>>>>>> Opening file";
	        call h5open_f(error) ! Initialize FORTRAN interface.
	        call h5fcreate_f(fileHDF5Name, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.

			write(*,*) ">>>>>>>>> Creating Coordinates dataset 'XYZ table'";

			allocate (grid_data(3, nPoints)) !3 lines to put X, Y and Z
	        grid_data = 0D0

	        dims = shape(grid_data)
	        write(*,*) "dims   = ", dims

	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
	        call h5dcreate_f(file_id, dsetXYZ, H5T_NATIVE_DOUBLE,  &
	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).



			do i = 1, nPoints
				call get_Permutation(i, xMax, xNStep, grid_data(1:nDim,i));
			end do

	 		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid_data, dims, error) ! Write the dataset.
	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.

			deallocate (grid_data)

			write(*,*) ">>>>>>>>> Creating Quantities dataset 'random field'";
	        dims = shape(randField)

			write(*,*) "dims   = ", dims;

	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
	        call h5dcreate_f(file_id, dsetRF, H5T_NATIVE_DOUBLE,  &
	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
	        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField, dims, error) ! Write the dataset.
	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.

			write(*,*) ">>>>>>>>> Closing file";
	        call h5fclose_f(file_id, error) ! Close the file.
	        call h5close_f(error) ! Close FORTRAN interface.

	        call writeXMF_RF(randField, fileHDF5Name, fileName)

        end if
! END - VERSION TESTED FOR Nmc = 1

        write(*,*) "";
        write(*,*) "------------END Writing result HDF5 file-----------------------";
        write(*,*) "";

    end subroutine writeResultHDF5

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine writeXMF_RF(randField, fileHDF5Name, fileName)
        implicit none

        !INPUTS
        double precision,  dimension(:,:), intent(in) :: randField;
        character (len=*),                 intent(in) :: filename;
        character (len=*),                 intent(in) :: fileHDF5Name;

		!LOCAL VARIABLES
        integer            :: Nmc, i, file, nElem
        character (len=35) :: fileXMFName;
        character (len=35) :: eventName;
        character (len=12) :: numberStr, nElemStr, NmcStr;

        write(*,*) "";
        write(*,*) "------------START Writing result XMF file-----------------------";
        write(*,*) "";

		fileXMFName = trim(fileName)//".xmf"
        Nmc     = size(randField, 2)
        nElem   = size(randField, 1)
        write(NmcStr,'(I)'  ) Nmc
		write(nElemStr,'(I)') nElem
		NmcStr   = adjustL(NmcStr)
		nElemStr = adjustL(nElemStr)

		write(*,*) "Nmc   = ", NmcStr
        write(*,*) "nElem = ", nElemStr
        write(*,*) "fileXMFName", fileXMFName

        file=21;

! START - TO BE MADE VERSION (for every Nmc > 1)
!        open (unit = file , file = trim(fileXMFName), action = 'write')
!
!			write (file,'(A)'    )'<?xml version="1.0" ?>'
!			write (file,'(A)'    )'<Xdmf Version="2.0">'
!			write (file,'(A)'    )' <Domain>'
!			write (file,'(A)'    )'  <Grid Name="meshRF" GridType="Uniform">'
!			write (file,'(A,A,A)')'   <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',trim(nElemStr),'">'
!			write (file,'(A)'    )'   </Topology>'
!			write (file,'(A)'    )'    <Geometry GeometryType="XYZ">'
!			write (file,'(A,A,A)')'     <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(nElemStr), ' 3">'
!			write (file,'(A,A,A)')'     	 ',trim(fileHDF5Name),':/XYZ'
!			write (file,'(A)'    )'     </DataItem>'
!			write (file,'(A)'    )'    </Geometry>'
!			write (file,'(A)'    )'    <Attribute Name="randomField" Center="Node" AttributeType="Scalar">'
!			write (file,'(A)'    )'     <DataItem Name="Random Field Tree" ItemType="Tree">'
!			do i = 1, Nmc
!				write(numberStr,'(I)'  ) i
!				numberStr = adjustL(numberStr)
!
!				write(eventName,'(A,A)') ":/RF_", numberStr
!				write (file,'(A,A,A)')'      <DataItem Dimensions="',trim(nElemStr),'" Format="HDF" DataType="Float" Precision="8">'
!				write (file,'(A,A,A)')'         ',trim(fileHDF5Name), eventName
!				write (file,'(A)'    )'      </DataItem>'
!			end do
!			write (file,'(A)'    )'     </DataItem>'
!
!			write (file,'(A)'    )'    </Attribute>'
!			write (file,'(A)'    )'  </Grid>'
!			write (file,'(A)'    )' </Domain>'
!			write (file,'(A)'    )'</Xdmf>'
!
!        close(file)
! END - TO BE MADE VERSION (for every Nmc > 1)

! START - VERSION TESTED FOR Nmc = 1
        open (unit = file , file = trim(fileXMFName), action = 'write')

			write (file,'(A)'    )'<?xml version="1.0" ?>'
			write (file,'(A)'    )'<Xdmf Version="2.0">'
			write (file,'(A)'    )' <Domain>'
			write (file,'(A)'    )'  <Grid Name="meshRF" GridType="Uniform">'
			write (file,'(A,I,A)')'   <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',nElem,'">'
			write (file,'(A)'    )'   </Topology>'
			write (file,'(A)'    )'    <Geometry GeometryType="XYZ">'
			write (file,'(A,I,A)')'     <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',nElem, ' 3">'
			write (file,'(A,A,A)')'     	 ',trim(fileHDF5Name),':/XYZ'
			write (file,'(A)'    )'     </DataItem>'
			write (file,'(A)'    )'    </Geometry>'
			write (file,'(A)'    )'    <Attribute Name="randomField" Center="Node" AttributeType="Scalar">'
			write (file,'(A,I,A)')'     <DataItem Dimensions="',nElem,'" Format="HDF" DataType="Float" Precision="8">'
			write (file,'(A,A,A)')'        ',trim(fileHDF5Name),':/RF'
			write (file,'(A)'    )'     </DataItem>'
			write (file,'(A)'    )'    </Attribute>'
			write (file,'(A)'    )'  </Grid>'
			write (file,'(A)'    )' </Domain>'
			write (file,'(A)'    )'</Xdmf>'

        close(file)
! END - VERSION TESTED FOR Nmc = 1

        write(*,*) "";
        write(*,*) "------------END Writing result XMF file-----------------------";
        write(*,*) "";

    end subroutine writeXMF_RF

end module writeResultFile_RF
