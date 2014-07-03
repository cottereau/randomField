module writeResultFile_RF

	use statistics_RF, only : calculateAverage

contains


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
    						   kMax, xNStep, kNStep, randField, &
    						   average, stdDeviation, averageCorrL, fileName, time)

        implicit none
        !INPUT
		character (len=*), dimension(:),   intent(in) :: corrMod;
		double precision,  dimension(:,:), intent(in) :: corrL, xMax, xPeriod, kMax;
		integer,           dimension(:,:), intent(in) :: xNStep, kNStep;
        double precision,  dimension(:,:), intent(in) :: randField;
		double precision,  dimension(:),   intent(in) :: average, stdDeviation, averageCorrL;
		character (len=*), optional,       intent(in) :: filename;
		double precision,  optional,       intent(in) :: time;

        !OUTPUT

        !LOCAL VARIABLES
        integer            :: i, j, file, nLines, nColumns, nDim;
        character (len=40) :: titleFmt, stringFmt, doubleFmt, intFmt, dimFmt;
        character (len=50) :: path
        character (len=20), dimension(:), allocatable :: eventLabel;

		write(*,*) "";
        write(*,*) "------------START Writing result file-----------------------";
		write(*,*) "";

		if(present(filename)) path = "./Results/"//trim(filename);
		if(.not.present(filename)) path = "./Results/Result0.txt";

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
        write (file,"("//titleFmt//","//doubleFmt//")") "corrL",   ((corrL(i,j),j=1,   size(corrL,2)),   i=1, size(corrL,1));;
        write (file,"("//titleFmt//","//doubleFmt//")") "xMax",    ((xMax(i,j),j=1,    size(xMax,2)),    i=1, size(xMax,1));
        write (file,"("//titleFmt//","//doubleFmt//")") "xPeriod", ((xPeriod(i,j),j=1, size(xPeriod,2)), i=1, size(xPeriod,1));

        write (file, *) "";
        write (file, *) "COMPUTED";
        write (file,"("//titleFmt//","//doubleFmt//")") "kMax",   ((kMax(i,j),j=1,   size(kMax,2)),  i=1, size(kMax,1));
        write (file,"("//titleFmt//","//intFmt   //")") "xNStep", ((xNStep(i,j),j=1, size(xNStep,2)),i=1, size(xNStep,1));
        write (file,"("//titleFmt//","//intFmt   //")") "kNStep", ((kNStep(i,j),j=1, size(kNStep,2)),i=1, size(kNStep,1));

        write (file, *) "";
        write (file, *) "RANDOM FIELD";
        write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)),i=1, size(randField,1));;

        write (file, *) "";
        write (file, *) "STATISTICS - POINT BY POINT";
        write (file,"("//titleFmt//","//stringFmt//")") "", "average", "stdDeviation";
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
        open (unit = file , file = path//"MATLAB", action = 'write')

		write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)), i=1, size(randField,1));;

        close(file)

		write(*,*) "";
        write(*,*) "------------END Writing result file-----------------------";
		write(*,*) "";

    end subroutine

end module writeResultFile_RF
