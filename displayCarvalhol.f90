module displayCarvalhol

    !All display routines
    interface DispCarvalhol
		module procedure Disp1Ddble,  &
		                 Disp2Ddble,  &
		                 Disp1Dint,   &
		                 Disp2Dint,   &
		                 Disp1Dchar,  &
		                 Disp2Dchar
	end interface DispCarvalhol

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp1Ddble(vector, title, format, nColumns)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        double precision,     dimension(:),          intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns

        !LOCAL VARIABLES
        double precision, dimension(:,:), allocatable :: matrix2d

		allocate(matrix2D(size(vector),1));
		matrix2d(:,1) = vector;
		call Disp2Ddble(matrix2D, title, format, nColumns);
		deallocate(matrix2D);

    end subroutine Disp1Ddble

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp1Dint(vector, title, format, nColumns)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        integer,            dimension(:),            intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns

        !LOCAL VARIABLES
        integer, dimension(:,:), allocatable :: matrix2D

		allocate(matrix2D(size(vector),1));
		matrix2d(:,1) = vector;
		call Disp2Dint(matrix2D, title, format, nColumns);
		deallocate(matrix2D);

    end subroutine Disp1Dint

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp1Dchar(vector, title, format, nColumns)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        character(len=*),   dimension(:),            intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns

        !LOCAL VARIABLES
        character(len=30), dimension(:,:), allocatable :: matrix2D

		allocate(matrix2D(size(vector),1));
		matrix2d(:,1) = vector;
		call Disp2Dchar(matrix2D, title, format, nColumns);
		deallocate(matrix2D);

    end subroutine Disp1Dchar

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp2Ddble(matrix2D, title, format, nColumns)
        ! Displays 2D Matrix, "div" columns at a time

		implicit none

		!INPUT
        double precision, dimension(:, :),           intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
		character (len=40) :: doubleFmt;
		character (len=10) :: tempFmt;

		if(present(format))                     tempFmt = format;
		if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
		if(.not.present(format))    tempFmt = "F10.5";
		if(.not.present(nColumns))  div = 5;
		write(doubleFmt, fmt="(I3, A)") div, tempFmt;

		write(*,*) ""
		if(present(title)) write(*,*) "/////// ", title, " ///////";

		tol = 1E-10;

        write(*,*) ""
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)
    end subroutine Disp2Ddble

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp2Dint(matrix2D, title, format, nColumns)
        ! Displays 2D Matrix, "div" columns at a time

		implicit none

		!INPUT
        integer,           dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
		character (len=40) :: doubleFmt;
		character (len=10) :: tempFmt;

		if(present(format))                     tempFmt = format;
		if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
		if(.not.present(format))    tempFmt = "i6";
		if(.not.present(nColumns))  div = 5;
		write(doubleFmt, fmt="(I3, A)") div, tempFmt;

		write(*,*) ""
		if(present(title)) write(*,*) "/////// ", title, " ///////";

		tol = 1E-10;

        write(*,*) ""
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)
    end subroutine Disp2Dint

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp2Dchar(matrix2D, title, format, nColumns)
        ! Displays 2D Matrix, "div" columns at a time

		implicit none

		!INPUT
        character (len=*), dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
		character (len=40) :: charFmt;
		character (len=10) :: tempFmt;

		if(present(format))                     tempFmt = format;
		if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
		if(.not.present(format))    tempFmt = "A";
		if(.not.present(nColumns))  div = 5;
		write(charFmt, fmt="(I3, A)") div, tempFmt;

		write(*,*) ""
		if(present(title)) write(*,*) "/////// ", title, " ///////";

		tol = 1E-10;

        write(*,*) ""
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//charFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//charFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)
    end subroutine Disp2Dchar

end module displayCarvalhol
